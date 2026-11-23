"""
QGIS port of the ValleyMorph Tool (river_transect.py, R. Dalumpines 2012-2014).

Each function mirrors the function of the same name in the original ArcGIS
script and applies the same steps in the same order, so the two can be read
side by side.  ArcGIS geoprocessing calls are replaced by numpy (rasters, see
grid.py / hydro.py) and QgsGeometry / GEOS (vectors).
"""

import math

import numpy as np
from qgis.core import (
    QgsCoordinateTransform,
    QgsDistanceArea,
    QgsGeometry,
    QgsPointXY,
    QgsProject,
)

from . import grid as G
from . import hydro as H


class VMTSkip(ValueError):
    """Equivalent of the IOError/ValueError raised in the original: skip this watershed."""


class Feat:
    __slots__ = ("geom", "a")

    def __init__(self, geom, attrs=None):
        self.geom = geom
        self.a = dict(attrs or {})

    def copy(self, geom=None):
        return Feat(self.geom if geom is None else geom, self.a)


class Context:
    """Run-wide settings (CRS, XY tolerance, messages)."""

    def __init__(self, crs, feedback=None, transform_context=None):
        self.crs = crs
        self.feedback = feedback
        self.tc = transform_context or QgsProject.instance().transformContext()
        # ArcGIS default XY tolerance: 0.001 m (or its equivalent in degrees).
        self.tol = 1e-8 if (crs.isValid() and crs.isGeographic()) else 0.001
        self.da = None
        self.to_geo = self.from_geo = None
        if crs.isValid() and crs.ellipsoidAcronym() and not crs.isGeographic():
            geo = crs.toGeographicCrs()
            self.to_geo = QgsCoordinateTransform(crs, geo, self.tc)
            self.from_geo = QgsCoordinateTransform(geo, crs, self.tc)
        if crs.isValid() and crs.ellipsoidAcronym():
            self.da = QgsDistanceArea()
            self.da.setSourceCrs(crs.toGeographicCrs() if not crs.isGeographic() else crs, self.tc)
            self.da.setEllipsoid(crs.ellipsoidAcronym())

    def msg(self, text):
        if self.feedback:
            self.feedback.pushInfo(text)

    def warn(self, text):
        if self.feedback:
            self.feedback.pushWarning(text)

    def label(self, text):
        if self.feedback:
            self.feedback.setProgressText(text)
            self.feedback.pushDebugInfo(text)


# ---------------------------------------------------------------------------
# Small geometry helpers
# ---------------------------------------------------------------------------

def _pt(x, y):
    return QgsGeometry.fromPointXY(QgsPointXY(x, y))


def _paths(geom):
    """List of vertex lists [(x, y), ...] for each part of a line geometry."""
    if geom is None or geom.isEmpty():
        return []
    if geom.isMultipart():
        return [[(p.x(), p.y()) for p in part] for part in geom.asMultiPolyline()]
    return [[(p.x(), p.y()) for p in geom.asPolyline()]]


def _line(pts):
    return QgsGeometry.fromPolylineXY([QgsPointXY(x, y) for x, y in pts])


def _single_path(geom):
    """Vertex list of a line; multipart lines are merged first."""
    paths = _paths(geom)
    if len(paths) > 1:
        paths = _paths(geom.mergeLines())
    out = []
    for p in paths:
        if out and out[-1] == p[0]:
            out.extend(p[1:])
        else:
            out.extend(p)
    return out


def _lines_only(geom):
    """Keep only the linear parts of an intersection result."""
    if geom is None or geom.isEmpty():
        return None
    parts = []
    for g in geom.asGeometryCollection() if geom.isMultipart() or geom.type() != 1 else [geom]:
        if g.type() == 1 and not g.isEmpty():
            parts.extend(_paths(g))
    parts = [p for p in parts if len(p) >= 2]
    if not parts:
        return None
    if len(parts) == 1:
        return _line(parts[0])
    return QgsGeometry.fromMultiPolylineXY([[QgsPointXY(x, y) for x, y in p] for p in parts])


def _explode(feats):
    """MultipartToSinglepart."""
    out = []
    for f in feats:
        for p in _paths(f.geom):
            out.append(f.copy(_line(p)))
    return out


def _near(ctx, a, b):
    return a.distance(b) <= ctx.tol


def _near_angle(x, y, tx, ty):
    """ArcGIS Near ANGLE: from input to near location, east = 0, counter-clockwise, -180..180."""
    return math.degrees(math.atan2(ty - y, tx - x))


def _endpoints(geom):
    paths = _paths(geom)
    return paths[0][0], paths[-1][-1]


def _split_line_at_points(ctx, f, points):
    """SplitLineAtPoint for one single-part line feature."""
    path = _paths(f.geom)[0]
    g = f.geom
    L = g.length()
    cuts = []
    for p in points:
        if g.distance(p) <= ctx.tol:
            d = g.lineLocatePoint(p)
            if ctx.tol < d < L - ctx.tol:
                cuts.append(d)
    if not cuts:
        return [f.copy()]
    cuts = sorted(set(round(c, 9) for c in cuts))
    ls = g.constGet()
    out = []
    prev = 0.0
    for c in cuts + [L]:
        if c - prev > ctx.tol:
            sub = ls.curveSubstring(prev, c)
            out.append(f.copy(QgsGeometry(sub)))
        prev = c
    return out


def _unsplit(feats, field="MEASURE"):
    """UnsplitLine by field: merge lines of the same value that share end points."""
    groups = {}
    order = []
    for f in feats:
        k = f.a.get(field)
        if k not in groups:
            groups[k] = []
            order.append(k)
        groups[k].append(f)
    out = []
    for k in order:
        merged = QgsGeometry.unaryUnion([f.geom for f in groups[k]]).mergeLines()
        for p in _paths(merged):
            out.append(Feat(_line(p), {field: k}))
    return out


def _spatial_join_first(ctx, targets, joins, keep_fields):
    """
    SpatialJoin JOIN_ONE_TO_ONE / KEEP_COMMON / INTERSECT with default field
    mapping: target attributes win, join attributes come from the first
    intersecting join feature.
    """
    out = []
    for t in targets:
        hit = None
        for j in joins:
            if t.geom.distance(j.geom) <= ctx.tol:
                hit = j
                break
        if hit is None:
            continue
        a = dict(hit.a)
        a.update(t.a)
        out.append(Feat(t.geom, {k: a.get(k) for k in keep_fields if k in a}))
    return out


# ---------------------------------------------------------------------------
# Ported functions
# ---------------------------------------------------------------------------

def get_coordpriority(in_raster, lowpoint):
    xcoord, ycoord = lowpoint
    key1 = "UPPER" if abs(in_raster.top - ycoord) < abs(in_raster.bottom - ycoord) else "LOWER"
    key2 = "LEFT" if abs(in_raster.left - xcoord) < abs(in_raster.right - xcoord) else "RIGHT"
    return "%s_%s" % (key1, key2)


def get_maxbound(in_raster):
    return int(max(abs(in_raster.top - in_raster.bottom), abs(in_raster.right - in_raster.left)))


class Hydro:
    """Caches FlowDirection / FlowAccumulation for the reduced DEM."""

    def __init__(self, dem, fill=False):
        self.dem = dem
        self.valid = dem.valid
        self.fdir = H.flow_direction(H.fill_sinks(dem.data) if fill else dem.data)
        self.acc, self.down, self.order = H.flow_accumulation(self.fdir, self.valid)


def get_hiflowaccpoint(ctx, hyd):
    acc = hyd.acc
    mx = np.nanmax(acc)
    pts = G.raster_to_points(hyd.dem, (acc == mx) & hyd.valid)
    if len(pts) > 1:
        ctx.warn("%i cells share the highest flow accumulation; all are kept as lowpoints" % len(pts))
    r, c = hyd.dem.cell_of(pts[0][0], pts[0][1])
    if hyd.fdir[r, c] == 0:
        ctx.warn("the highest flow accumulation is in a pit inside the watershed; if this watershed is "
                 "skipped or its river looks wrong, try the 'Fill sinks' option")
    return [(x, y) for x, y, _, _ in pts]


def get_farthestpoint(ctx, dem):
    minv = dem.min()
    src = (dem.data == minv) & dem.valid
    eucd = G.euclidean_distance(dem, src)
    mx = np.nanmax(eucd.data)
    pts = G.raster_to_points(dem, (eucd.data == mx) & dem.valid)
    if len(pts) > 1:
        ctx.warn("%i cells are equally far from the lowest cell; all are kept as farpoints" % len(pts))
    return [(x, y) for x, y, _, _ in pts]


def extract_mainstream(ctx, hyd, wid, order="LOW_HIGH", simplify_tol=None):
    """Main stream network from Strahler order computed on the flow accumulation raster."""
    # The flow accumulation raster is the stream raster; cells with value 0 are background.
    stream = (hyd.valid & (hyd.acc > 0)).ravel()
    so = H.stream_order_strahler(stream, hyd.down, hyd.order).reshape(hyd.fdir.shape)
    instream = stream.reshape(hyd.fdir.shape)
    vals = so[instream]
    if order == "HIGH_ONLY":
        rivergrid = instream & (so >= vals.max())
    elif order == "LOW_HIGH":
        rivergrid = instream & (so >= vals.min() + 1)
    elif order == "ALL":
        rivergrid = instream & (so >= vals.min())
    else:
        raise VMTSkip('stream order must be specified either "HIGH_ONLY", "LOW_HIGH", or "ALL"')
    if not rivergrid.any():
        raise VMTSkip("rivergrid is empty, uncheck option for edge effects then try again")
    if simplify_tol is None:
        simplify_tol = 0.5 * hyd.dem.cs
    arcs = H.stream_to_feature(hyd.dem, rivergrid, hyd.fdir, simplify_tol)
    feats = []
    for k, (pts, s, e) in enumerate(arcs):
        feats.append(Feat(_line(pts), {"arcid": k + 1, "from_node": s, "to_node": e, "WID": wid}))
    return feats


def _node_key(xy):
    return (round(xy[0], 6), round(xy[1], 6))


def _dangles(feats):
    deg = {}
    for f in feats:
        s, e = _endpoints(f.geom)
        for p in (s, e):
            deg[_node_key(p)] = deg.get(_node_key(p), 0) + 1
    out = []
    for i, f in enumerate(feats):
        s, e = _endpoints(f.geom)
        for p in (s, e):
            if deg[_node_key(p)] == 1:
                out.append(p)
    return out


def get_nearpoint(points_xy, targets_xy):
    """Nearest (first in order) vertex to any target point."""
    best, bd = None, None
    for p in points_xy:
        d = min(math.hypot(p[0] - t[0], p[1] - t[1]) for t in targets_xy)
        if bd is None or d < bd:
            best, bd = p, d
    return best


def get_longest_riverline(ctx, riverline, sinkpoints, sourcepoints, interval_profile, wid):
    """Trace the network path between the river vertex nearest the lowpoint and
    the dangling vertex nearest the farpoint (geometric network FIND_PATH)."""
    dangles = _dangles(riverline)
    if not dangles:
        raise VMTSkip("riverline is an empty feature, failed to generate vertices")
    flag1 = get_nearpoint(dangles, sinkpoints)

    adj = {}
    for i, f in enumerate(riverline):
        s, e = _endpoints(f.geom)
        a, b = _node_key(s), _node_key(e)
        adj.setdefault(a, []).append((b, i))
        adj.setdefault(b, []).append((a, i))

    start = _node_key(flag1)
    parent = {start: (None, None)}
    queue = [start]
    while queue:
        n = queue.pop(0)
        for m, i in adj.get(n, []):
            if m not in parent:
                parent[m] = (n, i)
                queue.append(m)
    connected_ids = sorted({i for n, (p, i) in parent.items() if i is not None})
    connected = [riverline[i] for i in connected_ids]
    if interval_profile > sum(f.geom.length() for f in connected):
        raise VMTSkip("riverline_connected too short for the river profile interval")

    flag2 = get_nearpoint(_dangles(connected), sourcepoints)
    node = _node_key(flag2)
    path_ids = []
    while parent[node][0] is not None:
        node, i = parent[node]
        path_ids.append(i)
    ctx.msg("flags consisting of source and sink points created ...")
    # path_ids runs from the source flag toward the sink flag (downstream).
    pts = []
    cur = _node_key(flag2)
    for i in path_ids:
        p = _paths(riverline[i].geom)[0]
        if _node_key(p[0]) != cur:
            p = p[::-1]
        pts.extend(p if not pts else p[1:])
        cur = _node_key(p[-1])
    longest = Feat(_line(pts), {"WID": wid})
    if interval_profile > longest.geom.length():
        raise VMTSkip("main (longest) river too short for the river profile interval")
    ctx.msg("longest_riverline created ...")
    return longest


def get_reducedinputs(localdem, watershedpolygon, adjust_edge, reduceby=2):
    if adjust_edge:
        outline = QgsGeometry(watershedpolygon.constGet().boundary())
        buff = outline.buffer(localdem.cs * reduceby, 32)
        adj = watershedpolygon.difference(buff)
        if adj.isEmpty():
            raise VMTSkip("watershed is too small to reduce for edge effects")
        return G.extract_by_mask(localdem, adj), adj
    return localdem, watershedpolygon


def extract_valleyfloor(ctx, in_raster, riverline, wid, pthreshold=12.5):
    """Flat cells (slope <= 12.5 %) plus river cells, kept where they touch the
    river, merged with a half-cell river buffer."""
    percentslope = G.slope_percent(in_raster)
    valid = in_raster.valid
    with np.errstate(invalid="ignore"):
        steep = percentslope.data > pthreshold           # vafloor1 = IsNull(flatslope)
    river = G.polyline_to_raster_mask(in_raster, [riverline.geom])  # riverlinegrid
    vafloor12 = steep & ~river                           # vafloor1 * vafloor2
    valleygrid = valid & ~vafloor12                      # SetNull(vafloor12 == 1)
    polys = G.raster_to_polygons(in_raster, valleygrid, simplify=True)
    # update_valleyfloor
    sel = [p for p in polys if p.intersects(riverline.geom) or p.distance(riverline.geom) <= ctx.tol]
    buf = riverline.geom.buffer(in_raster.cs / 2.0, 32)
    vf = QgsGeometry.unaryUnion(sel + [buf])
    if not vf.isMultipart():
        vf.convertToMultiType()
    return Feat(vf, {"WID": wid})


def get_longaxis(lowpoints, farpoints, wid):
    pts = list(lowpoints) + list(farpoints)
    if len(pts) < 2:
        raise VMTSkip("Not enough points to create a line for WID: %i" % wid)
    return Feat(_line(pts), {"WID": wid})


def create_intervaltable(length, interval, count_interval=None):
    tcount = int(math.ceil(length / interval))
    if count_interval and (count_interval <= tcount - 1):
        intervals = count_interval + 1
    elif count_interval and (count_interval > tcount - 1):
        raise VMTSkip("requested number of midline intervals exceeds the length of longaxis")
    elif tcount < 2:
        raise VMTSkip("midline interval exceeds the length of longaxis")
    else:
        intervals = tcount
    return [interval * x for x in range(1, intervals)]


def get_axisroute(in_raster, lowpoint, line_geom):
    """CreateRoutes (LENGTH) with the coordinate priority corner nearest the lowpoint:
    measures start at the route end closest to that corner of the route's extent."""
    pts = _single_path(line_geom)
    cp = get_coordpriority(in_raster, lowpoint)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    cx = min(xs) if cp.endswith("LEFT") else max(xs)
    cy = max(ys) if cp.startswith("UPPER") else min(ys)
    ds = math.hypot(pts[0][0] - cx, pts[0][1] - cy)
    de = math.hypot(pts[-1][0] - cx, pts[-1][1] - cy)
    if de < ds:
        pts = pts[::-1]
    return _line(pts)


def create_intervalpoints(in_raster, lowpoints, line, wid, interval, count_interval):
    route = get_axisroute(in_raster, lowpoints[0], line.geom)
    length = route.length()
    maxbound = get_maxbound(in_raster)
    out = []
    for m in create_intervaltable(length, interval, count_interval):
        p = route.interpolate(m).asPoint()
        out.append(Feat(_pt(p.x(), p.y()), {
            "WID": wid, "MEASURE": float(m), "DISTANCE": float(maxbound),
            "RBEARING": None, "LBEARING": None,
            "POINT_X": p.x(), "POINT_Y": p.y()}))
    return out


def calculate_bearings(points, anglefield="NEAR_ANGLE"):
    for f in points:
        angle = float(f.a.get(anglefield) or 0.0)
        f.a["RBEARING"] = (90 - angle) - 90  # right bearing
        f.a["LBEARING"] = (90 - angle) + 90  # left bearing


def create_axismidpoints(in_raster, lowpoints, longaxis, wid, interval, count_interval):
    pts = create_intervalpoints(in_raster, lowpoints, longaxis, wid, interval, count_interval)
    # append_axisbearings: Near to lowpoint (radius = (count + 1) * interval)
    radius = (len(pts) + 1) * interval
    for f in pts:
        x, y = f.a["POINT_X"], f.a["POINT_Y"]
        best = min(lowpoints, key=lambda t: math.hypot(t[0] - x, t[1] - y))
        if math.hypot(best[0] - x, best[1] - y) <= radius:
            f.a["NEAR_ANGLE"] = _near_angle(x, y, best[0], best[1])
        else:
            f.a["NEAR_ANGLE"] = 0.0
    calculate_bearings(pts)
    return pts


def get_segments_with_bearings(line, lowpoints):
    """Dissolve, split at vertices and number the segments from the sink end of
    the line (the end nearest the lowpoint); ANGLE is the Near angle from each
    vertex to the next."""
    pts = _single_path(line.geom)
    lx, ly = lowpoints[0]
    if math.hypot(pts[-1][0] - lx, pts[-1][1] - ly) < math.hypot(pts[0][0] - lx, pts[0][1] - ly):
        pts = pts[::-1]
    segs = []
    for k in range(len(pts) - 1):
        (x1, y1), (x2, y2) = pts[k], pts[k + 1]
        if (x1, y1) == (x2, y2):
            continue
        segs.append(Feat(_line([pts[k], pts[k + 1]]),
                         {"SEQUENCE": len(segs) + 1, "ANGLE": _near_angle(x1, y1, x2, y2)}))
    return segs


def create_csectionpoints(ctx, in_raster, lowpoints, line, wid, interval, count_interval):
    intervalpoints = create_intervalpoints(in_raster, lowpoints, line, wid, interval, count_interval)
    rsegments = get_segments_with_bearings(line, lowpoints)
    keep = ["MEASURE", "WID", "POINT_X", "POINT_Y", "DISTANCE", "RBEARING", "LBEARING", "ANGLE"]
    out = _spatial_join_first(ctx, intervalpoints, rsegments, keep)
    calculate_bearings(out, "ANGLE")
    return out


def get_bearingtoline(ctx, points, bearing_field):
    """BearingDistanceToLine (GEODESIC, METERS): north = 0, clockwise."""
    out = []
    for f in points:
        x, y = f.a["POINT_X"], f.a["POINT_Y"]
        b = f.a[bearing_field]
        d = f.a["DISTANCE"]
        ex, ey = _project_geodesic(ctx, x, y, b, d)
        out.append(Feat(_line([(x, y), (ex, ey)]), {
            "POINT_X": x, "POINT_Y": y, "DISTANCE": d, bearing_field: b,
            "MEASURE": f.a["MEASURE"]}))
    return out


def _project_geodesic(ctx, x, y, bearing, dist):
    az = math.radians(bearing % 360.0)
    if ctx.da is None:
        return x + dist * math.sin(az), y + dist * math.cos(az)
    p = QgsPointXY(x, y)
    if ctx.to_geo:
        p = ctx.to_geo.transform(p)
    q = ctx.da.computeSpheroidProject(p, dist, az)
    if ctx.from_geo:
        q = ctx.from_geo.transform(q)
    return q.x(), q.y()


def _clip(feats, polygon):
    out = []
    for f in feats:
        g = _lines_only(f.geom.intersection(polygon))
        if g is not None:
            out.append(f.copy(g))
    return out


def create_sidefield(in_raster, lowpoint, feats, sidefield="SIDE"):
    sidekey = get_coordpriority(in_raster, lowpoint).split("_")[1].upper()
    fieldname = "LBEARING" if sidekey == "RIGHT" else "RBEARING"
    for f in feats:
        f.a[sidefield] = "RIGHT" if f.a.get(fieldname) is not None else "LEFT"


def dissolve_axismidlines(feats):
    """Dissolve by MEASURE, SINGLE_PART, DISSOLVE_LINES."""
    groups = {}
    for f in feats:
        groups.setdefault(f.a["MEASURE"], []).append(f.geom)
    out = []
    for m in sorted(groups):
        merged = QgsGeometry.unaryUnion(groups[m]).mergeLines()
        for p in _paths(merged):
            out.append(Feat(_line(p), {"MEASURE": m}))
    return out


def remove_duplicate_side2(feats):
    out = []
    order = sorted(range(len(feats)), key=lambda i: (feats[i].a["MEASURE"], -feats[i].geom.length()))
    pid, plen = None, 0
    keep = set()
    for i in order:
        cid = feats[i].a["MEASURE"]
        clen = feats[i].geom.length()
        if pid == cid:
            if not plen > clen:
                raise VMTSkip("dissolved lines not sorted by length properly")
        else:
            pid, plen = cid, clen
            keep.add(i)
    return [f for i, f in enumerate(feats) if i in keep]


def dissolve_axismidlines2(feats):
    return remove_duplicate_side2(dissolve_axismidlines(feats))


def remove_danglinglines2(axismidlines):
    diss = dissolve_axismidlines2(axismidlines)
    keep = []
    for f in axismidlines:
        # SHARE_A_LINE_SEGMENT_WITH any kept dissolved line
        for d in diss:
            inter = f.geom.intersection(d.geom)
            if inter is not None and not inter.isEmpty() and inter.length() > 0:
                keep.append(f)
                break
    return keep


def get_axismidlines2(ctx, in_raster, axismidpoints, watershedpolygon, wid, lowpoints, crosstype):
    rightmidlines = get_bearingtoline(ctx, axismidpoints, "RBEARING")
    leftmidlines = get_bearingtoline(ctx, axismidpoints, "LBEARING")
    out_merge = rightmidlines + leftmidlines
    ctx.msg("bearing lines created and merged ...")
    axismidlines = remove_danglinglines2(_clip(out_merge, watershedpolygon))
    ctx.msg("dangling lines removed ...")
    for f in axismidlines:
        f.a["WID"] = wid
    create_sidefield(in_raster, lowpoints[0], axismidlines, "SIDE")
    for f in axismidlines:
        f.a["PERTO"] = crosstype
    return axismidlines


def remove_danglinglines(ctx, crosslines, line, sidefield):
    mts = _explode(crosslines)
    radius = line.geom.length()
    for f in mts:
        d = f.geom.distance(line.geom)
        f.a["NEAR_DIST"] = d if d <= radius else -1
    order = sorted(range(len(mts)), key=lambda i: (mts[i].a["MEASURE"], mts[i].a[sidefield], mts[i].a["NEAR_DIST"]))
    keep = set()
    m, p = 0, ""
    for i in order:
        f = mts[i]
        if f.a["NEAR_DIST"] == -1:
            raise VMTSkip("nearest distance for some crosslines cannot be calculated")
        n, q = f.a["MEASURE"], f.a[sidefield]
        if not (m == n and q == p):
            keep.add(i)
        m, p = n, q
    out = []
    for i, f in enumerate(mts):
        if i in keep:
            del f.a["NEAR_DIST"]
            out.append(f)
    return out


def get_axismidlines(ctx, in_raster, points, watershedpolygon, wid, line, lowpoints, crosstype, sidefield="SIDE"):
    rightmidlines = get_bearingtoline(ctx, points, "RBEARING")
    leftmidlines = get_bearingtoline(ctx, points, "LBEARING")
    out_clip = _clip(rightmidlines + leftmidlines, watershedpolygon)
    ctx.msg("bearing lines created and merged ...")
    for f in out_clip:
        f.a["WID"] = wid
    create_sidefield(in_raster, lowpoints[0], out_clip, sidefield)
    crosslines = remove_danglinglines(ctx, out_clip, line, sidefield)
    ctx.msg("dangling lines removed ...")
    for f in crosslines:
        f.a["PERTO"] = crosstype
    return crosslines


def get_watershedmidpoints(axisdissolve, lowpoints, farpoints):
    # FeatureToPoint INSIDE: a point on each line, half way along it.
    mids = [Feat(f.geom.interpolate(f.geom.length() / 2.0), {"MEASURE": f.a["MEASURE"]}) for f in axisdissolve]
    linecount = len(mids)
    rows = [Feat(_pt(*p), {"ENDPOINT": 1, "MEASURE": None}) for p in lowpoints] + mids + \
           [Feat(_pt(*p), {"ENDPOINT": 2, "MEASURE": None}) for p in farpoints]
    order = sorted(range(len(rows)), key=lambda i: (rows[i].a["MEASURE"] is not None, rows[i].a["MEASURE"] or 0))
    n = 1
    for i in order:
        r = rows[i]
        e = r.a.get("ENDPOINT")
        if e == 1:
            r.a["SORT_ID"] = n
        elif e == 2:
            r.a["SORT_ID"] = linecount + 2
        else:
            n += 1
            r.a["SORT_ID"] = n
    return rows


def get_watershedmidline(watershedmidpoints, wid):
    order = sorted(range(len(watershedmidpoints)), key=lambda i: watershedmidpoints[i].a["SORT_ID"])
    pts = [watershedmidpoints[i].geom.asPoint() for i in order]
    return Feat(_line([(p.x(), p.y()) for p in pts]), {"WID": wid})


def _retain(feats, fields):
    for f in feats:
        f.a = {k: f.a.get(k) for k in fields}


def measure_midlinetoriver(ctx, tlooklines, riverline, watershedmidpoints):
    for f in tlooklines:
        f.a["DIST_AB"] = 0
    # Intersect (POINT) crosslines with the river, then split the crosslines.
    nts = []
    for f in tlooklines:
        inter = f.geom.intersection(riverline.geom)
        if inter is None or inter.isEmpty():
            continue
        for v in inter.vertices():
            nts.append(_pt(v.x(), v.y()))
    vup = []
    for f in tlooklines:
        vup.extend(_split_line_at_points(ctx, f, nts))
    # Pieces from the midline to the river.
    sel = [f for f in vup if _near(ctx, f.geom, riverline.geom)
           and any(_near(ctx, f.geom, p.geom) for p in watershedmidpoints)]
    for f in sel:
        f.a["DIST_AB"] = 1
    # remove_duplicate_measure
    n = 0
    for f in sorted(sel, key=lambda f: (f.a["MEASURE"], f.geom.length())):
        m = f.a["MEASURE"]
        if n == m:
            f.a["DIST_AB"] = 0
        n = m
    ctx.msg("measured distance from midline to river ...")
    _retain(vup, ["MEASURE", "WID", "SIDE", "PERTO", "DIST_AB"])
    for f in vup:
        side, distab = f.a["SIDE"], f.a["DIST_AB"]
        f.a["DISTCLASS"] = {("RIGHT", 1): "RDIST_AB", ("LEFT", 1): "LDIST_AB",
                            ("RIGHT", 0): "RDIST_BC", ("LEFT", 0): "LDIST_BC"}.get((side, distab), "ERROR")
    return vup


def get_pivot_table(feats, pivotfield, valuefield, inputfields=("MEASURE", "WID", "PERTO")):
    """PivotTable: one row per input-field combination; a repeated pivot value
    within the same combination starts a new row.  Missing values are 0."""
    def val(f):
        return f.geom.length() if valuefield == "Shape_Length" else f.a.get(valuefield)

    recs = sorted(range(len(feats)), key=lambda i: (tuple(_sortable(feats[i].a.get(k)) for k in inputfields),
                                                    _sortable(feats[i].a.get(pivotfield)), i))
    columns = sorted({str(feats[i].a.get(pivotfield)) for i in recs})
    rows = []
    cur_key, cur = None, None
    for i in recs:
        f = feats[i]
        key = tuple(f.a.get(k) for k in inputfields)
        col = str(f.a.get(pivotfield))
        if key != cur_key or col in cur["_set"]:
            cur = {k: v for k, v in zip(inputfields, key)}
            cur["_set"] = set()
            rows.append(cur)
            cur_key = key
        cur[col] = val(f)
        cur["_set"].add(col)
    out = []
    for r in rows:
        r.pop("_set")
        for c in columns:
            r.setdefault(c, 0)
        out.append(r)
    return out, columns


def _sortable(v):
    return (v is not None, v if v is not None else 0)


def get_elevationpoints(ctx, in_raster, crosslines, crosslinedissolve, valleyfloor, csectionpoints, sidefield="SIDE_RIVER"):
    keep = ["MEASURE", "WID", sidefield, "PERTO"]
    # diss_endpoints -> elevendpoints (WSR / WSL)
    diss_endpoints = []
    for f in crosslinedissolve:
        s, e = _endpoints(f.geom)
        for p in (s, e):
            diss_endpoints.append(Feat(_pt(*p), {"MEASURE": f.a["MEASURE"]}))
    elevendpoints = _spatial_join_first(ctx, diss_endpoints, crosslines, keep)
    for f in elevendpoints:
        f.a["LOCATION"] = {"RIGHT": "WSR", "LEFT": "WSL"}.get(f.a.get(sidefield))
    ctx.msg("elevendpoints created ...")

    clipcross = _explode(_clip(crosslinedissolve, valleyfloor.geom))
    tsm = _unsplit(clipcross, "MEASURE")
    clipcrosslyr = [f for f in tsm if any(_near(ctx, f.geom, p.geom) for p in csectionpoints)]

    vfendpoints = []
    for f in clipcrosslyr:
        s, e = _endpoints(f.geom)
        for p in (s, e):
            vfendpoints.append(Feat(_pt(*p), {"MEASURE": f.a["MEASURE"]}))
    elevflrpoints = _spatial_join_first(ctx, vfendpoints, crosslines, keep)
    for f in elevflrpoints:
        f.a["LOCATION"] = {"RIGHT": "VFR", "LEFT": "VFL"}.get(f.a.get(sidefield))
    ctx.msg("elevflrpoints created ...")

    vfmidpoints = [Feat(f.geom.interpolate(f.geom.length() / 2.0), {"MEASURE": f.a["MEASURE"]}) for f in clipcrosslyr]
    elevmidpoints = _spatial_join_first(ctx, vfmidpoints, crosslines, keep)
    for f in elevmidpoints:
        f.a["LOCATION"] = "VFM"
    ctx.msg("elevmidpoints created ...")

    elevationpoints = elevendpoints + elevflrpoints + elevmidpoints
    for f in elevationpoints:
        p = f.geom.asPoint()
        f.a["ELEVATION"] = G.sample(in_raster, p.x(), p.y())
    if elevationpoints and elevationpoints[0].a["ELEVATION"] <= -9999:
        ctx.warn("some elevation points have nodata values")
    ctx.msg("extracted elevation values to points ...")
    _retain(elevationpoints, ["MEASURE", "WID", sidefield, "PERTO", "LOCATION", "ELEVATION"])
    return elevationpoints


def measure_valleyfloorwidth(ctx, crosslinedissolve, csectionpoints, valleyfloor, riverprofile):
    for f in riverprofile:
        f.a["WID_VF"] = 0
    nts = _explode(_clip(crosslinedissolve, valleyfloor.geom))
    tsm = _unsplit(nts, "MEASURE")
    sel = [f for f in tsm if any(_near(ctx, f.geom, p.geom) for p in csectionpoints)]
    endpoints = []
    for f in sel:
        s, e = _endpoints(f.geom)
        endpoints.extend([_pt(*s), _pt(*e)])
    vup = []
    for f in riverprofile:
        vup.extend(_split_line_at_points(ctx, f, endpoints))
    for f in vup:
        if any(_near(ctx, f.geom, p.geom) for p in csectionpoints):
            f.a["WID_VF"] = 1
    _retain(vup, ["MEASURE", "WID", "SIDE_RIVER", "PERTO", "WID_VF"])
    for f in vup:
        side, distab = f.a["SIDE_RIVER"], f.a["WID_VF"]
        f.a["VFSEGMENT"] = {("RIGHT", 1): "RVSEG_Y", ("LEFT", 1): "LVSEG_Y",
                            ("RIGHT", 0): "RVSEG_N", ("LEFT", 0): "LVSEG_N"}.get((side, distab), "ERROR")
    return vup


# ---------------------------------------------------------------------------
# Per-watershed driver
# ---------------------------------------------------------------------------

def get_riverprofile(ctx, in_raster, localdem, localwatershed, wid, interval, interval_profile,
                     count_intervalpr, adjust_edge, count_interval, options=None):
    options = options or {}
    reduced_raster, watershedpolygon = get_reducedinputs(localdem, localwatershed, adjust_edge, reduceby=2)
    hyd = Hydro(reduced_raster, fill=options.get("fill_sinks", False))

    lowpoints = get_hiflowaccpoint(ctx, hyd)
    ctx.label("lowpoint created")
    farpoints = get_farthestpoint(ctx, reduced_raster)
    ctx.label("farpoint created")

    riverline = extract_mainstream(ctx, hyd, wid, "LOW_HIGH", options.get("stream_simplify"))
    ctx.label("mainstream extracted")
    longest_riverline = get_longest_riverline(ctx, riverline, lowpoints, farpoints, interval_profile, wid)
    valleyfloor = extract_valleyfloor(ctx, reduced_raster, longest_riverline, wid, pthreshold=12.5)
    ctx.label("valleyfloor extracted")

    longaxis = get_longaxis(lowpoints, farpoints, wid)
    ctx.label("longaxis created")
    axismidpoints = create_axismidpoints(reduced_raster, lowpoints, longaxis, wid, interval, count_interval)
    ctx.label("axismidpoints created")
    axismidlines = get_axismidlines2(ctx, reduced_raster, axismidpoints, watershedpolygon, wid, lowpoints, "LONGAXIS")
    ctx.label("axismidlines created")
    axisdissolve = dissolve_axismidlines2(axismidlines)
    ctx.label("axisdissolve created")
    watershedmidpoints = get_watershedmidpoints(axisdissolve, lowpoints, farpoints)
    ctx.label("watershedmidpoints created")
    watershedmidline = get_watershedmidline(watershedmidpoints, wid)
    ctx.label("watershedmidline created")

    tecdem_mode = options.get("midline_method") == "TECDEM"
    if tecdem_mode and options.get("user_midline"):
        # TecDEM: the basin midline is the one drawn by the user.
        from . import tecdem
        watershedmidline = Feat(_line(tecdem.orient_from_outlet(options["user_midline"], lowpoints[0])),
                                {"WID": wid})
        ctx.msg("using the supplied midline for watershed ID %i ..." % wid)

    # midline profile
    midcrosspoints = create_csectionpoints(ctx, reduced_raster, lowpoints, watershedmidline, wid,
                                           interval_profile, count_intervalpr)
    midcrosslines = get_axismidlines(ctx, reduced_raster, midcrosspoints, watershedpolygon, wid,
                                     watershedmidline, lowpoints, "MIDLINE", "SIDE")
    midcrosslines = measure_midlinetoriver(ctx, midcrosslines, longest_riverline, midcrosspoints)
    midcrosslinestable = get_pivot_table(midcrosslines, "DISTCLASS", "Shape_Length")
    ctx.msg("midlinetable created ...")

    # river profile
    csectionpoints = create_csectionpoints(ctx, reduced_raster, lowpoints, longest_riverline, wid,
                                           interval_profile, count_intervalpr)
    ctx.label("csectionpoints created")
    rivercrosslines = get_axismidlines(ctx, reduced_raster, csectionpoints, watershedpolygon, wid,
                                       longest_riverline, lowpoints, "RIVER", "SIDE_RIVER")
    ctx.label("rivercrosslines created")
    crosslinedissolve = dissolve_axismidlines(rivercrosslines)
    ctx.label("crosslinedissolve created")

    elevationpoints = get_elevationpoints(ctx, in_raster, rivercrosslines, crosslinedissolve, valleyfloor,
                                          csectionpoints, "SIDE_RIVER")
    elevtable = get_pivot_table(elevationpoints, "LOCATION", "ELEVATION")
    ctx.msg("elevationpoints_sumtable created ...")

    rivercrosslines = measure_valleyfloorwidth(ctx, crosslinedissolve, csectionpoints, valleyfloor, rivercrosslines)
    ctx.msg("valley floor width measured ...")
    vfwtable = get_pivot_table(rivercrosslines, "VFSEGMENT", "Shape_Length")
    ctx.msg("vfw_sumtable created ...")

    out = {}
    if tecdem_mode:
        from . import tecdem
        vectors, perps = tecdem.tindex(ctx, _single_path(watershedmidline.geom), longest_riverline.geom,
                                       watershedpolygon, wid, options.get("tecdem_segment", 2000.0),
                                       options.get("tecdem_dd_rule", tecdem.RIVER_SIDE))
        out["rt_tindex"] = vectors
        out["rt_tindexperp"] = perps
        ctx.msg("TecDEM T-index calculated for %i segments ..." % len(vectors))

    out.update({
        "rt_rivercrosslines": rivercrosslines,
        "rt_midcrosslines": midcrosslines,
        "rt_watershedmidline": [watershedmidline],
        "rt_longest_riverline": [longest_riverline],
        "rt_valleyfloor": [valleyfloor],
        "rt_elevationpoints": elevationpoints,
        "rt_midcrosslinestable": midcrosslinestable,
        "rt_elevationpointstable": elevtable,
        "rt_valleyfloorwidthtable": vfwtable,
        # intermediate layers useful for checking
        "_lowpoints": lowpoints,
        "_farpoints": farpoints,
        "_riverline": riverline,
        "_reduced": reduced_raster,
    })
    return out


def parse_wid(value):
    """Watershed ID must be an integer (or text/float convertible to one)."""
    if isinstance(value, bool):
        raise ValueError("watershed ID must be an integer, not unknown type")
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if value == int(value):
            return int(value)
        raise ValueError("watershed ID must be an integer, not a float")
    if isinstance(value, str):
        try:
            v = float(value.strip())
        except ValueError:
            raise ValueError("watershed ID must be an integer, not a string")
        if v == int(v):
            return int(v)
        raise ValueError("watershed ID must be an integer, not a float")
    raise ValueError("watershed ID must be an integer, not unknown type")
