"""
TecDEM-style transverse topographic symmetry (T-index), after
Shahzad & Gloaguen (2011), TecDEM: A MATLAB based toolbox for tectonic
geomorphology, Part 2, Computers & Geosciences 37, 261-271.

Follows TecDEM's Methods/Tindex.m and Methods/Tindex_intersections.m:

* The basin midline is divided into segments of a fixed length (TecDEM's
  "segment length", default 2 km).  Only whole segments are used.
* For each segment, a line perpendicular to the segment's chord (the straight
  line between the segment's two ends) is drawn through the midline point half
  way along the segment.  The line is unbounded, so it reaches both sides.
* Da is the distance from that midline point to where the perpendicular
  crosses the main stream; Dd is the distance from the midline point to where
  it crosses the basin boundary.  T = Da / Dd.  When the perpendicular crosses
  the stream more than once, the crossing with the smallest T is used.
* ORIENT is the compass bearing (0-360, clockwise from north) from the
  midline point to that stream crossing, i.e. the direction the stream is
  shifted from the midline.

Differences from TecDEM, both inherent to running without its GUI:
* TecDEM measures segments by counting midline raster cells; here they are
  measured in map units along the midline.
* In TecDEM the user digitises the midline and picks the main stream by
  clicking.  Here the midline comes from a line layer (or the automatic
  ValleyMorph midline) and the main stream is the longest river found by the
  tool.

Dd rule: TecDEM's code takes the first boundary crossing returned by its
intersection routine, which is the westernmost crossing and is not always on
the stream's side of the midline.  The default "RIVER_SIDE" rule instead takes
the boundary crossing on the stream's side, beyond the stream (the divide on
that side, as in Cox 1994).  "TECDEM_FIRST" reproduces the TecDEM code.
"""

import math

from qgis.core import QgsGeometry, QgsPointXY

from .pipeline import Feat, _line, _single_path

RIVER_SIDE = "RIVER_SIDE"
TECDEM_FIRST = "TECDEM_FIRST"


def _points_of(geom):
    if geom is None or geom.isEmpty():
        return []
    return [(v.x(), v.y()) for v in geom.vertices()]


def orient_from_outlet(midline_pts, lowpoint):
    """Start the midline at the end nearest the watershed outlet (lowpoint)."""
    lx, ly = lowpoint
    if math.hypot(midline_pts[-1][0] - lx, midline_pts[-1][1] - ly) < \
            math.hypot(midline_pts[0][0] - lx, midline_pts[0][1] - ly):
        return midline_pts[::-1]
    return midline_pts


def tindex(ctx, midline_pts, river_geom, watershed_geom, wid, seglen, dd_rule=RIVER_SIDE):
    """
    Returns (vectors, perpendiculars):
      vectors: one line feature per segment with a stream crossing, from the
               midline point to the stream (TecDEM's asymmetry arrows), with
               WID, SEGMENT, MEASURE, DA, DD, T_INDEX, ORIENT.
      perpendiculars: the perpendicular for every segment, clipped to the
               watershed, for checking.
    """
    line = _line(midline_pts)
    length = line.length()
    nseg = int(math.floor(length / seglen + 1e-9))
    boundary = QgsGeometry(watershed_geom.constGet().boundary())
    bb = watershed_geom.boundingBox()
    reach = 2.0 * math.hypot(bb.width(), bb.height())

    vectors, perps = [], []
    for i in range(nseg):
        a = line.interpolate(i * seglen).asPoint()
        b = line.interpolate((i + 1) * seglen).asPoint()
        m = line.interpolate((i + 0.5) * seglen).asPoint()
        dx, dy = b.x() - a.x(), b.y() - a.y()
        n = math.hypot(dx, dy)
        if n == 0:
            continue
        # unit vector perpendicular to the chord
        ux, uy = -dy / n, dx / n
        p1 = QgsPointXY(m.x() - ux * reach, m.y() - uy * reach)
        p2 = QgsPointXY(m.x() + ux * reach, m.y() + uy * reach)
        perp = QgsGeometry.fromPolylineXY([p1, p2])
        clipped = perp.intersection(watershed_geom)
        measure = (i + 0.5) * seglen
        if clipped is not None and not clipped.isEmpty():
            perps.append(Feat(clipped, {"WID": wid, "SEGMENT": i + 1, "MEASURE": measure}))

        def t_of(p):
            return (p[0] - m.x()) * ux + (p[1] - m.y()) * uy

        rpts = _points_of(perp.intersection(river_geom))
        bpts = _points_of(perp.intersection(boundary))
        if not rpts or not bpts:
            continue
        best = None
        for rp in rpts:
            tr = t_of(rp)
            da = abs(tr)
            if dd_rule == TECDEM_FIRST:
                bp = min(bpts, key=lambda p: (p[0], p[1]))
                dd = math.hypot(bp[0] - m.x(), bp[1] - m.y())
            else:
                side = [t_of(p) for p in bpts if (t_of(p) >= 0) == (tr >= 0)]
                beyond = [abs(t) for t in side if abs(t) >= da - ctx.tol]
                if beyond:
                    dd = min(beyond)
                elif side:
                    dd = max(abs(t) for t in side)
                else:
                    continue
            if dd <= 0:
                continue
            t = da / dd
            if best is None or t < best[0]:
                best = (t, rp, da, dd)
        if best is None:
            continue
        t, rp, da, dd = best
        orient = math.degrees(math.atan2(rp[0] - m.x(), rp[1] - m.y())) % 360.0
        end = rp if da > 0 else (m.x(), m.y() + ctx.tol)  # stream on the midline: T = 0
        vectors.append(Feat(_line([(m.x(), m.y()), end]),
                            {"WID": wid, "SEGMENT": i + 1, "MEASURE": measure,
                             "DA": da, "DD": dd, "T_INDEX": t, "ORIENT": orient}))
    if nseg == 0:
        ctx.warn("watershed ID %i: midline (%.0f m) is shorter than one TecDEM segment (%.0f m)"
                 % (wid, length, seglen))
    return vectors, perps


def pick_midline(ctx, midline_feats, watershed_geom, wid, id_field=None):
    """
    The drawn midline for a watershed.  With id_field, the line whose ID equals
    the watershed ID.  Otherwise the line lying most completely inside the
    watershed (ties: the longer one).  Watersheds can be nested, so a sub-basin's
    midline also lies inside its parent; a warning is given when the choice is
    ambiguous, and an ID field should be used.
    """
    from .pipeline import parse_wid
    if id_field:
        hits = []
        for f in midline_feats:
            try:
                if parse_wid(f[id_field]) == wid:
                    hits.append(f.geometry())
            except ValueError:
                continue
        if len(hits) > 1:
            ctx.warn("watershed ID %i has %i drawn midlines; using the first" % (wid, len(hits)))
        g = hits[0] if hits else None
    else:
        cands = []
        for f in midline_feats:
            g = f.geometry()
            if g.isEmpty() or not g.boundingBox().intersects(watershed_geom.boundingBox()):
                continue
            inter = g.intersection(watershed_geom)
            inside = inter.length() if inter is not None and not inter.isEmpty() else 0.0
            frac = inside / g.length() if g.length() else 0.0
            if frac >= 0.8:
                cands.append((frac, g.length(), g))
        cands.sort(key=lambda c: (-round(c[0], 3), -c[1]))
        if len(cands) > 1:
            ctx.warn("watershed ID %i: %i drawn midlines lie inside it (nested watersheds?); using the "
                     "longest. Set the midline ID field to match them exactly." % (wid, len(cands)))
        g = cands[0][2] if cands else None
    if g is None:
        return None
    pts = _single_path(g)
    return pts if len(pts) >= 2 else None
