"""
Raster helpers that reproduce the ArcGIS Spatial Analyst behaviour used by the
original ValleyMorph Tool (river_transect.py).

A Grid is a north-up raster held in memory as a float64 numpy array where
NoData cells are NaN.  Cell (r, c) has its centre at
    x = x0 + (c + 0.5) * cs,  y = y0 - (r + 0.5) * cs
"""

import math

import numpy as np
from osgeo import gdal
from qgis.core import QgsGeometry, QgsPointXY, QgsRectangle

gdal.UseExceptions()


class Grid:
    def __init__(self, data, x0, y0, cs, is_int=False):
        self.data = data
        self.x0 = float(x0)
        self.y0 = float(y0)
        self.cs = float(cs)
        self.is_int = is_int

    # --- geometry -----------------------------------------------------------
    @property
    def nrows(self):
        return self.data.shape[0]

    @property
    def ncols(self):
        return self.data.shape[1]

    @property
    def left(self):
        return self.x0

    @property
    def top(self):
        return self.y0

    @property
    def right(self):
        return self.x0 + self.ncols * self.cs

    @property
    def bottom(self):
        return self.y0 - self.nrows * self.cs

    @property
    def valid(self):
        return ~np.isnan(self.data)

    def center(self, r, c):
        return (self.x0 + (c + 0.5) * self.cs, self.y0 - (r + 0.5) * self.cs)

    def cell_of(self, x, y):
        c = int(math.floor((x - self.x0) / self.cs))
        r = int(math.floor((self.y0 - y) / self.cs))
        return r, c

    def like(self, data, is_int=None):
        return Grid(data, self.x0, self.y0, self.cs, self.is_int if is_int is None else is_int)

    def min(self):
        return float(np.nanmin(self.data))

    def max(self):
        return float(np.nanmax(self.data))


def _ogr_mem():
    from osgeo import ogr
    ogr.UseExceptions()
    # GDAL >= 3.11 renamed the in-memory vector driver to "MEM".
    if int(gdal.VersionInfo()) >= 3110000:
        return ogr.GetDriverByName("MEM")
    return ogr.GetDriverByName("Memory")


def read_raster(path):
    """Reads band 1 of a GDAL raster into a Grid (NoData -> NaN)."""
    ds = gdal.Open(path)
    gt = ds.GetGeoTransform()
    if gt[2] != 0 or gt[4] != 0:
        raise ValueError("rotated rasters are not supported")
    if abs(abs(gt[1]) - abs(gt[5])) > 1e-6 * abs(gt[1]):
        raise ValueError("the DEM must have square cells")
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float64)
    nd = band.GetNoDataValue()
    if nd is not None:
        arr[arr == nd] = np.nan
    is_int = band.DataType in (gdal.GDT_Byte, gdal.GDT_Int16, gdal.GDT_UInt16,
                               gdal.GDT_Int32, gdal.GDT_UInt32)
    return Grid(arr, gt[0], gt[3], abs(gt[1]), is_int)


def polygon_cell_mask(grid_x0, grid_y0, cs, nrows, ncols, polygon):
    """True where a cell centre falls inside the polygon (ArcGIS feature-to-raster rule)."""
    mask = np.zeros((nrows, ncols), dtype=bool)
    if nrows == 0 or ncols == 0:
        return mask
    # Rasterise with GDAL for speed: burns cells whose centre is inside.
    drv = gdal.GetDriverByName("MEM")
    ds = drv.Create("", ncols, nrows, 1, gdal.GDT_Byte)
    ds.SetGeoTransform((grid_x0, cs, 0, grid_y0, 0, -cs))
    from osgeo import ogr
    ogr.UseExceptions()
    mds = _ogr_mem().CreateDataSource("")
    lyr = mds.CreateLayer("p", geom_type=ogr.wkbMultiPolygon)
    f = ogr.Feature(lyr.GetLayerDefn())
    f.SetGeometry(ogr.CreateGeometryFromWkb(bytes(polygon.asWkb())))
    lyr.CreateFeature(f)
    gdal.RasterizeLayer(ds, [1], lyr, burn_values=[1])
    mask[:] = ds.GetRasterBand(1).ReadAsArray() == 1
    return mask


def extract_by_mask(src, polygon):
    """
    ArcGIS 10.1 ExtractByMask with a feature mask and no snap raster: the output
    grid is anchored at the mask's lower-left corner (xmin, ymin) with the input
    cell size, values are taken from the input cell containing each output cell
    centre (nearest neighbour), and cells whose centre is outside the polygon are
    NoData.
    """
    bb = polygon.boundingBox()
    cs = src.cs
    ncols = max(1, int(math.ceil((bb.xMaximum() - bb.xMinimum()) / cs - 1e-9)))
    nrows = max(1, int(math.ceil((bb.yMaximum() - bb.yMinimum()) / cs - 1e-9)))
    x0 = bb.xMinimum()
    y0 = bb.yMinimum() + nrows * cs
    # Nearest-neighbour resample from src.
    xs = x0 + (np.arange(ncols) + 0.5) * cs
    ys = y0 - (np.arange(nrows) + 0.5) * cs
    sc = np.floor((xs - src.x0) / src.cs).astype(int)
    sr = np.floor((src.y0 - ys) / src.cs).astype(int)
    out = np.full((nrows, ncols), np.nan)
    okc = (sc >= 0) & (sc < src.ncols)
    okr = (sr >= 0) & (sr < src.nrows)
    rr, cc = np.meshgrid(np.where(okr)[0], np.where(okc)[0], indexing="ij")
    out[rr, cc] = src.data[sr[rr], sc[cc]]
    inside = polygon_cell_mask(x0, y0, cs, nrows, ncols, polygon)
    out[~inside] = np.nan
    return Grid(out, x0, y0, cs, src.is_int)


def sample(grid, x, y, nodata=-9999):
    """ExtractValuesToPoints (no interpolation): value of the cell containing the point."""
    r, c = grid.cell_of(x, y)
    if 0 <= r < grid.nrows and 0 <= c < grid.ncols:
        v = grid.data[r, c]
        if not np.isnan(v):
            return int(v) if grid.is_int else float(v)
    return nodata


def slope_percent(grid):
    """
    ArcGIS Slope (PERCENT_RISE, z-factor 1) using Horn's 3x3 method.  NoData
    neighbours (including those outside the raster) take the centre cell's value.
    """
    z = grid.data
    p = np.pad(z, 1, constant_values=np.nan)
    ctr = p[1:-1, 1:-1]

    def nb(dr, dc):
        a = p[1 + dr:p.shape[0] - 1 + dr, 1 + dc:p.shape[1] - 1 + dc]
        return np.where(np.isnan(a), ctr, a)

    a, b, c = nb(-1, -1), nb(-1, 0), nb(-1, 1)
    d, f = nb(0, -1), nb(0, 1)
    g, h, i = nb(1, -1), nb(1, 0), nb(1, 1)
    cs = grid.cs
    dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * cs)
    dzdy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * cs)
    s = np.sqrt(dzdx ** 2 + dzdy ** 2) * 100.0
    s[np.isnan(z)] = np.nan
    return grid.like(s, is_int=False)


def euclidean_distance(grid, source):
    """EucDistance from source cells (bool array) measured between cell centres."""
    from scipy import ndimage
    if not source.any():
        raise ValueError("no source cells for euclidean distance")
    d = ndimage.distance_transform_edt(~source) * grid.cs
    d = d.astype(np.float64)
    d[~grid.valid] = np.nan
    return grid.like(d, is_int=False)


def raster_to_points(grid, cellmask):
    """RasterToPoint: one point per selected cell, row-major order (pointid 1..n)."""
    rows, cols = np.nonzero(cellmask)
    pts = []
    for k, (r, c) in enumerate(zip(rows, cols)):
        x, y = grid.center(r, c)
        pts.append((x, y, k + 1, grid.data[r, c]))
    return pts


def polyline_to_raster_mask(grid, geoms):
    """
    PolylineToRaster (MAXIMUM_LENGTH) presence: True for every cell crossed by a
    line over a non-zero length.  Uses the grid's extent, cell size and alignment.
    """
    mask = np.zeros(grid.data.shape, dtype=bool)
    cs = grid.cs
    eps = 1e-9 * cs
    for g in geoms:
        for part in g.asMultiPolyline() if g.isMultipart() else [g.asPolyline()]:
            for p, q in zip(part[:-1], part[1:]):
                x1, y1, x2, y2 = p.x(), p.y(), q.x(), q.y()
                c0 = int(math.floor((min(x1, x2) - grid.x0) / cs)) - 1
                c1 = int(math.floor((max(x1, x2) - grid.x0) / cs)) + 1
                r0 = int(math.floor((grid.y0 - max(y1, y2)) / cs)) - 1
                r1 = int(math.floor((grid.y0 - min(y1, y2)) / cs)) + 1
                for r in range(max(r0, 0), min(r1, grid.nrows - 1) + 1):
                    ytop = grid.y0 - r * cs
                    ybot = ytop - cs
                    for c in range(max(c0, 0), min(c1, grid.ncols - 1) + 1):
                        xl = grid.x0 + c * cs
                        if _clip_len(x1, y1, x2, y2, xl, ybot, xl + cs, ytop) > eps:
                            mask[r, c] = True
    return mask


def _clip_len(x1, y1, x2, y2, xmin, ymin, xmax, ymax):
    """Length of segment inside an axis-aligned box (Liang-Barsky)."""
    dx, dy = x2 - x1, y2 - y1
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, x1 - xmin), (dx, xmax - x1), (-dy, y1 - ymin), (dy, ymax - y1)):
        if p == 0:
            if q < 0:
                return 0.0
        else:
            t = q / p
            if p < 0:
                if t > t1:
                    return 0.0
                t0 = max(t0, t)
            else:
                if t < t0:
                    return 0.0
                t1 = min(t1, t)
    return max(0.0, t1 - t0) * math.hypot(dx, dy)


def raster_to_polygons(grid, cellmask, simplify=True):
    """
    RasterToPolygon of a single-valued region raster.  Cells sharing an edge form
    one polygon.  With simplify=True the stair-stepped outline is replaced by a
    line through the midpoints of the cell edges (as ArcGIS's SIMPLIFY option
    does), with collinear vertices removed.
    Returns a list of QgsGeometry polygons.
    """
    from osgeo import ogr
    ogr.UseExceptions()
    nrows, ncols = cellmask.shape
    drv = gdal.GetDriverByName("MEM")
    ds = drv.Create("", ncols, nrows, 1, gdal.GDT_Byte)
    ds.SetGeoTransform((grid.x0, grid.cs, 0, grid.y0, 0, -grid.cs))
    band = ds.GetRasterBand(1)
    band.WriteArray(cellmask.astype(np.uint8))
    band.SetNoDataValue(0)
    mds = _ogr_mem().CreateDataSource("")
    lyr = mds.CreateLayer("poly", geom_type=ogr.wkbPolygon)
    lyr.CreateField(ogr.FieldDefn("v", ogr.OFTInteger))
    gdal.Polygonize(band, band, lyr, 0, [], callback=None)
    out = []
    for f in lyr:
        g = QgsGeometry.fromWkt(f.GetGeometryRef().ExportToWkt())
        if simplify:
            poly = g.asPolygon()
            rings = [_midpoint_ring(r, grid.cs) for r in poly]
            rings = [r for r in rings if len(r) >= 4]
            if not rings:
                continue
            g = QgsGeometry.fromPolygonXY(rings)
            if not g.isGeosValid():
                g = g.makeValid()
        out.append(g)
    return out


def _midpoint_ring(ring, cs):
    pts = [(p.x(), p.y()) for p in ring]
    if pts[0] == pts[-1]:
        pts = pts[:-1]
    # densify to unit cell edges
    unit = []
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        steps = max(1, int(round(math.hypot(x2 - x1, y2 - y1) / cs)))
        for k in range(steps):
            unit.append((x1 + (x2 - x1) * k / steps, y1 + (y2 - y1) * k / steps))
    m = len(unit)
    mids = []
    for i in range(m):
        x1, y1 = unit[i]
        x2, y2 = unit[(i + 1) % m]
        mids.append(((x1 + x2) / 2.0, (y1 + y2) / 2.0))
    # remove collinear vertices
    res = []
    k = len(mids)
    for i in range(k):
        ax, ay = mids[i - 1]
        bx, by = mids[i]
        cx, cy = mids[(i + 1) % k]
        cross = (bx - ax) * (cy - by) - (by - ay) * (cx - bx)
        if abs(cross) > 1e-9 * cs * cs:
            res.append(QgsPointXY(bx, by))
    if res:
        res.append(res[0])
    return res
