"""
main() of the original tool: loops over the watersheds, runs get_riverprofile
for each one, and writes the nine rt_* outputs.
"""

from qgis.core import (
    QgsFeature,
    QgsField,
    QgsFields,
    QgsVectorFileWriter,
    QgsWkbTypes,
)

try:  # QGIS >= 3.38
    from qgis.PyQt.QtCore import QMetaType
    _T = {"double": QMetaType.Type.Double, "int": QMetaType.Type.Int, "str": QMetaType.Type.QString}
except (ImportError, AttributeError):  # older QGIS
    from qgis.PyQt.QtCore import QVariant
    _T = {"double": QVariant.Double, "int": QVariant.Int, "str": QVariant.String}

from . import grid as G
from .pipeline import Context, VMTSkip, get_riverprofile, parse_wid

FEATURE_OUTPUTS = [
    ("rt_rivercrosslines", QgsWkbTypes.MultiLineString,
     [("MEASURE", "double"), ("WID", "int"), ("SIDE_RIVER", "str"), ("PERTO", "str"),
      ("WID_VF", "int"), ("VFSEGMENT", "str")]),
    ("rt_midcrosslines", QgsWkbTypes.MultiLineString,
     [("MEASURE", "double"), ("WID", "int"), ("SIDE", "str"), ("PERTO", "str"),
      ("DIST_AB", "int"), ("DISTCLASS", "str")]),
    ("rt_watershedmidline", QgsWkbTypes.MultiLineString, [("WID", "int")]),
    ("rt_longest_riverline", QgsWkbTypes.MultiLineString, [("WID", "int")]),
    ("rt_valleyfloor", QgsWkbTypes.MultiPolygon, [("WID", "int")]),
    ("rt_elevationpoints", QgsWkbTypes.Point,
     [("MEASURE", "double"), ("WID", "int"), ("SIDE_RIVER", "str"), ("PERTO", "str"),
      ("LOCATION", "str"), ("ELEVATION", "elev")]),
]
TABLE_OUTPUTS = [
    ("rt_midcrosslinestable", "double"),
    ("rt_elevationpointstable", "elev"),
    ("rt_valleyfloorwidthtable", "double"),
]


def list_watershedid(features, fieldname, interval_profile, ctx):
    widlist, feats = [], {}
    for f in features:
        try:
            # No field: use the feature ID, as the original did with a shapefile's FID.
            wid = parse_wid(f.id() if not fieldname else f[fieldname])
        except ValueError as err:
            raise RuntimeError(str(err))
        if wid in widlist or wid in feats:
            raise RuntimeError("some watershed IDs have duplicates")
        if interval_profile > f.geometry().constGet().perimeter():
            ctx.warn("skip watershed ID %i - too small for river profile interval" % wid)
            feats[wid] = None
            continue
        widlist.append(wid)
        feats[wid] = f.geometry()
    if not widlist:
        raise RuntimeError("no valid watersheds to process")
    return widlist, feats


def run(dem_path, features, crs, fieldname, interval, interval_profile, count_intervalpr, adjust_edge,
        output_path, feedback=None, transform_context=None, options=None):
    """
    dem_path: DEM file; features: watershed polygons (QgsFeature, in processing
    order); crs: their CRS (same as the DEM); fieldname: watershed ID field, or
    None to use feature IDs.
    """
    ctx = Context(crs, feedback, transform_context)
    in_raster = G.read_raster(dem_path)
    if adjust_edge:
        ctx.warn("reduced inputs by 2 cells to account for edge effects")
    widlist, geoms = list_watershedid(features, fieldname, interval_profile, ctx)

    results = {name: [] for name, _, _ in FEATURE_OUTPUTS}
    tables = {name: ([], []) for name, _ in TABLE_OUTPUTS}
    done = []
    for n, wid in enumerate(widlist):
        if feedback and feedback.isCanceled():
            break
        ctx.msg("processing watershed ID %i (%i of %i) ..." % (wid, n + 1, len(widlist)))
        local_watershed = geoms[wid]
        local_dem = G.extract_by_mask(in_raster, local_watershed)
        try:
            out = get_riverprofile(ctx, in_raster, local_dem, local_watershed, wid, interval,
                                   interval_profile, count_intervalpr, adjust_edge, None, options)
        except VMTSkip as err:
            ctx.warn("%s - skip watershed ID %i" % (err, wid))
            if feedback:
                feedback.setProgress(100.0 * (n + 1) / len(widlist))
            continue
        for name, _, _ in FEATURE_OUTPUTS:
            results[name].extend(out[name])
        for name, _ in TABLE_OUTPUTS:
            rows, cols = out[name]
            tables[name][0].extend(rows)
            for c in cols:
                if c not in tables[name][1]:
                    tables[name][1].append(c)
        done.append(wid)
        if feedback:
            feedback.setProgress(100.0 * (n + 1) / len(widlist))

    if not done:
        raise RuntimeError("no watershed could be processed")
    elev_int = in_raster.is_int
    write_outputs(output_path, crs, results, tables, elev_int, ctx.tc)
    ctx.msg("process successfully completed")
    return done


def _qfield(name, kind, elev_int):
    if kind == "elev":
        kind = "int" if elev_int else "double"
    return QgsField(name, _T[kind], len=254 if kind == "str" else 0)


def write_outputs(path, crs, results, tables, elev_int, tc):
    is_gdb = path.lower().endswith(".gdb")
    driver = "OpenFileGDB" if is_gdb else "GPKG"
    first = True
    for name, wkb, schema in FEATURE_OUTPUTS:
        fields = QgsFields()
        for fname, kind in schema:
            fields.append(_qfield(fname, kind, elev_int))
        is_line = QgsWkbTypes.geometryType(wkb) == QgsWkbTypes.LineGeometry
        is_poly = QgsWkbTypes.geometryType(wkb) == QgsWkbTypes.PolygonGeometry
        if not is_gdb and (is_line or is_poly):
            fields.append(QgsField("Shape_Length", _T["double"]))
            if is_poly:
                fields.append(QgsField("Shape_Area", _T["double"]))
        w = _writer(path, name, fields, wkb, crs, tc, driver, first, is_gdb)
        first = False
        for f in results[name]:
            qf = QgsFeature(fields)
            g = f.geom
            if QgsWkbTypes.isMultiType(wkb) and not g.isMultipart():
                g = type(g)(g)
                g.convertToMultiType()
            qf.setGeometry(g)
            for fname, _ in schema:
                qf[fname] = f.a.get(fname)
            if not is_gdb and (is_line or is_poly):
                qf["Shape_Length"] = g.length() if is_line else g.constGet().perimeter()
                if is_poly:
                    qf["Shape_Area"] = g.area()
            w.addFeature(qf)
        del w
    for name, kind in TABLE_OUTPUTS:
        rows, cols = tables[name]
        fields = QgsFields()
        fields.append(_qfield("MEASURE", "double", elev_int))
        fields.append(_qfield("WID", "int", elev_int))
        fields.append(_qfield("PERTO", "str", elev_int))
        for c in sorted(cols):
            fields.append(_qfield(c, kind, elev_int))
        w = _writer(path, name, fields, QgsWkbTypes.NoGeometry, crs, tc, driver, False, is_gdb)
        for r in rows:
            qf = QgsFeature(fields)
            for i in range(fields.count()):
                name = fields.at(i).name()
                # pivot columns missing for a watershed are 0, as in PivotTable
                qf.setAttribute(i, r.get(name, 0 if i >= 3 else None))
            w.addFeature(qf)
        del w


def _writer(path, layer_name, fields, wkb, crs, tc, driver, first, is_gdb):
    opts = QgsVectorFileWriter.SaveVectorOptions()
    opts.driverName = driver
    opts.layerName = layer_name
    opts.fileEncoding = "UTF-8"
    if first:
        opts.actionOnExistingFile = QgsVectorFileWriter.CreateOrOverwriteFile
    else:
        opts.actionOnExistingFile = QgsVectorFileWriter.CreateOrOverwriteLayer
    if is_gdb and wkb != QgsWkbTypes.NoGeometry:
        opts.layerOptions = ["CREATE_SHAPE_AREA_AND_LENGTH_FIELDS=YES"]
    w = QgsVectorFileWriter.create(path, fields, wkb, crs, tc, opts)
    if w.hasError() != QgsVectorFileWriter.NoError:
        raise RuntimeError("cannot write %s: %s" % (layer_name, w.errorMessage()))
    return w
