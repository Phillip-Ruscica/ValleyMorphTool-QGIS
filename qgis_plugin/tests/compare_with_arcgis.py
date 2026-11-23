"""
Runs the QGIS ValleyMorph Tool on the sample data in test/ and compares the
results with the original ArcGIS output in test/output.gdb.

Run with QGIS's Python, for example on macOS:
    /Applications/QGIS.app/Contents/MacOS/python3 qgis_plugin/tests/compare_with_arcgis.py
or from the QGIS Python console:
    exec(open('/path/to/repo/qgis_plugin/tests/compare_with_arcgis.py').read())
"""

import os
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "qgis_plugin", "valleymorph_tool"))

from osgeo import ogr  # noqa: E402
from qgis.core import QgsApplication, QgsProcessingFeedback, QgsVectorLayer  # noqa: E402

ogr.UseExceptions()
_app = None
if QgsApplication.instance() is None:
    _app = QgsApplication([], False)
    _app.initQgis()

from vmt import runner  # noqa: E402


class Feedback(QgsProcessingFeedback):
    def pushWarning(self, text):
        print("  warning:", text)


def run(out):
    lyr = QgsVectorLayer(os.path.join(REPO, "test", "wshedbound.shp"), "watersheds", "ogr")
    # Same settings as test/output.gdb: FID as watershed ID, 50 m, 500 m, no edge adjustment.
    runner.run(os.path.join(REPO, "test", "wsheddem"), list(lyr.getFeatures()), lyr.crs(), None,
               50, 500, None, False, out, Feedback())


def lengths(ds, name):
    out = {}
    for f in ds.GetLayerByName(name):
        g = f.GetGeometryRef().GetLinearGeometry()
        d = out.setdefault(f.GetField("WID"), [0, 0.0, 0.0])
        d[0] += 1
        d[1] += g.Length() if g.GetDimension() == 1 else g.Boundary().Length()
        d[2] += g.GetArea()
    return out


def table_agreement(ref, new, name):
    rows = {}
    for f in new.GetLayerByName(name):
        r = f.items()
        rows.setdefault((r["WID"], r["MEASURE"]), r)
    diffs = []
    for f in ref.GetLayerByName(name):
        r = f.items()
        o = rows.get((r["WID"], r["MEASURE"]))
        if not o:
            continue
        for k, a in r.items():
            if k in ("MEASURE", "WID", "PERTO"):
                continue
            b = o.get(k)
            if a is None or b is None or a <= -9999 or b <= -9999:
                continue
            diffs.append(abs(a - b))
    diffs.sort()
    within = sum(d < 1 for d in diffs)
    return len(diffs), within, diffs[len(diffs) // 2] if diffs else 0


def main():
    out = os.path.join(tempfile.mkdtemp(), "valleymorph_sample.gpkg")
    run(out)
    ref = ogr.Open(os.path.join(REPO, "test", "output.gdb"))
    new = ogr.Open(out)
    print("\nOutput:", out)
    print("\n%-22s %4s %18s %18s %20s" % ("layer", "WID", "features ArcGIS/QGIS", "length ArcGIS/QGIS", "area ArcGIS/QGIS"))
    for name in ["rt_longest_riverline", "rt_watershedmidline", "rt_valleyfloor",
                 "rt_rivercrosslines", "rt_midcrosslines", "rt_elevationpoints"]:
        a, b = lengths(ref, name), lengths(new, name)
        for wid in sorted(set(a) | set(b)):
            ra, rb = a.get(wid, [0, 0, 0]), b.get(wid, [0, 0, 0])
            print("%-22s %4s %9i/%-8i %9.0f/%-8.0f %10.0f/%-10.0f" % (name, wid, ra[0], rb[0], ra[1], rb[1], ra[2], rb[2]))
    print("\n%-26s %7s %16s %12s" % ("summary table", "values", "within 1 m/unit", "median diff"))
    for name in ["rt_midcrosslinestable", "rt_elevationpointstable", "rt_valleyfloorwidthtable"]:
        n, w, med = table_agreement(ref, new, name)
        print("%-26s %7i %10i (%3.0f%%) %12.1f" % (name, n, w, 100.0 * w / max(n, 1), med))


if __name__ == "__main__" or _app is None:
    main()
