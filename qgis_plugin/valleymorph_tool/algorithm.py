"""
QGIS Processing algorithm: ValleyMorph Tool.

Same inputs as the ArcGIS toolbox version (river_transect.py / ValleyMorph Tool.tbx):
output workspace, watershed DEM, watershed boundary, watershed ID field,
longaxis interval, midline and river profile interval (MRPI), number of MRPI
and adjust for edge effects.  Writes the same nine rt_* outputs.
"""

import os

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingContext,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterDefinition,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterLayer,
    QgsProcessingOutputMultipleLayers,
)
from qgis.PyQt.QtCore import QCoreApplication
from qgis.PyQt.QtGui import QIcon

from .vmt import runner

OUTPUT_LAYERS = [name for name, _, _ in runner.FEATURE_OUTPUTS] + [name for name, _ in runner.TABLE_OUTPUTS]
TECDEM_LAYERS = [name for name, _, _ in runner.TECDEM_OUTPUTS]
MIDLINE_METHODS = ["VALLEYMORPH", "TECDEM"]
DD_RULES = ["RIVER_SIDE", "TECDEM_FIRST"]


def tr(text):
    return QCoreApplication.translate("ValleyMorphTool", text)


class ValleyMorphAlgorithm(QgsProcessingAlgorithm):
    OUTPUT = "OUTPUT"
    DEM = "DEM"
    WATERSHEDS = "WATERSHEDS"
    WATERSHED_ID = "WATERSHED_ID"
    INTERVAL = "LONGAXIS_INTERVAL"
    MRPI = "MRPI"
    NUMBER_MRPI = "NUMBER_MRPI"
    ADJUST_EDGE = "ADJUST_EDGE"
    FILL_SINKS = "FILL_SINKS"
    MIDLINE_METHOD = "MIDLINE_METHOD"
    MIDLINES = "MIDLINES"
    MIDLINES_ID = "MIDLINES_ID"
    SEGMENT = "TECDEM_SEGMENT"
    DD_RULE = "TECDEM_DD_RULE"
    LAYERS = "LAYERS"

    def name(self):
        return "valleymorph"

    def displayName(self):
        return tr("ValleyMorph Tool")

    def group(self):
        return ""

    def groupId(self):
        return ""

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(__file__), "icon.svg"))

    def shortHelpString(self):
        return tr(
            "<p>Extracts the transverse topographic symmetry (T-) factor and valley width to valley "
            "height (Vf) ratio inputs for each watershed: crosslines perpendicular to the watershed "
            "midline and to the main river, distances from the midline to the river and from the river "
            "to the watershed boundary, elevations where the crosslines meet the watershed boundary and "
            "the valley floor, and the valley floor width.</p>"
            "<p>QGIS port of the ArcGIS ValleyMorph Tool (Daxberger, Dalumpines, Scott &amp; Riller 2014, "
            "Computers &amp; Geosciences, doi:10.1016/j.cageo.2014.05.015). Outputs use the same names, "
            "fields and summary tables as the original: rt_rivercrosslines, rt_midcrosslines, "
            "rt_watershedmidline, rt_longest_riverline, rt_valleyfloor, rt_elevationpoints, "
            "rt_midcrosslinestable, rt_elevationpointstable and rt_valleyfloorwidthtable.</p>"
            "<p><b>Watershed DEM</b>: should extend at least two cells beyond the watersheds, otherwise "
            "elevations on the watershed edge are NoData (-9999). Use <b>Adjust for edge effects</b> to "
            "shrink each watershed by two cells instead.</p>"
            "<p><b>Watershed ID</b>: integer field (or text convertible to integer). Leave empty to use "
            "the feature ID, like the FID of a shapefile in the original.</p>"
            "<p><b>Longaxis interval</b>: spacing (m) of the crosslines used to build the watershed "
            "midline, measured from the point of highest flow accumulation. Lower it if part of the "
            "midline falls outside the watershed.</p>"
            "<p><b>Fill sinks</b>: off by default, as in the original. Turn it on if the DEM has not been "
            "hydrologically filled and watersheds are skipped because the river drains into a pit inside "
            "the watershed. It only affects flow routing (lowpoint and river); elevations and slopes still "
            "come from the original DEM.</p>"
            "<p><b>Midline method</b>: <i>ValleyMorph</i> is the original tool. <i>TecDEM</i> also "
            "computes the T-index the way TecDEM does (Shahzad &amp; Gloaguen 2011): the midline is cut "
            "into segments of the given length, a perpendicular is drawn through the middle of each "
            "segment, and T = Da/Dd, where Da is the distance from the midline to the main river and Dd "
            "the distance from the midline to the basin boundary. Results go to rt_tindex (one arrow per "
            "segment from the midline to the river, with DA, DD, T_INDEX and ORIENT, the bearing of the "
            "river from the midline) and rt_tindexperp (the perpendiculars). In TecDEM the midline is "
            "drawn by hand: supply a line layer of drawn midlines to use them, ideally with a field holding "
            "the watershed ID (without one, each watershed takes the line lying inside it, which is "
            "ambiguous for nested watersheds). The rest of the tool then uses that midline too. "
            "Without one, the automatic ValleyMorph midline is used.</p>"
            "<p><b>Midline and river profile interval (MRPI)</b>: spacing (m) of the midcrosslines and "
            "rivercrosslines. <b>Number of MRPI</b> limits how many are created (0 = no limit).</p>"
        )

    def createInstance(self):
        return ValleyMorphAlgorithm()

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.DEM, tr("Watershed DEM")))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.WATERSHEDS, tr("Watershed boundary"), [QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterField(
            self.WATERSHED_ID, tr("Watershed ID (empty = feature ID)"), None, self.WATERSHEDS,
            QgsProcessingParameterField.Any, optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            self.INTERVAL, tr("Longaxis interval in meters"), QgsProcessingParameterNumber.Double,
            50.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.MRPI, tr("Midline and river profile interval (MRPI) in meters"),
            QgsProcessingParameterNumber.Double, 500.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.NUMBER_MRPI, tr("Number of MRPI (0 = no limit)"), QgsProcessingParameterNumber.Integer,
            0, minValue=0))
        self.addParameter(QgsProcessingParameterBoolean(
            self.ADJUST_EDGE, tr("Adjust for edge effects"), False))
        self.addParameter(QgsProcessingParameterBoolean(
            self.FILL_SINKS, tr("Fill sinks in the DEM before tracing the river (not in the original tool)"),
            False))
        self.addParameter(QgsProcessingParameterEnum(
            self.MIDLINE_METHOD, tr("Midline method"),
            [tr("ValleyMorph: automatic long-axis midline (original)"),
             tr("TecDEM: midline segments and T-index (Shahzad & Gloaguen 2011)")], defaultValue=0))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.MIDLINES, tr("TecDEM: drawn basin midlines (optional; one line per watershed)"),
            [QgsProcessing.TypeVectorLine], optional=True))
        self.addParameter(QgsProcessingParameterField(
            self.MIDLINES_ID, tr("TecDEM: watershed ID field of the drawn midlines (optional)"), None,
            self.MIDLINES, QgsProcessingParameterField.Any, optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            self.SEGMENT, tr("TecDEM: segment length in meters"), QgsProcessingParameterNumber.Double,
            2000.0, minValue=1.0))
        dd = QgsProcessingParameterEnum(
            self.DD_RULE, tr("TecDEM: Dd measured to"),
            [tr("Basin boundary on the river's side (Cox 1994)"),
             tr("First boundary crossing, as in the TecDEM code")], defaultValue=0)
        dd.setFlags(dd.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(dd)
        self.addParameter(QgsProcessingParameterFileDestination(
            self.OUTPUT, tr("Output workspace"),
            tr("GeoPackage (*.gpkg);;ESRI File Geodatabase (*.gdb)")))
        self.addOutput(QgsProcessingOutputMultipleLayers(self.LAYERS, tr("Output layers")))

    def processAlgorithm(self, parameters, context, feedback):
        dem = self.parameterAsRasterLayer(parameters, self.DEM, context)
        source = self.parameterAsSource(parameters, self.WATERSHEDS, context)
        if dem is None or source is None:
            raise QgsProcessingException(tr("Invalid inputs"))
        if dem.providerType() != "gdal":
            raise QgsProcessingException(tr("The DEM must be a file-based (GDAL) raster"))
        if dem.crs().isValid() and source.sourceCrs().isValid() and dem.crs() != source.sourceCrs():
            raise QgsProcessingException(
                tr("The watershed boundaries must use the same coordinate system as the DEM"))
        if source.sourceCrs().isGeographic():
            feedback.pushWarning(tr("The data use geographic coordinates; intervals are still read as "
                                    "map units. A projected coordinate system in meters is recommended."))
        field = self.parameterAsString(parameters, self.WATERSHED_ID, context) or None

        interval = self.parameterAsDouble(parameters, self.INTERVAL, context)
        if not interval:
            interval = 500.0  # default value in meters (as in the original script)
            feedback.pushWarning(tr("midline interval set to default value: 500 meters"))
        interval_profile = self.parameterAsDouble(parameters, self.MRPI, context)
        if not interval_profile:
            interval_profile = 500.0
            feedback.pushWarning(tr("river profile interval set to default value: 500 meters"))
        count = self.parameterAsInt(parameters, self.NUMBER_MRPI, context) or None
        adjust = self.parameterAsBoolean(parameters, self.ADJUST_EDGE, context)
        fill = self.parameterAsBoolean(parameters, self.FILL_SINKS, context)
        method = MIDLINE_METHODS[self.parameterAsEnum(parameters, self.MIDLINE_METHOD, context)]
        options = {"fill_sinks": fill, "midline_method": method}
        if method == "TECDEM":
            options["tecdem_segment"] = self.parameterAsDouble(parameters, self.SEGMENT, context)
            options["tecdem_dd_rule"] = DD_RULES[self.parameterAsEnum(parameters, self.DD_RULE, context)]
            midlines = self.parameterAsSource(parameters, self.MIDLINES, context)
            if midlines is not None:
                if midlines.sourceCrs().isValid() and midlines.sourceCrs() != source.sourceCrs():
                    raise QgsProcessingException(
                        tr("The midlines must use the same coordinate system as the watersheds"))
                options["user_midlines"] = list(midlines.getFeatures())
                options["midline_id_field"] = self.parameterAsString(parameters, self.MIDLINES_ID, context) or None
        out = self.parameterAsFileOutput(parameters, self.OUTPUT, context)
        if not out.lower().endswith((".gpkg", ".gdb")):
            out += ".gpkg"

        try:
            runner.run(dem.source().split("|")[0], list(source.getFeatures()), source.sourceCrs(), field,
                       interval, interval_profile, count, adjust, out, feedback, context.transformContext(),
                       options=options)
        except RuntimeError as err:
            raise QgsProcessingException(str(err))

        uris = []
        for name in OUTPUT_LAYERS + (TECDEM_LAYERS if method == "TECDEM" else []):
            uri = "%s|layername=%s" % (out, name)
            uris.append(uri)
            if context.willLoadLayerOnCompletion(uri) is False:
                context.addLayerToLoadOnCompletion(
                    uri, QgsProcessingContext.LayerDetails(name, context.project(), self.LAYERS))
        return {self.OUTPUT: out, self.LAYERS: uris}
