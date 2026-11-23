import os

from qgis.core import QgsApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction

from .provider import ValleyMorphProvider


class ValleyMorphPlugin:
    """Registers the ValleyMorph Tool in the Processing Toolbox and adds a menu/toolbar entry."""

    def __init__(self, iface):
        self.iface = iface
        self.provider = None
        self.action = None

    def initProcessing(self):
        self.provider = ValleyMorphProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):
        self.initProcessing()
        icon = QIcon(os.path.join(os.path.dirname(__file__), "icon.svg"))
        self.action = QAction(icon, "ValleyMorph Tool…", self.iface.mainWindow())
        self.action.triggered.connect(self.run)
        self.iface.addPluginToRasterMenu("ValleyMorph Tool", self.action)
        self.iface.addRasterToolBarIcon(self.action)

    def unload(self):
        if self.action:
            self.iface.removePluginRasterMenu("ValleyMorph Tool", self.action)
            self.iface.removeRasterToolBarIcon(self.action)
        if self.provider:
            QgsApplication.processingRegistry().removeProvider(self.provider)

    def run(self):
        import processing
        processing.execAlgorithmDialog("valleymorph:valleymorph")
