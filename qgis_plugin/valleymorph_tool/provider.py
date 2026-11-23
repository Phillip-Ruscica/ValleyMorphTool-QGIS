import os

from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon

from .algorithm import ValleyMorphAlgorithm


class ValleyMorphProvider(QgsProcessingProvider):
    def id(self):
        return "valleymorph"

    def name(self):
        return "ValleyMorph Tool"

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(__file__), "icon.svg"))

    def loadAlgorithms(self):
        self.addAlgorithm(ValleyMorphAlgorithm())
