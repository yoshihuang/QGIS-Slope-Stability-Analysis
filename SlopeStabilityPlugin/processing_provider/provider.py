# -*- coding: utf-8 -*-
import os

from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon

from ..core.i18n import tr


class SlopeStabilityProvider(QgsProcessingProvider):
    def loadAlgorithms(self):
        from .algorithm import SlopeStabilityAlgorithm
        self.addAlgorithm(SlopeStabilityAlgorithm())

    def id(self):
        return "slopestability"

    def name(self):
        return tr("Slope Stability")

    def longName(self):
        return tr("Slope Stability (limit equilibrium)")

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(os.path.dirname(__file__)), "icon.png"))
