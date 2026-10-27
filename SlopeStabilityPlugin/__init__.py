# -*- coding: utf-8 -*-
"""QGIS Slope Stability Analysis plugin."""


def classFactory(iface):
    try:
        import numpy  # noqa: F401
        import shapely  # noqa: F401
    except ImportError as e:
        from qgis.PyQt.QtWidgets import QMessageBox
        from qgis.PyQt.QtCore import QCoreApplication

        def tr(s):
            return QCoreApplication.translate("SlopeStability", s)
        QMessageBox.critical(
            None, tr("Missing dependencies"),
            tr("The Slope Stability plugin needs the Python libraries numpy and shapely, "
               "which are included in the standard QGIS installers (OSGeo4W, macOS, Flatpak). "
               "On Linux install the package python3-shapely.") + "\n\n" + str(e))
        raise

    from .plugin import SlopeStabilityPlugin
    return SlopeStabilityPlugin(iface)
