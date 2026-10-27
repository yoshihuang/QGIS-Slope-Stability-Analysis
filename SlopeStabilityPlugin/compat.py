# -*- coding: utf-8 -*-
"""Small compatibility layer: QGIS 3.28+ (Qt5) and QGIS 4.x (Qt6)."""

from qgis.PyQt.QtCore import QVariant
from qgis.core import QgsField, Qgis, QgsWkbTypes

try:  # Qt6 / QGIS >= 3.38
    from qgis.PyQt.QtCore import QMetaType
    _HAS_METATYPE = hasattr(QMetaType, "Type")
except ImportError:  # pragma: no cover
    QMetaType = None
    _HAS_METATYPE = False


def _qgis_int():
    try:
        return Qgis.versionInt()
    except Exception:
        return Qgis.QGIS_VERSION_INT


QGIS_VERSION_INT = _qgis_int()


def make_field(name, kind):
    """kind: 'double' | 'int' | 'string'."""
    if _HAS_METATYPE and QGIS_VERSION_INT >= 33800:
        t = {"double": QMetaType.Type.Double, "int": QMetaType.Type.Int,
             "string": QMetaType.Type.QString}[kind]
        return QgsField(name, t)
    t = {"double": QVariant.Double, "int": QVariant.Int, "string": QVariant.String}[kind]
    return QgsField(name, t)


def layer_filter(kind):
    """Filter flags for QgsMapLayerComboBox: 'polygon' | 'line' | 'vector'."""
    names = {"polygon": "PolygonLayer", "line": "LineLayer", "vector": "VectorLayer",
             "point": "PointLayer"}
    name = names[kind]
    lf = getattr(Qgis, "LayerFilter", None)
    if lf is not None and hasattr(lf, name):
        return getattr(lf, name)
    from qgis.core import QgsMapLayerProxyModel
    return getattr(QgsMapLayerProxyModel, name)


def geometry_kind(layer):
    """'polygon' | 'line' | 'point' | 'unknown'."""
    gt = layer.geometryType()
    try:
        gtypes = Qgis.GeometryType
        if gt == gtypes.Polygon:
            return "polygon"
        if gt == gtypes.Line:
            return "line"
        if gt == gtypes.Point:
            return "point"
    except AttributeError:
        pass
    try:
        if gt == QgsWkbTypes.GeometryType.PolygonGeometry:
            return "polygon"
        if gt == QgsWkbTypes.GeometryType.LineGeometry:
            return "line"
        if gt == QgsWkbTypes.GeometryType.PointGeometry:
            return "point"
    except AttributeError:
        pass
    return "unknown"


def exec_dialog(dlg):
    f = getattr(dlg, "exec", None) or getattr(dlg, "exec_")
    return f()


def enum(owner, scoped, name):
    """Return owner.scoped.name (Qt6) falling back to owner.name (Qt5)."""
    sc = getattr(owner, scoped, None)
    if sc is not None and hasattr(sc, name):
        return getattr(sc, name)
    return getattr(owner, name)
