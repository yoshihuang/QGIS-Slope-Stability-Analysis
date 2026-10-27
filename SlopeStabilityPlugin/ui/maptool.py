# -*- coding: utf-8 -*-
"""Map tool: pick a crack / scarp point on the ground surface of the section."""

from qgis.PyQt.QtCore import pyqtSignal, Qt
from qgis.PyQt.QtGui import QColor
from qgis.gui import QgsMapToolEmitPoint, QgsVertexMarker

from ..compat import enum


class CrackPointTool(QgsMapToolEmitPoint):
    """Emits ``picked(QgsPointXY)`` (canvas CRS) on a left click; Esc / right
    click cancels (``cancelled``)."""

    picked = pyqtSignal(object)
    cancelled = pyqtSignal()

    def __init__(self, canvas):
        super().__init__(canvas)
        self.canvas = canvas
        self.marker = None
        self.setCursor(enum(Qt, "CursorShape", "CrossCursor"))

    def canvasReleaseEvent(self, e):
        if e.button() == enum(Qt, "MouseButton", "RightButton"):
            self.cancelled.emit()
            return
        pt = self.toMapCoordinates(e.pos())
        self.show_marker(pt)
        self.picked.emit(pt)

    def keyPressEvent(self, e):
        if e.key() == enum(Qt, "Key", "Key_Escape"):
            self.cancelled.emit()

    def show_marker(self, pt):
        if self.marker is None:
            self.marker = QgsVertexMarker(self.canvas)
            self.marker.setColor(QColor("#c0392b"))
            self.marker.setIconSize(14)
            self.marker.setPenWidth(3)
            self.marker.setIconType(QgsVertexMarker.IconType.ICON_TRIANGLE)
        self.marker.setCenter(pt)
        self.marker.show()

    def clear_marker(self):
        if self.marker is not None:
            self.canvas.scene().removeItem(self.marker)
            self.marker = None
