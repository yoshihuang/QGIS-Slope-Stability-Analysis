# -*- coding: utf-8 -*-
"""QGIS plugin controller."""

import os
import tempfile
import datetime

from qgis.PyQt.QtCore import Qt, QUrl
from qgis.PyQt.QtGui import QIcon, QDesktopServices
from qgis.PyQt.QtWidgets import QMessageBox, QProgressDialog
try:  # Qt6: QAction lives in QtGui
    from qgis.PyQt.QtGui import QAction
except ImportError:  # Qt5
    from qgis.PyQt.QtWidgets import QAction
from qgis.core import (QgsApplication, QgsProject, QgsVectorLayer, QgsMessageLog, Qgis)

from . import plugin_i18n
from .compat import enum
from .core.i18n import tr
from .core.search import SearchSettings
from .core.analysis import Scenario
from .core.report import build_report

PLUGIN_DIR = os.path.dirname(__file__)
LOG_TAG = "Slope Stability"


def plugin_version():
    try:
        with open(os.path.join(PLUGIN_DIR, "metadata.txt"), encoding="utf-8") as f:
            for line in f:
                if line.startswith("version="):
                    return line.split("=", 1)[1].strip()
    except OSError:
        return ""
    return ""


class SlopeStabilityPlugin(object):
    def __init__(self, iface):
        self.iface = iface
        self.actions = []
        self.toolbar = None
        self.provider = None
        self.task = None
        self.progress = None
        self.last_config = None
        self.dlg = None
        self.crack_tool = None
        self.crack_action = None
        plugin_i18n.install()

    # ----------------------------------------------------------------- GUI
    def initProcessing(self):
        from .processing_provider.provider import SlopeStabilityProvider
        self.provider = SlopeStabilityProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):
        self.initProcessing()
        self.menu = tr("&Slope Stability")
        icon = QIcon(os.path.join(PLUGIN_DIR, "icon.png"))
        self.toolbar = self.iface.addToolBar(tr("Slope Stability"))
        self.toolbar.setObjectName("SlopeStabilityToolBar")
        a = QAction(icon, tr("Slope stability analysis…"), self.iface.mainWindow())
        a.triggered.connect(self.run)
        self._add(a, toolbar=True)
        a = QAction(QIcon(os.path.join(PLUGIN_DIR, "icon_crack.svg")), tr("Pick a crack point and analyse…"), self.iface.mainWindow())
        a.setCheckable(True)
        a.triggered.connect(self.pick_crack)
        self.crack_action = a
        self._add(a, toolbar=True)
        a = QAction(tr("Load sample data"), self.iface.mainWindow())
        a.triggered.connect(self.load_sample_data)
        self._add(a)
        a = QAction(tr("User guide"), self.iface.mainWindow())
        a.triggered.connect(self.open_help)
        self._add(a)

    def _add(self, action, toolbar=False):
        self.iface.addPluginToMenu(self.menu, action)
        if toolbar:
            self.toolbar.addAction(action)
        self.actions.append(action)

    def unload(self):
        self._release_plugin_files()
        if self.crack_tool is not None:
            self.iface.mapCanvas().unsetMapTool(self.crack_tool)
            self.crack_tool.clear_marker()
            self.crack_tool = None
        if self.dlg is not None:
            self.dlg.close()
            self.dlg = None
        for a in self.actions:
            self.iface.removePluginMenu(self.menu, a)
            if self.toolbar is not None:
                self.toolbar.removeAction(a)
        self.actions = []
        if self.toolbar is not None:
            self.toolbar.deleteLater()
            self.toolbar = None
        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None
        plugin_i18n.uninstall()

    # ----------------------------------------------------------------- actions
    def open_help(self):
        lang = plugin_i18n.active_language()
        for name in ("user_guide_%s.html" % lang, "user_guide_en.html"):
            p = os.path.join(PLUGIN_DIR, "docs", name)
            if os.path.exists(p):
                QDesktopServices.openUrl(QUrl.fromLocalFile(p))
                return

    def _sample_folder(self):
        """Copy the sample data out of the plugin folder and return that folder.

        Layers opened directly from the plugin folder keep their files open
        (GDAL), which on Windows prevents QGIS from deleting the folder when
        the plugin is upgraded or uninstalled."""
        import shutil
        src = os.path.join(PLUGIN_DIR, "testdata")
        dst = os.path.join(QgsApplication.qgisSettingsDirPath(), "slope_stability_sample_data")
        try:
            os.makedirs(dst, exist_ok=True)
            for fn in os.listdir(src):
                a, b = os.path.join(src, fn), os.path.join(dst, fn)
                if os.path.isfile(a) and (not os.path.exists(b) or os.path.getmtime(b) < os.path.getmtime(a)):
                    shutil.copy2(a, b)
            return dst
        except Exception:
            return src

    def _release_plugin_files(self):
        """Remove project layers that read files inside the plugin folder."""
        root = os.path.normcase(os.path.abspath(PLUGIN_DIR))
        prj = QgsProject.instance()
        ids = []
        for lid, lyr in prj.mapLayers().items():
            path = lyr.source().split("|")[0]
            if path and os.path.normcase(os.path.abspath(path)).startswith(root):
                ids.append(lid)
        if ids:
            prj.removeMapLayers(ids)

    def load_sample_data(self):
        from qgis.core import QgsLineSymbol, QgsSingleSymbolRenderer
        d = self._sample_folder()
        group = QgsProject.instance().layerTreeRoot().insertGroup(0, tr("Slope stability sample data"))
        # (file, layer name, line style) - listed top to bottom in the layer tree
        samples = (
            ("testdata_user_slip_surfaces.geojson", tr("User-defined slip surfaces"),
             {"color": "#1f6fd1", "width": "0.6", "line_style": "dash"}),
            ("testdata_reinforcement.geojson", tr("Reinforcement"), {"color": "#2c3e50", "width": "0.8"}),
            ("testdata_distributed_loads.geojson", tr("Distributed loads"), {"color": "#8e44ad", "width": "1.4"}),
            ("retaining_wall.geojson", tr("Retaining wall"), {"color": "#7f8c8d", "width": "1.6"}),
            ("testdata_water_level_event.geojson", tr("Piezometric line"),
             {"color": "#1f6fd1", "width": "0.5", "line_style": "dot"}),
            ("testdata_geomodel.geojson", tr("Geology"), None),
        )
        for fname, label, style in samples:
            p = os.path.join(d, fname)
            if not os.path.exists(p):
                continue
            vl = QgsVectorLayer(p, label, "ogr")
            if not vl.isValid():
                continue
            if style:
                vl.setRenderer(QgsSingleSymbolRenderer(QgsLineSymbol.createSimple(style)))
            QgsProject.instance().addMapLayer(vl, False)
            group.addLayer(vl)
        self.iface.mapCanvas().refreshAllLayers()

    def run(self, crack_point=None, crack_crs=None):
        if self.task is not None:
            QMessageBox.information(self.iface.mainWindow(), tr("Slope Stability"), tr("An analysis is already running."))
            return
        if self.dlg is not None and self.dlg.isVisible():
            dlg = self.dlg
        else:
            from .ui.dialog import SlopeStabilityDialog
            dlg = SlopeStabilityDialog(self.iface.mainWindow(), self.iface)
            dlg.setModal(False)
            dlg.accepted.connect(lambda d=dlg: self.start(d.get_config()))
            self.dlg = dlg
        if crack_point is not None:
            dlg.set_crack_point(crack_point, crack_crs)
        dlg.show()
        dlg.raise_()
        dlg.activateWindow()

    def pick_crack(self, checked=True):
        """Toolbar tool: click the crack point, then the dialog opens."""
        from .ui.maptool import CrackPointTool
        canvas = self.iface.mapCanvas()
        if not checked:
            if self.crack_tool is not None:
                canvas.unsetMapTool(self.crack_tool)
            return
        self.crack_tool = CrackPointTool(canvas)
        self.crack_tool.setAction(self.crack_action)

        def picked(pt):
            canvas.unsetMapTool(self.crack_tool)
            self.crack_tool.clear_marker()
            self.run(pt, canvas.mapSettings().destinationCrs())
        self.crack_tool.picked.connect(picked)
        self.crack_tool.cancelled.connect(lambda: canvas.unsetMapTool(self.crack_tool))
        canvas.setMapTool(self.crack_tool)
        self.iface.messageBar().pushMessage(tr("Slope Stability"), tr("Click on the ground surface where the crack or scarp is. Right click or Esc cancels."), duration=6)

    # ----------------------------------------------------------------- run
    def start(self, cfg):
        from .qgis_io import build_model, InputError
        from .task import AnalysisTask
        try:
            model, user, msgs = build_model(cfg["layers"], cfg["opts"], cfg["selected_only"])
        except InputError as e:
            QMessageBox.warning(self.iface.mainWindow(), tr("Error"), u"%s" % e)
            return
        for level, m in msgs:
            QgsMessageLog.logMessage(m, LOG_TAG, Qgis.MessageLevel.Warning)
        search = dict(cfg["search"])
        search["user_surfaces"] = user
        if search["surface_type"] == "user" and not user:
            QMessageBox.warning(self.iface.mainWindow(), tr("Error"), tr("Select a layer with user-defined slip surfaces."))
            return
        settings = SearchSettings(**search)
        scenarios = []
        if cfg["main"]:
            scenarios.append(Scenario("Main", label=tr("Main analysis")))
        scan = cfg.get("scan")
        if scan:
            for kh in scan["kh"]:
                for d in scan["depths"]:
                    if d is None:
                        name = "Scan_Kh%g_Dry" % kh
                    else:
                        name = "Scan_Kh%g_D%gm" % (kh, d)
                    scenarios.append(Scenario(name, kh=kh, water_depth=d, density=scan["density"]))
        if not scenarios:
            return
        self.last_config = cfg
        self.last_config["model_msgs"] = msgs

        task = AnalysisTask(model, settings, scenarios, cfg["methods"], cfg["reliability"], self._done)
        self.progress = QProgressDialog(tr("Preparing the analysis…"), tr("Cancel"), 0, 100, self.iface.mainWindow())
        self.progress.setWindowTitle(tr("Slope stability analysis"))
        self.progress.setWindowModality(enum(Qt, "WindowModality", "NonModal"))
        self.progress.setMinimumDuration(0)
        self.progress.setMinimumWidth(420)
        self.progress.setValue(0)
        self.progress.canceled.connect(task.cancel)
        task.message.connect(self._on_progress)
        self.task = task
        QgsApplication.taskManager().addTask(task)

    def _on_progress(self, p, msg):
        if self.progress is not None:
            self.progress.setValue(p)
            if msg:
                self.progress.setLabelText(msg)

    def _done(self, task, ok):
        self.task = None
        if self.progress is not None:
            self.progress.close()
            self.progress = None
        win = self.iface.mainWindow()
        if not ok:
            if task.cancelled_by_user or task.isCanceled():
                self.iface.messageBar().pushMessage(tr("Slope Stability"), tr("Analysis cancelled."), level=_lvl("Info"), duration=5)
            else:
                QMessageBox.critical(win, tr("Analysis error"), task.error or tr("Unknown error"))
            return
        res = task.result
        errors = [m for lvl, m in res.preflight if lvl == "error"]
        if errors:
            QMessageBox.warning(win, tr("Input errors"), u"\n".join(errors))
            return
        try:
            self._publish(res)
        except Exception as e:
            import traceback
            QMessageBox.critical(win, tr("Analysis error"), u"%s\n\n%s" % (e, traceback.format_exc()))

    def _publish(self, res):
        from .qgis_io import result_rows, memory_layers, save_gpkg, style_layers
        cfg = self.last_config
        geo = cfg["layers"]["geology"]
        crs = geo.crs() if geo is not None else None
        rows = result_rows(res)
        layers = memory_layers(rows, crs)
        out = cfg.get("output")
        if out:
            layers = save_gpkg(layers, out)
        style_layers(layers)
        stamp = datetime.datetime.now().strftime("%H:%M:%S")
        root = QgsProject.instance().layerTreeRoot()
        group = root.insertGroup(0, u"%s %s" % (tr("Slope stability"), stamp))
        names = {"critical": tr("Critical slip surface"), "surfaces": tr("Slip surfaces"),
                 "slices": tr("Slices (critical surface)"), "centers": tr("Search grid centres"),
                 "crack": tr("Crack point")}
        for key in ("crack", "critical", "surfaces", "slices", "centers"):
            vl = layers.get(key)
            if vl is None:
                continue
            vl.setName(names[key])
            QgsProject.instance().addMapLayer(vl, False)
            node = group.addLayer(vl)
            if key == "centers" and node is not None:
                node.setItemVisibilityChecked(False)

        report_path = None
        if cfg.get("report"):
            if out:
                report_path = os.path.splitext(out)[0] + "_report.html"
            else:
                report_path = os.path.join(tempfile.gettempdir(), "slope_stability_report_%s.html" % datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            res.preflight = list(cfg.get("model_msgs", [])) + res.preflight
            html = build_report(res, {"plugin": "Slope Stability Analysis", "version": plugin_version(),
                                      "project": geo.name() if geo is not None else ""})
            with open(report_path, "w", encoding="utf-8") as f:
                f.write(html)

        from .ui.result_dialog import ResultDialog
        dlg = ResultDialog(self.iface.mainWindow(), res, report_path)
        dlg.show()
        self._result_dialog = dlg
        if report_path and cfg.get("open_report"):
            QDesktopServices.openUrl(QUrl.fromLocalFile(report_path))


def _lvl(name):
    ml = getattr(Qgis, "MessageLevel", None)
    if ml is not None and hasattr(ml, name):
        return getattr(ml, name)
    return getattr(Qgis, name)
