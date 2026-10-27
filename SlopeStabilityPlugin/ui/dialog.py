# -*- coding: utf-8 -*-
"""Main analysis dialog."""

import os

from qgis.PyQt.QtCore import Qt
from qgis.PyQt.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
                                 QDoubleSpinBox, QSpinBox, QPushButton, QFileDialog,
                                 QGroupBox, QFormLayout, QLineEdit, QCheckBox, QTabWidget,
                                 QWidget, QTextBrowser, QGridLayout, QMessageBox,
                                 QDialogButtonBox, QScrollArea)
from qgis.core import QgsSettings
from qgis.gui import QgsMapLayerComboBox

from ..compat import layer_filter
from ..core.methods import METHODS, method_label
from ..core.i18n import tr
from .. import plugin_i18n

S_PREFIX = "SlopeStability/ui/"


def _spin(lo, hi, val, step=0.1, dec=3, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(dec)
    s.setSingleStep(step)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    return s


def _ispin(lo, hi, val, step=1, suffix=""):
    s = QSpinBox()
    s.setRange(lo, hi)
    s.setSingleStep(step)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    return s


def _layer_combo(kind, optional):
    cb = QgsMapLayerComboBox()
    cb.setFilters(layer_filter(kind))
    if optional:
        try:
            cb.setAllowEmptyLayer(True, tr("-- not used --"))
        except TypeError:
            cb.setAllowEmptyLayer(True)
        cb.setLayer(None)
    return cb


def parse_list(text, allow_dry=False):
    out = []
    for tok in text.replace(";", ",").split(","):
        t = tok.strip().lower()
        if not t:
            continue
        if allow_dry and t in ("dry", "none", "-", "n"):
            out.append(None)
            continue
        out.append(float(t))
    return out


class SlopeStabilityDialog(QDialog):
    def __init__(self, parent=None, iface=None):
        super().__init__(parent)
        self.iface = iface
        self.setWindowTitle(tr("Slope Stability Analysis (Limit Equilibrium)"))
        self.resize(640, 680)
        self.settings = QgsSettings()
        self._build()
        self._load_settings()
        self._update_state()
        self.finished.connect(self._clear_marker)

    def _clear_marker(self, *args):
        if getattr(self, "_tool", None) is not None:
            self._tool.clear_marker()

    # ------------------------------------------------------------------ UI
    def _build(self):
        main = QVBoxLayout(self)
        self.tabs = QTabWidget()
        main.addWidget(self.tabs)
        self.tabs.addTab(self._tab_input(), tr("Input"))
        self.tabs.addTab(self._tab_loads(), tr("Loads && supports"))
        self.tabs.addTab(self._tab_analysis(), tr("Analysis"))
        self.tabs.addTab(self._tab_scenarios(), tr("Scenarios"))
        self.tabs.addTab(self._tab_reliability(), tr("Reliability"))
        self.tabs.addTab(self._tab_output(), tr("Output"))

        btns = QHBoxLayout()
        self.btn_check = QPushButton(tr("Check inputs"))
        self.btn_check.clicked.connect(self.check_inputs)
        btns.addWidget(self.btn_check)
        btns.addStretch()
        self.btn_close = QPushButton(tr("Close"))
        self.btn_run = QPushButton(tr("Run analysis"))
        self.btn_run.setDefault(True)
        self.btn_run.setStyleSheet("background-color:#2b5c8f;color:white;font-weight:bold;padding:6px 14px;")
        self.btn_close.clicked.connect(self.reject)
        self.btn_run.clicked.connect(self._accept)
        btns.addWidget(self.btn_close)
        btns.addWidget(self.btn_run)
        main.addLayout(btns)

    def _tab_input(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        g = QGroupBox(tr("Section layers (x = distance, y = elevation)"))
        f = QFormLayout(g)
        self.cb_geo = _layer_combo("polygon", False)
        self.cb_water = _layer_combo("line", True)
        self.cb_loads = _layer_combo("line", True)
        self.cb_reinf = _layer_combo("line", True)
        self.cb_piles = _layer_combo("vector", True)
        self.cb_surf = _layer_combo("line", True)
        self.chk_selected = QCheckBox(tr("Use selected geology features only"))
        f.addRow(tr("Geology / materials (polygon) *"), self.cb_geo)
        f.addRow(tr("Piezometric line (line)"), self.cb_water)
        f.addRow(tr("Distributed loads (line)"), self.cb_loads)
        f.addRow(tr("Reinforcement (line)"), self.cb_reinf)
        f.addRow(tr("Piles / retaining walls"), self.cb_piles)
        f.addRow(tr("User-defined slip surfaces (line)"), self.cb_surf)
        f.addRow("", self.chk_selected)
        lay.addWidget(g)
        self.info = QTextBrowser()
        self.info.setOpenExternalLinks(True)
        self.info.setHtml(self._field_help())
        lay.addWidget(self.info, 1)
        return w

    def _field_help(self):
        rows = [
            (tr("Geology"), "unit_weight, cohesion, phi", "name / layer_name, unit_weight_sat, model (mc, undrained, hb), sigci, gsi, mi, d, ru, piezo, c_std, phi_std, gamma_std"),
            (tr("Distributed loads"), "-", "q (kPa)"),
            (tr("Reinforcement"), "-", "T (kN/m), pullout (kN/m/m), mode (active / passive)"),
            (tr("Piles / retaining walls"), "-", "V (kN/m), mode (active / passive)"),
        ]
        h = [u"<p><b>%s</b></p>" % tr("Attribute fields recognised on the input layers (names are not case-sensitive):")]
        h.append(u"<table border='1' cellspacing='0' cellpadding='3'><tr><th>%s</th><th>%s</th><th>%s</th></tr>" % (tr("Layer"), tr("Required"), tr("Optional")))
        for a, b, c in rows:
            h.append(u"<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (a, b, c))
        h.append(u"</table><p>%s</p>" % tr("Units: kPa, kN/m3, degrees, metres. Fields missing on load, reinforcement or pile layers use the default values on the Loads & supports tab."))
        return u"".join(h)

    def _tab_loads(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        g = QGroupBox(tr("Seismic and water"))
        f = QFormLayout(g)
        self.sp_kh = _spin(0, 1, 0.0, 0.05, 3)
        self.sp_kh.setToolTip(tr("Horizontal pseudo-static seismic coefficient, acting in the direction of sliding."))
        self.sp_kv = _spin(-1, 1, 0.0, 0.05, 3)
        self.sp_kv.setToolTip(tr("Vertical seismic coefficient (positive = downward)."))
        self.sp_gw = _spin(1, 100, 9.81, 0.01, 3, " kN/m3")
        f.addRow(tr("Horizontal seismic coefficient kh"), self.sp_kh)
        f.addRow(tr("Vertical seismic coefficient kv"), self.sp_kv)
        f.addRow(tr("Unit weight of water"), self.sp_gw)
        lay.addWidget(g)

        g = QGroupBox(tr("Tension crack"))
        f = QFormLayout(g)
        self.sp_crack = _spin(0, 1000, 0.0, 0.5, 2, " m")
        self.sp_fill = _ispin(0, 100, 0, 10, " %")
        f.addRow(tr("Crack depth"), self.sp_crack)
        f.addRow(tr("Water in crack"), self.sp_fill)
        lay.addWidget(g)

        g = QGroupBox(tr("Default values (used when the layer has no attribute)"))
        f = QFormLayout(g)
        self.sp_q = _spin(0, 1e5, 0.0, 5, 2, " kPa")
        self.sp_T = _spin(0, 1e6, 100.0, 10, 1, " kN/m")
        self.sp_po = _spin(0, 1e6, 0.0, 10, 1, " kN/m/m")
        self.sp_po.setToolTip(tr("Pull-out resistance per metre of reinforcement behind the slip surface. 0 = tensile capacity governs."))
        self.cb_rmode = QComboBox()
        self.cb_rmode.addItem(tr("Passive (factored by FS)"), True)
        self.cb_rmode.addItem(tr("Active (reduces driving force)"), False)
        self.sp_V = _spin(0, 1e6, 0.0, 10, 1, " kN/m")
        self.sp_V.setToolTip(tr("Shear resistance of a pile or retaining wall where the slip surface cuts it (v2.x 'wall force')."))
        self.cb_pmode = QComboBox()
        self.cb_pmode.addItem(tr("Passive (factored by FS)"), True)
        self.cb_pmode.addItem(tr("Active (reduces driving force)"), False)
        self.chk_global = QCheckBox(tr("No pile layer: apply the shear force at the toe of every surface (v2.x mode)"))
        f.addRow(tr("Surcharge q"), self.sp_q)
        f.addRow(tr("Reinforcement capacity T"), self.sp_T)
        f.addRow(tr("Reinforcement pull-out resistance"), self.sp_po)
        f.addRow(tr("Reinforcement force mode"), self.cb_rmode)
        f.addRow(tr("Pile / wall shear force V"), self.sp_V)
        f.addRow(tr("Pile / wall force mode"), self.cb_pmode)
        f.addRow("", self.chk_global)
        lay.addWidget(g)
        lay.addStretch()
        return w

    def _tab_analysis(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        g = QGroupBox(tr("Slip surface search"))
        f = QFormLayout(g)
        self.cb_type = QComboBox()
        self.cb_type.addItem(tr("Circular - grid search"), "grid")
        self.cb_type.addItem(tr("Circular - entry and exit"), "entry_exit")
        self.cb_type.addItem(tr("Non-circular - path search"), "path")
        self.cb_type.addItem(tr("User-defined surfaces only"), "user")
        self.cb_type.addItem(tr("From a crack point picked on the map"), "point")
        self.cb_type.currentIndexChanged.connect(self._update_state)
        self.cb_method = QComboBox()
        for m in METHODS:
            self.cb_method.addItem(method_label(m), m)
        self.cb_method.setCurrentIndex(METHODS.index("bishop"))
        self.cb_method.setToolTip(tr("Method used to evaluate every trial surface during the search."))
        self.sp_density = _ispin(25, 400, 100, 25, " %")
        self.sp_density.setToolTip(tr("100 % = standard density. Increase for a finer search, decrease for a quick test."))
        self.sp_slices = _ispin(5, 200, 30, 5)
        self.chk_opt = QCheckBox(tr("Optimise the critical surface"))
        self.chk_opt.setChecked(True)
        self.chk_bound = QCheckBox(tr("Ignore surfaces ending at the left/right model boundary"))
        self.chk_bound.setChecked(True)
        self.cb_dir = QComboBox()
        self.cb_dir.addItem(tr("Automatic"), 0)
        self.cb_dir.addItem(tr("Slides towards +x (left to right)"), 1)
        self.cb_dir.addItem(tr("Slides towards -x (right to left)"), -1)
        self.sp_vert = _ispin(3, 12, 6)
        self.sp_mindepth = _spin(0, 1000, 0.0, 0.5, 2, " m")
        self.sp_mindepth.setToolTip(tr("Surfaces shallower than this are ignored (0 = no limit)."))
        f.addRow(tr("Surface type"), self.cb_type)
        f.addRow(tr("Search method"), self.cb_method)
        f.addRow(tr("Search density"), self.sp_density)
        f.addRow(tr("Number of slices"), self.sp_slices)
        f.addRow(tr("Path vertices (non-circular)"), self.sp_vert)
        f.addRow(tr("Minimum slide depth"), self.sp_mindepth)
        f.addRow(tr("Slope direction"), self.cb_dir)
        f.addRow("", self.chk_opt)
        f.addRow("", self.chk_bound)
        lay.addWidget(g)

        g = QGroupBox(tr("Crack point (slip surfaces start here and extend downslope)"))
        self.grp_crack = g
        f = QGridLayout(g)
        self.sp_cx = _spin(-1e9, 1e9, 0, 1, 2)
        self.sp_cy = _spin(-1e9, 1e9, 0, 1, 2)
        self.sp_cy.setToolTip(tr("For information only: the point is projected vertically onto the ground surface."))
        self.btn_pick = QPushButton(tr("Pick on map"))
        self.btn_pick.clicked.connect(self.pick_crack_point)
        self.lbl_pick = QLabel(tr("Click on the ground surface where the crack or scarp is. Right click or Esc cancels."))
        self.lbl_pick.setWordWrap(True)
        self._crack_picked = False
        f.addWidget(QLabel("x"), 0, 0)
        f.addWidget(self.sp_cx, 0, 1)
        f.addWidget(QLabel("y"), 0, 2)
        f.addWidget(self.sp_cy, 0, 3)
        f.addWidget(self.btn_pick, 0, 4)
        f.addWidget(self.lbl_pick, 1, 0, 1, 5)
        lay.addWidget(g)

        g = QGroupBox(tr("Methods computed on the critical surface"))
        grid = QGridLayout(g)
        self.chk_methods = {}
        for i, m in enumerate(METHODS):
            c = QCheckBox(method_label(m))
            c.setChecked(m in ("ordinary", "bishop", "janbu", "spencer", "morgenstern_price"))
            self.chk_methods[m] = c
            grid.addWidget(c, i // 2, i % 2)
        lay.addWidget(g)

        g = QGroupBox(tr("Search limits (world coordinates)"))
        f = QGridLayout(g)
        self.chk_auto_ee = QCheckBox(tr("Automatic entry / exit ranges"))
        self.chk_auto_ee.setChecked(True)
        self.chk_auto_ee.toggled.connect(self._update_state)
        self.sp_e0 = _spin(-1e9, 1e9, 0, 1, 2)
        self.sp_e1 = _spin(-1e9, 1e9, 0, 1, 2)
        self.sp_x0 = _spin(-1e9, 1e9, 0, 1, 2)
        self.sp_x1 = _spin(-1e9, 1e9, 0, 1, 2)
        f.addWidget(self.chk_auto_ee, 0, 0, 1, 4)
        f.addWidget(QLabel(tr("Entry x from")), 1, 0)
        f.addWidget(self.sp_e0, 1, 1)
        f.addWidget(QLabel(tr("to")), 1, 2)
        f.addWidget(self.sp_e1, 1, 3)
        f.addWidget(QLabel(tr("Exit x from")), 2, 0)
        f.addWidget(self.sp_x0, 2, 1)
        f.addWidget(QLabel(tr("to")), 2, 2)
        f.addWidget(self.sp_x1, 2, 3)
        self.chk_auto_grid = QCheckBox(tr("Automatic centre grid"))
        self.chk_auto_grid.setChecked(True)
        self.chk_auto_grid.toggled.connect(self._update_state)
        self.sp_gx0 = _spin(-1e9, 1e9, 0, 1, 2)
        self.sp_gx1 = _spin(-1e9, 1e9, 0, 1, 2)
        self.sp_gy0 = _spin(-1e9, 1e9, 0, 1, 2)
        self.sp_gy1 = _spin(-1e9, 1e9, 0, 1, 2)
        f.addWidget(self.chk_auto_grid, 3, 0, 1, 4)
        f.addWidget(QLabel(tr("Centre x from")), 4, 0)
        f.addWidget(self.sp_gx0, 4, 1)
        f.addWidget(QLabel(tr("to")), 4, 2)
        f.addWidget(self.sp_gx1, 4, 3)
        f.addWidget(QLabel(tr("Centre y from")), 5, 0)
        f.addWidget(self.sp_gy0, 5, 1)
        f.addWidget(QLabel(tr("to")), 5, 2)
        f.addWidget(self.sp_gy1, 5, 3)
        self.btn_fill = QPushButton(tr("Fill from geology extent"))
        self.btn_fill.clicked.connect(self._fill_limits)
        f.addWidget(self.btn_fill, 6, 0, 1, 4)
        lay.addWidget(g)
        lay.addStretch()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(w)
        return scroll

    def _tab_scenarios(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        self.chk_main = QCheckBox(tr("Main analysis (kh and piezometric line from the previous tabs)"))
        self.chk_main.setChecked(True)
        lay.addWidget(self.chk_main)
        g = QGroupBox(tr("Parametric scan (kh x water table depth)"))
        g.setCheckable(True)
        g.setChecked(False)
        self.grp_scan = g
        f = QFormLayout(g)
        self.ed_kh = QLineEdit("0, 0.1, 0.2")
        self.ed_depth = QLineEdit("0, 15, 30")
        self.ed_depth.setToolTip(tr("Water table parallel to the ground surface at these depths. Use 'dry' for no water."))
        self.sp_scan_density = _ispin(25, 400, 75, 25, " %")
        f.addRow(tr("kh values"), self.ed_kh)
        f.addRow(tr("Water table depths (m)"), self.ed_depth)
        f.addRow(tr("Search density"), self.sp_scan_density)
        lab = QLabel(tr("Equivalent to the 'General Scan' mode of version 2.x: every combination is searched and stored with its scenario name."))
        lab.setWordWrap(True)
        f.addRow(lab)
        lay.addWidget(g)
        lay.addStretch()
        return w

    def _tab_reliability(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        g = QGroupBox(tr("Probabilistic analysis of the critical surface"))
        g.setCheckable(True)
        g.setChecked(False)
        self.grp_rel = g
        f = QFormLayout(g)
        self.sp_mc = _ispin(0, 100000, 1000, 100)
        self.chk_tspm = QCheckBox(tr("Taylor series method (TSPM)"))
        self.chk_tspm.setChecked(True)
        self.sp_seed = _ispin(0, 999999, 1)
        f.addRow(tr("Monte Carlo samples"), self.sp_mc)
        f.addRow("", self.chk_tspm)
        f.addRow(tr("Random seed"), self.sp_seed)
        lab = QLabel(tr("Random variables are read from the geology layer fields c_std, phi_std and gamma_std (standard deviations, normal distribution)."))
        lab.setWordWrap(True)
        f.addRow(lab)
        lay.addWidget(g)
        lay.addStretch()
        return w

    def _tab_output(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        g = QGroupBox(tr("Results"))
        f = QFormLayout(g)
        row = QHBoxLayout()
        self.ed_out = QLineEdit()
        self.ed_out.setPlaceholderText(tr("Empty = temporary layers"))
        b = QPushButton(u"…")
        b.clicked.connect(self._browse)
        row.addWidget(self.ed_out)
        row.addWidget(b)
        f.addRow(tr("GeoPackage file"), row)
        self.sp_keep = _ispin(1, 5000, 300, 50)
        self.sp_keep.setToolTip(tr("Number of lowest-FS surfaces stored per scenario."))
        f.addRow(tr("Surfaces to keep"), self.sp_keep)
        self.chk_report = QCheckBox(tr("Create HTML report"))
        self.chk_report.setChecked(True)
        self.chk_open = QCheckBox(tr("Open the report when finished"))
        self.chk_open.setChecked(True)
        f.addRow("", self.chk_report)
        f.addRow("", self.chk_open)
        lay.addWidget(g)
        g = QGroupBox(tr("Language"))
        f = QFormLayout(g)
        self.cb_lang = QComboBox()
        for code, name in plugin_i18n.LANGUAGES:
            self.cb_lang.addItem(tr("Automatic (QGIS language)") if code == "auto" else name, code)
        cur = plugin_i18n.configured_language()
        i = self.cb_lang.findData(cur)
        self.cb_lang.setCurrentIndex(max(i, 0))
        self.cb_lang.currentIndexChanged.connect(self._language_changed)
        f.addRow(tr("Interface and report language"), self.cb_lang)
        lay.addWidget(g)
        lay.addStretch()
        return w

    # ------------------------------------------------------------- behaviour
    def _language_changed(self):
        code = self.cb_lang.currentData()
        plugin_i18n.set_configured_language(code)
        plugin_i18n.install(code)
        QMessageBox.information(self, tr("Language"), tr("The new language is used the next time the dialog is opened."))

    def _update_state(self):
        t = self.cb_type.currentData()
        auto_ee = self.chk_auto_ee.isChecked()
        for s in (self.sp_e0, self.sp_e1, self.sp_x0, self.sp_x1):
            s.setEnabled(not auto_ee and t in ("entry_exit", "path"))
        self.chk_auto_ee.setEnabled(t in ("entry_exit", "path"))
        auto_g = self.chk_auto_grid.isChecked()
        for s in (self.sp_gx0, self.sp_gx1, self.sp_gy0, self.sp_gy1):
            s.setEnabled(not auto_g and t == "grid")
        self.chk_auto_grid.setEnabled(t == "grid")
        self.sp_vert.setEnabled(t in ("path", "point"))
        self.grp_crack.setEnabled(t == "point")

    def _fill_limits(self):
        lyr = self.cb_geo.currentLayer()
        if lyr is None:
            return
        e = lyr.extent()
        x0, x1, y0, y1 = e.xMinimum(), e.xMaximum(), e.yMinimum(), e.yMaximum()
        w = x1 - x0
        h = y1 - y0
        self.sp_e0.setValue(x0); self.sp_e1.setValue(x0 + 0.5 * w)
        self.sp_x0.setValue(x0 + 0.5 * w); self.sp_x1.setValue(x1)
        self.sp_gx0.setValue(x0 + 0.1 * w); self.sp_gx1.setValue(x1 - 0.1 * w)
        self.sp_gy0.setValue(y1 + 0.1 * h); self.sp_gy1.setValue(y1 + 1.5 * h)

    # ------------------------------------------------------------- crack point
    def set_crack_point(self, pt, crs=None):
        """Set the crack point from a map point (``crs`` = its CRS)."""
        lyr = self.cb_geo.currentLayer()
        if crs is not None and lyr is not None and crs.isValid() and lyr.crs().isValid() and crs != lyr.crs():
            from qgis.core import QgsCoordinateTransform, QgsProject
            try:
                pt = QgsCoordinateTransform(crs, lyr.crs(), QgsProject.instance()).transform(pt)
            except Exception:
                pass
        self.sp_cx.setValue(pt.x())
        self.sp_cy.setValue(pt.y())
        self._crack_picked = True
        i = self.cb_type.findData("point")
        self.cb_type.setCurrentIndex(i)
        self.tabs.setCurrentIndex(2)
        self._update_state()

    def pick_crack_point(self):
        if self.iface is None:
            return
        from .maptool import CrackPointTool
        if getattr(self, "_tool", None) is not None:
            self._tool.clear_marker()
        canvas = self.iface.mapCanvas()
        self._prev_tool = canvas.mapTool()
        self._tool = CrackPointTool(canvas)

        def done(pt=None):
            canvas.unsetMapTool(self._tool)
            if self._prev_tool is not None:
                canvas.setMapTool(self._prev_tool)
            if pt is not None:
                self.set_crack_point(pt, canvas.mapSettings().destinationCrs())
            self.show()
            self.raise_()
            self.activateWindow()
        self._tool.picked.connect(lambda pt: done(pt))
        self._tool.cancelled.connect(lambda: done(None))
        self.hide()
        canvas.setMapTool(self._tool)
        if self.iface is not None:
            self.iface.messageBar().pushMessage(tr("Slope Stability"), tr("Click on the ground surface where the crack or scarp is. Right click or Esc cancels."), duration=6)

    def _browse(self):
        path, _ = QFileDialog.getSaveFileName(self, tr("Save results"), self.ed_out.text(), "GeoPackage (*.gpkg)")
        if path:
            if not path.lower().endswith(".gpkg"):
                path += ".gpkg"
            self.ed_out.setText(path)

    def _accept(self):
        if self.cb_geo.currentLayer() is None:
            QMessageBox.warning(self, tr("Error"), tr("Please select the geology layer."))
            return
        try:
            self.get_config()
        except ValueError as e:
            QMessageBox.warning(self, tr("Error"), u"%s\n%s" % (tr("Invalid list of values in the Scenarios tab."), e))
            return
        if self.cb_type.currentData() == "point" and not self._crack_picked:
            QMessageBox.warning(self, tr("Error"), tr("Pick the crack point on the map first (Analysis tab)."))
            return
        if not self.chk_main.isChecked() and not self.grp_scan.isChecked():
            QMessageBox.warning(self, tr("Error"), tr("Select the main analysis and/or the parametric scan."))
            return
        self._save_settings()
        self.accept()

    def check_inputs(self):
        from ..qgis_io import build_model, InputError
        from ..core.methods import METHODS
        cfg = self.get_config()
        try:
            model, user, msgs = build_model(cfg["layers"], cfg["opts"], cfg["selected_only"])
            model.prepare()
            msgs = msgs + model.preflight()
        except InputError as e:
            self.info.setHtml(u"<p style='color:#c0392b'><b>%s</b></p>" % e)
            return
        except Exception as e:
            self.info.setHtml(u"<p style='color:#c0392b'><b>%s</b></p>" % e)
            return
        h = [u"<h3>%s</h3>" % tr("Input checks")]
        if not msgs:
            h.append(u"<p style='color:#27ae60'>%s</p>" % tr("No problems found."))
        for level, m in msgs:
            col = {"error": "#c0392b", "warning": "#b9770e"}.get(level, "#555")
            h.append(u"<p style='color:%s'>&bull; %s</p>" % (col, m))
        h.append(u"<h3>%s</h3><table border='1' cellspacing='0' cellpadding='3'>" % tr("Materials"))
        h.append(u"<tr><th>%s</th><th>%s</th><th>&gamma;</th><th>&gamma;sat</th><th>c/Su</th><th>&phi;</th></tr>" % (tr("Material"), tr("Strength model")))
        for mt in model.materials:
            h.append(u"<tr><td>%s</td><td>%s</td><td>%g</td><td>%g</td><td>%g</td><td>%g</td></tr>" % (mt.name, mt.model, mt.gamma, mt.gamma_sat, mt.cohesion, mt.phi))
        h.append(u"</table>")
        h.append(u"<p>%s: %s &nbsp; %s: %d &nbsp; %s: %d &nbsp; %s: %d &nbsp; %s: %d</p>" % (
            tr("Piezometric line"), tr("Yes") if model.piezo is not None else tr("No"),
            tr("Distributed loads"), len(model.loads), tr("Reinforcements"), len(model.reinforcements),
            tr("Piles / walls"), len(model.piles), tr("User-defined slip surfaces"), len(user)))
        h.append(u"<p>%s: %s</p>" % (tr("Slope direction"), tr("Slides towards +x") if model.sign > 0 else tr("Slides towards -x")))
        self.info.setHtml(u"".join(h))

    # ------------------------------------------------------------- config
    def get_config(self):
        t = self.cb_type.currentData()
        layers = {
            "geology": self.cb_geo.currentLayer(),
            "water": self.cb_water.currentLayer(),
            "loads": self.cb_loads.currentLayer(),
            "reinforcement": self.cb_reinf.currentLayer(),
            "piles": self.cb_piles.currentLayer(),
            "surfaces": self.cb_surf.currentLayer(),
        }
        opts = {
            "kh": self.sp_kh.value(), "kv": self.sp_kv.value(), "gamma_w": self.sp_gw.value(),
            "crack_depth": self.sp_crack.value(), "crack_fill": self.sp_fill.value() / 100.0,
            "default_q": self.sp_q.value(), "default_T": self.sp_T.value(),
            "default_pullout": self.sp_po.value(), "reinf_passive": self.cb_rmode.currentData(),
            "default_V": self.sp_V.value(), "pile_passive": self.cb_pmode.currentData(),
            "global_support": self.sp_V.value() if self.chk_global.isChecked() else 0.0,
            "direction": self.cb_dir.currentData(), "min_depth": self.sp_mindepth.value(),
        }
        search = {
            "surface_type": t, "method": self.cb_method.currentData(),
            "density": self.sp_density.value() / 100.0, "n_slices": self.sp_slices.value(),
            "optimize": self.chk_opt.isChecked(), "exclude_boundary": self.chk_bound.isChecked(),
            "n_keep": self.sp_keep.value(), "path_vertices": self.sp_vert.value(),
        }
        if not self.chk_auto_ee.isChecked() and t in ("entry_exit", "path"):
            search["entry_range"] = (self.sp_e0.value(), self.sp_e1.value())
            search["exit_range"] = (self.sp_x0.value(), self.sp_x1.value())
        if t == "point":
            search["crack_point"] = (self.sp_cx.value(), self.sp_cy.value())
        if not self.chk_auto_grid.isChecked() and t == "grid":
            search["grid"] = (self.sp_gx0.value(), self.sp_gx1.value(), self.sp_gy0.value(), self.sp_gy1.value())
        scan = None
        if self.grp_scan.isChecked():
            scan = {"kh": parse_list(self.ed_kh.text()),
                    "depths": parse_list(self.ed_depth.text(), allow_dry=True),
                    "density": self.sp_scan_density.value() / 100.0}
        rel = None
        if self.grp_rel.isChecked():
            rel = {"mc_samples": self.sp_mc.value(), "tspm": self.chk_tspm.isChecked(), "seed": self.sp_seed.value()}
        return {
            "layers": layers, "opts": opts, "search": search,
            "methods": [m for m, c in self.chk_methods.items() if c.isChecked()],
            "main": self.chk_main.isChecked(), "scan": scan, "reliability": rel,
            "output": self.ed_out.text().strip(), "report": self.chk_report.isChecked(),
            "open_report": self.chk_open.isChecked(), "selected_only": self.chk_selected.isChecked(),
        }

    # ------------------------------------------------------------- settings
    def _widgets(self):
        return {
            "kh": self.sp_kh, "kv": self.sp_kv, "gw": self.sp_gw, "crack": self.sp_crack, "fill": self.sp_fill,
            "q": self.sp_q, "T": self.sp_T, "po": self.sp_po, "V": self.sp_V, "density": self.sp_density,
            "slices": self.sp_slices, "vert": self.sp_vert, "mindepth": self.sp_mindepth,
            "e0": self.sp_e0, "e1": self.sp_e1, "x0": self.sp_x0, "x1": self.sp_x1,
            "gx0": self.sp_gx0, "gx1": self.sp_gx1, "gy0": self.sp_gy0, "gy1": self.sp_gy1,
            "scan_density": self.sp_scan_density, "mc": self.sp_mc, "seed": self.sp_seed, "keep": self.sp_keep,
        }

    def _save_settings(self):
        s = self.settings
        for k, wdg in self._widgets().items():
            s.setValue(S_PREFIX + k, wdg.value())
        for k, cb in (("type", self.cb_type), ("method", self.cb_method), ("dir", self.cb_dir),
                      ("rmode", self.cb_rmode), ("pmode", self.cb_pmode)):
            s.setValue(S_PREFIX + k, cb.currentIndex())
        for k, c in (("opt", self.chk_opt), ("bound", self.chk_bound), ("auto_ee", self.chk_auto_ee),
                     ("auto_grid", self.chk_auto_grid), ("main", self.chk_main), ("report", self.chk_report),
                     ("open", self.chk_open), ("tspm", self.chk_tspm), ("global", self.chk_global),
                     ("scan", self.grp_scan), ("rel", self.grp_rel)):
            s.setValue(S_PREFIX + k, c.isChecked())
        for m, c in self.chk_methods.items():
            s.setValue(S_PREFIX + "m_" + m, c.isChecked())
        s.setValue(S_PREFIX + "ed_kh", self.ed_kh.text())
        s.setValue(S_PREFIX + "ed_depth", self.ed_depth.text())
        s.setValue(S_PREFIX + "out", self.ed_out.text())

    def _load_settings(self):
        s = self.settings
        for k, wdg in self._widgets().items():
            v = s.value(S_PREFIX + k, None)
            if v is not None:
                try:
                    wdg.setValue(type(wdg.value())(float(v)))
                except (TypeError, ValueError):
                    pass
        for k, cb in (("type", self.cb_type), ("method", self.cb_method), ("dir", self.cb_dir),
                      ("rmode", self.cb_rmode), ("pmode", self.cb_pmode)):
            v = s.value(S_PREFIX + k, None)
            if v is not None:
                try:
                    cb.setCurrentIndex(int(v))
                except (TypeError, ValueError):
                    pass
        for k, c in (("opt", self.chk_opt), ("bound", self.chk_bound), ("auto_ee", self.chk_auto_ee),
                     ("auto_grid", self.chk_auto_grid), ("main", self.chk_main), ("report", self.chk_report),
                     ("open", self.chk_open), ("tspm", self.chk_tspm), ("global", self.chk_global),
                     ("scan", self.grp_scan), ("rel", self.grp_rel)):
            v = s.value(S_PREFIX + k, None)
            if v is not None:
                c.setChecked(str(v).lower() in ("true", "1"))
        for m, c in self.chk_methods.items():
            v = s.value(S_PREFIX + "m_" + m, None)
            if v is not None:
                c.setChecked(str(v).lower() in ("true", "1"))
        for k, ed in (("ed_kh", self.ed_kh), ("ed_depth", self.ed_depth), ("out", self.ed_out)):
            v = s.value(S_PREFIX + k, None)
            if v:
                ed.setText(str(v))
