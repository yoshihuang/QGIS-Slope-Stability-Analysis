# -*- coding: utf-8 -*-
"""Processing algorithm: slope stability analysis (limit equilibrium)."""

import os

from qgis.core import (Qgis, QgsProcessing, QgsProcessingAlgorithm, QgsProcessingException,
                       QgsProcessingParameterVectorLayer, QgsProcessingParameterNumber,
                       QgsProcessingParameterEnum, QgsProcessingParameterBoolean,
                       QgsProcessingParameterFeatureSink, QgsProcessingParameterFileDestination,
                       QgsProcessingOutputNumber, QgsWkbTypes)
from qgis.PyQt.QtGui import QIcon

from ..core.i18n import tr
from ..core.methods import METHODS, method_label
from ..core.search import SearchSettings, Cancelled
from ..core.analysis import run_analysis, Scenario
from ..core.report import build_report

SURFACE_TYPES = ["grid", "entry_exit", "path", "user", "point"]


def _src(name):
    st = getattr(Qgis, "ProcessingSourceType", None)
    if st is not None and hasattr(st, name):
        return getattr(st, name)
    return getattr(QgsProcessing, "Type" + name)


def _num_double():
    nt = getattr(Qgis, "ProcessingNumberParameterType", None)
    if nt is not None:
        return nt.Double
    return QgsProcessingParameterNumber.Double


def _num_int():
    nt = getattr(Qgis, "ProcessingNumberParameterType", None)
    if nt is not None:
        return nt.Integer
    return QgsProcessingParameterNumber.Integer


def _wkb(name):
    wt = getattr(Qgis, "WkbType", None)
    if wt is not None and hasattr(wt, name):
        return getattr(wt, name)
    return getattr(QgsWkbTypes, name)


class SlopeStabilityAlgorithm(QgsProcessingAlgorithm):
    GEOLOGY = "GEOLOGY"
    WATER = "WATER"
    LOADS = "LOADS"
    REINF = "REINFORCEMENT"
    PILES = "PILES"
    SURFACES = "SURFACES"
    KH = "KH"
    KV = "KV"
    CRACK = "CRACK_DEPTH"
    CRACK_FILL = "CRACK_FILL"
    TYPE = "SURFACE_TYPE"
    METHOD = "METHOD"
    DENSITY = "DENSITY"
    SLICES = "SLICES"
    OPTIMIZE = "OPTIMIZE"
    DEFAULT_V = "DEFAULT_V"
    CRACK_X = "CRACK_X"
    OUT_SURF = "OUTPUT_SURFACES"
    OUT_CRIT = "OUTPUT_CRITICAL"
    OUT_SLICES = "OUTPUT_SLICES"
    OUT_REPORT = "OUTPUT_REPORT"
    MIN_FS = "MIN_FS"

    def createInstance(self):
        return SlopeStabilityAlgorithm()

    def name(self):
        return "slopestabilitylem"

    def displayName(self):
        return tr("Slope stability analysis (limit equilibrium)")

    def group(self):
        return tr("Slope stability")

    def groupId(self):
        return "slopestability"

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(os.path.dirname(__file__)), "icon.png"))

    def shortHelpString(self):
        return tr("Two-dimensional slope stability analysis by limit equilibrium methods of slices "
                  "(Ordinary, Bishop, Janbu, Corps of Engineers, Lowe-Karafiath, Spencer, "
                  "Morgenstern-Price) with circular and non-circular critical surface search.\n\n"
                  "The geology layer must be drawn in section coordinates (x = distance, y = elevation) "
                  "and contain the fields unit_weight (kN/m3), cohesion (kPa) and phi (degrees).")

    def initAlgorithm(self, config=None):
        P = QgsProcessingParameterVectorLayer
        self.addParameter(P(self.GEOLOGY, tr("Geology / materials (polygon)"), [_src("VectorPolygon")]))
        self.addParameter(P(self.WATER, tr("Piezometric line (line)"), [_src("VectorLine")], optional=True))
        self.addParameter(P(self.LOADS, tr("Distributed loads (line)"), [_src("VectorLine")], optional=True))
        self.addParameter(P(self.REINF, tr("Reinforcement (line)"), [_src("VectorLine")], optional=True))
        self.addParameter(P(self.PILES, tr("Piles / retaining walls"), [_src("VectorAnyGeometry")], optional=True))
        self.addParameter(P(self.SURFACES, tr("User-defined slip surfaces (line)"), [_src("VectorLine")], optional=True))
        N = QgsProcessingParameterNumber
        self.addParameter(N(self.KH, tr("Horizontal seismic coefficient kh"), _num_double(), 0.0, False, 0.0, 1.0))
        self.addParameter(N(self.KV, tr("Vertical seismic coefficient kv"), _num_double(), 0.0, False, -1.0, 1.0))
        self.addParameter(N(self.CRACK, tr("Crack depth"), _num_double(), 0.0, False, 0.0))
        self.addParameter(N(self.CRACK_FILL, tr("Water in crack") + " (0-1)", _num_double(), 0.0, False, 0.0, 1.0))
        self.addParameter(N(self.DEFAULT_V, tr("Pile / wall shear force V"), _num_double(), 0.0, False, 0.0))
        self.addParameter(QgsProcessingParameterEnum(self.TYPE, tr("Surface type"), [
            tr("Circular - grid search"), tr("Circular - entry and exit"),
            tr("Non-circular - path search"), tr("User-defined surfaces only"),
            tr("From a crack point picked on the map")], defaultValue=0))
        p = QgsProcessingParameterNumber(self.CRACK_X, tr("Crack point x (for the crack point search)"), _num_double(), None, True)
        self.addParameter(p)
        self.addParameter(QgsProcessingParameterEnum(self.METHOD, tr("Search method"),
                                                     [method_label(m) for m in METHODS], defaultValue=1))
        self.addParameter(N(self.DENSITY, tr("Search density") + " (%)", _num_int(), 100, False, 25, 400))
        self.addParameter(N(self.SLICES, tr("Number of slices"), _num_int(), 30, False, 5, 200))
        self.addParameter(QgsProcessingParameterBoolean(self.OPTIMIZE, tr("Optimise the critical surface"), True))
        self.addParameter(QgsProcessingParameterFeatureSink(self.OUT_SURF, tr("Slip surfaces"), _src("VectorLine")))
        self.addParameter(QgsProcessingParameterFeatureSink(self.OUT_CRIT, tr("Critical slip surface"), _src("VectorLine")))
        self.addParameter(QgsProcessingParameterFeatureSink(self.OUT_SLICES, tr("Slices (critical surface)"), _src("VectorPolygon"), optional=True))
        self.addParameter(QgsProcessingParameterFileDestination(self.OUT_REPORT, tr("HTML report"), "HTML (*.html)", optional=True))
        self.addOutput(QgsProcessingOutputNumber(self.MIN_FS, tr("Minimum factor of safety")))

    def processAlgorithm(self, parameters, context, feedback):
        from ..qgis_io import (build_model, InputError, result_rows, make_fields, make_feature,
                               OUTPUTS)
        lay = lambda k: self.parameterAsVectorLayer(parameters, k, context)
        layers = {"geology": lay(self.GEOLOGY), "water": lay(self.WATER), "loads": lay(self.LOADS),
                  "reinforcement": lay(self.REINF), "piles": lay(self.PILES), "surfaces": lay(self.SURFACES)}
        opts = {
            "kh": self.parameterAsDouble(parameters, self.KH, context),
            "kv": self.parameterAsDouble(parameters, self.KV, context),
            "crack_depth": self.parameterAsDouble(parameters, self.CRACK, context),
            "crack_fill": self.parameterAsDouble(parameters, self.CRACK_FILL, context),
            "default_V": self.parameterAsDouble(parameters, self.DEFAULT_V, context),
        }
        try:
            model, user, msgs = build_model(layers, opts)
        except InputError as e:
            raise QgsProcessingException(str(e))
        for lvl, m in msgs:
            feedback.pushInfo(m)
        stype = SURFACE_TYPES[self.parameterAsEnum(parameters, self.TYPE, context)]
        method = METHODS[self.parameterAsEnum(parameters, self.METHOD, context)]
        settings = SearchSettings(surface_type=stype, method=method,
                                  density=self.parameterAsInt(parameters, self.DENSITY, context) / 100.0,
                                  n_slices=self.parameterAsInt(parameters, self.SLICES, context),
                                  optimize=self.parameterAsBool(parameters, self.OPTIMIZE, context),
                                  user_surfaces=user)
        crack = None
        if stype == "point":
            if parameters.get(self.CRACK_X) in (None, ""):
                raise QgsProcessingException(tr("Crack point x is required for the crack point search."))
            crack = (self.parameterAsDouble(parameters, self.CRACK_X, context), 0.0)
            settings.crack_point = crack
        if stype == "user" and not user:
            raise QgsProcessingException(tr("Select a layer with user-defined slip surfaces."))

        def progress(f, msg=""):
            feedback.setProgress(100.0 * f)
        try:
            res = run_analysis(model, settings, [Scenario("Main", label=tr("Main analysis"))], list(METHODS),
                               None, progress, feedback.isCanceled, feedback.pushInfo)
        except Cancelled:
            return {}
        for lvl, m in res.preflight:
            (feedback.reportError if lvl == "error" else feedback.pushWarning if hasattr(feedback, "pushWarning") else feedback.pushInfo)(m)
        if any(lvl == "error" for lvl, _ in res.preflight):
            raise QgsProcessingException(tr("Input errors"))

        rows = result_rows(res)
        crs = layers["geology"].crs()
        results = {}
        sinks = {"surfaces": (self.OUT_SURF, "LineString"), "critical": (self.OUT_CRIT, "LineString"),
                 "slices": (self.OUT_SLICES, "Polygon")}
        schemas = {k: (kind, schema) for k, _n, kind, _u, schema in OUTPUTS}
        for key, (pname, wkb) in sinks.items():
            kind, schema = schemas[key]
            fields = make_fields(schema)
            sink, dest = self.parameterAsSink(parameters, pname, context, fields, _wkb(wkb), crs)
            if sink is None:
                continue
            for coords, attrs in rows.get(key, []):
                sink.addFeature(make_feature(fields, schema, kind, coords, attrs))
            results[pname] = dest
        sr = res.scenarios[0]
        results[self.MIN_FS] = sr.critical.fs if sr.critical is not None else None
        if sr.critical is not None:
            feedback.pushInfo(u"%s = %.3f" % (tr("Minimum factor of safety"), sr.critical.fs))
        rpt = self.parameterAsFileOutput(parameters, self.OUT_REPORT, context)
        if rpt:
            with open(rpt, "w", encoding="utf-8") as f:
                f.write(build_report(res, {"plugin": "Slope Stability Analysis", "project": layers["geology"].name()}))
            results[self.OUT_REPORT] = rpt
        return results
