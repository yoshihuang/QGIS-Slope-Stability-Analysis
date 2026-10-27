# -*- coding: utf-8 -*-
"""High level analysis driver (QGIS independent)."""

import time

from .search import Searcher, Cancelled
from .slices import build_slices
from .methods import solve, slice_table, CIRCULAR_ONLY, METHODS
from .reliability import run_reliability
from .i18n import tr


class Scenario(object):
    def __init__(self, name, kh=None, water_depth="keep", use_piezo=True, label="", density=None):
        self.name = name
        self.density = density              # optional search density override
        self.kh = kh
        self.water_depth = water_depth      # "keep" | None (dry) | float
        self.use_piezo = use_piezo
        self.label = label or name


class ScenarioResult(object):
    def __init__(self, scenario):
        self.scenario = scenario
        self.records = []
        self.critical = None
        self.method_results = {}
        self.slices = None
        self.slice_rows = []
        self.center_grid = []
        self.n_eval = 0
        self.n_valid = 0
        self.model = None
        self.elapsed = 0.0
        self.message = ""


class AnalysisResult(object):
    def __init__(self):
        self.scenarios = []
        self.reliability = None
        self.preflight = []
        self.model = None
        self.settings = None
        self.methods = []
        self.elapsed = 0.0


def scenario_model(model, sc):
    piezo = "keep"
    wd = "keep"
    if sc.water_depth != "keep":
        wd = sc.water_depth
        piezo = None
    elif not sc.use_piezo:
        piezo = None
    return model.with_conditions(kh=sc.kh, water_depth=wd, piezo=piezo)


def evaluate_all_methods(model, surface, methods, n_slices):
    S = build_slices(model, surface, n_slices)
    out = {}
    if not S.valid:
        return S, out
    for m in methods:
        if m in CIRCULAR_ONLY and surface.kind != "circular":
            continue
        out[m] = solve(S, model, m)
    return S, out


def run_analysis(model, settings, scenarios, methods, reliability=None,
                 progress=None, cancelled=None, log=None):
    """Run every scenario and return an :class:`AnalysisResult`.

    ``settings``    SearchSettings
    ``scenarios``   list[Scenario] (first one is the main analysis)
    ``methods``     methods computed on each critical surface
    ``reliability`` dict(mc_samples=int, tspm=bool, seed=int) or None
    """
    progress = progress or (lambda f, m="": None)
    cancelled = cancelled or (lambda: False)
    log = log or (lambda m: None)
    t0 = time.time()
    res = AnalysisResult()
    res.settings = settings
    res.methods = [m for m in METHODS if m in methods]
    if not model._prepared:
        model.prepare()
    res.model = model
    res.preflight = model.preflight()
    if any(level == "error" for level, _ in res.preflight):
        return res

    n_sc = len(scenarios)
    rel_share = 0.15 if reliability else 0.0
    for k, sc in enumerate(scenarios):
        if cancelled():
            raise Cancelled()
        t1 = time.time()
        mdl = scenario_model(model, sc)
        base = (1.0 - rel_share) * k / n_sc
        span = (1.0 - rel_share) / n_sc

        def prog(f, msg="", _b=base, _s=span, _n=sc.label):
            progress(_b + _s * min(max(f, 0.0), 1.0), u"[{0}] {1}".format(_n, msg))

        st = settings
        if sc.density is not None:
            import copy
            st = copy.copy(settings)
            st.density = sc.density
        searcher = Searcher(mdl, st, prog, cancelled)
        searcher.run()
        sr = ScenarioResult(sc)
        sr.model = mdl
        sr.records = searcher.records
        for r in sr.records:
            r.scenario = sc.name
        sr.center_grid = searcher.center_grid
        sr.n_eval, sr.n_valid = searcher.n_eval, searcher.n_valid
        if sr.records:
            sr.critical = sr.records[0]
            prog(0.97, tr("Evaluating all methods on the critical surface"))
            need = set(res.methods) | {sr.critical.method or searcher.method}
            S, mres = evaluate_all_methods(mdl, sr.critical.surface,
                                           [m for m in METHODS if m in need], settings.n_slices)
            sr.slices = S
            sr.method_results = mres
            pm = sr.critical.method or searcher.method
            primary = mres.get(pm)
            if primary is None or not primary.ok:
                for m in res.methods:
                    if m in mres and mres[m].ok:
                        primary = mres[m]
                        break
            sr.primary_method = primary.method if (primary is not None and primary.ok) else pm
            if primary is not None and primary.ok:
                sr.slice_rows = slice_table(S, mdl, primary)
        else:
            sr.message = tr("No valid slip surface was found. Check the geometry, the search ranges and the material parameters.")
        sr.elapsed = time.time() - t1
        log(u"{0}: {1} surfaces evaluated, {2} valid, min FS = {3}".format(
            sc.label, sr.n_eval, sr.n_valid,
            "-" if sr.critical is None else "%.3f" % sr.critical.fs))
        res.scenarios.append(sr)

    if reliability and res.scenarios and res.scenarios[0].critical is not None:
        sr = res.scenarios[0]
        base = 1.0 - rel_share

        def rprog(f, msg=""):
            progress(base + rel_share * f, msg)
        res.reliability = run_reliability(
            sr.model, sr.critical.surface, sr.primary_method, settings.n_slices,
            reliability.get("mc_samples", 1000), reliability.get("tspm", True),
            reliability.get("seed", 1), rprog, cancelled)
    res.elapsed = time.time() - t0
    progress(1.0, tr("Finished"))
    return res
