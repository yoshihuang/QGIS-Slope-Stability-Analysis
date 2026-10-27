# -*- coding: utf-8 -*-
"""Verification of the computation core (run from the repository root: python tests/test_benchmarks.py).

Needs numpy and shapely only (no QGIS)."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from shapely.geometry import Polygon, LineString
from SlopeStabilityPlugin.core.model import SlopeModel
from SlopeStabilityPlugin.core.materials import Material
from SlopeStabilityPlugin.core.slices import build_slices, PolylineSurface
from SlopeStabilityPlugin.core.methods import solve
from SlopeStabilityPlugin.core.search import Searcher, SearchSettings


def acads_1a():
    m = SlopeModel()
    i = m.add_material(Material("soil", 20.0, 3.0, 19.6))
    m.add_region(Polygon([(20, 25), (30, 25), (50, 35), (70, 35), (70, 20), (20, 20)]), i)
    return m.prepare()


def test_acads_1a_circular():
    m = acads_1a()
    r = Searcher(m, SearchSettings(surface_type="grid", method="bishop")).run()[0]
    assert abs(r.fs - 0.987) < 0.01          # Slide: Bishop 0.987, referee ~1.00
    S = build_slices(m, r.surface, 30)
    assert abs(solve(S, m, "spencer").fs - 0.986) < 0.01
    assert abs(solve(S, m, "morgenstern_price").fs - 0.986) < 0.01


def _infinite(**kw):
    m = SlopeModel()
    i = m.add_material(Material("s", 18.0, 0.0, 30.0))
    m.add_region(Polygon([(0, 120), (20, 120), (220, 20), (260, 20), (260, 0), (0, 0)]), i)
    for k, v in kw.items():
        setattr(m, k, v)
    m.prepare()
    beta = math.atan(0.5)
    x = np.linspace(21, 219, 40)
    y = 120 - (x - 20) * 0.5 - 2.0 / math.cos(beta)
    s = PolylineSurface(np.r_[20.5, x, 219.5], np.r_[119.75, y, 20.25])
    return m, build_slices(m, s, 80), beta


def test_infinite_slope_dry():
    m, S, beta = _infinite()
    exact = math.tan(math.radians(30)) / math.tan(beta)
    for meth in ("spencer", "morgenstern_price", "corps1", "lowe_karafiath"):
        assert abs(solve(S, m, meth).fs - exact) / exact < 0.025, meth


def test_infinite_slope_seismic():
    m, S, beta = _infinite(kh=0.1)
    kh = 0.1
    exact = (math.cos(beta) - kh * math.sin(beta)) * math.tan(math.radians(30)) / (math.sin(beta) + kh * math.cos(beta))
    assert abs(solve(S, m, "spencer").fs - exact) / exact < 0.02


def test_taylor_vertical_cut_undrained():
    m = SlopeModel()
    i = m.add_material(Material("clay", 20.0, 50.0, 0.0, model="undrained"))
    m.add_region(Polygon([(0, 20), (30, 20), (30, 10), (80, 10), (80, 0), (0, 0)]), i)
    m.prepare()
    r = Searcher(m, SearchSettings(surface_type="grid", method="bishop", density=1.5)).run()[0]
    assert abs(r.fs - 3.83 * 50 / (20 * 10)) < 0.02   # Taylor Ns = 3.83


def test_crack_point_search():
    """Every surface of a crack point search starts at the picked point."""
    m = acads_1a()                       # slides towards -x, crest at x = 50..70
    st = SearchSettings(surface_type="point", method="bishop", crack_point=(55.0, 35.0), density=0.6)
    s = Searcher(m, st)
    recs = s.run()
    assert recs
    for r in recs:
        assert abs(m.to_world_x(r.x_in) - 55.0) < 0.01 * m.width + 1e-6
    # constrained search can never beat the global minimum (0.987)
    assert recs[0].fs > 0.97


if __name__ == "__main__":
    for name, f in sorted(globals().items()):
        if name.startswith("test_"):
            f(); print("OK", name)
