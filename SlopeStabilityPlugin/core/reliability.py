# -*- coding: utf-8 -*-
"""Probabilistic analysis on a given slip surface.

* Monte Carlo simulation (normal distributions truncated to physical limits)
* Taylor Series Probabilistic Method (TSPM, Duncan 2000) with a lognormal FS

Random variables: for every material with a non-zero standard deviation
``c_std`` (cohesion / Su), ``phi_std`` and ``gamma_std``.
"""

import math
import numpy as np

from .slices import build_slices
from .methods import solve
from .i18n import tr


def _phi_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _lognormal_beta(mean, std):
    if mean <= 0 or std <= 0:
        return None
    v = std / mean
    s2 = math.log(1.0 + v * v)
    return math.log(mean / math.sqrt(1.0 + v * v)) / math.sqrt(s2)


def random_variables(model):
    out = []
    for i, m in enumerate(model.materials):
        if m.c_std > 0:
            out.append((i, "cohesion", m.c_std))
        if m.phi_std > 0 and m.model == "mc":
            out.append((i, "phi", m.phi_std))
        if m.gamma_std > 0:
            out.append((i, "gamma", m.gamma_std))
    return out


def _apply(model, values):
    mats = []
    for i, m in enumerate(model.materials):
        ch = {}
        for (j, prop), v in values.items():
            if j != i:
                continue
            if prop == "cohesion":
                ch["cohesion"] = max(v, 0.0)
            elif prop == "phi":
                ch["phi"] = min(max(v, 0.0), 89.0)
            elif prop == "gamma":
                g = max(v, 0.1)
                ch["gamma_sat"] = m.gamma_sat + (g - m.gamma)
                ch["gamma"] = g
        mats.append(m.copy(**ch) if ch else m)
    return model.with_materials(mats)


def run_reliability(model, surface, method, n_slices=30, mc_samples=1000,
                    tspm=True, seed=1, progress=None, cancelled=None):
    progress = progress or (lambda f, m="": None)
    cancelled = cancelled or (lambda: False)
    rvs = random_variables(model)
    if not rvs:
        return {"error": tr("No material has a standard deviation (c_std, phi_std, gamma_std); reliability analysis skipped.")}

    def fs_of(values):
        mdl = _apply(model, values) if values else model
        S = build_slices(mdl, surface, n_slices)
        if not S.valid:
            return None
        r = solve(S, mdl, method)
        return r.fs if r.ok else None

    mean_vals = {}
    for (i, prop, sd) in rvs:
        m = model.materials[i]
        mean_vals[(i, prop)] = {"cohesion": m.cohesion, "phi": m.phi, "gamma": m.gamma}[prop]
    F_det = fs_of({})
    out = {"variables": [(model.materials[i].name, prop, mean_vals[(i, prop)], sd) for (i, prop, sd) in rvs],
           "fs_deterministic": F_det, "method": method}

    # ---- TSPM -----------------------------------------------------------------
    if tspm and F_det is not None:
        var = 0.0
        contrib = []
        for (i, prop, sd) in rvs:
            if cancelled():
                return None
            mu = mean_vals[(i, prop)]
            fp = fs_of({(i, prop): mu + sd})
            fm = fs_of({(i, prop): mu - sd})
            if fp is None or fm is None:
                continue
            d = (fp - fm) / 2.0
            var += d * d
            contrib.append((model.materials[i].name, prop, fp, fm, d * d))
        sf = math.sqrt(var)
        b = _lognormal_beta(F_det, sf)
        tot = sum(c[4] for c in contrib) or 1.0
        out["tspm"] = {
            "sigma": sf,
            "cov": sf / F_det if F_det else None,
            "beta_ln": b,
            "pf": _phi_cdf(-b) if b is not None else None,
            "contributions": [(n, p, fp, fm, v / tot) for (n, p, fp, fm, v) in contrib],
        }

    # ---- Monte Carlo ----------------------------------------------------------------
    if mc_samples and mc_samples > 0:
        rng = np.random.default_rng(seed)
        fs = []
        failed = 0
        n = int(mc_samples)
        for k in range(n):
            if k % 25 == 0:
                if cancelled():
                    return None
                progress(k / float(n), tr("Monte Carlo simulation"))
            vals = {}
            for (i, prop, sd) in rvs:
                vals[(i, prop)] = float(rng.normal(mean_vals[(i, prop)], sd))
            f = fs_of(vals)
            if f is None:
                failed += 1
                continue
            fs.append(f)
        if fs:
            a = np.asarray(fs)
            mu, sd = float(a.mean()), float(a.std(ddof=1)) if len(a) > 1 else 0.0
            hist, edges = np.histogram(a, bins=30)
            out["mc"] = {
                "n": len(a), "n_failed_solutions": failed,
                "mean": mu, "std": sd, "min": float(a.min()), "max": float(a.max()),
                "pf": float(np.mean(a < 1.0)),
                "beta_normal": (mu - 1.0) / sd if sd > 0 else None,
                "beta_ln": _lognormal_beta(mu, sd),
                "hist": hist.tolist(), "edges": edges.tolist(),
            }
    return out
