# -*- coding: utf-8 -*-
"""Limit equilibrium methods of slices.

One consistent formulation is used for every method (model coordinates, the
sliding mass moves towards +x, see ``slices.py`` for sign conventions).

Per slice i (unit width), unknowns N_i (base normal force) and the interslice
forces (E, X) acting on the right side.  The force exerted by slice i-1 on
slice i is (E_{i-1}, X_{i-1}); slice i receives (-E_i, -X_i) from slice i+1.

    x:  E_{i-1} - E_i + N sin(a) - S cos(a) + Fx = 0
    y:  X_{i-1} - X_i + N cos(a) + S sin(a) + Fy = 0
        S = [c l + (N - u l) tan(phi)] / F

Fx, Fy are the external forces (weight, seismic, surcharge, ponded water and
point forces).  "Passive" point forces (reinforcement, piles) are factored by
1/F like the soil strength; "active" ones act at full value.

Whole-mass moment equilibrium about a point O and horizontal force
equilibrium then give F_m and F_f.  Methods differ only in the interslice
force assumption X_i = k_i E_i:

    Ordinary (Fellenius)       interslice forces ignored, moment
    Bishop simplified          X = 0, moment
    Janbu simplified           X = 0, force (+ Janbu f0 correction)
    Corps of Engineers #1      X/E = tan(theta), theta = entry-exit chord
    Corps of Engineers #2      X/E = tan(theta), theta = ground slope
    Lowe-Karafiath             X/E = tan(avg(ground slope, base slope)), force
    Spencer                    X = lambda E, moment + force
    Morgenstern-Price          X = lambda f(x) E, f = half-sine, moment + force
"""

import math
import numpy as np

from .i18n import tr

METHODS = ("ordinary", "bishop", "janbu", "corps1", "corps2",
           "lowe_karafiath", "spencer", "morgenstern_price")

CIRCULAR_ONLY = ("ordinary", "bishop")


def method_label(key):
    labels = {
        "ordinary": tr("Ordinary (Fellenius)"),
        "bishop": tr("Bishop simplified"),
        "janbu": tr("Janbu simplified"),
        "corps1": tr("Corps of Engineers #1"),
        "corps2": tr("Corps of Engineers #2"),
        "lowe_karafiath": tr("Lowe-Karafiath"),
        "spencer": tr("Spencer"),
        "morgenstern_price": tr("Morgenstern-Price"),
    }
    return labels.get(key, key)


def method_short(key):
    return {
        "ordinary": "OMS", "bishop": "Bishop", "janbu": "Janbu",
        "corps1": "Corps1", "corps2": "Corps2", "lowe_karafiath": "LK",
        "spencer": "Spencer", "morgenstern_price": "MP",
    }.get(key, key)


class LEMResult(object):
    def __init__(self, method, fs=None, converged=False, message=""):
        self.method = method
        self.fs = fs
        self.converged = converged
        self.message = message
        self.lam = None
        self.theta = None
        self.fs_uncorrected = None
        self.N = None
        self.E = None
        self.X = None
        self.malpha_warning = False

    @property
    def ok(self):
        return self.fs is not None and self.converged and math.isfinite(self.fs)

    def __repr__(self):
        return "LEMResult(%s, FS=%s)" % (self.method, self.fs)


# --------------------------------------------------------------------------
#  external forces
# --------------------------------------------------------------------------
class _Ext(object):
    """External forces resolved per slice + their moments about O."""

    def __init__(self, S, model, O):
        n = S.n
        xo, yo = O
        kh, kv = model.kh, model.kv
        W = S.W
        # active forces per slice
        self.Fx_a = kh * W + S.pond_fx
        self.Fy_a = -(W + S.Q + S.pond) - kv * W
        self.Fx_p = np.zeros(n)
        self.Fy_p = np.zeros(n)
        # moments (counter-clockwise positive) about O
        Ma = np.sum((S.xm - xo) * (-(W + kv * W)))          # weight
        Ma += np.sum((S.xm - xo) * (-(S.Q + S.pond)))       # surcharge + water
        Ma += np.sum(-(S.yc - yo) * (kh * W))               # seismic at centroid
        Ma += np.sum(-(S.ground_m - yo) * S.pond_fx)        # ponded water (horizontal)
        Mp = 0.0
        for (k, px, py, fx, fy, passive, _lbl) in S.point_forces:
            m = (px - xo) * fy - (py - yo) * fx
            if passive:
                self.Fx_p[k] += fx
                self.Fy_p[k] += fy
                Mp += m
            else:
                self.Fx_a[k] += fx
                self.Fy_a[k] += fy
                Ma += m
        self.Ma, self.Mp = Ma, Mp
        sa, ca = np.sin(S.alpha), np.cos(S.alpha)
        dx, dy = S.xm - xo, S.yb - yo
        self.mn = dx * ca - dy * sa          # moment arm of unit normal
        self.ms = dx * sa + dy * ca          # moment arm of unit (resisting) shear
        self.sa, self.ca = sa, ca


def moment_point(S, model):
    """Moment centre: circle centre or a point above the sliding mass."""
    surf = S.surface
    if surf.kind == "circular":
        return (surf.cx, surf.cy)
    xo = 0.5 * (S.x_in + S.x_out)
    yo = float(np.max(S.ground_m)) + 0.5 * (S.x_out - S.x_in)
    return (xo, yo)


# --------------------------------------------------------------------------
#  solvers
# --------------------------------------------------------------------------
_MAX_IT = 100
_TOL = 1e-5
_MALPHA_MIN = 0.2


def _strength(S, N):
    return S.c * S.l + (N - S.u * S.l) * S.tphi


def _N_vertical(S, ext, F):
    """Base normal force from vertical equilibrium with X = 0 (Bishop/Janbu)."""
    ma = ext.ca + ext.sa * S.tphi / F
    warn = bool(np.any(ma < _MALPHA_MIN))
    ma = np.maximum(ma, _MALPHA_MIN)
    N = (-ext.Fy_a - ext.Fy_p / F - (S.c * S.l - S.u * S.l * S.tphi) * ext.sa / F) / ma
    return N, warn


def _fs_moment(S, ext, N):
    den = ext.Ma + float(np.sum(N * ext.mn))
    num = -(float(np.sum(_strength(S, N) * ext.ms)) + ext.Mp)
    if den <= 1e-9:
        return None
    return num / den


def _fs_force(S, ext, N):
    den = float(np.sum(ext.Fx_a)) + float(np.sum(N * ext.sa))
    num = float(np.sum(_strength(S, N) * ext.ca)) - float(np.sum(ext.Fx_p))
    if den <= 1e-9:
        return None
    return num / den


def solve_ordinary(S, ext):
    r = LEMResult("ordinary")
    Fx = ext.Fx_a + ext.Fx_p
    Fy = ext.Fy_a + ext.Fy_p
    N = -(Fx * ext.sa + Fy * ext.ca)
    Neff = np.maximum(N - S.u * S.l, 0.0)
    N = Neff + S.u * S.l
    den = ext.Ma + ext.Mp + float(np.sum(N * ext.mn))
    num = -float(np.sum((S.c * S.l + Neff * S.tphi) * ext.ms))
    if den <= 1e-9:
        r.message = tr("No driving moment.")
        return r
    r.fs = num / den
    r.converged = True
    r.N = N
    return r


def solve_bishop(S, ext, F0=1.5):
    r = LEMResult("bishop")
    F = F0
    for _ in range(_MAX_IT):
        N, warn = _N_vertical(S, ext, F)
        Fn = _fs_moment(S, ext, N)
        if Fn is None or not math.isfinite(Fn) or Fn <= 0:
            r.message = tr("No driving moment.")
            return r
        if abs(Fn - F) < _TOL:
            F = Fn
            r.converged = True
            break
        F = Fn
    r.fs = F
    r.N = _N_vertical(S, ext, F)[0]
    r.malpha_warning = warn
    if not r.converged:
        r.message = tr("Did not converge.")
    return r


def _janbu_f0(S):
    L = math.hypot(S.x_out - S.x_in, S.y_out - S.y_in)
    if L <= 0:
        return 1.0
    # perpendicular distance from chord to slip surface
    nx, ny = (S.y_out - S.y_in) / L, -(S.x_out - S.x_in) / L
    d = np.abs((S.xm - S.x_in) * nx + (S.yb - S.y_in) * ny)
    dl = float(np.max(d)) / L
    has_c = bool(np.any(S.c > 0))
    has_phi = bool(np.any(S.tphi > 0))
    if has_c and has_phi:
        b1 = 0.50
    elif has_c:
        b1 = 0.69
    else:
        b1 = 0.31
    return 1.0 + b1 * (dl - 1.4 * dl * dl)


def solve_janbu(S, ext, F0=1.5, correct=True):
    r = LEMResult("janbu")
    F = F0
    for _ in range(_MAX_IT):
        N, warn = _N_vertical(S, ext, F)
        Fn = _fs_force(S, ext, N)
        if Fn is None or not math.isfinite(Fn) or Fn <= 0:
            r.message = tr("No driving force.")
            return r
        if abs(Fn - F) < _TOL:
            F = Fn
            r.converged = True
            break
        F = Fn
    r.fs_uncorrected = F
    r.fs = F * (_janbu_f0(S) if correct else 1.0)
    r.N = _N_vertical(S, ext, F)[0]
    r.malpha_warning = warn
    if not r.converged:
        r.message = tr("Did not converge.")
    return r


# ---------------------------------------------------------------- GLE march
def _prep_lists(S, ext):
    return dict(
        sa=ext.sa.tolist(), ca=ext.ca.tolist(),
        cl=(S.c * S.l).tolist(), ult=(S.u * S.l * S.tphi).tolist(),
        tp=S.tphi.tolist(),
        fxa=ext.Fx_a.tolist(), fya=ext.Fy_a.tolist(),
        fxp=ext.Fx_p.tolist(), fyp=ext.Fy_p.tolist(),
    )


def _march(L, F, k):
    """Return (N list, E list, X list).  k = interslice factor per boundary
    (len n+1; k[0] unused)."""
    sa, ca, cl, ult, tp = L["sa"], L["ca"], L["cl"], L["ult"], L["tp"]
    fxa, fya, fxp, fyp = L["fxa"], L["fya"], L["fxp"], L["fyp"]
    n = len(sa)
    Ep = 0.0
    Xp = 0.0
    Ns, Es, Xs = [0.0] * n, [0.0] * n, [0.0] * n
    invF = 1.0 / F
    for i in range(n):
        a = (cl[i] - ult[i]) * invF
        bb = tp[i] * invF
        Fx = fxa[i] + fxp[i] * invF
        Fy = fya[i] + fyp[i] * invF
        ki = k[i + 1]
        s, c = sa[i], ca[i]
        den = c + bb * s - ki * (s - bb * c)
        if abs(den) < 1e-12:
            den = 1e-12 if den >= 0 else -1e-12
        N = (ki * (Ep - a * c + Fx) - Xp - a * s - Fy) / den
        E = Ep + N * (s - bb * c) - a * c + Fx
        X = ki * E
        Ns[i], Es[i], Xs[i] = N, E, X
        Ep, Xp = E, X
    return Ns, Es, Xs


def _fs_from_march(S, ext, N, which):
    N = np.asarray(N)
    return _fs_moment(S, ext, N) if which == "m" else _fs_force(S, ext, N)


def _residual_fun(S, ext, L, k, which):
    Wsum = max(float(np.sum(S.W + S.Q + S.pond)), 1e-9)
    Lsc = max(S.x_out - S.x_in, 1e-9)
    cl, ult = S.c * S.l, S.u * S.l * S.tphi

    def r(F):
        N, E, X = _march(L, F, k)
        if which == "f":
            return E[-1] / Wsum
        Na = np.asarray(N)
        strength = cl + Na * S.tphi - ult
        M = ext.Ma + float(np.dot(Na, ext.mn)) + (float(np.dot(strength, ext.ms)) + ext.Mp) / F
        return M / (Wsum * Lsc)
    return r


def _solve_F(S, ext, L, k, which, F0):
    """Root of the force (E_n = 0) or moment residual for fixed interslice
    factors k.  When several roots exist the one closest to F0 is returned
    (fixed-point iteration is avoided: it can be unstable at the physical
    root)."""
    r = _residual_fun(S, ext, L, k, which)
    F0 = min(max(F0, 0.05), 50.0)
    # secant from the starting estimate
    x0, x1 = F0, F0 * 1.05
    try:
        g0, g1 = r(x0), r(x1)
        for _ in range(40):
            if abs(g1 - g0) < 1e-15:
                break
            x2 = x1 - g1 * (x1 - x0) / (g1 - g0)
            if not (0.02 < x2 < 100.0) or abs(x2 - F0) > 0.6 * F0 + 0.5:
                break
            x0, g0, x1, g1 = x1, g1, x2, r(x2)
            if abs(g1) < 1e-10 or abs(x1 - x0) < 1e-7:
                return x1
    except (ZeroDivisionError, OverflowError, ValueError):
        pass
    # bracket scan, root nearest to F0
    grid = np.geomspace(0.03, 30.0, 90)
    vals = []
    for F in grid:
        try:
            vals.append(r(F))
        except (ZeroDivisionError, OverflowError, ValueError):
            vals.append(None)
    roots = []
    for i in range(len(grid) - 1):
        a, b = vals[i], vals[i + 1]
        if a is None or b is None or (a > 0) == (b > 0):
            continue
        lo, hi, glo = grid[i], grid[i + 1], a
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            gm = r(mid)
            if (gm > 0) == (glo > 0):
                lo, glo = mid, gm
            else:
                hi = mid
        roots.append(0.5 * (lo + hi))
    if not roots:
        return None
    return min(roots, key=lambda x: abs(math.log(x / F0)))


def _secant_root(g, x0, x1, tol=1e-5, it=40):
    g0, g1 = g(x0), g(x1)
    if g0 is None or g1 is None:
        return None
    for _ in range(it):
        if abs(g1 - g0) < 1e-14:
            break
        x2 = x1 - g1 * (x1 - x0) / (g1 - g0)
        x2 = max(min(x2, 3.0), -3.0)
        g2 = g(x2)
        if g2 is None:
            return None
        x0, g0, x1, g1 = x1, g1, x2, g2
        if abs(g1) < tol:
            return x1
    return x1 if abs(g1) < 10 * tol else None


def _bracket_root(g, lo=-2.0, hi=2.0, steps=40, tol=1e-5):
    xs = np.linspace(lo, hi, steps + 1)
    prev_x, prev_g = None, None
    for x in xs:
        gx = g(x)
        if gx is None:
            prev_x, prev_g = None, None
            continue
        if prev_g is not None and (gx == 0 or (gx > 0) != (prev_g > 0)):
            a, ga, b = prev_x, prev_g, x
            for _ in range(60):
                m = 0.5 * (a + b)
                gm = g(m)
                if gm is None:
                    return None
                if abs(gm) < tol:
                    return m
                if (gm > 0) == (ga > 0):
                    a, ga = m, gm
                else:
                    b = m
            return 0.5 * (a + b)
        prev_x, prev_g = x, gx
    return None


def _interslice_f(S, shape):
    xb = S.boundaries
    if shape == "constant":
        f = np.ones(len(xb))
    else:  # half-sine
        t = (xb - xb[0]) / max(xb[-1] - xb[0], 1e-12)
        f = np.sin(math.pi * t)
    f[0] = 0.0
    return f


def _newton_gle(S, ext, L, f, F0, lam0=0.0):
    """Solve E_n(F, lam) = 0 and M(F, lam) = 0 simultaneously (2-D Newton)."""
    Wsum = max(float(np.sum(S.W + S.Q + S.pond)), 1e-9)
    Lsc = max(S.x_out - S.x_in, 1e-9)
    mn, ms = ext.mn, ext.ms
    cl, ult = S.c * S.l, S.u * S.l * S.tphi
    fl = f.tolist()

    def res(F, lam):
        if F <= 0.01:
            return None
        k = [lam * v for v in fl]
        N, E, X = _march(L, F, k)
        Na = np.asarray(N)
        strength = cl + Na * S.tphi - ult
        M = ext.Ma + float(np.dot(Na, mn)) + (float(np.dot(strength, ms)) + ext.Mp) / F
        return np.array([E[-1] / Wsum, M / (Wsum * Lsc)])

    x = np.array([F0, lam0], dtype=float)
    r0 = res(*x)
    if r0 is None:
        return None
    n0 = float(np.hypot(*r0))
    for _ in range(40):
        if n0 < 1e-9:
            return x
        J = np.empty((2, 2))
        for j, h in enumerate((1e-6 * max(1.0, x[0]), 1e-6)):
            xp = x.copy()
            xp[j] += h
            rp = res(*xp)
            if rp is None:
                return None
            J[:, j] = (rp - r0) / h
        try:
            dx = -np.linalg.solve(J, r0)
        except np.linalg.LinAlgError:
            return None
        # limit the step
        sc = 1.0
        if abs(dx[0]) > 0.5 * x[0]:
            sc = min(sc, 0.5 * x[0] / abs(dx[0]))
        if abs(dx[1]) > 0.5:
            sc = min(sc, 0.5 / abs(dx[1]))
        dx *= sc
        t = 1.0
        while t > 1e-3:
            xn = x + t * dx
            rn = res(*xn)
            if rn is not None and float(np.hypot(*rn)) < n0:
                break
            t *= 0.5
        else:
            return None
        x, r0, n0 = xn, rn, float(np.hypot(*rn))
        if np.max(np.abs(t * dx)) < 1e-9:
            break
    return x if n0 < 1e-7 else None


def _admissible(S, res):
    """Reject numerically spurious full-equilibrium solutions."""
    if res.N is None:
        return True
    N = np.asarray(res.N)
    Neff = N - S.u * S.l
    W = S.W + S.Q + S.pond
    scale = float(np.sum(W)) + 1e-9
    # large negative effective normal forces or interslice tension indicate a
    # non-physical root of the equations
    if float(np.sum(np.clip(-Neff, 0, None))) > 0.10 * scale:
        return False
    if res.E is not None:
        E = np.asarray(res.E)
        if float(np.min(E)) < -0.25 * scale:
            return False
    if res.lam is not None and abs(res.lam) > 3.0:
        return False
    return True


def solve_gle(S, ext, method, F0=1.5):
    """Spencer (f = 1) or Morgenstern-Price (f = half-sine)."""
    r = LEMResult(method)
    L = _prep_lists(S, ext)
    shape = "constant" if method == "spencer" else "half-sine"
    f = _interslice_f(S, shape)
    for lam0 in (-0.3, 0.0, -0.8):
        sol = _newton_gle(S, ext, L, f, F0, lam0)
        if sol is not None and 0.01 < sol[0] < 100:
            F, lam = float(sol[0]), float(sol[1])
            N, E, X = _march(L, F, (lam * f).tolist())
            r.fs, r.lam = F, lam
            r.theta = math.degrees(math.atan(-lam))
            r.N, r.E, r.X = np.asarray(N), np.asarray(E), np.asarray(X)
            r.converged = True
            if _admissible(S, r):
                return r
            r = LEMResult(method)
    return _solve_gle_slow(S, ext, method, F0, L, f)


def _solve_gle_slow(S, ext, method, F0, L, f):
    r = LEMResult(method)
    cache = {}

    def fm_ff(lam):
        key = round(lam, 10)
        if key in cache:
            return cache[key]
        k = (lam * f).tolist()
        Fm = _solve_F(S, ext, L, k, "m", F0)
        Ff = _solve_F(S, ext, L, k, "f", F0)
        cache[key] = (Fm, Ff)
        return Fm, Ff

    def g(lam):
        Fm, Ff = fm_ff(lam)
        if Fm is None or Ff is None:
            return None
        return Fm - Ff

    lam = _secant_root(g, 0.0, -0.15)
    if lam is None:
        lam = _bracket_root(g)
    if lam is None:
        r.message = tr("Did not converge.")
        return r
    Fm, Ff = fm_ff(lam)
    if Fm is None or Ff is None or abs(Fm - Ff) > 1e-3 * max(1.0, abs(Fm)):
        r.message = tr("Did not converge.")
        return r
    F = 0.5 * (Fm + Ff)
    k = (lam * f).tolist()
    N, E, X = _march(L, F, k)
    r.fs = F
    r.lam = lam
    r.theta = math.degrees(math.atan(-lam))
    r.N, r.E, r.X = np.asarray(N), np.asarray(E), np.asarray(X)
    r.converged = True
    if not _admissible(S, r):
        bad = LEMResult(method)
        bad.message = tr("Solution rejected: non-physical interslice or base forces.")
        return bad
    return r


def _ground_slope(model, x):
    h = max(model.width * 1e-4, 1e-4)
    g = model.ground(np.array([x - h, x + h]))
    return math.atan2(g[0] - g[1], 2 * h)       # descending positive


def solve_force_method(S, ext, model, method, F0=1.5):
    """Corps of Engineers #1/#2 and Lowe-Karafiath (force equilibrium only)."""
    r = LEMResult(method)
    xb = S.boundaries
    n = S.n
    if method == "corps1":
        th = math.atan2(S.y_in - S.y_out, S.x_out - S.x_in)
        theta = np.full(n + 1, th)
    elif method == "corps2":
        theta = np.array([_ground_slope(model, x) for x in xb])
    else:  # Lowe-Karafiath
        gs = np.array([_ground_slope(model, x) for x in xb])
        a = S.alpha
        base = np.empty(n + 1)
        base[0], base[-1] = a[0], a[-1]
        base[1:-1] = 0.5 * (a[:-1] + a[1:])
        theta = 0.5 * (gs + base)
    k = (-np.tan(theta)).tolist()
    k[0] = 0.0
    L = _prep_lists(S, ext)
    F = _solve_F(S, ext, L, k, "f", F0)
    if F is None:
        r.message = tr("Did not converge.")
        return r
    N, E, X = _march(L, F, k)
    r.fs, r.converged = F, True
    r.N, r.E, r.X = np.asarray(N), np.asarray(E), np.asarray(X)
    if not _admissible(S, r):
        bad = LEMResult(method)
        bad.message = tr("Solution rejected: non-physical interslice or base forces.")
        return bad
    return r


# --------------------------------------------------------------------------
#  public entry point
# --------------------------------------------------------------------------
def _dispatch(S, model, method, ext, F0):
    if method == "ordinary":
        return solve_ordinary(S, ext)
    if method == "bishop":
        return solve_bishop(S, ext, F0)
    if method == "janbu":
        return solve_janbu(S, ext, F0)
    if method in ("spencer", "morgenstern_price"):
        return solve_gle(S, ext, method, F0)
    if method in ("corps1", "corps2", "lowe_karafiath"):
        return solve_force_method(S, ext, model, method, F0)
    raise ValueError("Unknown method %s" % method)


def solve(S, model, method, F0=None):
    """Compute the factor of safety of slice set ``S`` with ``method``."""
    if not S.valid:
        return LEMResult(method, message=S.reason)
    if method in CIRCULAR_ONLY and S.surface.kind != "circular":
        return LEMResult(method, message=tr("Method requires a circular slip surface."))
    O = moment_point(S, model)
    ext = _Ext(S, model, O)
    if F0 is None:
        F0 = 1.5
        if method not in ("ordinary", "bishop"):
            rb = solve_bishop(S, ext) if S.surface.kind == "circular" else solve_janbu(S, ext, correct=False)
            if rb.ok and 0.05 < rb.fs < 50:
                F0 = rb.fs
    if not np.any(S.nonlinear):
        return _dispatch(S, model, method, ext, F0)
    # non-linear strength: update tangent c, phi from the base normal stress
    res = None
    for _ in range(25):
        res = _dispatch(S, model, method, ext, F0)
        if not res.ok or res.N is None:
            return res
        sig = np.maximum((res.N - S.u * S.l) / S.l, 0.0)
        changed = 0.0
        for i in np.nonzero(S.nonlinear)[0]:
            c, t = model.materials[S.mat[i]].linear_parameters(sig[i])
            changed = max(changed, abs(c - S.c[i]) / max(1.0, abs(c)), abs(t - S.tphi[i]))
            S.c[i], S.tphi[i] = c, t
        F0 = res.fs
        if changed < 1e-4:
            break
    return res


def slice_table(S, model, res):
    """Per-slice output values for reports and the slice layer."""
    rows = []
    if res is None or res.N is None:
        return rows
    F = res.fs
    N = np.asarray(res.N)
    Neff = N - S.u * S.l
    Sm = (S.c * S.l + Neff * S.tphi) / F
    for i in range(S.n):
        rows.append({
            "slice": i + 1,
            "x": float(S.xm[i]),
            "width": float(S.b[i]),
            "height": float(S.ground_m[i] - S.yb[i]),
            "alpha": math.degrees(float(S.alpha[i])),
            "base_len": float(S.l[i]),
            "weight": float(S.W[i]),
            "surcharge": float(S.Q[i]),
            "material": model.materials[int(S.mat[i])].name,
            "c": float(S.c[i]),
            "phi": math.degrees(math.atan(float(S.tphi[i]))),
            "u": float(S.u[i]),
            "N_eff": float(Neff[i]),
            "sigma_eff": float(Neff[i] / S.l[i]),
            "shear_mob": float(Sm[i]),
            "E": float(res.E[i]) if res.E is not None else None,
            "X": float(res.X[i]) if res.X is not None else None,
        })
    return rows
