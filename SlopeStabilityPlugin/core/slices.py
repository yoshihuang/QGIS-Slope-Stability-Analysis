# -*- coding: utf-8 -*-
"""Slip surfaces and discretisation of the sliding mass into vertical slices.

Sign conventions (model coordinates, mass moves towards +x):

* ``alpha``  base inclination, positive when the base descends towards +x
  (driving part of the surface), negative under the toe.
* Normal force N acts on the slice base along n = (sin a, cos a).
* Mobilised shear S acts on the base against the motion: s = (-cos a, sin a).
"""

import math
import numpy as np
from shapely.geometry import LineString, Point, Polygon

from .i18n import tr


# --------------------------------------------------------------------------
#  slip surfaces
# --------------------------------------------------------------------------
class CircularSurface(object):
    kind = "circular"

    def __init__(self, cx, cy, r):
        self.cx, self.cy, self.r = float(cx), float(cy), float(r)

    def x_range(self):
        return self.cx - self.r, self.cx + self.r

    def y(self, xs):
        d = np.clip(self.r ** 2 - (np.asarray(xs) - self.cx) ** 2, 0.0, None)
        return self.cy - np.sqrt(d)

    def params(self):
        return (self.cx, self.cy, self.r)

    def vertices_x(self):
        return []

    def __repr__(self):
        return "Circle(cx=%.3f, cy=%.3f, r=%.3f)" % (self.cx, self.cy, self.r)


class PolylineSurface(object):
    kind = "noncircular"

    def __init__(self, xs, ys):
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        order = np.argsort(xs, kind="stable")
        self.xs, self.ys = xs[order], ys[order]

    def x_range(self):
        return float(self.xs[0]), float(self.xs[-1])

    def y(self, xs):
        return np.interp(xs, self.xs, self.ys)

    def params(self):
        return tuple(np.column_stack([self.xs, self.ys]).ravel())

    def vertices_x(self):
        return self.xs[1:-1].tolist()

    def __repr__(self):
        return "Polyline(%d vertices)" % len(self.xs)


# --------------------------------------------------------------------------
#  slice set
# --------------------------------------------------------------------------
class SliceSet(object):
    """Discretised sliding mass (all arrays are per slice, model coordinates)."""

    def __init__(self):
        self.valid = False
        self.reason = ""

    @property
    def n(self):
        return len(self.b)


def _refine_crossing(fun, a, b, it=40):
    fa = fun(a)
    for _ in range(it):
        m = 0.5 * (a + b)
        fm = fun(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def _below_ground_interval(model, surface, x_lo, x_hi, samples=240):
    """Find the (entry, exit) interval where the surface is below the ground."""
    xs = np.linspace(x_lo, x_hi, samples)
    d = model.ground(xs) - surface.y(xs)
    tol = model.relief * 1e-6
    below = d > tol
    if not below.any():
        return None, tr("Slip surface does not cut the slope.")
    # contiguous runs of "below"
    idx = np.nonzero(below)[0]
    splits = np.nonzero(np.diff(idx) > 1)[0]
    starts = np.concatenate([[idx[0]], idx[splits + 1]])
    ends = np.concatenate([idx[splits], [idx[-1]]])
    best, best_area = None, -1.0
    for s, e in zip(starts, ends):
        area = float(np.sum(d[s:e + 1]))
        if area > best_area:
            best, best_area = (s, e), area
    s, e = best
    if s == 0 or e == len(xs) - 1:
        return None, tr("Slip surface leaves the model through its side boundary.")
    f = lambda x: float(model.ground(np.array([x]))[0] - surface.y(np.array([x]))[0])
    x_in = _refine_crossing(f, xs[s - 1], xs[s])
    x_out = _refine_crossing(f, xs[e], xs[e + 1])
    return (x_in, x_out), ""


def build_slices(model, surface, n_slices=30, check_bottom=True):
    """Discretise the sliding mass above ``surface`` into vertical slices."""
    S = SliceSet()
    S.surface = surface
    gx0, gx1 = float(model.gx[0]), float(model.gx[-1])
    sx0, sx1 = surface.x_range()
    x_lo = max(gx0, sx0)
    x_hi = min(gx1, sx1)
    if x_hi - x_lo <= model.width * 1e-4:
        S.reason = tr("Slip surface does not cut the slope.")
        return S
    if surface.kind == "circular":
        # the side boundary of the model must not be crossed below ground
        interval, reason = _below_ground_interval(model, surface, x_lo, x_hi)
        if interval is not None:
            if (sx0 >= gx0 and abs(interval[0] - sx0) < 1e-9) or (sx1 <= gx1 and abs(interval[1] - sx1) < 1e-9):
                interval, reason = None, tr("Slip surface is steeper than vertical at its end.")
    else:
        interval, reason = _below_ground_interval(model, surface, x_lo, x_hi)
    if interval is None:
        S.reason = reason
        return S
    x_in, x_out = interval
    if x_out - x_in < model.width * 1e-3:
        S.reason = tr("Sliding mass is too small.")
        return S

    # ---- tension crack ----------------------------------------------------------
    crack_h = 0.0
    if model.crack_depth > 0:
        dc = model.crack_depth
        depth = lambda x: float(model.ground(np.array([x]))[0] - surface.y(np.array([x]))[0])
        xs = np.linspace(x_in, x_out, 200)
        dd = model.ground(xs) - surface.y(xs)
        k = np.nonzero(dd >= dc)[0]
        if len(k) == 0:
            S.reason = tr("Sliding mass is shallower than the tension crack.")
            return S
        k = k[0]
        if k > 0:
            x_in = _refine_crossing(lambda x: depth(x) - dc, xs[k - 1], xs[k])
            crack_h = dc
        else:
            crack_h = max(dd[0], 0.0)

    # ---- slice boundaries -------------------------------------------------------
    nb = max(int(n_slices), 3)
    xb = np.linspace(x_in, x_out, nb + 1)
    extra = [x for x in model.gvx if x_in < x < x_out]
    extra += [x for x in surface.vertices_x() if x_in < x < x_out]
    if model._px is not None:
        extra += [x for x in model._px if x_in < x < x_out]
    if extra:
        xb = np.unique(np.concatenate([xb, np.asarray(extra, dtype=float)]))
        min_w = (x_out - x_in) * 2e-3
        keep = [xb[0]]
        for x in xb[1:-1]:
            if x - keep[-1] >= min_w and x_out - x >= min_w:
                keep.append(x)
        keep.append(xb[-1])
        xb = np.asarray(keep)
    xl, xr = xb[:-1], xb[1:]
    b = xr - xl
    xm = 0.5 * (xl + xr)

    # ---- base geometry -------------------------------------------------------------
    if surface.kind == "circular":
        cx, cy, r = surface.cx, surface.cy, surface.r
        def ang(x):
            return np.arcsin(np.clip((x - cx) / r, -1.0, 1.0))
        tl, tr_, tm = ang(xl), ang(xr), ang(xm)
        alpha = -tm                       # descending towards +x -> positive
        lb = r * (tr_ - tl)
        yb = cy - r * np.cos(tm)
        ybl = cy - r * np.cos(tl)
        ybr = cy - r * np.cos(tr_)
    else:
        ybl = surface.y(xl)
        ybr = surface.y(xr)
        alpha = np.arctan2(ybl - ybr, b)
        lb = np.hypot(b, ybr - ybl)
        yb = 0.5 * (ybl + ybr)

    if np.any(np.abs(alpha) > math.radians(89.0)):
        S.reason = tr("Slip surface is steeper than vertical at its end.")
        return S

    col = model.columns(xm, yb)
    mat = col["mat"]
    if np.any(mat < 0):
        S.reason = tr("Part of the slip surface lies outside the geology polygons.")
        return S
    ground_m = col["ground"]
    if np.any(~np.isfinite(ground_m)):
        S.reason = tr("Part of the slip surface lies outside the geology polygons.")
        return S
    if check_bottom:
        bot = col["bottom"]
        tol = model.relief * 1e-3
        if np.any(yb < bot - tol):
            S.reason = tr("Slip surface passes below the bottom of the model.")
            return S

    height = ground_m - yb
    depth_lim = model.min_depth if model.min_depth > 0 else 0.0
    if depth_lim > 0 and np.max(height) < depth_lim:
        S.reason = tr("Sliding mass is shallower than the minimum depth.")
        return S

    gl = model.ground(xl)
    gr = model.ground(xr)
    gw = model.gamma_w

    # ---- weights -------------------------------------------------------------------
    W = col["wpa"] * b

    # ---- pore pressure at base -------------------------------------------------------
    hw = col["hw"]
    u = np.zeros_like(xm)
    ru = np.array([model.materials[m].ru for m in mat])
    use_pz = np.array([model.materials[m].use_piezo for m in mat])
    if model.has_water:
        u = np.where(use_pz & (hw > yb), gw * (hw - yb), 0.0)
    u = np.where(ru > 0, np.maximum(u, ru * col["wpa"]), u)

    # ---- surface water (ponded) as pressure on the ground surface -------------------
    pond = np.zeros_like(xm)
    pond_fx = np.zeros_like(xm)
    if model.has_water:
        hm = hw - 0.5 * (gl + gr)
        hm = np.where(hm > 0, hm, 0.0)
        pond = gw * hm * b                    # vertical (downward) component
        pond_fx = gw * hm * (gr - gl)         # horizontal component (into slope if gr<gl)

    # ---- surcharge ---------------------------------------------------------------------
    Q = np.zeros_like(xm)
    for (a0, a1, q) in model._loads:
        ov = np.clip(np.minimum(xr, a1) - np.maximum(xl, a0), 0.0, None)
        Q += q * ov

    # ---- assemble ----------------------------------------------------------------------
    S.valid = True
    S.xl, S.xr, S.xm, S.b = xl, xr, xm, b
    S.alpha = alpha
    S.l = lb
    S.yb, S.ybl, S.ybr = yb, ybl, ybr
    S.ground_m, S.gl, S.gr = ground_m, gl, gr
    S.yc = 0.5 * (ground_m + yb)
    S.W = W
    S.Q = Q
    S.pond = pond
    S.pond_fx = pond_fx
    S.u = u
    S.mat = mat
    S.x_in, S.x_out = float(x_in), float(x_out)
    S.y_in = float(surface.y(np.array([x_in]))[0])
    S.y_out = float(surface.y(np.array([x_out]))[0])
    S.crack_h = crack_h
    S.max_depth = float(np.max(height))
    S.area = float(np.sum(height * b))

    # strength parameters (tangent parameters for non-linear materials are
    # updated by the solver)
    S.c = np.zeros_like(xm)
    S.tphi = np.zeros_like(xm)
    S.nonlinear = np.zeros(len(xm), dtype=bool)
    for i, m in enumerate(mat):
        M = model.materials[m]
        sv = max(col["wpa"][i] - u[i], 0.0)
        c, t = M.linear_parameters(sv)
        S.c[i], S.tphi[i] = c, t
        S.nonlinear[i] = M.is_nonlinear

    # ---- point forces: (slice, px, py, fx, fy, passive, label) ------------------------
    pf = []
    if crack_h > 0 and model.crack_fill > 0:
        hwc = crack_h * min(model.crack_fill, 1.0)
        Pw = 0.5 * gw * hwc ** 2
        pf.append((0, x_in, S.y_in + hwc / 3.0, Pw, 0.0, False, "crack water"))
    if model._reinf or model._piles:
        path = LineString(list(zip(np.concatenate([[x_in], xm, [x_out]]),
                                   np.concatenate([[S.y_in], yb, [S.y_out]]))))
        mass = _mass_polygon(model, S)
        for geom, rf in model._reinf:
            for P in _intersections(path, geom):
                force, d = _reinforcement_force(rf, geom, P, mass)
                if force <= 0:
                    continue
                k = _slice_of(xb, P.x)
                pf.append((k, P.x, P.y, force * d[0], force * d[1], rf.passive, "reinforcement"))
        for geom, pl in model._piles:
            if pl.shear <= 0:
                continue
            hit = path.intersection(geom)
            if hit.is_empty:
                continue
            P = hit.centroid if hit.geom_type != "Point" else hit
            k = _slice_of(xb, P.x)
            pf.append((k, P.x, float(surface.y(np.array([P.x]))[0]), -pl.shear, 0.0, pl.passive, "pile"))
    if model.global_support > 0:
        k = len(xm) - 1
        a = alpha[k]
        pf.append((k, xm[k], yb[k], -model.global_support * math.cos(a),
                   model.global_support * math.sin(a), model.global_support_passive, "support"))
    S.point_forces = pf
    S.boundaries = xb
    return S


def _slice_of(xb, x):
    k = int(np.searchsorted(xb, x) - 1)
    return min(max(k, 0), len(xb) - 2)


def _intersections(path, geom):
    hit = path.intersection(geom)
    if hit.is_empty:
        return []
    if hit.geom_type == "Point":
        return [hit]
    if hit.geom_type == "MultiPoint":
        return list(hit.geoms)
    pts = []
    for g in getattr(hit, "geoms", [hit]):
        if g.geom_type == "Point":
            pts.append(g)
        else:
            pts.append(g.centroid)
    return pts


def _mass_polygon(model, S):
    xs = np.concatenate([[S.x_in], S.xm, [S.x_out]])
    yb = np.concatenate([[S.y_in], S.yb, [S.y_out]])
    gx = np.concatenate([[S.x_in], S.xm, [S.x_out]])
    gy = model.ground(gx)
    if S.crack_h > 0:
        gy[0] = S.y_in + S.crack_h
    ring = list(zip(xs, yb)) + list(zip(gx[::-1], gy[::-1]))
    try:
        poly = Polygon(ring)
        if not poly.is_valid:
            poly = poly.buffer(0)
        return poly
    except Exception:
        return None


def _reinforcement_force(rf, geom, P, mass):
    """Force magnitude and unit direction (pointing out of the sliding mass)."""
    coords = list(geom.coords) if geom.geom_type == "LineString" else []
    if len(coords) < 2:
        return 0.0, (0.0, 0.0)
    a = np.array(coords[0], dtype=float)
    b = np.array(coords[-1], dtype=float)
    p = np.array([P.x, P.y])
    # which end lies outside the sliding mass (the anchored end)?
    def outside(pt):
        if mass is None:
            return True
        return not mass.contains(Point(pt[0], pt[1]))
    # use the parts of the line on each side of P
    proj = geom.project(P)
    L = geom.length
    side_a_len, side_b_len = proj, L - proj
    pa = np.array(geom.interpolate(max(proj - min(0.5 * side_a_len, 1.0), 0.0)).coords[0]) if side_a_len > 0 else a
    pb = np.array(geom.interpolate(min(proj + min(0.5 * side_b_len, 1.0), L)).coords[0]) if side_b_len > 0 else b
    if outside(pb) and not outside(pa):
        tgt, embed = pb, side_b_len
    elif outside(pa) and not outside(pb):
        tgt, embed = pa, side_a_len
    else:
        # ambiguous: assume the end with the lower x (into the slope) is anchored
        tgt, embed = (pa, side_a_len) if pa[0] < pb[0] else (pb, side_b_len)
    d = tgt - p
    nrm = float(np.hypot(d[0], d[1]))
    if nrm == 0:
        return 0.0, (0.0, 0.0)
    d = d / nrm
    force = rf.capacity
    if rf.pullout > 0:
        force = min(force, rf.pullout * embed)
    return force, (float(d[0]), float(d[1]))
