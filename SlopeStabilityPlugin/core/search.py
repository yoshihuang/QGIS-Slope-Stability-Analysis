# -*- coding: utf-8 -*-
"""Critical slip surface search algorithms.

* ``grid``        circular, grid of centres x tangent (radius) levels
* ``entry_exit``  circular, entry range x exit range x arc depth
* ``path``        non-circular, random convex paths between entry and exit
* ``user``        user-defined surfaces only

Circular searches can be followed by Nelder-Mead optimisation of
(cx, cy, R); non-circular searches by a vertex-moving pattern search
(surface altering optimisation).
"""

import math
import numpy as np

from .slices import CircularSurface, PolylineSurface, build_slices
from .methods import solve, CIRCULAR_ONLY
from .i18n import tr


class Cancelled(Exception):
    pass


class SearchSettings(object):
    def __init__(self, **kw):
        self.surface_type = kw.get("surface_type", "grid")
        self.method = kw.get("method", "bishop")
        self.density = float(kw.get("density", 1.0))
        self.n_slices = int(kw.get("n_slices", 30))
        self.optimize = bool(kw.get("optimize", True))
        self.entry_range = kw.get("entry_range")      # world (x0, x1) or None
        self.exit_range = kw.get("exit_range")
        self.grid = kw.get("grid")                    # world (x0, x1, y0, y1) or None
        self.n_keep = int(kw.get("n_keep", 300))
        self.path_vertices = int(kw.get("path_vertices", 6))
        self.seed = kw.get("seed", 12345)
        self.user_surfaces = kw.get("user_surfaces", [])  # list of world coord lists
        self.include_user = bool(kw.get("include_user", True))
        # reject surfaces whose entry/exit lies on the lateral model boundary
        self.exclude_boundary = bool(kw.get("exclude_boundary", True))
        # "point" search: world (x, y) of a crack / scarp picked on the ground
        self.crack_point = kw.get("crack_point")


class SurfaceRecord(object):
    def __init__(self, fs, surface, S):
        self.fs = fs
        self.surface = surface
        self.x_in, self.x_out = S.x_in, S.x_out
        self.y_in, self.y_out = S.y_in, S.y_out
        self.scenario = ""
        self.source = ""
        self.method = ""

    def polyline(self, model, npts=72):
        """Slip surface between entry and exit, world coordinates."""
        s = self.surface
        if s.kind == "circular":
            xs = np.linspace(self.x_in, self.x_out, npts)
            ys = s.y(xs)
        else:
            inner = [(x, y) for x, y in zip(s.xs, s.ys) if self.x_in < x < self.x_out]
            xs = np.array([self.x_in] + [p[0] for p in inner] + [self.x_out])
            ys = np.array([self.y_in] + [p[1] for p in inner] + [self.y_out])
        wx, wy = model.to_world_xy(xs, ys)
        return list(zip(wx.tolist(), wy.tolist()))

    def center_world(self, model):
        if self.surface.kind != "circular":
            return None
        return (model.to_world_x(self.surface.cx), self.surface.cy)


class Searcher(object):
    def __init__(self, model, settings, progress=None, cancelled=None):
        self.model = model
        self.st = settings
        self.progress = progress or (lambda f, m="": None)
        self.cancelled = cancelled or (lambda: False)
        self.records = []
        self.n_eval = 0
        self.n_valid = 0
        self.center_grid = []          # (cx_world, cy, min FS)
        # numpy generator: reproducible random trial paths (not used for security)
        self.rng = np.random.default_rng(settings.seed)
        self._xp = None                # crack point (model x) for "point" searches
        method = settings.method
        if settings.surface_type == "path" and method in CIRCULAR_ONLY:
            method = "spencer"
        self.method = method

    # ------------------------------------------------------------- helpers
    def _check(self):
        if self.cancelled():
            raise Cancelled()

    def evaluate(self, surface, keep=True, source=""):
        self.n_eval += 1
        S = build_slices(self.model, surface, self.st.n_slices)
        if not S.valid:
            return None
        if self._xp is not None and source != "user":
            # every surface must start at the picked crack point
            w = self.model.width
            lo = self._xp - 0.01 * w
            hi = self._xp + (0.05 * w if self.model.crack_depth > 0 else 0.01 * w)
            if not (lo <= S.x_in <= hi):
                return None
        if self.st.exclude_boundary and source != "user":
            tol = 0.01 * self.model.width
            if S.x_in - self.model.gx[0] < tol or self.model.gx[-1] - S.x_out < tol:
                return None
        method = self.method
        if surface.kind != "circular" and method in CIRCULAR_ONLY:
            method = "spencer"
        r = solve(S, self.model, method)
        if not r.ok or not (0.01 < r.fs < 100.0):
            return None
        self.n_valid += 1
        if keep:
            rec = SurfaceRecord(r.fs, surface, S)
            rec.source = source
            rec.method = method
            self.records.append(rec)
            if len(self.records) > 4 * max(self.st.n_keep, 50):
                self._trim()
        return r.fs

    def _trim(self):
        self.records.sort(key=lambda r: r.fs)
        self.records = self.records[:max(self.st.n_keep, 50) * 2]

    def ranges(self):
        m = self.model
        gx, gy = m.gx, m.gy
        mid = 0.5 * (gy.max() + gy.min())
        k = np.nonzero(gy <= mid)[0]
        x_mid = gx[k[0]] if len(k) else 0.5 * (gx[0] + gx[-1])
        entry = (gx[0], x_mid)
        exit_ = (x_mid, gx[-1])
        if self.st.entry_range:
            a, b = [m.to_model_x(v) for v in self.st.entry_range]
            entry = (min(a, b), max(a, b))
        if self.st.exit_range:
            a, b = [m.to_model_x(v) for v in self.st.exit_range]
            exit_ = (min(a, b), max(a, b))
        return entry, exit_

    # ------------------------------------------------------------- run
    def run(self):
        t = self.st.surface_type
        if t == "grid":
            self.grid_search()
        elif t == "entry_exit":
            self.entry_exit_search()
        elif t == "path":
            self.path_search()
        elif t == "point":
            self.point_search()
        if self.st.include_user and self.st.user_surfaces:
            self.user_search()
        self._trim()
        self.records.sort(key=lambda r: r.fs)
        # remove near-duplicates
        uniq = []
        seen = set()
        for r in self.records:
            key = (round(r.x_in, 2), round(r.x_out, 2), round(r.fs, 4))
            if key in seen:
                continue
            seen.add(key)
            uniq.append(r)
        self.records = uniq[:self.st.n_keep]
        return self.records

    # --------------------------------------------------- circular: grid
    def grid_search(self):
        m = self.model
        d = self.st.density
        nx = max(5, int(round(18 * d)))
        ny = max(5, int(round(14 * d)))
        nr = max(5, int(round(14 * d)))
        h = m.relief
        gx0, gx1 = float(m.gx[0]), float(m.gx[-1])
        w = gx1 - gx0
        if self.st.grid:
            a, b, y0, y1 = self.st.grid
            a, b = m.to_model_x(a), m.to_model_x(b)
            x0, x1 = min(a, b), max(a, b)
            y0, y1 = min(y0, y1), max(y0, y1)
        else:
            entry, exit_ = self.ranges()
            x0 = entry[0] + 0.1 * (exit_[1] - entry[0])
            x1 = exit_[1] - 0.1 * (exit_[1] - entry[0])
            gmax = float(m.gy.max())
            y0 = gmax + 0.1 * h
            y1 = gmax + max(1.5 * h, 0.4 * w)
        xs = np.linspace(x0, x1, nx)
        ys = np.linspace(y0, y1, ny)
        ylo = m.ymin + 0.01 * h
        yhi = float(m.gy.max()) - 0.01 * h
        total = nx * ny
        done = 0
        for cx in xs:
            for cy in ys:
                self._check()
                best = None
                for yb in np.linspace(ylo, yhi, nr):
                    r = cy - yb
                    if r <= 0:
                        continue
                    fs = self.evaluate(CircularSurface(cx, cy, r), source="grid")
                    if fs is not None and (best is None or fs < best):
                        best = fs
                if best is not None:
                    self.center_grid.append((m.to_world_x(cx), cy, best))
                done += 1
                if done % max(1, total // 50) == 0:
                    self.progress(0.75 * done / total, tr("Grid search"))
        if self.st.optimize:
            self.optimize_circles()

    # ------------------------------------------- circular: entry & exit
    def _circle_through(self, p1, p2, ratio):
        L = math.hypot(p2[0] - p1[0], p2[1] - p1[1])
        if L <= 0:
            return None
        s = ratio * L
        R = (L * L / 4.0 + s * s) / (2.0 * s)
        tx, ty = (p2[0] - p1[0]) / L, (p2[1] - p1[1]) / L
        nx_, ny_ = -ty, tx
        if ny_ < 0:
            nx_, ny_ = -nx_, -ny_
        mx, my = 0.5 * (p1[0] + p2[0]), 0.5 * (p1[1] + p2[1])
        return CircularSurface(mx + nx_ * (R - s), my + ny_ * (R - s), R)

    def entry_exit_search(self):
        m = self.model
        d = self.st.density
        ne = max(4, int(round(14 * d)))
        nx = max(4, int(round(14 * d)))
        nr = max(4, int(round(10 * d)))
        entry, exit_ = self.ranges()
        xe = np.linspace(entry[0], entry[1], ne + 2)[1:-1]
        xx = np.linspace(exit_[0], exit_[1], nx + 2)[1:-1]
        ratios = np.linspace(0.04, 0.5, nr)
        total = ne * nx
        done = 0
        for a in xe:
            ya = float(m.ground(np.array([a]))[0])
            for b in xx:
                self._check()
                done += 1
                if b - a < m.width * 0.01:
                    continue
                yb = float(m.ground(np.array([b]))[0])
                for q in ratios:
                    c = self._circle_through((a, ya), (b, yb), q)
                    if c is not None:
                        self.evaluate(c, source="entry_exit")
                if done % max(1, total // 50) == 0:
                    self.progress(0.75 * done / total, tr("Entry-exit search"))
        if self.st.optimize:
            self.optimize_circles()

    # -------------------------------------------- circular: Nelder-Mead
    def _distinct_best(self, k, kind):
        recs = sorted([r for r in self.records if r.surface.kind == kind], key=lambda r: r.fs)
        out = []
        for r in recs:
            if all(abs(r.x_in - o.x_in) + abs(r.x_out - o.x_out) > 0.05 * self.model.width for o in out):
                out.append(r)
            if len(out) >= k:
                break
        return out

    def optimize_circles(self, starts=3):
        seeds = self._distinct_best(starts, "circular")
        w = self.model.width
        for i, rec in enumerate(seeds):
            self._check()
            self.progress(0.75 + 0.2 * i / max(len(seeds), 1), tr("Optimising slip surface"))
            s = rec.surface

            def f(p):
                if p[2] <= 0:
                    return 1e6
                fs = self.evaluate(CircularSurface(*p), source="optimized")
                return 1e6 if fs is None else fs
            nelder_mead(f, [s.cx, s.cy, s.r], [0.04 * w, 0.04 * w, 0.03 * w],
                        max_eval=220, check=self._check)

    # ------------------------------------------------ non-circular paths
    def _admissible(self, xs, ys, tol_deg=5.0):
        """Kinematic admissibility of a non-circular surface: (almost) concave
        upwards, passive exit not steeper than 45 deg, active part < 85 deg."""
        if np.any(np.diff(xs) <= 1e-9):
            return False
        a = np.arctan2(ys[:-1] - ys[1:], np.diff(xs))      # descending positive
        if np.min(a) < -math.radians(45.0) or np.max(a) > math.radians(85.0):
            return False
        return bool(np.all(np.diff(a) <= math.radians(tol_deg)))

    def _random_path(self, entry, exit_, nv):
        m = self.model
        a = entry[0] if entry[0] == entry[1] else self.rng.uniform(*entry)
        b = self.rng.uniform(*exit_)
        if b - a < 0.02 * m.width:
            return None
        ya = float(m.ground(np.array([a]))[0])
        yb = float(m.ground(np.array([b]))[0])
        L = math.hypot(b - a, yb - ya)
        ratio = self.rng.uniform(0.04, 0.45)
        c = self._circle_through((a, ya), (b, yb), ratio)
        xs = np.linspace(a, b, nv + 2)
        ys = c.y(xs)
        ys[1:-1] += np.array([self.rng.uniform(-0.25, 0.25) * ratio * L for _ in range(nv)])
        ys[0], ys[-1] = ya, yb
        # keep inside the model
        g = m.ground(xs)
        bot = m.bottom(xs)
        bot = np.where(np.isfinite(bot), bot, m.ymin)
        ys[1:-1] = np.minimum(ys[1:-1], g[1:-1] - 1e-3 * m.relief)
        ys[1:-1] = np.maximum(ys[1:-1], bot[1:-1] + 1e-3 * m.relief)
        # lower convex hull -> kinematically admissible
        xs, ys = _lower_hull(xs, ys)
        if len(xs) < 3:
            return None
        return PolylineSurface(xs, ys)

    def path_search(self):
        d = self.st.density
        n = max(100, int(round(1500 * d)))
        entry, exit_ = self.ranges()
        nv = max(3, self.st.path_vertices)
        for i in range(n):
            if i % 20 == 0:
                self._check()
                self.progress(0.6 * i / n, tr("Path search"))
            s = self._random_path(entry, exit_, nv)
            if s is not None and self._admissible(s.xs, s.ys):
                self.evaluate(s, source="path")
        if self.st.optimize:
            seeds = self._distinct_best(3, "noncircular")
            for j, rec in enumerate(seeds):
                self.progress(0.6 + 0.35 * j / max(len(seeds), 1), tr("Optimising slip surface"))
                self.optimize_path(rec.surface, rec.fs)

    def optimize_path(self, surface, fs, max_eval=900, fixed_entry=False):
        m = self.model
        xs, ys = surface.xs.copy(), surface.ys.copy()
        best = fs
        step = 0.08 * m.relief
        evals = 0
        refinements = 0

        def trial(nx_, ny_):
            nonlocal evals
            evals += 1
            if not self._admissible(nx_, ny_):
                return None
            return self.evaluate(PolylineSurface(nx_, ny_), keep=True, source="optimized")

        while evals < max_eval:
            self._check()
            improved = False
            n = len(xs)
            for j in range(n):
                for dxy in ((0, -1), (0, 1), (-1, 0), (1, 0)):
                    nx_, ny_ = xs.copy(), ys.copy()
                    if j == 0 and fixed_entry:
                        continue
                    if j in (0, n - 1):
                        if dxy[1] != 0:
                            continue
                        nx_[j] += dxy[0] * step
                        ny_[j] = float(m.ground(np.array([nx_[j]]))[0])
                    else:
                        nx_[j] += dxy[0] * step * 0.5
                        ny_[j] += dxy[1] * step
                    v = trial(nx_, ny_)
                    if v is not None and v < best - 1e-5:
                        best, xs, ys = v, nx_, ny_
                        improved = True
            if not improved:
                step *= 0.5
                if step < 0.004 * m.relief:
                    if refinements < 2 and len(xs) < 14:
                        xs, ys = _insert_midpoints(xs, ys)
                        refinements += 1
                        step = 0.03 * m.relief
                    else:
                        break
        return best

    # ------------------------------------------------- crack point
    def point_search(self):
        """Slip surfaces starting at a crack / scarp point picked on the ground
        and extending downslope: circles through the point, random
        non-circular paths from the point, then optimisation with the entry
        kept fixed."""
        m = self.model
        if not self.st.crack_point:
            raise ValueError(tr("No crack point has been picked."))
        xp = m.to_model_x(float(self.st.crack_point[0]))
        gx0, gx1 = float(m.gx[0]), float(m.gx[-1])
        if not (gx0 < xp < gx1):
            raise ValueError(tr("The crack point is outside the geology layer."))
        if gx1 - xp < 0.03 * m.width:
            raise ValueError(tr("The crack point is too close to the toe; there is no slope below it."))
        self._xp = xp
        yp = float(m.ground(np.array([xp]))[0])
        d = self.st.density
        nx = max(6, int(round(30 * d)))
        nr = max(5, int(round(14 * d)))
        exits = np.linspace(xp + 0.02 * m.width, gx1, nx + 1)
        ratios = np.linspace(0.03, 0.5, nr)
        total = len(exits)
        for k, b in enumerate(exits):
            self._check()
            yb = float(m.ground(np.array([b]))[0])
            for q in ratios:
                c = self._circle_through((xp, yp), (b, yb), q)
                if c is not None:
                    self.evaluate(c, source="crack_circle")
            self.progress(0.35 * (k + 1) / total, tr("Circular surfaces from the crack point"))
        # non-circular paths from the crack point
        n = max(100, int(round(800 * d)))
        nv = max(3, self.st.path_vertices)
        for i in range(n):
            if i % 20 == 0:
                self._check()
                self.progress(0.35 + 0.3 * i / n, tr("Non-circular surfaces from the crack point"))
            s = self._random_path((xp, xp), (xp + 0.02 * m.width, gx1), nv)
            if s is not None and self._admissible(s.xs, s.ys):
                self.evaluate(s, source="crack_path")
        if self.st.optimize:
            # circles: optimise (exit x, arc depth ratio) with the entry fixed
            seeds = [r for r in self._distinct_best(2, "circular")]
            for j, rec in enumerate(seeds):
                self._check()
                self.progress(0.65 + 0.15 * j / max(len(seeds), 1), tr("Optimising slip surface"))
                L = math.hypot(rec.x_out - xp, rec.y_out - yp)
                s0 = rec.surface
                sag = s0.r - math.hypot((xp + rec.x_out) / 2 - s0.cx, (yp + rec.y_out) / 2 - s0.cy)

                def f(p):
                    b, q = p
                    if not (xp + 0.01 * m.width < b < gx1) or not (0.01 < q < 0.6):
                        return 1e6
                    yb = float(m.ground(np.array([b]))[0])
                    c = self._circle_through((xp, yp), (b, yb), q)
                    fs = None if c is None else self.evaluate(c, source="optimized")
                    return 1e6 if fs is None else fs
                q0 = min(max(sag / L if L > 0 else 0.2, 0.03), 0.5)
                nelder_mead(f, [rec.x_out, q0], [0.03 * m.width, 0.05], max_eval=120, check=self._check)
            seeds = self._distinct_best(2, "noncircular")
            for j, rec in enumerate(seeds):
                self.progress(0.8 + 0.15 * j / max(len(seeds), 1), tr("Optimising slip surface"))
                self.optimize_path(rec.surface, rec.fs, max_eval=600, fixed_entry=True)

    # ------------------------------------------------- user surfaces
    def user_search(self):
        m = self.model
        for coords in self.st.user_surfaces:
            self._check()
            if len(coords) < 2:
                continue
            xs = np.array([m.to_model_x(c[0]) for c in coords])
            ys = np.array([c[1] for c in coords])
            surf = circle_from_points(xs, ys)
            if surf is None:
                # extend the end segments so that the surface reaches the ground
                surf = PolylineSurface(*_extend_ends(xs, ys, m))
            method = self.method
            if surf.kind != "circular" and method in CIRCULAR_ONLY:
                self.method = "spencer"
            self.evaluate(surf, source="user")
            self.method = method


# --------------------------------------------------------------------------
#  utilities
# --------------------------------------------------------------------------
def _lower_hull(xs, ys):
    pts = sorted(zip(xs.tolist(), ys.tolist()))
    hull = []
    for p in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) <= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    return np.array([p[0] for p in hull]), np.array([p[1] for p in hull])


def _insert_midpoints(xs, ys):
    nx_ = [xs[0]]
    ny_ = [ys[0]]
    for i in range(1, len(xs)):
        nx_.append(0.5 * (xs[i - 1] + xs[i]))
        ny_.append(0.5 * (ys[i - 1] + ys[i]))
        nx_.append(xs[i])
        ny_.append(ys[i])
    return np.array(nx_), np.array(ny_)


def _extend_ends(xs, ys, m):
    order = np.argsort(xs)
    xs, ys = xs[order], ys[order]
    ext = m.width
    for end in (0, -1):
        i, j = (0, 1) if end == 0 else (-1, -2)
        dx, dy = xs[i] - xs[j], ys[i] - ys[j]
        L = math.hypot(dx, dy)
        if L == 0:
            continue
        g = float(m.ground(np.array([xs[i]]))[0])
        if ys[i] < g - 1e-6:
            nx_ = xs[i] + dx / L * ext
            ny_ = ys[i] + dy / L * ext
            if end == 0:
                xs = np.concatenate([[nx_], xs])
                ys = np.concatenate([[ny_], ys])
            else:
                xs = np.concatenate([xs, [nx_]])
                ys = np.concatenate([ys, [ny_]])
    return xs, ys


def circle_from_points(xs, ys, tol=0.002):
    """Fit a circle to a polyline; return CircularSurface if it is an arc."""
    if len(xs) < 5:
        return None
    A = np.column_stack([xs, ys, np.ones_like(xs)])
    bvec = xs ** 2 + ys ** 2
    try:
        sol, *_ = np.linalg.lstsq(A, bvec, rcond=None)
    except Exception:
        return None
    cx, cy = sol[0] / 2.0, sol[1] / 2.0
    r = math.sqrt(max(sol[2] + cx * cx + cy * cy, 0.0))
    if r <= 0:
        return None
    res = np.abs(np.hypot(xs - cx, ys - cy) - r)
    if np.max(res) > tol * r or cy < np.max(ys):
        return None
    return CircularSurface(cx, cy, r)


def nelder_mead(f, x0, steps, max_eval=200, tol=1e-5, check=None):
    n = len(x0)
    pts = [np.array(x0, dtype=float)]
    for i in range(n):
        p = np.array(x0, dtype=float)
        p[i] += steps[i]
        pts.append(p)
    vals = [f(p) for p in pts]
    evals = n + 1
    while evals < max_eval:
        if check:
            check()
        order = np.argsort(vals)
        pts = [pts[i] for i in order]
        vals = [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol and max(np.max(np.abs(p - pts[0])) for p in pts) < 1e-3 * max(steps):
            break
        c = np.mean(pts[:-1], axis=0)
        xr = c + (c - pts[-1])
        fr = f(xr)
        evals += 1
        if fr < vals[0]:
            xe = c + 2.0 * (c - pts[-1])
            fe = f(xe)
            evals += 1
            if fe < fr:
                pts[-1], vals[-1] = xe, fe
            else:
                pts[-1], vals[-1] = xr, fr
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = c + 0.5 * (pts[-1] - c)
            fc = f(xc)
            evals += 1
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                for i in range(1, len(pts)):
                    pts[i] = pts[0] + 0.5 * (pts[i] - pts[0])
                    vals[i] = f(pts[i])
                    evals += 1
    i = int(np.argmin(vals))
    return pts[i], vals[i]
