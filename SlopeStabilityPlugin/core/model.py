# -*- coding: utf-8 -*-
"""2-D slope model: geometry, materials, water, loads and supports.

All geometry is given in section coordinates: x = horizontal distance,
y = elevation (any consistent length unit, normally metres).

Internally the model is normalised so that the sliding mass always moves in
the +x direction (the slope "faces right").  If the input slope faces left the
geometry is mirrored (x -> -x) on :meth:`SlopeModel.prepare` and every result
is mirrored back by :meth:`SlopeModel.to_world`.
"""

import numpy as np
from shapely.geometry import (LineString, Polygon, MultiPolygon,
                              MultiLineString, GeometryCollection)
from shapely.ops import unary_union
from shapely import affinity

from .i18n import tr

GAMMA_W = 9.81


# --------------------------------------------------------------------------
#  helpers
# --------------------------------------------------------------------------
def _polygons(geom):
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, (MultiPolygon, GeometryCollection)):
        out = []
        for g in geom.geoms:
            out.extend(_polygons(g))
        return out
    return []


def _lines(geom):
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, LineString):
        return [geom]
    if isinstance(geom, (MultiLineString, GeometryCollection)):
        out = []
        for g in geom.geoms:
            out.extend(_lines(g))
        return out
    if isinstance(geom, Polygon):
        return [LineString(geom.exterior.coords)]
    return []


def _mirror(geom):
    return affinity.scale(geom, xfact=-1.0, yfact=1.0, origin=(0, 0))


class DistributedLoad(object):
    """Uniform vertical surcharge q (kPa = kN/m per m of horizontal length)
    applied on the ground surface over the horizontal extent of ``geom``."""

    def __init__(self, geom, q, name=""):
        self.geom = geom
        self.q = float(q)
        self.name = name

    def interval(self):
        x0, _, x1, _ = self.geom.bounds
        return min(x0, x1), max(x0, x1)


class Reinforcement(object):
    """Linear reinforcement (geogrid, nail, tie-back, anchor...).

    ``capacity``  tensile capacity T (kN per m of slope width).
    ``pullout``   pull-out resistance per unit length behind the slip surface
                  (kN/m/m).  0 = unlimited (capacity governs).
    ``passive``   True: force is treated as a resisting force and factored by FS
                  (passive), False: subtracted from the driving forces (active).
    """

    def __init__(self, geom, capacity, pullout=0.0, passive=True, name=""):
        self.geom = geom
        self.capacity = float(capacity)
        self.pullout = float(pullout or 0.0)
        self.passive = bool(passive)
        self.name = name


class Pile(object):
    """Pile / retaining wall / shear key providing a shear resistance V
    (kN per m of slope width) where the slip surface cuts it."""

    def __init__(self, geom, shear, passive=True, name=""):
        self.geom = geom
        self.shear = float(shear)
        self.passive = bool(passive)
        self.name = name


# --------------------------------------------------------------------------
#  model
# --------------------------------------------------------------------------
class SlopeModel(object):
    """Container for every analysis input."""

    def __init__(self):
        self.materials = []          # list[Material]
        self.regions = []            # list[(Polygon, material index)]
        self.piezo = None            # LineString or None
        self.water_depth = None      # water table parallel to ground at depth d
        self.gamma_w = GAMMA_W
        self.loads = []              # list[DistributedLoad]
        self.reinforcements = []     # list[Reinforcement]
        self.piles = []              # list[Pile]
        self.global_support = 0.0    # v2.x compatible "wall force" w/o geometry
        self.global_support_passive = True
        self.kh = 0.0
        self.kv = 0.0
        self.crack_depth = 0.0
        self.crack_fill = 0.0        # 0..1 fraction of crack filled with water
        self.direction = 0           # 0 auto, +1 slides to +x, -1 slides to -x
        self.min_depth = 0.0         # minimum slide depth filter (0 = auto)
        self._prepared = False

    # ------------------------------------------------------------ building
    def add_material(self, material):
        self.materials.append(material)
        return len(self.materials) - 1

    def add_region(self, geom, material_index):
        for p in _polygons(geom):
            if p.area > 0:
                self.regions.append((p, int(material_index)))

    # --------------------------------------------------------- validation
    def preflight(self):
        """Return a list of (level, message) describing input problems."""
        msgs = []
        if not self.regions:
            msgs.append(("error", tr("The geology layer contains no valid polygons.")))
            return msgs
        for i, m in enumerate(self.materials):
            for level, text in m.validate():
                msgs.append((level, u"{0}: {1}".format(m.name, text)))
        union = unary_union([r[0] for r in self.regions])
        area_sum = sum(r[0].area for r in self.regions)
        if area_sum > union.area * 1.001:
            msgs.append(("warning", tr("Geology polygons overlap; the first polygon found is used in overlapping areas.")))
        if isinstance(union, MultiPolygon):
            msgs.append(("warning", tr("Geology polygons do not form one continuous body (gaps between polygons).")))
        x0, y0, x1, y1 = union.bounds
        if (x1 - x0) > 1e5 or (y1 - y0) > 1e5:
            msgs.append(("warning", tr("The model is very large. Make sure the layer uses section coordinates (x = distance, y = elevation) in metres.")))
        if self.kh < 0 or self.kh > 1:
            msgs.append(("error", tr("Horizontal seismic coefficient kh must be between 0 and 1.")))
        if self.crack_depth < 0:
            msgs.append(("error", tr("Tension crack depth cannot be negative.")))
        if self.piezo is not None:
            xs = [c[0] for c in self.piezo.coords]
            if any(b < a for a, b in zip(xs, xs[1:])) and any(b > a for a, b in zip(xs, xs[1:])):
                msgs.append(("warning", tr("The piezometric line is not monotonic in x; it is resampled as a function of x.")))
        return msgs

    # ------------------------------------------------------------ prepare
    def prepare(self):
        if not self.regions:
            raise ValueError(tr("The geology layer contains no valid polygons."))
        self._union_world = unary_union([r[0] for r in self.regions])
        self._build_edges([r[0] for r in self.regions])
        sign = self.direction
        if sign == 0:
            sign = self._detect_direction()
        self.sign = 1 if sign >= 0 else -1
        if self.sign < 0:
            self._regions = [(_mirror(p), m) for p, m in self.regions]
        else:
            self._regions = list(self.regions)
        polys = [r[0] for r in self._regions]
        self.union = unary_union(polys)
        self._build_edges(polys)
        self._build_ground()
        self._build_piezo()
        self._loads = []
        for ld in self.loads:
            g = self.to_model_geom(ld.geom)
            x0, _, x1, _ = g.bounds
            self._loads.append((min(x0, x1), max(x0, x1), ld.q))
        self._reinf = [(self.to_model_geom(r.geom), r) for r in self.reinforcements]
        self._piles = [(self.to_model_geom(p.geom), p) for p in self.piles]
        self.relief = max(self.gy.max() - self.gy.min(), 1e-6)
        self.width = max(self.xmax - self.xmin, 1e-6)
        self._prepared = True
        return self

    # coordinate transforms -------------------------------------------------
    def to_model_geom(self, geom):
        return _mirror(geom) if getattr(self, "sign", 1) < 0 else geom

    def to_world_xy(self, xs, ys):
        xs = np.asarray(xs, dtype=float)
        return (xs * self.sign, np.asarray(ys, dtype=float))

    def to_world_x(self, x):
        return x * self.sign

    def to_model_x(self, x):
        return x * self.sign

    # edges -------------------------------------------------------------------
    def _build_edges(self, polys):
        x1, y1, x2, y2, rid = [], [], [], [], []
        for i, p in enumerate(polys):
            rings = [p.exterior] + list(p.interiors)
            for ring in rings:
                c = np.asarray(ring.coords)
                x1.append(c[:-1, 0])
                y1.append(c[:-1, 1])
                x2.append(c[1:, 0])
                y2.append(c[1:, 1])
                rid.append(np.full(len(c) - 1, i))
        self._ex1 = np.concatenate(x1)
        self._ey1 = np.concatenate(y1)
        self._ex2 = np.concatenate(x2)
        self._ey2 = np.concatenate(y2)
        self._erid = np.concatenate(rid)
        self._elo = np.minimum(self._ex1, self._ex2)
        self._ehi = np.maximum(self._ex1, self._ex2)
        dx = self._ex2 - self._ex1
        self._eslope = np.where(np.abs(dx) > 0, (self._ey2 - self._ey1) / np.where(dx == 0, 1, dx), 0.0)
        self._region_cols = [np.nonzero(self._erid == i)[0] for i in range(len(polys))]
        allx = np.concatenate([self._ex1, self._ex2])
        self.xmin, self.xmax = float(allx.min()), float(allx.max())
        ally = np.concatenate([self._ey1, self._ey2])
        self.ymin, self.ymax = float(ally.min()), float(ally.max())

    def _crossings(self, xs):
        """Matrix (len(xs) x n_edges) of edge elevations at xs (NaN = no crossing)."""
        xs = np.asarray(xs, dtype=float)[:, None]
        inside = (xs >= self._elo) & (xs < self._ehi)
        y = self._ey1 + (xs - self._ex1) * self._eslope
        return np.where(inside, y, np.nan)

    def _detect_direction(self):
        xs = np.linspace(self.xmin, self.xmax, 201)[1:-1]
        Y = self._crossings(xs)
        with np.errstate(all="ignore"):
            top = np.nanmax(np.where(np.isnan(Y), -np.inf, Y), axis=1)
        n = max(3, len(top) // 5)
        return 1 if np.mean(top[:n]) >= np.mean(top[-n:]) else -1

    # ground surface -----------------------------------------------------------
    def _build_ground(self):
        """Ground surface = upper envelope of the union polygon."""
        xs = set()
        for p, _ in self._regions:
            xs.update(np.asarray(p.exterior.coords)[:, 0].tolist())
        xs = np.array(sorted(xs))
        # sample slightly inside each vertex pair to obtain the envelope
        eps = max(self.xmax - self.xmin, 1.0) * 1e-9
        pts_x = []
        for i, x in enumerate(xs):
            for xx in (x - eps, x + eps) if 0 < i < len(xs) - 1 else ((x + eps,) if i == 0 else (x - eps,)):
                pts_x.append(xx)
        pts_x = np.array(pts_x)
        top = self._top(pts_x)
        ok = np.isfinite(top)
        self.gx = pts_x[ok]
        self.gy = top[ok]
        self.gvx = xs  # vertex x positions of the surface

    def _top(self, xs):
        Y = self._crossings(xs)
        with np.errstate(all="ignore"):
            t = np.nanmax(np.where(np.isnan(Y), -np.inf, Y), axis=1)
        t[~np.isfinite(t)] = np.nan
        return t

    def ground(self, xs):
        return np.interp(xs, self.gx, self.gy)

    def ground_line(self):
        return LineString(list(zip(self.gx, self.gy)))

    def bottom(self, xs):
        Y = self._crossings(np.atleast_1d(xs))
        with np.errstate(all="ignore"):
            b = np.nanmin(np.where(np.isnan(Y), np.inf, Y), axis=1)
        b[~np.isfinite(b)] = np.nan
        return b

    # water ----------------------------------------------------------------------
    def _build_piezo(self):
        self._px = None
        if self.piezo is not None and not self.piezo.is_empty:
            lines = _lines(self.to_model_geom(self.piezo))
            pts = []
            for ln in lines:
                pts.extend(list(ln.coords))
            pts.sort(key=lambda p: p[0])
            if len(pts) >= 2:
                self._px = np.array([p[0] for p in pts])
                self._py = np.array([p[1] for p in pts])

    @property
    def has_water(self):
        return self._px is not None or self.water_depth is not None

    def water_level(self, xs):
        xs = np.asarray(xs, dtype=float)
        if self._px is not None:
            return np.interp(xs, self._px, self._py)
        if self.water_depth is not None:
            return self.ground(xs) - self.water_depth
        return np.full(xs.shape, -np.inf)

    def water_line_world(self):
        """Polyline of the active water table in world coordinates (for plots)."""
        if self._px is not None:
            x, y = self._px, self._py
        elif self.water_depth is not None:
            x = self.gx
            y = self.gy - self.water_depth
        else:
            return None
        wx, wy = self.to_world_xy(x, y)
        return list(zip(wx.tolist(), wy.tolist()))

    # columns ---------------------------------------------------------------------
    def columns(self, xs, ybase):
        """Vertical column properties at abscissae ``xs`` above elevations ``ybase``.

        Returns dict with arrays: ground, bottom, wpa (weight per unit width of
        column = sum gamma*thickness, i.e. total vertical stress at base),
        mat (material index at base, -1 if none), hw (water level).
        """
        xs = np.asarray(xs, dtype=float)
        yb = np.asarray(ybase, dtype=float)
        Y = self._crossings(xs)
        hw = self.water_level(xs)
        n = len(xs)
        wpa = np.zeros(n)
        mat = np.full(n, -1, dtype=int)
        with np.errstate(all="ignore"):
            ground = np.nanmax(np.where(np.isnan(Y), -np.inf, Y), axis=1)
            bottom = np.nanmin(np.where(np.isnan(Y), np.inf, Y), axis=1)
        gap = np.full(n, np.inf)
        near = np.full(n, -1, dtype=int)
        for r, cols in enumerate(self._region_cols):
            if len(cols) == 0:
                continue
            m = self.materials[self._regions[r][1]]
            Yr = np.sort(Y[:, cols], axis=1)  # NaN sorted to the end
            k = Yr.shape[1]
            for j in range(0, k - 1, 2):
                lo = Yr[:, j]
                hi = Yr[:, j + 1]
                valid = ~(np.isnan(lo) | np.isnan(hi))
                if not valid.any():
                    break
                lo_v = np.where(valid, lo, 0.0)
                hi_v = np.where(valid, hi, 0.0)
                bot = np.maximum(lo_v, yb)
                t = np.where(valid, np.maximum(hi_v - bot, 0.0), 0.0)
                tw = np.where(valid, np.clip(np.minimum(hi_v, hw) - bot, 0.0, None), 0.0)
                tw = np.minimum(tw, t)
                wpa += m.gamma * (t - tw) + m.gamma_sat * tw
                at_base = valid & (lo_v <= yb) & (yb < hi_v) & (mat < 0)
                mat[at_base] = self._regions[r][1]
                # nearest material above the base (for slip surfaces that run
                # exactly along a boundary or through tiny gaps)
                d = np.where(valid & (lo_v >= yb), lo_v - yb, np.inf)
                upd = d < gap
                gap = np.where(upd, d, gap)
                near = np.where(upd, self._regions[r][1], near)
        tol = max(getattr(self, "relief", 1.0) * 1e-3, 1e-3)
        fix = (mat < 0) & (gap < tol)
        mat[fix] = near[fix]
        ground[~np.isfinite(ground)] = np.nan
        bottom[~np.isfinite(bottom)] = np.nan
        return {"ground": ground, "bottom": bottom, "wpa": wpa, "mat": mat, "hw": hw}

    # copies for reliability analysis ----------------------------------------------
    def with_materials(self, materials):
        """Shallow copy of a prepared model with a different material list."""
        m = SlopeModel.__new__(SlopeModel)
        m.__dict__.update(self.__dict__)
        m.materials = materials
        return m

    def with_conditions(self, kh=None, water_depth="keep", piezo="keep"):
        m = SlopeModel.__new__(SlopeModel)
        m.__dict__.update(self.__dict__)
        if kh is not None:
            m.kh = kh
        if water_depth != "keep":
            m.water_depth = water_depth
        if piezo != "keep":
            m.piezo = piezo
            m._build_piezo()
        return m
