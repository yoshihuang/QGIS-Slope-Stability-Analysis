# -*- coding: utf-8 -*-
"""Self-contained HTML analysis report with inline SVG figures."""

import datetime
import html
import math

import numpy as np

from .methods import method_label, METHODS
from .i18n import tr

MAT_COLORS = ["#c9b18f", "#a7c4a0", "#d9c27e", "#b7a8c9", "#9fc1d0", "#d6a59a",
              "#c2c98f", "#e0b5d0", "#a8b8a0", "#d8cfb0", "#b0c4de", "#e3c08d"]


def fs_color(fs):
    if fs is None:
        return "#888888"
    if fs < 1.0:
        return "#c0392b"
    if fs < 1.2:
        return "#e67e22"
    if fs < 1.5:
        return "#f1c40f"
    if fs < 2.0:
        return "#7fb800"
    return "#2e8b57"


def _e(s):
    return html.escape(u"%s" % s)


def _f(v, nd=3):
    if v is None:
        return "-"
    try:
        if not math.isfinite(v):
            return "-"
    except TypeError:
        return _e(v)
    return ("%." + str(nd) + "f") % v


# --------------------------------------------------------------------------
#  SVG cross section
# --------------------------------------------------------------------------
class _Frame(object):
    def __init__(self, x0, y0, x1, y1, width=900, pad=46):
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1
        sx = (width - 2 * pad) / max(x1 - x0, 1e-9)
        self.s = sx
        self.pad = pad
        self.w = width
        self.h = int((y1 - y0) * sx + 2 * pad)
        if self.h > width * 0.9:
            self.s = (width * 0.9 - 2 * pad) / max(y1 - y0, 1e-9)
            self.h = int(width * 0.9)

    def X(self, x):
        return self.pad + (x - self.x0) * self.s

    def Y(self, y):
        return self.h - self.pad - (y - self.y0) * self.s

    def pts(self, coords):
        return " ".join("%.1f,%.1f" % (self.X(x), self.Y(y)) for x, y in coords)


def _nice_ticks(a, b, n=8):
    span = max(b - a, 1e-9)
    raw = span / n
    mag = 10 ** math.floor(math.log10(raw))
    step = min((s * mag for s in (1, 2, 2.5, 5, 10) if s * mag >= raw), default=mag * 10)
    t = math.ceil(a / step) * step
    out = []
    while t <= b + 1e-9:
        out.append(t)
        t += step
    return out


def section_svg(model, sr, n_surfaces=40, show_slices=True, width=900, crack_point=None):
    """Cross-section with materials, water, surfaces and the critical slip."""
    regions = model.regions
    xs, ys = [], []
    for p, _ in regions:
        x0, y0, x1, y1 = p.bounds
        xs += [x0, x1]
        ys += [y0, y1]
    crit = sr.critical if sr is not None else None
    if crit is not None and crit.surface.kind == "circular":
        c = crit.center_world(model)
        xs.append(c[0])
        ys.append(c[1])
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    dx, dy = x1 - x0, y1 - y0
    fr = _Frame(x0 - 0.02 * dx, y0 - 0.02 * dy, x1 + 0.02 * dx, y1 + 0.05 * dy, width)
    out = [u'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%" style="max-width:%dpx;background:#fff;border:1px solid #ccc">' % (fr.w, fr.h, fr.w)]
    # axes
    for t in _nice_ticks(fr.x0, fr.x1):
        X = fr.X(t)
        out.append(u'<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#eee"/>' % (X, fr.pad / 2, X, fr.h - fr.pad))
        out.append(u'<text x="%.1f" y="%d" font-size="10" text-anchor="middle" fill="#666">%g</text>' % (X, fr.h - fr.pad + 14, round(t, 6)))
    for t in _nice_ticks(fr.y0, fr.y1, 6):
        Y = fr.Y(t)
        out.append(u'<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#eee"/>' % (fr.pad, Y, fr.w - fr.pad / 2, Y))
        out.append(u'<text x="%d" y="%.1f" font-size="10" text-anchor="end" fill="#666">%g</text>' % (fr.pad - 4, Y + 3, round(t, 6)))
    # materials
    for p, mi in regions:
        col = MAT_COLORS[mi % len(MAT_COLORS)]
        out.append(u'<polygon points="%s" fill="%s" stroke="#555" stroke-width="0.8"/>' % (fr.pts(p.exterior.coords), col))
    # water
    if sr is not None:
        wl = sr.model.water_line_world()
        if wl:
            out.append(u'<polyline points="%s" fill="none" stroke="#1f6fd1" stroke-width="1.6" stroke-dasharray="6,3"/>' % fr.pts(wl))
    # supports
    for ld in model.loads:
        a, b = ld.interval()
        g = model.ground(np.array([model.to_model_x(a), model.to_model_x(b)]))
        for xx in np.linspace(a, b, 6):
            gy = float(model.ground(np.array([model.to_model_x(xx)]))[0])
            out.append(u'<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#8e44ad" stroke-width="1.2" marker-end="url(#arr)"/>' % (fr.X(xx), fr.Y(gy) - 18, fr.X(xx), fr.Y(gy) - 2))
    for rf in model.reinforcements:
        for ln in _iter_lines(rf.geom):
            out.append(u'<polyline points="%s" fill="none" stroke="#2c3e50" stroke-width="1.6"/>' % fr.pts(ln.coords))
    for pl in model.piles:
        g = pl.geom
        if g.geom_type in ("Polygon", "MultiPolygon"):
            for p in getattr(g, "geoms", [g]):
                out.append(u'<polygon points="%s" fill="#7f8c8d" stroke="#2c3e50"/>' % fr.pts(p.exterior.coords))
        else:
            for ln in _iter_lines(g):
                out.append(u'<polyline points="%s" fill="none" stroke="#7f8c8d" stroke-width="3"/>' % fr.pts(ln.coords))
    out.append(u'<defs><marker id="arr" markerWidth="6" markerHeight="6" refX="3" refY="5" orient="auto"><path d="M0,0 L6,0 L3,6 z" fill="#8e44ad"/></marker></defs>')
    # surfaces
    if sr is not None and sr.records:
        recs = sr.records[:n_surfaces][::-1]
        for r in recs:
            out.append(u'<polyline points="%s" fill="none" stroke="%s" stroke-width="0.8" opacity="0.7"/>' % (fr.pts(r.polyline(sr.model)), fs_color(r.fs)))
    if crit is not None:
        S = sr.slices
        if show_slices and S is not None and S.valid:
            m = sr.model
            for i in range(S.n):
                xa, ya = m.to_world_xy([S.xl[i], S.xr[i], S.xr[i], S.xl[i]], [S.ybl[i], S.ybr[i], S.gr[i], S.gl[i]])
                out.append(u'<polygon points="%s" fill="none" stroke="#c0392b" stroke-width="0.4" opacity="0.8"/>' % fr.pts(zip(xa, ya)))
        out.append(u'<polyline points="%s" fill="none" stroke="#c0392b" stroke-width="2.6"/>' % fr.pts(crit.polyline(sr.model)))
        if crit.surface.kind == "circular":
            c = crit.center_world(sr.model)
            pl = crit.polyline(sr.model)
            for p in (pl[0], pl[-1]):
                out.append(u'<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#c0392b" stroke-width="0.6" stroke-dasharray="3,3"/>' % (fr.X(c[0]), fr.Y(c[1]), fr.X(p[0]), fr.Y(p[1])))
            out.append(u'<circle cx="%.1f" cy="%.1f" r="3.5" fill="#c0392b"/>' % (fr.X(c[0]), fr.Y(c[1])))
            lx, ly = fr.X(c[0]) + 6, fr.Y(c[1]) - 6
        else:
            pl = crit.polyline(sr.model)
            lx, ly = fr.X(pl[0][0]), fr.Y(pl[0][1]) - 10
        out.append(u'<text x="%.1f" y="%.1f" font-size="15" font-weight="bold" fill="#c0392b">FS = %.3f</text>' % (lx, ly, crit.fs))
    # crack point
    if crack_point is not None:
        cx_ = crack_point[0]
        gy = float(model.ground(np.array([model.to_model_x(cx_)]))[0])
        X, Y = fr.X(cx_), fr.Y(gy)
        out.append(u'<path d="M%.1f,%.1f L%.1f,%.1f L%.1f,%.1f Z" fill="#c0392b" stroke="#fff" stroke-width="1"/>' % (X, Y, X - 7, Y - 13, X + 7, Y - 13))
    # legend
    ly = 16
    seen = set()
    for p, mi in regions:
        if mi in seen:
            continue
        seen.add(mi)
        out.append(u'<rect x="%d" y="%d" width="12" height="10" fill="%s" stroke="#555"/>' % (fr.w - 190, ly, MAT_COLORS[mi % len(MAT_COLORS)]))
        out.append(u'<text x="%d" y="%d" font-size="11">%s</text>' % (fr.w - 172, ly + 9, _e(model.materials[mi].name)))
        ly += 15
    out.append(u"</svg>")
    return u"\n".join(out)


def _iter_lines(g):
    if g is None or g.is_empty:
        return []
    if g.geom_type == "LineString":
        return [g]
    return [x for x in getattr(g, "geoms", []) if x.geom_type == "LineString"]


def histogram_svg(rel, width=640, height=240):
    mc = rel.get("mc")
    if not mc:
        return u""
    hist, edges = mc["hist"], mc["edges"]
    pad = 36
    hmax = max(hist) or 1
    bw = (width - 2 * pad) / float(len(hist))
    out = [u'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%" style="max-width:%dpx;background:#fff;border:1px solid #ccc">' % (width, height, width)]
    for i, h in enumerate(hist):
        x = pad + i * bw
        hh = (height - 2 * pad) * h / hmax
        col = "#c0392b" if edges[i + 1] <= 1.0 else "#3d7ab8"
        out.append(u'<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>' % (x + 1, height - pad - hh, bw - 2, hh, col))
    for i in range(0, len(edges), 5):
        out.append(u'<text x="%.1f" y="%d" font-size="10" text-anchor="middle">%.2f</text>' % (pad + i * bw, height - pad + 14, edges[i]))
    if edges[0] < 1.0 < edges[-1]:
        x1 = pad + (1.0 - edges[0]) / (edges[-1] - edges[0]) * (width - 2 * pad)
        out.append(u'<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#000" stroke-dasharray="4,3"/>' % (x1, pad / 2, x1, height - pad))
        out.append(u'<text x="%.1f" y="%d" font-size="11">FS = 1</text>' % (x1 + 4, pad / 2 + 10))
    out.append(u'<text x="%d" y="%d" font-size="11" text-anchor="middle">%s</text>' % (width / 2, height - 6, _e(tr("Factor of safety"))))
    out.append(u"</svg>")
    return u"\n".join(out)


# --------------------------------------------------------------------------
#  HTML
# --------------------------------------------------------------------------
CSS = u"""
body{font-family:-apple-system,'Segoe UI','Noto Sans','Noto Sans CJK TC','Microsoft JhengHei',sans-serif;margin:24px auto;max-width:1000px;color:#222;line-height:1.45;padding:0 16px}
h1{font-size:24px;margin-bottom:4px}h2{font-size:19px;border-bottom:2px solid #2b5c8f;padding-bottom:3px;margin-top:30px}
h3{font-size:16px;margin-top:22px}
table{border-collapse:collapse;margin:8px 0;font-size:13px}th,td{border:1px solid #ccc;padding:3px 7px;text-align:right}
th{background:#eef3f8}td.l,th.l{text-align:left}.muted{color:#777;font-size:12px}
.fs{font-weight:bold;font-size:20px}.warn{color:#b9770e}.err{color:#c0392b}
.badge{display:inline-block;padding:1px 8px;border-radius:10px;color:#fff;font-weight:bold}
details{margin:8px 0}summary{cursor:pointer;color:#2b5c8f}
.scroll{overflow-x:auto}
"""


def build_report(result, info=None):
    info = info or {}
    model = result.model
    o = []
    o.append(u'<!DOCTYPE html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">')
    o.append(u"<title>%s</title><style>%s</style></head><body>" % (_e(tr("Slope stability analysis report")), CSS))
    o.append(u"<h1>%s</h1>" % _e(tr("Slope stability analysis report")))
    o.append(u'<p class="muted">%s %s &middot; %s &middot; %s %.1f s</p>' % (
        _e(info.get("plugin", "Slope Stability")), _e(info.get("version", "")),
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        _e(tr("Run time")), result.elapsed))
    if info.get("project"):
        o.append(u"<p>%s: %s</p>" % (_e(tr("Geology layer")), _e(info["project"])))

    # ---- pre-flight ------------------------------------------------------------
    if result.preflight:
        o.append(u"<h2>%s</h2><ul>" % _e(tr("Input checks")))
        for level, msg in result.preflight:
            cls = "err" if level == "error" else "warn"
            o.append(u'<li class="%s">%s</li>' % (cls, _e(msg)))
        o.append(u"</ul>")

    # ---- summary ---------------------------------------------------------------
    o.append(u"<h2>%s</h2>" % _e(tr("Summary")))
    o.append(u'<table><tr><th class="l">%s</th><th>kh</th><th class="l">%s</th><th>%s</th><th class="l">%s</th><th>%s</th><th>%s</th></tr>' % (
        _e(tr("Scenario")), _e(tr("Water")), _e(tr("Minimum FS")), _e(tr("Method")),
        _e(tr("Entry x")), _e(tr("Exit x"))))
    for sr in result.scenarios:
        c = sr.critical
        m = sr.model
        wtxt = _water_text(m)
        if c is None:
            o.append(u'<tr><td class="l">%s</td><td>%s</td><td class="l">%s</td><td colspan="4" class="l err">%s</td></tr>' % (
                _e(sr.scenario.label), _f(m.kh, 2), _e(wtxt), _e(sr.message)))
            continue
        ex = sorted([m.to_world_x(c.x_in), m.to_world_x(c.x_out)])
        o.append(u'<tr><td class="l">%s</td><td>%s</td><td class="l">%s</td><td><span class="badge" style="background:%s">%s</span></td><td class="l">%s</td><td>%s</td><td>%s</td></tr>' % (
            _e(sr.scenario.label), _f(m.kh, 2), _e(wtxt), fs_color(c.fs), _f(c.fs),
            _e(method_label(getattr(sr, "primary_method", ""))), _f(ex[0], 2), _f(ex[1], 2)))
    o.append(u"</table>")

    # ---- inputs ------------------------------------------------------------------
    o.append(u"<h2>%s</h2>" % _e(tr("Input data")))
    o.append(u'<table><tr><th class="l">%s</th><th class="l">%s</th><th>&gamma;</th><th>&gamma;<sub>sat</sub></th><th>c / Su</th><th>&phi;</th><th class="l">%s</th><th>ru</th></tr>' % (
        _e(tr("Material")), _e(tr("Strength model")), _e(tr("Other parameters"))))
    for m in model.materials:
        other = u""
        if m.model == "hb":
            other = u"&sigma;ci=%g, GSI=%g, mi=%g, D=%g" % (m.sigci, m.gsi, m.mi, m.disturbance)
        if m.c_std or m.phi_std or m.gamma_std:
            other += u" &sigma;(c,&phi;,&gamma;)=(%g, %g, %g)" % (m.c_std, m.phi_std, m.gamma_std)
        o.append(u'<tr><td class="l">%s</td><td class="l">%s</td><td>%g</td><td>%g</td><td>%g</td><td>%g</td><td class="l">%s</td><td>%g</td></tr>' % (
            _e(m.name), _e(_model_name(m.model)), m.gamma, m.gamma_sat, m.cohesion, m.phi, other, m.ru))
    o.append(u"</table>")
    st = result.settings
    rows = [
        (tr("Slip surface type"), _surface_type_name(st.surface_type)),
        (tr("Search method"), method_label(st.method)),
        (tr("Number of slices"), st.n_slices),
        (tr("Search density"), u"%d %%" % round(st.density * 100)),
        (tr("Optimisation"), tr("Yes") if st.optimize else tr("No")),
        (tr("Seismic coefficient kh / kv"), u"%g / %g" % (model.kh, model.kv)),
        (tr("Tension crack depth"), u"%g (%s %d %%)" % (model.crack_depth, tr("water filled"), round(model.crack_fill * 100))),
        (tr("Unit weight of water"), model.gamma_w),
        (tr("Distributed loads"), len(model.loads)),
        (tr("Reinforcements"), len(model.reinforcements)),
        (tr("Piles / walls"), len(model.piles)),
        (tr("Slope direction"), tr("Slides towards +x") if model.sign > 0 else tr("Slides towards -x")),
    ]
    if st.surface_type == "point" and getattr(st, "crack_point", None):
        cpx = st.crack_point[0]
        cpy = float(model.ground(np.array([model.to_model_x(cpx)]))[0])
        rows.insert(1, (tr("Crack point (x, ground y)"), u"(%s, %s)" % (_f(cpx, 2), _f(cpy, 2))))
    o.append(u"<table>")
    for k, v in rows:
        o.append(u'<tr><th class="l">%s</th><td class="l">%s</td></tr>' % (_e(k), _e(v)))
    o.append(u"</table>")

    # ---- scenario details -------------------------------------------------------------
    for k, sr in enumerate(result.scenarios):
        o.append(u"<h2>%s: %s</h2>" % (_e(tr("Scenario")), _e(sr.scenario.label)))
        if sr.critical is None:
            o.append(u'<p class="err">%s</p>' % _e(sr.message))
            continue
        c = sr.critical
        o.append(u'<p>%s: <span class="fs" style="color:%s">%s</span> &nbsp; <span class="muted">%s: %d, %s: %d, %.1f s</span></p>' % (
            _e(tr("Minimum factor of safety")), fs_color(c.fs), _f(c.fs), _e(tr("surfaces evaluated")),
            sr.n_eval, _e(tr("valid")), sr.n_valid, sr.elapsed))
        if k == 0 or len(result.scenarios) <= 4:
            o.append(u'<div class="scroll">%s</div>' % section_svg(model, sr, crack_point=getattr(st, "crack_point", None) if st.surface_type == "point" else None))
        geo = []
        if c.surface.kind == "circular":
            cw = c.center_world(sr.model)
            geo.append(u"%s (x, y) = (%s, %s), R = %s" % (_e(tr("Centre")), _f(cw[0], 2), _f(cw[1], 2), _f(c.surface.r, 2)))
        pl = c.polyline(sr.model)
        geo.append(u"%s (%s, %s) &rarr; %s (%s, %s)" % (_e(tr("Entry")), _f(pl[0][0], 2), _f(pl[0][1], 2), _e(tr("Exit")), _f(pl[-1][0], 2), _f(pl[-1][1], 2)))
        o.append(u"<p>%s</p>" % u"<br>".join(geo))
        o.append(u'<table><tr><th class="l">%s</th><th>FS</th><th>&lambda;</th><th>&theta; (&deg;)</th><th class="l">%s</th></tr>' % (_e(tr("Method")), _e(tr("Note"))))
        for mk in METHODS:
            r = sr.method_results.get(mk)
            if r is None:
                continue
            note = r.message
            if r.fs_uncorrected is not None:
                note = u"%s %s" % (tr("uncorrected FS ="), _f(r.fs_uncorrected))
            if r.malpha_warning:
                note += u" " + tr("(m-alpha limited)")
            o.append(u'<tr><td class="l">%s</td><td>%s</td><td>%s</td><td>%s</td><td class="l">%s</td></tr>' % (
                _e(method_label(mk)), _f(r.fs) if r.ok else "-", _f(r.lam), _f(r.theta, 1), _e(note)))
        o.append(u"</table>")
        S = sr.slices
        if S is not None and S.valid and S.point_forces:
            o.append(u'<table><tr><th class="l">%s</th><th>x</th><th>y</th><th>Fx</th><th>Fy</th><th class="l">%s</th></tr>' % (_e(tr("Point force")), _e(tr("Mode"))))
            for (i, px, py, fx, fy, passive, lbl) in S.point_forces:
                wx = sr.model.to_world_x(px)
                o.append(u'<tr><td class="l">%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td class="l">%s</td></tr>' % (
                    _e(_pf_label(lbl)), _f(wx, 2), _f(py, 2), _f(fx * sr.model.sign, 1), _f(fy, 1), _e(tr("passive") if passive else tr("active"))))
            o.append(u"</table>")
        if sr.slice_rows:
            o.append(u"<details><summary>%s</summary><div class='scroll'><table>" % _e(tr("Slice data (critical surface)")))
            o.append(u"<tr><th>#</th><th>x</th><th>b</th><th>h</th><th>&alpha;</th><th>W</th><th>Q</th><th class='l'>%s</th><th>c</th><th>&phi;</th><th>u</th><th>N'</th><th>&sigma;'n</th><th>S<sub>m</sub></th></tr>" % _e(tr("Material")))
            for r in sr.slice_rows:
                o.append(u"<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td class='l'>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
                    r["slice"], _f(sr.model.to_world_x(r["x"]), 2), _f(r["width"], 2), _f(r["height"], 2), _f(r["alpha"] * sr.model.sign, 1),
                    _f(r["weight"], 1), _f(r["surcharge"], 1), _e(r["material"]), _f(r["c"], 1), _f(r["phi"], 1),
                    _f(r["u"], 1), _f(r["N_eff"], 1), _f(r["sigma_eff"], 1), _f(r["shear_mob"], 1)))
            o.append(u"</table></div></details>")

    # ---- reliability ------------------------------------------------------------------
    rel = result.reliability
    if rel:
        o.append(u"<h2>%s</h2>" % _e(tr("Probabilistic analysis")))
        if rel.get("error"):
            o.append(u'<p class="warn">%s</p>' % _e(rel["error"]))
        else:
            o.append(u"<p>%s: %s, %s = %s</p>" % (_e(tr("Method")), _e(method_label(rel["method"])), _e(tr("deterministic FS")), _f(rel["fs_deterministic"])))
            o.append(u'<table><tr><th class="l">%s</th><th class="l">%s</th><th>%s</th><th>%s</th></tr>' % (_e(tr("Material")), _e(tr("Variable")), _e(tr("Mean")), _e(tr("Std. dev."))))
            for (n, p, mu, sd) in rel["variables"]:
                o.append(u'<tr><td class="l">%s</td><td class="l">%s</td><td>%g</td><td>%g</td></tr>' % (_e(n), _e(p), mu, sd))
            o.append(u"</table>")
            t = rel.get("tspm")
            mc = rel.get("mc")
            o.append(u'<table><tr><th class="l"></th><th>%s</th><th>%s</th><th>&beta;<sub>normal</sub></th><th>&beta;<sub>lognormal</sub></th><th>P<sub>f</sub></th></tr>' % (_e(tr("Mean FS")), _e(tr("Std. dev."))))
            if mc:
                o.append(u'<tr><td class="l">%s (n = %d)</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (
                    _e(tr("Monte Carlo")), mc["n"], _f(mc["mean"]), _f(mc["std"]), _f(mc["beta_normal"], 2), _f(mc["beta_ln"], 2), _pct(mc["pf"])))
            if t:
                o.append(u'<tr><td class="l">%s</td><td>%s</td><td>%s</td><td>-</td><td>%s</td><td>%s</td></tr>' % (
                    _e(tr("Taylor series (TSPM)")), _f(rel["fs_deterministic"]), _f(t["sigma"]), _f(t["beta_ln"], 2), _pct(t["pf"])))
            o.append(u"</table>")
            if mc:
                o.append(histogram_svg(rel))
            if t and t.get("contributions"):
                o.append(u'<table><tr><th class="l">%s</th><th class="l">%s</th><th>FS(+&sigma;)</th><th>FS(-&sigma;)</th><th>%s</th></tr>' % (_e(tr("Material")), _e(tr("Variable")), _e(tr("Share of variance"))))
                for (n, p, fp, fm, sh) in sorted(t["contributions"], key=lambda c: -c[4]):
                    o.append(u'<tr><td class="l">%s</td><td class="l">%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (_e(n), _e(p), _f(fp), _f(fm), _pct(sh)))
                o.append(u"</table>")

    o.append(u'<h2>%s</h2><ul class="muted">' % _e(tr("Notes")))
    for note in (
            tr("Two-dimensional limit equilibrium analysis (plane strain), unit width of slope."),
            tr("Pore pressure from the piezometric line is u = gamma_w x vertical distance below the line."),
            tr("Janbu simplified results include the Janbu correction factor f0."),
            tr("Passive support forces are divided by the factor of safety; active support forces reduce the driving forces."),
            tr("Results must be checked by a qualified engineer.")):
        o.append(u"<li>%s</li>" % _e(note))
    o.append(u"</ul></body></html>")
    return u"\n".join(o)


def _pct(v):
    if v is None:
        return "-"
    if v < 1e-4 and v > 0:
        return "%.2e" % v
    return "%.2f %%" % (100.0 * v)


def _water_text(m):
    if m._px is not None:
        return tr("Piezometric line")
    if m.water_depth is not None:
        return tr("Depth {0} m").format(_f(m.water_depth, 1))
    return tr("Dry")


def _model_name(key):
    return {"mc": tr("Mohr-Coulomb"), "undrained": tr("Undrained (Su)"),
            "hb": tr("Hoek-Brown"), "none": tr("No strength")}.get(key, key)


def _surface_type_name(key):
    return {"grid": tr("Circular - grid search"), "entry_exit": tr("Circular - entry and exit"),
            "path": tr("Non-circular - path search"), "user": tr("User-defined surfaces"),
            "point": tr("From a crack point picked on the map")}.get(key, key)


def _pf_label(lbl):
    return {"crack water": tr("Water in tension crack"), "reinforcement": tr("Reinforcement"),
            "pile": tr("Pile / wall"), "support": tr("Global support force")}.get(lbl, lbl)
