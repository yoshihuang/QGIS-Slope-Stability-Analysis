# -*- coding: utf-8 -*-
"""Bridge between QGIS layers and the computation core.

Reading (main thread): QGIS layers -> :class:`core.model.SlopeModel`.
Writing: analysis results -> rows -> memory / GeoPackage layers with styles.
"""

import math
import os

from qgis.core import (QgsFeature, QgsGeometry, QgsPointXY, QgsVectorLayer,
                       QgsProject, QgsVectorFileWriter, QgsFields,
                       QgsCoordinateTransformContext, QgsSymbol,
                       QgsGraduatedSymbolRenderer, QgsRendererRange,
                       QgsSingleSymbolRenderer, QgsFeatureRequest,
                       QgsPalLayerSettings, QgsVectorLayerSimpleLabeling,
                       QgsTextFormat, QgsTextBufferSettings, QgsLineSymbol,
                       QgsFillSymbol, QgsMarkerSymbol)
from qgis.PyQt.QtGui import QColor

from shapely import wkt as shp_wkt

from .compat import make_field, enum, geometry_kind
from .core.model import SlopeModel, DistributedLoad, Reinforcement, Pile
from .core.materials import Material
from .core.methods import METHODS, method_short
from .core.i18n import tr

# --------------------------------------------------------------------------
#  attribute names accepted on input layers (case-insensitive)
# --------------------------------------------------------------------------
FIELD_ALIASES = {
    "name": ("name", "layer_name", "material", "mat", "soil", "unit"),
    "gamma": ("unit_weight", "gamma", "uw", "unitweight", "unit_wt"),
    "gamma_sat": ("unit_weight_sat", "gamma_sat", "uw_sat", "sat_weight", "gammasat"),
    "cohesion": ("cohesion", "c", "su", "cu", "c_eff"),
    "phi": ("phi", "friction", "friction_angle", "phi_eff"),
    "model": ("model", "strength_model", "strength", "type"),
    "sigci": ("sigci", "ucs", "sigma_ci"),
    "gsi": ("gsi",),
    "mi": ("mi",),
    "disturbance": ("d", "disturbance", "hb_d"),
    "ru": ("ru",),
    "use_piezo": ("piezo", "use_piezo", "use_water"),
    "c_std": ("c_std", "cohesion_std", "sd_c", "c_sd"),
    "phi_std": ("phi_std", "sd_phi", "phi_sd"),
    "gamma_std": ("gamma_std", "unit_weight_std", "sd_gamma", "gamma_sd"),
    "q": ("q", "load", "surcharge", "kpa", "pressure"),
    "capacity": ("t", "capacity", "tensile", "force", "tmax"),
    "pullout": ("pullout", "bond", "pull_out"),
    "mode": ("mode", "application", "force_mode"),
    "shear": ("v", "shear", "fw", "force", "resistance", "wall_force"),
}

REQUIRED_GEOLOGY = ("gamma", "cohesion", "phi")


def _field_map(layer):
    names = {f.name().lower(): f.name() for f in layer.fields()}
    out = {}
    for key, aliases in FIELD_ALIASES.items():
        for a in aliases:
            if a in names:
                out[key] = names[a]
                break
    return out


def _is_null(v):
    if v is None:
        return True
    try:
        if hasattr(v, "isNull") and v.isNull():
            return True
    except Exception:
        pass
    return isinstance(v, str) and v.strip() in ("", "NULL", "null", "None")


def _num(v, default=None):
    if _is_null(v):
        return default
    try:
        return float(str(v).replace(",", "."))
    except (TypeError, ValueError):
        return default


def _attr(feat, fmap, key, default=None):
    name = fmap.get(key)
    if name is None:
        return default
    v = feat[name]
    return default if _is_null(v) else v


def _shapely(geom):
    g = QgsGeometry(geom)
    try:
        ab = g.get()
        ab.dropZValue()
        ab.dropMValue()
    except Exception:
        pass
    return shp_wkt.loads(g.asWkt())


def _features(layer, selected_only=False):
    if selected_only and layer.selectedFeatureCount() > 0:
        return layer.selectedFeatures()
    return layer.getFeatures()


def _passive(value, default=True):
    if _is_null(value):
        return default
    s = str(value).strip().lower()
    if s in ("active", "a", "0", "false", "act"):
        return False
    if s in ("passive", "p", "1", "true", "pas"):
        return True
    return default


# --------------------------------------------------------------------------
#  input
# --------------------------------------------------------------------------
class InputError(Exception):
    pass


def build_model(layers, opts, selected_only=False):
    """Create a SlopeModel from QGIS layers.

    ``layers``  dict(geology, water, loads, reinforcement, piles, surfaces)
    ``opts``    dict of scalar options (see dialog / processing algorithm)
    Returns (model, user_surfaces, messages)
    """
    msgs = []
    geo = layers.get("geology")
    if geo is None:
        raise InputError(tr("Please select the geology layer."))
    if geometry_kind(geo) != "polygon":
        raise InputError(tr("The geology layer must be a polygon layer."))
    fmap = _field_map(geo)
    missing = [k for k in REQUIRED_GEOLOGY if k not in fmap]
    if missing:
        raise InputError(tr("The geology layer is missing required fields: {0}. Required: unit_weight, cohesion, phi.").format(", ".join(missing)))
    if geo.crs().isValid() and geo.crs().isGeographic():
        msgs.append(("warning", tr("The geology layer uses a geographic CRS. Coordinates are used as section coordinates (x = distance, y = elevation); make sure they are in metres.")))

    m = SlopeModel()
    m.gamma_w = float(opts.get("gamma_w", 9.81))
    mat_index = {}
    for feat in _features(geo, selected_only):
        geom = feat.geometry()
        if geom is None or geom.isEmpty():
            continue
        gamma = _num(_attr(feat, fmap, "gamma"))
        c = _num(_attr(feat, fmap, "cohesion"), 0.0)
        phi = _num(_attr(feat, fmap, "phi"), 0.0)
        if gamma is None:
            raise InputError(tr("Feature {0}: unit weight is empty or not a number.").format(feat.id()))
        name = _attr(feat, fmap, "name", None)
        name = str(name) if name is not None else u"%s %d" % (tr("Material"), feat.id())
        params = dict(
            name=name, gamma=gamma, cohesion=c, phi=phi,
            gamma_sat=_num(_attr(feat, fmap, "gamma_sat"), None),
            model=_attr(feat, fmap, "model", "mc"),
            sigci=_num(_attr(feat, fmap, "sigci"), 0.0),
            gsi=_num(_attr(feat, fmap, "gsi"), 0.0),
            mi=_num(_attr(feat, fmap, "mi"), 0.0),
            disturbance=_num(_attr(feat, fmap, "disturbance"), 0.0),
            ru=_num(_attr(feat, fmap, "ru"), 0.0),
            use_piezo=_passive(_attr(feat, fmap, "use_piezo", None), True),
            c_std=_num(_attr(feat, fmap, "c_std"), 0.0),
            phi_std=_num(_attr(feat, fmap, "phi_std"), 0.0),
            gamma_std=_num(_attr(feat, fmap, "gamma_std"), 0.0),
        )
        key = tuple(sorted((k, str(v)) for k, v in params.items()))
        if key not in mat_index:
            mat_index[key] = m.add_material(Material(**params))
        m.add_region(_shapely(geom), mat_index[key])
    if not m.regions:
        raise InputError(tr("The geology layer contains no valid polygons."))

    # water ----------------------------------------------------------------
    wl = layers.get("water")
    if wl is not None:
        from shapely.ops import unary_union, linemerge
        lines = [_shapely(f.geometry()) for f in wl.getFeatures() if not f.geometry().isEmpty()]
        if lines:
            u = unary_union(lines)
            try:
                u = linemerge(u) if u.geom_type == "MultiLineString" else u
            except Exception:
                pass
            m.piezo = u
            if u.geom_type == "MultiLineString":
                msgs.append(("warning", tr("The water layer contains several lines; they are merged into one piezometric line.")))

    # loads ------------------------------------------------------------------
    ll = layers.get("loads")
    if ll is not None:
        lm = _field_map(ll)
        for f in ll.getFeatures():
            if f.geometry().isEmpty():
                continue
            q = _num(_attr(f, lm, "q"), opts.get("default_q", 0.0))
            if q:
                m.loads.append(DistributedLoad(_shapely(f.geometry()), q))

    # reinforcement ------------------------------------------------------------
    rl = layers.get("reinforcement")
    if rl is not None:
        rm = _field_map(rl)
        for f in rl.getFeatures():
            if f.geometry().isEmpty():
                continue
            g = _shapely(f.geometry())
            T = _num(_attr(f, rm, "capacity"), opts.get("default_T", 0.0))
            po = _num(_attr(f, rm, "pullout"), opts.get("default_pullout", 0.0))
            pas = _passive(_attr(f, rm, "mode", None), opts.get("reinf_passive", True))
            for ln in getattr(g, "geoms", [g]):
                if ln.geom_type == "LineString" and T > 0:
                    m.reinforcements.append(Reinforcement(ln, T, po, pas))

    # piles / walls ---------------------------------------------------------------
    pl = layers.get("piles")
    if pl is not None:
        pm = _field_map(pl)
        for f in pl.getFeatures():
            if f.geometry().isEmpty():
                continue
            V = _num(_attr(f, pm, "shear"), opts.get("default_V", 0.0))
            pas = _passive(_attr(f, pm, "mode", None), opts.get("pile_passive", True))
            if V and V > 0:
                m.piles.append(Pile(_shapely(f.geometry()), V, pas))
        if not m.piles:
            msgs.append(("warning", tr("Pile / wall layer selected but no shear force given (field V or default value).")))
    elif opts.get("global_support", 0.0) > 0:
        m.global_support = float(opts["global_support"])
        m.global_support_passive = bool(opts.get("pile_passive", True))
        msgs.append(("info", tr("A global support force is applied at the toe of every slip surface (v2.x compatible mode).")))

    # user surfaces -------------------------------------------------------------------
    user = []
    sl = layers.get("surfaces")
    if sl is not None:
        for f in sl.getFeatures():
            if f.geometry().isEmpty():
                continue
            g = _shapely(f.geometry())
            for ln in getattr(g, "geoms", [g]):
                if ln.geom_type == "LineString":
                    user.append([(c[0], c[1]) for c in ln.coords])

    m.kh = float(opts.get("kh", 0.0))
    m.kv = float(opts.get("kv", 0.0))
    m.crack_depth = float(opts.get("crack_depth", 0.0))
    m.crack_fill = float(opts.get("crack_fill", 0.0))
    m.direction = int(opts.get("direction", 0))
    m.min_depth = float(opts.get("min_depth", 0.0))
    return m, user, msgs


# --------------------------------------------------------------------------
#  output rows (QGIS independent structures)
# --------------------------------------------------------------------------
SURFACE_FIELDS = [("scenario", "string"), ("rank", "int"), ("fs", "double"), ("method", "string"),
                  ("critical", "int"), ("kind", "string"), ("source", "string"),
                  ("center_x", "double"), ("center_y", "double"), ("radius", "double"),
                  ("entry_x", "double"), ("exit_x", "double"), ("kh", "double"), ("water", "string")]
CRITICAL_FIELDS = ([("scenario", "string"), ("fs", "double"), ("method", "string"), ("kind", "string"),
                    ("center_x", "double"), ("center_y", "double"), ("radius", "double"),
                    ("entry_x", "double"), ("entry_y", "double"), ("exit_x", "double"), ("exit_y", "double"),
                    ("kh", "double"), ("water", "string")]
                   + [("fs_" + method_short(k).lower(), "double") for k in METHODS]
                   + [("lambda_spencer", "double"), ("lambda_mp", "double")])
SLICE_FIELDS = [("scenario", "string"), ("slice", "int"), ("x", "double"), ("width", "double"),
                ("height", "double"), ("alpha", "double"), ("weight", "double"), ("surcharge", "double"),
                ("material", "string"), ("c", "double"), ("phi", "double"), ("u", "double"),
                ("n_eff", "double"), ("sigma_eff", "double"), ("shear_mob", "double")]
CENTER_FIELDS = [("scenario", "string"), ("fs", "double")]
CRACK_FIELDS = [("x", "double"), ("y", "double"), ("min_fs", "double")]


def _water_desc(mdl):
    if mdl._px is not None:
        return "piezometric line"
    if mdl.water_depth is not None:
        return "depth %g m" % mdl.water_depth
    return "dry"


def result_rows(result):
    out = {"surfaces": [], "critical": [], "slices": [], "centers": [], "crack": []}
    st = result.settings
    if st is not None and st.surface_type == "point" and getattr(st, "crack_point", None) and result.model is not None:
        import numpy as np
        mdl0 = result.model
        cx_ = float(st.crack_point[0])
        cy_ = float(mdl0.ground(np.array([mdl0.to_model_x(cx_)]))[0])
        fs_ = [sr.critical.fs for sr in result.scenarios if sr.critical is not None]
        out["crack"].append(((cx_, cy_), {"x": cx_, "y": cy_, "min_fs": round(min(fs_), 4) if fs_ else None}))
    for sr in result.scenarios:
        mdl = sr.model
        name = sr.scenario.name
        water = _water_desc(mdl)
        method = method_short(getattr(sr, "primary_method", ""))
        for rank, r in enumerate(sr.records, 1):
            cw = r.center_world(mdl)
            xe = sorted([mdl.to_world_x(r.x_in), mdl.to_world_x(r.x_out)])
            out["surfaces"].append((r.polyline(mdl), {
                "scenario": name, "rank": rank, "fs": round(r.fs, 4),
                "method": method_short(r.method) if r.method else method,
                "critical": 1 if rank == 1 else 0, "kind": r.surface.kind, "source": r.source,
                "center_x": cw[0] if cw else None, "center_y": cw[1] if cw else None,
                "radius": r.surface.r if r.surface.kind == "circular" else None,
                "entry_x": xe[0], "exit_x": xe[1], "kh": mdl.kh, "water": water}))
        c = sr.critical
        if c is not None:
            pl = c.polyline(mdl)
            cw = c.center_world(mdl)
            attrs = {"scenario": name, "fs": round(c.fs, 4), "method": method, "kind": c.surface.kind,
                     "center_x": cw[0] if cw else None, "center_y": cw[1] if cw else None,
                     "radius": c.surface.r if c.surface.kind == "circular" else None,
                     "entry_x": pl[0][0], "entry_y": pl[0][1], "exit_x": pl[-1][0], "exit_y": pl[-1][1],
                     "kh": mdl.kh, "water": water}
            for k in METHODS:
                r = sr.method_results.get(k)
                attrs["fs_" + method_short(k).lower()] = round(r.fs, 4) if (r is not None and r.ok) else None
            sp = sr.method_results.get("spencer")
            mp = sr.method_results.get("morgenstern_price")
            attrs["lambda_spencer"] = sp.lam if sp is not None and sp.ok else None
            attrs["lambda_mp"] = mp.lam if mp is not None and mp.ok else None
            out["critical"].append((pl, attrs))
            S = sr.slices
            if S is not None and S.valid and sr.slice_rows:
                for i, row in enumerate(sr.slice_rows):
                    xa, ya = mdl.to_world_xy([S.xl[i], S.xr[i], S.xr[i], S.xl[i], S.xl[i]],
                                             [S.ybl[i], S.ybr[i], S.gr[i], S.gl[i], S.ybl[i]])
                    attrs = {"scenario": name}
                    for k, _t in SLICE_FIELDS[1:]:
                        v = row.get(k if k != "n_eff" else "N_eff")
                        attrs[k] = v
                    attrs["x"] = mdl.to_world_x(row["x"])
                    attrs["alpha"] = row["alpha"] * mdl.sign
                    out["slices"].append((list(zip(xa.tolist(), ya.tolist())), attrs))
        for (x, y, fs) in sr.center_grid:
            out["centers"].append(((x, y), {"scenario": name, "fs": round(fs, 4)}))
    return out


def _qgs_geometry(kind, coords):
    if kind == "line":
        return QgsGeometry.fromPolylineXY([QgsPointXY(x, y) for x, y in coords])
    if kind == "polygon":
        return QgsGeometry.fromPolygonXY([[QgsPointXY(x, y) for x, y in coords]])
    return QgsGeometry.fromPointXY(QgsPointXY(coords[0], coords[1]))


def make_fields(schema):
    fields = QgsFields()
    for name, kind in schema:
        fields.append(make_field(name, kind))
    return fields


def _py(v):
    """numpy scalars -> Python, NaN/inf -> NULL (QGIS cannot store them)."""
    if v is None:
        return None
    if hasattr(v, "item"):
        v = v.item()
    if isinstance(v, float) and not math.isfinite(v):
        return None
    return v


def make_feature(fields, schema, kind, coords, attrs):
    f = QgsFeature(fields)
    coords = [(float(x), float(y)) for x, y in coords] if kind != "point" else (float(coords[0]), float(coords[1]))
    f.setGeometry(_qgs_geometry(kind, coords))
    f.setAttributes([_py(attrs.get(n)) for n, _ in schema])
    return f


OUTPUTS = [
    # key, layer name, geometry kind, uri type, schema
    ("surfaces", "slip_surfaces", "line", "LineString", SURFACE_FIELDS),
    ("critical", "critical_surface", "line", "LineString", CRITICAL_FIELDS),
    ("slices", "slices", "polygon", "Polygon", SLICE_FIELDS),
    ("centers", "search_centers", "point", "Point", CENTER_FIELDS),
    ("crack", "crack_point", "point", "Point", CRACK_FIELDS),
]


def memory_layers(rows, crs):
    layers = {}
    for key, lname, kind, uri, schema in OUTPUTS:
        if not rows.get(key):
            continue
        authid = crs.authid() if crs is not None and crs.isValid() else ""
        vl = QgsVectorLayer("%s?crs=%s" % (uri, authid) if authid else uri, lname, "memory")
        if crs is not None and crs.isValid() and not authid:
            vl.setCrs(crs)
        pr = vl.dataProvider()
        fields = make_fields(schema)
        pr.addAttributes(list(fields))
        vl.updateFields()
        feats = [make_feature(vl.fields(), schema, kind, c, a) for c, a in rows[key]]
        ok = pr.addFeatures(feats)
        if isinstance(ok, tuple):
            ok = ok[0]
        if not ok:
            from qgis.core import QgsMessageLog
            err = pr.lastError() if hasattr(pr, "lastError") else ""
            QgsMessageLog.logMessage(u"%s: %s" % (lname, err), "Slope Stability")
        vl.updateExtents()
        layers[key] = vl
    return layers


def save_gpkg(layers, path):
    """Write memory layers to one GeoPackage; return the re-loaded layers."""
    ctx = QgsProject.instance().transformContext() if QgsProject.instance() else QgsCoordinateTransformContext()
    first = True
    out = {}
    for key, lname, kind, uri, schema in OUTPUTS:
        vl = layers.get(key)
        if vl is None:
            continue
        o = QgsVectorFileWriter.SaveVectorOptions()
        o.driverName = "GPKG"
        o.layerName = lname
        o.fileEncoding = "UTF-8"
        if first and os.path.exists(path):
            o.actionOnExistingFile = enum(QgsVectorFileWriter, "ActionOnExistingFile", "CreateOrOverwriteFile")
        elif not first:
            o.actionOnExistingFile = enum(QgsVectorFileWriter, "ActionOnExistingFile", "CreateOrOverwriteLayer")
        res = QgsVectorFileWriter.writeAsVectorFormatV3(vl, path, ctx, o)
        err = res[0] if isinstance(res, tuple) else res
        if err != enum(QgsVectorFileWriter, "WriterError", "NoError"):
            raise IOError(tr("Could not write {0}: {1}").format(path, res[1] if isinstance(res, tuple) else err))
        first = False
        nl = QgsVectorLayer("%s|layername=%s" % (path, lname), lname, "ogr")
        out[key] = nl
    return out


# --------------------------------------------------------------------------
#  styling
# --------------------------------------------------------------------------
FS_CLASSES = [
    (0.0, 1.0, "#c0392b", "FS < 1.0"),
    (1.0, 1.2, "#e67e22", "1.0 - 1.2"),
    (1.2, 1.5, "#f1c40f", "1.2 - 1.5"),
    (1.5, 2.0, "#7fb800", "1.5 - 2.0"),
    (2.0, 1000.0, "#2e8b57", "FS > 2.0"),
]


def _geom_symbol(layer):
    try:
        return QgsSymbol.defaultSymbol(layer.geometryType())
    except Exception:
        return None


def _graduated(layer, width=0.3, point_size=None):
    ranges = []
    for lo, hi, col, label in FS_CLASSES:
        sym = _geom_symbol(layer)
        sym.setColor(QColor(col))
        if point_size is not None:
            sym.setSize(point_size)
        else:
            try:
                sym.setWidth(width)
            except Exception:
                pass
        ranges.append(QgsRendererRange(lo, hi, sym, label))
    r = QgsGraduatedSymbolRenderer("fs", ranges)
    try:
        r.setOrderBy(QgsFeatureRequest.OrderBy([QgsFeatureRequest.OrderByClause("fs", False)]))
        r.setOrderByEnabled(True)
    except Exception:
        pass
    layer.setRenderer(r)


def _fs_labels(layer, expr, size=11, color="#c0392b", line=True):
    s = QgsPalLayerSettings()
    s.fieldName = expr
    s.isExpression = True
    fmt = QgsTextFormat()
    fmt.setSize(size)
    fmt.setColor(QColor(color))
    buf = QgsTextBufferSettings()
    buf.setEnabled(True)
    buf.setSize(1.0)
    buf.setColor(QColor("white"))
    fmt.setBuffer(buf)
    s.setFormat(fmt)
    if line:
        try:
            s.placement = enum(Qgis, "LabelPlacement", "Line")
        except Exception:
            pass
    layer.setLabeling(QgsVectorLayerSimpleLabeling(s))
    layer.setLabelsEnabled(True)


def style_layers(layers):
    if "surfaces" in layers:
        _graduated(layers["surfaces"], 0.3)
    if "critical" in layers:
        vl = layers["critical"]
        sym = QgsLineSymbol.createSimple({"color": "#c0392b", "width": "1.2"})
        vl.setRenderer(QgsSingleSymbolRenderer(sym))
        _fs_labels(vl, "'FS = ' || format_number(\"fs\", 3) || ' (' || \"scenario\" || ')'")
    if "slices" in layers:
        vl = layers["slices"]
        sym = QgsFillSymbol.createSimple({"color": "255,0,0,0", "outline_color": "#c0392b",
                                          "outline_width": "0.15", "style": "no"})
        vl.setRenderer(QgsSingleSymbolRenderer(sym))
    if "crack" in layers:
        vl = layers["crack"]
        sym = QgsMarkerSymbol.createSimple({"name": "triangle", "color": "#c0392b", "size": "4",
                                            "outline_color": "white", "angle": "180"})
        vl.setRenderer(QgsSingleSymbolRenderer(sym))
    if "centers" in layers:
        _graduated(layers["centers"], point_size=1.6)
    for vl in layers.values():
        vl.triggerRepaint()
