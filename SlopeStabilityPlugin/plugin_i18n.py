# -*- coding: utf-8 -*-
"""Multi-language support.

The interface language follows the QGIS user interface language
(Settings > Options > General > User interface translation).  It can be
overridden in the plugin dialog (Output tab > Language).

Translations are standard Qt catalogues in ``i18n/SlopeStability_<lang>.qm``
(sources: ``i18n/*.ts``).  All strings use the Qt context ``SlopeStability``.
"""

import os

from qgis.PyQt.QtCore import QCoreApplication, QTranslator, QLocale
from qgis.core import QgsSettings, QgsApplication

from .core import i18n as core_i18n

CONTEXT = "SlopeStability"
SETTINGS_KEY = "SlopeStability/language"
PLUGIN_DIR = os.path.dirname(__file__)

# (code, native name) - "auto" follows QGIS
LANGUAGES = [
    ("auto", None),
    ("en", "English"),
    ("zh_TW", u"繁體中文"),
    ("zh_CN", u"简体中文"),
    ("ja", u"日本語"),
    ("ko", u"한국어"),
    ("es", u"Español"),
    ("fr", u"Français"),
    ("de", u"Deutsch"),
    ("pt", u"Português"),
]
AVAILABLE = [c for c, _ in LANGUAGES if c not in ("auto", "en")]

_translator = None
_active = "en"


def tr(text):
    return core_i18n.tr(text)


def _qgis_locale():
    loc = ""
    try:
        loc = QgsApplication.locale() or ""
    except Exception:
        loc = ""
    if not loc:
        s = QgsSettings()
        override = s.value("locale/overrideFlag", False, type=bool)
        loc = s.value("locale/userLocale", "", type=str) if override else ""
        if not loc:
            loc = s.value("locale/userLocale", "", type=str) or QLocale.system().name()
    return loc or "en"


def resolve_language(code):
    """Map a QGIS / Qt locale name to one of the shipped catalogues."""
    if not code:
        return "en"
    c = code.replace("-", "_")
    low = c.lower()
    if low.startswith("zh"):
        if any(t in low for t in ("hant", "_tw", "_hk", "_mo")):
            return "zh_TW"
        return "zh_CN"
    if c in AVAILABLE:
        return c
    short = low.split("_")[0]
    for a in AVAILABLE:
        if a.lower() == short:
            return a
    return "en"


def configured_language():
    return QgsSettings().value(SETTINGS_KEY, "auto", type=str) or "auto"


def set_configured_language(code):
    QgsSettings().setValue(SETTINGS_KEY, code or "auto")


def active_language():
    return _active


def install(code=None):
    """(Re)install the translator. ``code`` None -> configured / QGIS language."""
    global _translator, _active
    if code is None:
        code = configured_language()
    lang = resolve_language(_qgis_locale() if code == "auto" else code)
    app = QCoreApplication.instance()
    if _translator is not None and app is not None:
        app.removeTranslator(_translator)
        _translator = None
    _active = lang
    if lang != "en":
        qm = os.path.join(PLUGIN_DIR, "i18n", "SlopeStability_%s.qm" % lang)
        if os.path.exists(qm):
            t = QTranslator()
            if t.load(qm) and app is not None:
                app.installTranslator(t)
                _translator = t
    core_i18n.set_translator(lambda s: QCoreApplication.translate(CONTEXT, s))
    return lang


def uninstall():
    global _translator
    app = QCoreApplication.instance()
    if _translator is not None and app is not None:
        app.removeTranslator(_translator)
    _translator = None
    core_i18n.set_translator(None)
