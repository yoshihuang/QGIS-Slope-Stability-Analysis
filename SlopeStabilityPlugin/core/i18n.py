# -*- coding: utf-8 -*-
"""Translation hook for the computation core.

The core package is independent of QGIS/Qt so it can be unit-tested and reused.
The QGIS plugin installs a translator function at start-up (see ``plugin_i18n``),
after which every ``tr("...")`` call in the core returns the localised string.

All translatable strings must be passed to ``tr`` as *literal* strings so that
``tools/i18n_build.py`` can extract them into the Qt ``.ts`` catalogues.
"""

_translate = None


def set_translator(func):
    """Install a callable ``func(text) -> str`` used by :func:`tr`."""
    global _translate
    _translate = func


def tr(text):
    if _translate is None:
        return text
    try:
        out = _translate(text)
        return out if out else text
    except Exception:  # pragma: no cover - never break a run on i18n
        return text
