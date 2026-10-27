# QGIS Slope Stability Analysis 3.1

2-D limit equilibrium slope stability analysis for QGIS 3.28+ and QGIS 4, with a feature set modelled on [XSLOPE](https://xslope.org) (LEM part).

* **Methods**: Ordinary (Fellenius), Bishop simplified, Janbu simplified (f0), Corps of Engineers #1/#2, Lowe-Karafiath, Spencer, Morgenstern-Price
* **Search**: circular grid, circular entry-exit, non-circular path search, Nelder-Mead / vertex optimisation, user-defined surfaces
* **Crack point analysis**: click a crack / scarp on the ground in the map; circular and non-circular surfaces starting there are analysed
* **Materials**: Mohr-Coulomb, undrained Su, Hoek-Brown; saturated unit weight; ru
* **Loads and supports**: piezometric line, ponded water, surcharge, tension crack with water, kh/kv, reinforcement with pull-out, piles/walls (active or passive)
* **Scenarios**: parametric kh × water-depth scan (the v2 *General Scan*)
* **Reliability**: Monte Carlo and Taylor series (TSPM)
* **Outputs**: styled GeoPackage layers, HTML report with cross-section, Processing algorithm
* **Languages**: English, 繁體中文, 简体中文, 日本語, 한국어, Español, Français, Deutsch, Português – follows the QGIS locale

See `docs/user_guide_en.html` (English) or `docs/user_guide_zh_TW.html` (繁體中文).

## Layout

```
core/                 QGIS-independent engine (numpy + shapely)
  model.py            geometry, materials, water, loads, supports
  slices.py           slip surfaces and slicing
  methods.py          the 8 limit-equilibrium methods (one formulation)
  search.py           critical surface searches and optimisation
  reliability.py      Monte Carlo / TSPM
  analysis.py         scenario driver
  report.py           HTML report with SVG figures
qgis_io.py            QGIS layers <-> core, output layers and styles
ui/                   dialogs
processing_provider/  Processing algorithm
plugin_i18n.py        language selection (QTranslator)
i18n/                 .ts sources and compiled .qm catalogues
tools/                i18n_build.py + translations/*.json
tests/                benchmark tests (python tests/test_benchmarks.py)
```

## Translations

Every user-visible string is written as `tr("English text")`. To add or update a language:

1. edit `tools/translations/<lang>.json` (`{"English source": "translation"}`),
2. run `python tools/i18n_build.py` (needs `lrelease` from Qt Linguist or `pyside6-lrelease`),
3. `python tools/i18n_build.py --check` lists missing strings.

To add a new language, add its code to `LANGS` in `tools/i18n_build.py` and to `LANGUAGES` in `plugin_i18n.py`.

## Verification

ACADS 1(a): Bishop 0.985, Spencer 0.984 (Slide 0.987 / 0.986, referee ≈ 1.00); infinite slope dry / seismic / seepage within 1–2 % (end effects); Taylor undrained vertical cut 0.963 vs 0.958.

## License

GNU GPL v2 or later (see LICENSE).
