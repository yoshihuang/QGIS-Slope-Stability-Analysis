# QGIS-Slope-Stability-Analysis

QGIS plugin for 2-D limit equilibrium slope stability analysis (version 3.1). The plugin source is in [`SlopeStabilityPlugin/`](SlopeStabilityPlugin/); see [`SlopeStabilityPlugin/README.md`](SlopeStabilityPlugin/README.md) for the developer notes.

QGIS 二維極限平衡邊坡穩定分析外掛（3.1 版）。外掛原始碼位於 [`SlopeStabilityPlugin/`](SlopeStabilityPlugin/)。

## Features / 功能

* 8 methods of slices / 8 種切片法: Ordinary (Fellenius), Bishop, Janbu (f0), Corps of Engineers #1/#2, Lowe-Karafiath, Spencer, Morgenstern-Price
* Critical surface search / 臨界滑動面搜尋: circular grid, entry-exit, non-circular path search, optimisation, user-defined surfaces
* **Crack point analysis / 裂隙點分析**: pick a crack or scarp on the ground in the map; circular and non-circular slip surfaces starting there are analysed / 在地圖上點選地表裂隙，由該點向下延伸滑動面並計算安全係數
* Materials / 材料: Mohr-Coulomb, undrained Su, Hoek-Brown
* Water, surcharge, tension crack, seismic kh/kv, reinforcement, piles / retaining walls / 地下水、載重、張力裂縫、地震、加勁材、樁與擋土牆
* Parametric kh × water depth scan (v2 "General Scan"), Monte Carlo and Taylor series reliability / 參數掃描、可靠度分析
* Styled GeoPackage output, HTML report, Processing algorithm / 成果圖層、HTML 報告、處理工具箱
* 9 languages following the QGIS locale / 9 種語言，依 QGIS 介面語言切換: English, 繁體中文, 简体中文, 日本語, 한국어, Español, Français, Deutsch, Português
* QGIS 3.28 – 4.x

## Quick start / 快速開始

1. Install from the QGIS plugin manager (*QGIS Slope Stability Analysis*) or zip the `SlopeStabilityPlugin` folder and use *Install from ZIP*.
2. *Plugins ▸ Slope Stability ▸ Load sample data* / 外掛 ▸ 邊坡穩定分析 ▸ 載入範例資料
3. *Slope stability analysis…* or the toolbar button *Pick a crack point and analyse…*

User guide / 使用說明: `SlopeStabilityPlugin/docs/user_guide_en.html`, `SlopeStabilityPlugin/docs/user_guide_zh_TW.html`.
The PDF at the repository root is the manual of version 1.x / 根目錄的 PDF 為 1.x 版手冊。

## History / 開發紀錄

Versions 1.x–2.x were developed with the assistance of Gemini; version 3 was developed with the assistance of Claude.
1.x–2.x 版由 Gemini 協助開發；3 版由 Claude 協助開發。

License: GPL v2 or later.
