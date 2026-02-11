import os
from qgis.PyQt import uic
from qgis.PyQt import QtWidgets
from qgis.core import QgsMapLayerProxyModel, QgsProject

class SlopeStabilityDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("邊坡穩定分析系統 (QGIS Pro)")
        self.resize(600, 500)
        
        layout = QtWidgets.QVBoxLayout(self)
        self.tab_widget = QtWidgets.QTabWidget()
        layout.addWidget(self.tab_widget)
        
        # --- Tab 1: 分析設定 ---
        self.tab_analysis = QtWidgets.QWidget()
        self.init_analysis_tab()
        self.tab_widget.addTab(self.tab_analysis, "分析設定")
        
        # --- Tab 2: 說明手冊 ---
        self.tab_help = QtWidgets.QWidget()
        self.init_help_tab()
        self.tab_widget.addTab(self.tab_help, "使用說明與限制")

        # --- 底部按鈕 ---
        btn_layout = QtWidgets.QHBoxLayout()
        self.run_button = QtWidgets.QPushButton("開始分析")
        self.close_button = QtWidgets.QPushButton("關閉")
        self.close_button.clicked.connect(self.close)
        
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setValue(0)
        
        btn_layout.addWidget(self.progress_bar)
        btn_layout.addWidget(self.run_button)
        btn_layout.addWidget(self.close_button)
        layout.addLayout(btn_layout)

    def init_analysis_tab(self):
        layout = QtWidgets.QFormLayout(self.tab_analysis)
        
        # 1. 地質模型圖層
        self.cb_geology = QtWidgets.QComboBox()
        self.populate_layers(self.cb_geology, layer_type=0) # Polygon
        layout.addRow("地質模型圖層 (Polygon):", self.cb_geology)
        
        # 2. 地下水位圖層
        self.cb_water = QtWidgets.QComboBox()
        self.cb_water.addItem("無 (使用深度參數)", None)
        self.populate_layers(self.cb_water, layer_type=1) # Line
        layout.addRow("地下水位圖層 (Line) [事件用]:", self.cb_water)
        
        # 3. 地震係數
        self.sb_kh = QtWidgets.QDoubleSpinBox()
        self.sb_kh.setRange(0.0, 1.0)
        self.sb_kh.setSingleStep(0.01)
        self.sb_kh.setValue(0.15)
        layout.addRow("事件地震係數 (Kh):", self.sb_kh)
        
        # 4. 分析模式
        self.rb_general = QtWidgets.QRadioButton("僅廣域掃描 (General Scan)")
        self.rb_event = QtWidgets.QRadioButton("僅事件分析 (Event Specific)")
        self.rb_both = QtWidgets.QRadioButton("全部執行 (Both)")
        self.rb_both.setChecked(True)
        
        mode_layout = QtWidgets.QVBoxLayout()
        mode_layout.addWidget(self.rb_general)
        mode_layout.addWidget(self.rb_event)
        mode_layout.addWidget(self.rb_both)
        layout.addRow("分析模式:", mode_layout)
        
        # 5. 輸出路徑
        self.le_output = QtWidgets.QLineEdit()
        self.btn_browse = QtWidgets.QPushButton("瀏覽...")
        self.btn_browse.clicked.connect(self.select_output_file)
        
        out_layout = QtWidgets.QHBoxLayout()
        out_layout.addWidget(self.le_output)
        out_layout.addWidget(self.btn_browse)
        layout.addRow("輸出破壞面路徑 (.geojson):", out_layout)

    def init_help_tab(self):
        layout = QtWidgets.QVBoxLayout(self.tab_help)
        text_edit = QtWidgets.QTextEdit()
        text_edit.setReadOnly(True)
        manual_text = """
        <h3>邊坡穩定分析系統 - 使用說明</h3>
        <hr>
        <b>1. 參數與單位需求</b><br>
        地質圖層 (Polygon) 必須包含以下屬性欄位：
        <ul>
            <li><b>cohesion</b> (凝聚力): kPa (kN/m²)</li>
            <li><b>phi</b> (內摩擦角): 度 (Degree)</li>
            <li><b>unit_weight</b> (單位重): kN/m³</li>
        </ul>
        
        <b>2. 功能模式說明</b><br>
        <ul>
            <li><b>廣域掃描</b>: 自動對 Kh=0.0~0.2 及深度 0~30m 進行矩陣分析。</li>
            <li><b>事件分析</b>: 針對指定的 Kh 與特定的「地下水位圖層」進行精確分析。</li>
        </ul>
        
        <b>3. 使用限制</b><br>
        <ul>
            <li>座標系統建議使用投影座標 (如 TWD97 / EPSG:3826)，單位為公尺。</li>
            <li>Bishop 方法假設滑動面為圓弧形。</li>
            <li>若未選擇水位圖層，事件模式將無法執行或改用深度參數。</li>
        </ul>
        """
        text_edit.setHtml(manual_text)
        layout.addWidget(text_edit)

    def populate_layers(self, combobox, layer_type):
        # 0=Polygon, 1=Line
        layers = QgsProject.instance().mapLayers().values()
        for layer in layers:
            if layer_type == 0 and layer.geometryType() == 2: # Polygon
                combobox.addItem(layer.name(), layer)
            elif layer_type == 1 and layer.geometryType() == 1: # Line
                combobox.addItem(layer.name(), layer)

    def select_output_file(self):
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "選擇輸出檔案", "", "GeoJSON Files (*.geojson)")
        if filename:
            self.le_output.setText(filename)