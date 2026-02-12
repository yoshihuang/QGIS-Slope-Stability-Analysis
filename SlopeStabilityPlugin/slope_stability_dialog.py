import os
from qgis.PyQt import uic
from qgis.PyQt import QtWidgets
from qgis.core import QgsMapLayerProxyModel, QgsProject

class SlopeStabilityDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Slope Stability Analysis (QGIS Pro)")
        self.resize(600, 500)
        
        layout = QtWidgets.QVBoxLayout(self)
        self.tab_widget = QtWidgets.QTabWidget()
        layout.addWidget(self.tab_widget)
        
        # --- Tab 1: Analysis Settings ---
        self.tab_analysis = QtWidgets.QWidget()
        self.init_analysis_tab()
        self.tab_widget.addTab(self.tab_analysis, "Analysis Settings")
        
        # --- Tab 2: Help & Manual ---
        self.tab_help = QtWidgets.QWidget()
        self.init_help_tab()
        self.tab_widget.addTab(self.tab_help, "Help & Limitations")

        # --- Bottom Buttons ---
        btn_layout = QtWidgets.QHBoxLayout()
        self.run_button = QtWidgets.QPushButton("Run Analysis")
        self.close_button = QtWidgets.QPushButton("Close")
        self.close_button.clicked.connect(self.close)
        
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setValue(0)
        
        btn_layout.addWidget(self.progress_bar)
        btn_layout.addWidget(self.run_button)
        btn_layout.addWidget(self.close_button)
        layout.addLayout(btn_layout)

    def init_analysis_tab(self):
        layout = QtWidgets.QFormLayout(self.tab_analysis)
        
        # 1. Geology Model Layer
        self.cb_geology = QtWidgets.QComboBox()
        self.populate_layers(self.cb_geology, layer_type=0) # Polygon
        layout.addRow("Geology Model Layer (Polygon):", self.cb_geology)
        
        # 2. Water Table Layer
        self.cb_water = QtWidgets.QComboBox()
        self.cb_water.addItem("None (Use Depth Parameter)", None)
        self.populate_layers(self.cb_water, layer_type=1) # Line
        layout.addRow("Water Table Layer (Line) [For Event Mode]:", self.cb_water)
        
        # 3. Seismic Coefficient
        self.sb_kh = QtWidgets.QDoubleSpinBox()
        self.sb_kh.setRange(0.0, 1.0)
        self.sb_kh.setSingleStep(0.01)
        self.sb_kh.setValue(0.15)
        layout.addRow("Seismic Coefficient (Kh):", self.sb_kh)
        
        # 4. Analysis Mode
        self.rb_general = QtWidgets.QRadioButton("General Scan Only")
        self.rb_event = QtWidgets.QRadioButton("Event Specific Only")
        self.rb_both = QtWidgets.QRadioButton("Run Both Modes")
        self.rb_both.setChecked(True)
        
        mode_layout = QtWidgets.QVBoxLayout()
        mode_layout.addWidget(self.rb_general)
        mode_layout.addWidget(self.rb_event)
        mode_layout.addWidget(self.rb_both)
        layout.addRow("Analysis Mode:", mode_layout)
        
        # 5. Output Path
        self.le_output = QtWidgets.QLineEdit()
        self.btn_browse = QtWidgets.QPushButton("Browse...")
        self.btn_browse.clicked.connect(self.select_output_file)
        
        out_layout = QtWidgets.QHBoxLayout()
        out_layout.addWidget(self.le_output)
        out_layout.addWidget(self.btn_browse)
        layout.addRow("Output Path (.geojson):", out_layout)

    def init_help_tab(self):
        layout = QtWidgets.QVBoxLayout(self.tab_help)
        text_edit = QtWidgets.QTextEdit()
        text_edit.setReadOnly(True)
        manual_text = """
        <h3>Slope Stability Analysis System - User Manual</h3>
        <hr>
        <b>1. Parameters and Units Requirements</b><br>
        The Geology Layer (Polygon) must contain the following attribute fields:
        <ul>
            <li><b>cohesion</b> (Cohesion): kPa (kN/m²)</li>
            <li><b>phi</b> (Friction Angle): Degree (°)</li>
            <li><b>unit_weight</b> (Unit Weight): kN/m³</li>
        </ul>
        
        <b>2. Function Modes</b><br>
        <ul>
            <li><b>General Scan</b>: Automatically performs a matrix analysis for Kh=0.0~0.2 and Depth 0~30m.</li>
            <li><b>Event Specific</b>: Performs a precise analysis for a specific Kh and a specific "Water Table Layer".</li>
        </ul>
        
        <b>3. Limitations</b><br>
        <ul>
            <li>Coordinate System: Projected CRS (e.g., EPSG:3826) with units in Meters is recommended.</li>
            <li>Method: Bishop's Simplified Method (Circular Slip Surface).</li>
            <li>If no water table layer is selected, Event Mode will fail or default to using depth parameters.</li>
            <li>Source Code,test data,detail information:https://github.com/yoshihuang/QGIS-Slope-Stability-Analysis</li>
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
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Select Output File", "", "GeoJSON Files (*.geojson)")
        if filename:
            self.le_output.setText(filename)