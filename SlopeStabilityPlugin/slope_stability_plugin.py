import os
from qgis.PyQt.QtWidgets import QAction, QMessageBox
from qgis.PyQt.QtCore import QThread, pyqtSignal
from qgis.core import QgsVectorLayer, QgsProject
from .slope_stability_dialog import SlopeStabilityDialog
from .slope_stability_logic import SlopeStabilityWorker

class SlopeStabilityPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.dlg = None
        self.worker = None

    def initGui(self):
        self.action = QAction("Slope Stability Analysis (Bishop)", self.iface.mainWindow())
        self.action.triggered.connect(self.run)
        self.iface.addPluginToMenu("&Slope Stability", self.action)
        self.iface.addToolBarIcon(self.action)

    def unload(self):
        self.iface.removePluginMenu("&Slope Stability", self.action)
        self.iface.removeToolBarIcon(self.action)

    def run(self):
        if not self.dlg:
            self.dlg = SlopeStabilityDialog()
        
        # Refresh layer list
        self.dlg.cb_geology.clear()
        self.dlg.populate_layers(self.dlg.cb_geology, 0)
        self.dlg.cb_water.clear()
        self.dlg.cb_water.addItem("None (Use Depth Parameter)", None)
        self.dlg.populate_layers(self.dlg.cb_water, 1)
        
        self.dlg.run_button.clicked.connect(self.start_analysis)
        self.dlg.show()
        self.dlg.exec_()

    def start_analysis(self):
        # 1. Get parameters
        geo_layer = self.dlg.cb_geology.currentData()
        water_layer = self.dlg.cb_water.currentData()
        kh = self.dlg.sb_kh.value()
        output_path = self.dlg.le_output.text()
        
        mode = "both"
        if self.dlg.rb_general.isChecked(): mode = "general"
        elif self.dlg.rb_event.isChecked(): mode = "event"

        # Validation
        if not geo_layer:
            self.iface.messageBar().pushMessage("Error", "Please select a geological model layer.", level=3)
            return
        if not output_path:
            self.iface.messageBar().pushMessage("Error", "Please specify the output file path.", level=3)
            return

        # 2. Start background thread (avoid freezing QGIS)
        self.dlg.run_button.setEnabled(False)
        self.dlg.progress_bar.setValue(10)
        
        self.worker = SlopeStabilityWorker(geo_layer, water_layer, kh, mode, output_path)
        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.analysis_finished)
        self.worker.error.connect(self.analysis_error)
        self.worker.start()

    def update_progress(self, val, msg):
        self.dlg.progress_bar.setValue(val)
        self.iface.messageBar().pushMessage("Analyzing", msg, level=0)

    def analysis_finished(self, output_path):
        self.dlg.progress_bar.setValue(100)
        self.dlg.run_button.setEnabled(True)
        QMessageBox.information(self.dlg, "Finished", f"Analysis completed!\nFile saved to: {output_path}")
        
        # Automatically load result into QGIS
        layer = QgsVectorLayer(output_path, "Analysis Results (Slip Surfaces)", "ogr")
        if layer.isValid():
            QgsProject.instance().addMapLayer(layer)

    def analysis_error(self, msg):
        self.dlg.run_button.setEnabled(True)
        self.dlg.progress_bar.setValue(0)
        QMessageBox.critical(self.dlg, "Error", str(msg))