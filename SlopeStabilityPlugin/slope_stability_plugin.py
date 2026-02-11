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
        self.action = QAction("邊坡穩定分析 (Bishop)", self.iface.mainWindow())
        self.action.triggered.connect(self.run)
        self.iface.addPluginToMenu("&Slope Stability", self.action)
        self.iface.addToolBarIcon(self.action)

    def unload(self):
        self.iface.removePluginMenu("&Slope Stability", self.action)
        self.iface.removeToolBarIcon(self.action)

    def run(self):
        if not self.dlg:
            self.dlg = SlopeStabilityDialog()
        
        # 重新整理圖層清單
        self.dlg.cb_geology.clear()
        self.dlg.populate_layers(self.dlg.cb_geology, 0)
        self.dlg.cb_water.clear()
        self.dlg.cb_water.addItem("無 (使用深度參數)", None)
        self.dlg.populate_layers(self.dlg.cb_water, 1)
        
        self.dlg.run_button.clicked.connect(self.start_analysis)
        self.dlg.show()
        self.dlg.exec_()

    def start_analysis(self):
        # 1. 獲取參數
        geo_layer = self.dlg.cb_geology.currentData()
        water_layer = self.dlg.cb_water.currentData()
        kh = self.dlg.sb_kh.value()
        output_path = self.dlg.le_output.text()
        
        mode = "both"
        if self.dlg.rb_general.isChecked(): mode = "general"
        elif self.dlg.rb_event.isChecked(): mode = "event"

        # 驗證
        if not geo_layer:
            self.iface.messageBar().pushMessage("錯誤", "請選擇地質模型圖層", level=3)
            return
        if not output_path:
            self.iface.messageBar().pushMessage("錯誤", "請指定輸出檔案路徑", level=3)
            return

        # 2. 啟動後台執行緒 (避免卡死 QGIS)
        self.dlg.run_button.setEnabled(False)
        self.dlg.progress_bar.setValue(10)
        
        self.worker = SlopeStabilityWorker(geo_layer, water_layer, kh, mode, output_path)
        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.analysis_finished)
        self.worker.error.connect(self.analysis_error)
        self.worker.start()

    def update_progress(self, val, msg):
        self.dlg.progress_bar.setValue(val)
        self.iface.messageBar().pushMessage("分析中", msg, level=0)

    def analysis_finished(self, output_path):
        self.dlg.progress_bar.setValue(100)
        self.dlg.run_button.setEnabled(True)
        QMessageBox.information(self.dlg, "完成", f"分析完成！\n檔案已儲存至: {output_path}")
        
        # 自動載入結果到 QGIS
        layer = QgsVectorLayer(output_path, "分析結果 (破壞面)", "ogr")
        if layer.isValid():
            QgsProject.instance().addMapLayer(layer)

    def analysis_error(self, msg):
        self.dlg.run_button.setEnabled(True)
        self.dlg.progress_bar.setValue(0)
        QMessageBox.critical(self.dlg, "錯誤", str(msg))