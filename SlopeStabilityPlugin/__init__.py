    # 在 slope_stability_plugin.py 或 __init__.py 開頭加入

def classFactory(iface):
    try:
        import numpy
        import shapely
    except ImportError as e:
        from qgis.PyQt.QtWidgets import QMessageBox
        QMessageBox.critical(
            None, 
            "Missing Dependencies",
            f"Your QGIS Python environment is missing required libraries to run this plugin.\n\nError message: {str(e)}\n\nPlease ensure you are using a standard QGIS installation."
        )
        # 這裡可以回傳一個空物件或 Dummy，防止 QGIS 載入錯誤
        # 但通常直接拋出錯誤讓使用者知道會比較好
        raise e

    from .slope_stability_plugin import SlopeStabilityPlugin
    return SlopeStabilityPlugin(iface)