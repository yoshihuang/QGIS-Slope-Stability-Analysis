# -*- coding: utf-8 -*-
"""Background task running the analysis (QgsTask)."""

import traceback

from qgis.core import QgsTask, QgsMessageLog
from qgis.PyQt.QtCore import pyqtSignal

from .core.analysis import run_analysis
from .core.search import Cancelled
from .core.i18n import tr

LOG_TAG = "Slope Stability"


class AnalysisTask(QgsTask):
    message = pyqtSignal(int, str)

    def __init__(self, model, settings, scenarios, methods, reliability, on_done):
        super().__init__(tr("Slope stability analysis"))  # default flags allow cancelling
        self.model = model
        self.settings = settings
        self.scenarios = scenarios
        self.methods = methods
        self.reliability = reliability
        self.on_done = on_done
        self.result = None
        self.error = None
        self.cancelled_by_user = False

    def run(self):
        try:
            def progress(f, msg=""):
                p = int(round(100 * f))
                self.setProgress(p)
                self.message.emit(p, msg)

            def log(m):
                QgsMessageLog.logMessage(m, LOG_TAG)

            self.result = run_analysis(self.model, self.settings, self.scenarios, self.methods,
                                       self.reliability, progress, self.isCanceled, log)
            return True
        except Cancelled:
            self.cancelled_by_user = True
            return False
        except Exception as e:  # report every error to the user
            self.error = u"%s\n\n%s" % (e, traceback.format_exc())
            return False

    def finished(self, ok):
        try:
            self.on_done(self, ok)
        except Exception:
            QgsMessageLog.logMessage(traceback.format_exc(), LOG_TAG)
