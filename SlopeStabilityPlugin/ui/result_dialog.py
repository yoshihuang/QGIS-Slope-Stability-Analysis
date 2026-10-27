# -*- coding: utf-8 -*-
"""Summary of the analysis results."""

from qgis.PyQt.QtCore import QUrl
from qgis.PyQt.QtGui import QDesktopServices
from qgis.PyQt.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QTextBrowser, QPushButton

from ..core.methods import METHODS, method_label
from ..core.report import fs_color
from ..core.i18n import tr


class ResultDialog(QDialog):
    def __init__(self, parent, result, report_path=None):
        super().__init__(parent)
        self.setWindowTitle(tr("Slope stability results"))
        self.resize(720, 420)
        self.report_path = report_path
        lay = QVBoxLayout(self)
        tb = QTextBrowser()
        tb.setHtml(self._html(result))
        lay.addWidget(tb)
        row = QHBoxLayout()
        row.addStretch()
        if report_path:
            b = QPushButton(tr("Open report"))
            b.clicked.connect(lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(report_path)))
            row.addWidget(b)
        c = QPushButton(tr("Close"))
        c.clicked.connect(self.close)
        row.addWidget(c)
        lay.addLayout(row)

    def _html(self, res):
        methods = [m for m in METHODS if m in res.methods]
        h = [u"<h3>%s</h3>" % tr("Minimum factor of safety")]
        h.append(u"<table border='1' cellspacing='0' cellpadding='4'><tr><th>%s</th>" % tr("Scenario"))
        for m in methods:
            h.append(u"<th>%s</th>" % method_label(m))
        h.append(u"</tr>")
        for sr in res.scenarios:
            h.append(u"<tr><td>%s</td>" % sr.scenario.label)
            if sr.critical is None:
                h.append(u"<td colspan='%d' style='color:#c0392b'>%s</td></tr>" % (len(methods), sr.message))
                continue
            for m in methods:
                r = sr.method_results.get(m)
                if r is not None and r.ok:
                    h.append(u"<td align='right' style='color:%s'><b>%.3f</b></td>" % (fs_color(r.fs), r.fs))
                else:
                    h.append(u"<td align='center'>-</td>")
            h.append(u"</tr>")
        h.append(u"</table>")
        rel = res.reliability
        if rel and not rel.get("error"):
            h.append(u"<h3>%s</h3>" % tr("Probabilistic analysis"))
            mc = rel.get("mc")
            if mc:
                h.append(u"<p>%s: %s = %.3f, &sigma; = %.3f, P<sub>f</sub> = %.3g %%</p>" % (
                    tr("Monte Carlo"), tr("Mean FS"), mc["mean"], mc["std"], 100 * mc["pf"]))
            t = rel.get("tspm")
            if t and t.get("pf") is not None:
                h.append(u"<p>%s: &beta; = %.2f, P<sub>f</sub> = %.3g %%</p>" % (tr("Taylor series (TSPM)"), t["beta_ln"], 100 * t["pf"]))
        elif rel and rel.get("error"):
            h.append(u"<p style='color:#b9770e'>%s</p>" % rel["error"])
        warn = [m for lvl, m in res.preflight if lvl == "warning"]
        if warn:
            h.append(u"<h3>%s</h3><ul>" % tr("Warnings"))
            for m in warn:
                h.append(u"<li>%s</li>" % m)
            h.append(u"</ul>")
        h.append(u"<p style='color:#777'>%s %.1f s</p>" % (tr("Run time"), res.elapsed))
        return u"".join(h)
