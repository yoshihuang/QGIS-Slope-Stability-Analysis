# -*- coding: utf-8 -*-
"""Build the Qt translation catalogues of the plugin.

Usage (from the plugin folder)::

    python tools/i18n_build.py            # update i18n/*.ts and compile *.qm
    python tools/i18n_build.py --check    # list untranslated strings

* Source strings are collected from every ``tr("...")`` call (and
  ``QCoreApplication.translate("SlopeStability", "...")``) using the Python AST,
  so only literal strings are extracted.
* Translations are kept in ``tools/translations/<lang>.json`` ({source: text}).
* ``.qm`` files are compiled with ``lrelease`` (Qt Linguist tools), or
  ``pyside6-lrelease`` / ``lrelease-qt5`` if available.
"""

import ast
import json
import os
import shutil
import subprocess
import sys
from xml.sax.saxutils import escape

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CONTEXT = "SlopeStability"
LANGS = ["zh_TW", "zh_CN", "ja", "ko", "es", "fr", "de", "pt"]
TS_LANG = {"zh_TW": "zh_TW", "zh_CN": "zh_CN", "ja": "ja_JP", "ko": "ko_KR", "es": "es_ES",
           "fr": "fr_FR", "de": "de_DE", "pt": "pt_PT"}


def extract():
    found = {}
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in ("tools", "__pycache__", "i18n", "testdata")]
        for fn in sorted(filenames):
            if not fn.endswith(".py"):
                continue
            path = os.path.join(dirpath, fn)
            with open(path, encoding="utf-8") as f:
                tree = ast.parse(f.read(), path)
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                func = node.func
                name = func.id if isinstance(func, ast.Name) else (func.attr if isinstance(func, ast.Attribute) else "")
                args = node.args
                text = None
                if name == "tr" and args and isinstance(args[0], ast.Constant) and isinstance(args[0].value, str):
                    text = args[0].value
                elif name == "translate" and len(args) >= 2 and isinstance(args[1], ast.Constant) \
                        and isinstance(args[0], ast.Constant) and args[0].value == CONTEXT:
                    text = args[1].value
                if text:
                    rel = os.path.relpath(path, ROOT).replace(os.sep, "/")
                    found.setdefault(text, (rel, node.lineno))
    return found


def write_ts(lang, strings, trans):
    lines = ['<?xml version="1.0" encoding="utf-8"?>', "<!DOCTYPE TS>",
             '<TS version="2.1" language="%s">' % TS_LANG.get(lang, lang), "<context>",
             "    <name>%s</name>" % CONTEXT]
    missing = []
    for src in sorted(strings):
        loc, line = strings[src]
        t = trans.get(src)
        lines.append("    <message>")
        lines.append('        <location filename="../%s" line="%d"/>' % (loc, line))
        lines.append("        <source>%s</source>" % escape(src))
        if t:
            lines.append("        <translation>%s</translation>" % escape(t))
        else:
            missing.append(src)
            lines.append('        <translation type="unfinished"></translation>')
        lines.append("    </message>")
    lines += ["</context>", "</TS>", ""]
    path = os.path.join(ROOT, "i18n", "SlopeStability_%s.ts" % lang)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return path, missing


def lrelease():
    for exe in ("lrelease", "lrelease-qt5", "lrelease-qt6", "pyside6-lrelease"):
        p = shutil.which(exe)
        if p:
            return p
    return None


def main():
    check = "--check" in sys.argv
    strings = extract()
    print("%d translatable strings" % len(strings))
    os.makedirs(os.path.join(ROOT, "i18n"), exist_ok=True)
    tool = lrelease()
    bad = 0
    for lang in LANGS:
        jp = os.path.join(HERE, "translations", "%s.json" % lang)
        trans = {}
        if os.path.exists(jp):
            with open(jp, encoding="utf-8") as f:
                trans = json.load(f)
        unused = [k for k in trans if k not in strings]
        ts, missing = write_ts(lang, strings, trans)
        print("%-6s translated %d/%d%s" % (lang, len(strings) - len(missing), len(strings),
                                          ("  unused: %d" % len(unused)) if unused else ""))
        if check:
            for m in missing:
                print("    MISSING:", m)
            for u in unused:
                print("    UNUSED:", u)
        bad += len(missing)
        if tool and not check:
            qm = ts[:-3] + ".qm"
            subprocess.check_call([tool, "-silent", ts, "-qm", qm])
    if not tool and not check:
        print("lrelease not found: .qm files were not compiled")
    return 1 if (check and bad) else 0


if __name__ == "__main__":
    sys.exit(main())
