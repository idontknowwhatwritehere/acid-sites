#!/usr/bin/env python3
"""Prepare GitHub Pages copies from the static exports of the two Sites projects."""

import argparse
from pathlib import Path
import re
import shutil


def replace_once(content, old, new):
    if content.count(old) != 1:
        raise ValueError(f"Expected one occurrence of {old!r}")
    return content.replace(old, new, 1)


def rename_page(html, name):
    for tag in ("title", "h1"):
        html, count = re.subn(fr"<{tag}>[^<]*</{tag}>", f"<{tag}>{name}</{tag}>", html)
        if count != 1:
            raise ValueError(f"Expected one {tag} on {name}")
    return html


def write_page(path, source):
    path.write_text(source, encoding="utf-8")


def prepare_calibration(source, target):
    shutil.copytree(source, target)
    index = target / "index.html"
    html = rename_page(index.read_text(encoding="utf-8"), "Калибровка кислот")
    html = replace_once(
        html,
        "'https://hexbn-acid-permeability.nocanbewant.chatgpt.site'",
        "location.origin",
    )
    for tag, filename, include in (
        ("style", "calibration.css", '<link rel="stylesheet" href="./calibration.css">'),
        ("script", "calibration.js", '<script src="./calibration.js"></script>'),
    ):
        matches = list(re.finditer(fr"<{tag}>(.*?)</{tag}>", html, re.DOTALL))
        if len(matches) != 1:
            raise ValueError(f"Expected one inline {tag} block")
        match = matches[0]
        write_page(target / filename, match.group(1))
        html = html[:match.start()] + include + html[match.end():]
    write_page(index, html)


def prepare_permeability(source, target):
    shutil.copytree(source, target)
    index = target / "index.html"
    write_page(index, rename_page(index.read_text(encoding="utf-8"), "Проницаемость кислот"))
    app = target / "app.js"
    js = replace_once(
        app.read_text(encoding="utf-8"),
        "const origin='https://ibu-bn-acid-calibration.nocanbewant.chatgpt.site';",
        "const calibrationUrl=new URL('../acid-calibration/#hexbn-calibrations',location.href);\nconst origin=calibrationUrl.origin;",
    )
    js = replace_once(js, "origin+'/#hexbn-calibrations'", "calibrationUrl.href")
    write_page(app, js)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("calibration", type=Path, help="Calibration dist directory")
    parser.add_argument("permeability", type=Path, help="Permeability dist directory")
    parser.add_argument("output", type=Path, help="A new directory for the two builds")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output must be a new directory")
    prepare_calibration(args.calibration, args.output / "acid-calibration")
    prepare_permeability(args.permeability, args.output / "acid-permeability")


if __name__ == "__main__":
    main()
