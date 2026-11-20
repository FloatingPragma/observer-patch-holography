"""The replay app as a static site, with its own address for every experiment and OpenWave row.

    python -m simulator app [RUN ...] [--out DIR]

The site holds the code once (app.js, app.css), the list of pages with the OpenWave rows
(site.js) and each run once (data/NAME.js). Every page is a small HTML file that loads them:

    index.html                  the first experiment that is not an OpenWave row's replay
    experiments/                every experiment
    experiments/NAME/           one experiment; NAME is the run file's name without its number
    openwave/                   the 31 OpenWave rows
    openwave/row-N/             one row: the problem, an animated diagram, OPH's result
    openwave/row-N/replay/      the replay of a row's run

The pages open from disk as well as from a web server, and request nothing outside the site.
"""
from __future__ import annotations

import base64
import gzip
import html
import json
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = HERE / "app"
ROWS = HERE / "openwave_rows.json"
DEFAULT_OUT = HERE / "runs" / "site"
NAME = "OPH Simulator"
# Copied as they are: the code, the style and the Pragma Research icons.
STATIC = ("app.js", "app.css", "favicon.svg", "favicon.ico", "apple-touch-icon.png")


def slug(path: Path) -> str:
    """01-continuum.run.json.gz -> continuum."""
    name = path.name
    for suffix in (".run.json.gz", ".run.json", ".json.gz", ".json", ".gz"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    name = re.sub(r"^\d+-", "", name)
    return re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-") or "run"


def first_sentence(text: str) -> str:
    match = re.match(r"(.+?[.!?])(\s|$)", text.strip())
    return match.group(1) if match else text.strip()


def row_path(number: int) -> str:
    return f"openwave/row-{number}/"


def _read(path: Path) -> tuple[bytes, dict]:
    raw = path.read_bytes()
    text = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
    return raw, json.loads(text)


def _script_json(value) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def build(run_paths: list[Path], out: Path = DEFAULT_OUT) -> list[str]:
    """Write the site to ``out`` and return the address of every page, relative to the site."""
    template = (APP / "index.html").read_text(encoding="utf-8")
    for marker in ("{{TITLE}}", "{{DESCRIPTION}}", "{{ROOT}}", "{{PAGE}}", "{{DATA}}"):
        if marker not in template:
            raise SystemExit(f"the app template lacks {marker}")
    openwave = json.loads(ROWS.read_text(encoding="utf-8"))
    for row in openwave.get("rows", []):
        row["path"] = row_path(row["number"])

    runs, data, taken = [], {}, set()
    for path in map(Path, run_paths):
        raw, doc = _read(path)
        name = slug(path)
        while name in taken:
            name += "-2"
        taken.add(name)
        number = doc.get("openwave_row")
        if number:
            page = row_path(number) + "replay/"
            while any(r["path"] == page for r in runs):
                page = page[:-1] + "-2/"
        else:
            page = f"experiments/{name}/"
        step = doc.get("step") or ""
        runs.append({"slug": name, "path": page, "step": step, "title": doc["title"],
                     "label": (step + " · " if step else "") + doc["title"],
                     "kind": doc["kind"], "openwave_row": number,
                     "summary": first_sentence(doc.get("summary") or ""), "data": f"data/{name}.js"})
        data[name] = base64.b64encode(raw).decode()

    if out.exists():
        if any(out.iterdir()) and not (out / "site.js").is_file():
            raise SystemExit(f"{out} holds files that are not an app build; choose an empty folder")
        shutil.rmtree(out)
    (out / "data").mkdir(parents=True)
    for name in STATIC:
        shutil.copyfile(APP / name, out / name)
    (out / "site.js").write_text(
        "window.OPH_SITE = " + _script_json({"runs": runs, "openwave": openwave}) + ";\n", encoding="utf-8")
    for name, b64 in data.items():
        (out / "data" / f"{name}.js").write_text(f'window.OPH_RUNS.push("{b64}");\n', encoding="utf-8")

    pages = []

    def page(address: str, title: str, description: str, state: dict, run: dict | None = None):
        root = "../" * address.count("/") or "./"
        scripts = f'<script src="{root}{run["data"]}"></script>' if run else ""
        text = (template.replace("{{TITLE}}", html.escape(title))
                .replace("{{DESCRIPTION}}", html.escape(description))
                .replace("{{PAGE}}", _script_json(dict(state, root=root)))
                .replace("{{DATA}}", scripts)
                .replace("{{ROOT}}", root))
        target = out / address / "index.html"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")
        pages.append(address)

    about = "Experiments run on the OPH law, replayed from their records, and the 31 rows of the OpenWave comparison."
    if runs:
        first = next((r for r in runs if not r["openwave_row"]), runs[0])
        page("", f"{NAME} · Pragma Research", about, {"view": "run", "run": first["slug"], "home": True}, first)
    else:
        page("", f"{NAME} · Pragma Research", about, {"view": "empty", "home": True})
    page("experiments/", f"Experiments · {NAME}", "Every experiment of the OPH simulator, each on its own page.",
         {"view": "contents", "section": "experiments"})
    for run in runs:
        page(run["path"], f"{run['label']} · {NAME}", run["summary"] or run["label"],
             {"view": "run", "run": run["slug"], "row": run["openwave_row"]}, run)
    if openwave.get("rows"):
        page("openwave/", f"OpenWave rows · {NAME}",
             "The 31 rows of the OpenWave comparison, each with the problem, an animated diagram and OPH's result.",
             {"view": "contents", "section": "openwave"})
        for row in openwave["rows"]:
            page(row["path"], f"OpenWave Row {row['number']}: {row['name']} · {NAME}",
                 f"OpenWave's test: {row['test']}. OPH: {row['status_label']}.",
                 {"view": "row", "row": row["number"]})
    return pages
