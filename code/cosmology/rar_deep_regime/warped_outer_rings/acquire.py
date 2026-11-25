"""Recover the published Table C.1 columns directly from the authors' TeX source."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import tarfile

HERE = Path(__file__).resolve().parent
SOURCE = "https://arxiv.org/src/0810.3025v1"
KEYS = (
    "radius_arcsec",
    "NHI_1e19_cm2",
    "error_NHI",
    "speed_kms",
    "error_speed",
    "inclination_deg",
    "error_inclination",
)


def extract(archive):
    with tarfile.open(archive) as bundle:
        src = bundle.extractfile("xx.tex").read()
    text = (
        src.decode("latin1")
        .split("Radially dependent best-fit parameters.")[1]
        .split(r"\end{tabular}")[0]
    )
    rows = []
    for line in text.splitlines():
        if not line.strip().startswith("$") or "&" not in line:
            continue
        data = re.findall(r"\$\s*(-?\d+(?:\.\d+)?)\s*\$", line)
        if len(data) < 11:
            continue
        values = [
            float(data[i]) if i < len(data) else None for i in (0, 5, 6, 9, 10, 11, 12)
        ]
        rows.append(dict(zip(KEYS, values, strict=True)))
    assert len(rows) == 21 and sorted(set(r["radius_arcsec"] for r in rows)) == [
        r["radius_arcsec"] for r in rows
    ]
    return dict(
        source=SOURCE, table="C.1", tex_sha256=sha256(src).hexdigest(), rows=rows
    )


def check_pdf(pdf, rows):
    """Check the TeX extraction against the independently typeset table."""
    text = subprocess.check_output(["pdftotext", "-layout", str(pdf), "-"], text=True)
    text = text.split("Table C.1. Radially dependent best-fit parameters.")[-1]
    matched = []
    for line in text.splitlines():
        fields = line.replace("−", "-").split()
        try:
            data = list(map(float, fields))
        except ValueError:
            continue
        if len(data) not in (11, 21):
            continue
        matched.append(
            [data[i] if i < len(data) else None for i in (0, 5, 6, 9, 10, 11, 12)]
        )
    if matched != [[row[key] for key in KEYS] for row in rows]:
        raise ValueError(
            "Typeset Table C.1 disagrees with the extracted scientific columns"
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    p.add_argument("--pdf", type=Path)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    record = extract(args.archive)
    if args.pdf:
        check_pdf(args.pdf, record["rows"])
    data = json.dumps(record, indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "rings.json").write_text(data)
    else:
        print(data, end="")
