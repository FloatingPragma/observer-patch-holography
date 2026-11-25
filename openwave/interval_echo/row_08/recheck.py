"""Check the free blank-memory annihilation channel and retain both wall denominators."""

import hashlib
import json
from pathlib import Path

import sympy as sp

from model.candidate_laws.interval_echo import candidate

ROOT = Path(__file__).resolve().parents[3]
MODES = (
    ROOT
    / "model/candidate_laws/interval_echo/openwave/row_08_annihilation/daughter_modes/results.json"
)


def run():
    law = candidate.native_law
    checks = law.exact_checks()
    if not all(checks["checks"].values()):
        raise ArithmeticError("Native symbolic reduction failed")
    swap, contrast = (sp.Matrix(m.tolist()) for m in law._load_ops())
    h = (sp.eye(9) - swap) / 2 + sp.Symbol("kappa") * contrast
    counts = [sp.diag(*(word.count(n) for word in law.SRC)) for n in range(3)]
    if any(h * c != c * h for c in counts):
        raise ArithmeticError("Free blank contact changes its load histogram")
    pair = law.SRC.index((2, 0))
    vacuum = law.SRC.index((1, 1))
    # Positive instrument control: a conversion conserving total load can change the histogram.
    conversion = sp.zeros(9)
    conversion[vacuum, pair] = conversion[pair, vacuum] = 1
    if any(conversion * c == c * conversion for c in counts):
        raise ArithmeticError("Histogram detector misses the planted conversion")
    data = json.loads(MODES.read_bytes())
    (selected,) = [
        r
        for r in data["energy_counts"]
        if r["contrast"] == 8 and not r["cut"] and r["wall"] and r["exposure"] == 32
    ]
    ordered = sum(map(sum, selected["mode_number_counts"]))
    total = ordered + selected["ring_miss"]
    both = selected["mode_number_counts"][1][1]
    if total != 2**26 or sum(map(sum, selected["one_each_mode_counts"])) != both:
        raise ArithmeticError("Incomplete wall spectrum receipt")
    return dict(
        native_checks=len(checks["checks"]),
        free_blank_histogram_commutators=[0] * 3,
        free_20_to_11_matrix_element=str(h[vacuum, pair]),
        planted_20_to_11_matrix_element=str(conversion[vacuum, pair]),
        wall_archive_sha256=hashlib.sha256(MODES.read_bytes()).hexdigest(),
        wall_launches=total,
        ordered_ring=ordered,
        one_material_mode_each=both,
        all_launch_probability=both / total,
        conditional_on_ordered=both / ordered,
        native_guide_gap=data["native_uniform_bulk_gap"],
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
