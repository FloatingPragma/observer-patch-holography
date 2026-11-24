"""Complete cuts and tapered aperture interference for both written-wall preparations."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import currents
from . import field_probe, polarization as pol


def run():
    path, _, _ = pol.layout()
    old = pol.common_heads(path)
    base = pol.potential(0, 0)
    # The reference and sample can each point either way on each boundary.
    cut = np.flatnonzero(np.diff(base))
    assert len(cut) == 2
    menu = []
    for flags in product(range(2), repeat=4):
        before, after = old.copy(), old.copy()
        before[cut] = flags[:2]
        after[cut] = flags[2:]
        q = int(np.diff(base) @ ((before - after) * np.sign(np.diff(path))))
        reads = field_probe.interferometer(before, after, path, base, np.pi, 1)
        y = reads["quadratures"]["y"]
        delta = y[0] - y[1]
        assert abs(delta) < 1e-12
        menu.append(
            dict(
                heads=flags,
                net=q,
                imaginary_contrast=delta,
                visibility=reads["visibility"],
            )
        )
    sources = []
    rng = np.random.default_rng(120131)
    for shift in (1, -2):
        bulk, n, flux, _ = pol.translated(path, shift)
        new = old ^ (flux != 0).astype(int)
        for width in (0, 8, 32):
            f = pol.potential(0, width)
            for sign in (1, -1):
                read = field_probe.interferometer(old, new, path, f, np.pi, sign)
                y = np.array(read["quadratures"]["y"])
                shots = 2**20
                counts = rng.multinomial(shots, y).tolist()
                sources.append(
                    dict(
                        shift=shift,
                        width=width,
                        sign=sign,
                        applied_cut_charge=float(f @ (n - bulk)),
                        y_probabilities=y,
                        y_counts=counts,
                        measured_y=(counts[0] - counts[1]) / shots,
                    )
                )
    return currents._portable(
        dict(
            cut_sites=cut.tolist(),
            all_two_end_assignments=menu,
            panels=sources,
            shots_per_panel=shots,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = run()
    with args.output.open("x") as f:
        f.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
