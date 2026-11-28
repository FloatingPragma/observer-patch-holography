"""Transfer the independent nearest-gateway-line instrument to one and two winding charges."""

import argparse
import json
from hashlib import sha256

import numpy as np
from scipy.linalg import eigh

from . import copies, geometry

HERE = copies.HERE
PROTOCOL = "5c2d78cfc8d6292736bda80935b310af90468919cc3121fea46bbc91175b6636"
SHOTS = 2**22


def read_filter(weights, spectrum, selected, stages, rng):
    dif = spectrum - selected
    target = abs(dif) < 1e-10
    f = np.ones_like(dif)
    for j in range(stages):
        f *= np.cos((2**j / 8) * dif / 2) ** 2
    filtered = weights * f
    accepted = filtered.sum()
    line = filtered[target].sum()
    diagnostic = rng.multinomial(
        SHOTS, [line, max(0.0, accepted - line), max(0.0, 1 - accepted)]
    )
    out = []
    for t in [0, 1, 17, 257, 8192]:
        z = np.sum(filtered * np.exp(-1j * t * dif))
        raw = []
        for val in (z.real, z.imag):
            p = np.maximum(
                [(accepted + val) / 2, (accepted - val) / 2, 1 - accepted], 0
            )
            raw.append(rng.multinomial(SHOTS, p / p.sum()))
        raw = np.array(raw)
        sample = (raw[:, 0] - raw[:, 1]) / (raw[:, 0] + raw[:, 1])
        phase = float(np.angle(sample[0] + 1j * sample[1]))
        out.append(
            dict(
                time=t,
                counts=raw.tolist(),
                phase=phase,
                exact_phase=float(np.angle(z)),
                exact_conditional_contrast=abs(z) / accepted,
            )
        )
    return dict(
        stages=stages,
        accepted=accepted,
        conditional_target_weight=line / accepted,
        ideal_diagnostic_counts=diagnostic.tolist(),
        phase_reads=out,
        sampled_max_phase=max(abs(r["phase"]) for r in out),
        pass_both=line / accepted > 0.99 and max(abs(r["phase"]) for r in out) < 0.02,
    )


def run():
    assert sha256((HERE / "COPIES_SOURCE.md").read_bytes()).hexdigest() == PROTOCOL
    data = json.loads((HERE / "copies.json").read_text())
    h, _ = geometry.microscopic(data["stock"]["routes"][0])
    h = h / 8
    energy, modes = eigh(h)
    prob = modes**2
    pair = energy[:, None] + energy
    single_records, pairs = [], []
    for sharp in np.unique(h.diagonal()):
        x = np.flatnonzero(h.diagonal() == sharp)
        line = energy[np.argmin(abs(energy - sharp))]
        sel = abs(energy - line) < 1e-10
        w = prob[:, sel].sum(axis=1)
        single_records.append(
            dict(
                diagonal=sharp,
                positions=len(x),
                selected=line,
                rank=int(sel.sum()),
                min_input_weight=min(w[x]),
                max_input_weight=max(w[x]),
            )
        )
    for a in np.unique(h.diagonal()):
        for b in np.unique(h.diagonal()):
            target = a + b
            selected = pair.flat[np.argmin(abs(pair - target))]
            sel = abs(pair - selected) < 1e-10
            w = prob @ sel.astype(float) @ prob.T
            rows = np.flatnonzero(h.diagonal() == a)
            cols = np.flatnonzero(h.diagonal() == b)
            sub = w[np.ix_(rows, cols)]
            pairs.append(
                dict(
                    diagonals=[a, b],
                    sharp_placements=sub.size,
                    selected=selected,
                    rank=int(sel.sum()),
                    nearest_distinct_gap=min(abs(pair[~sel] - selected)),
                    min_input_weight=sub.min(),
                    mean_input_weight=sub.mean(),
                    max_input_weight=sub.max(),
                )
            )
    rng = np.random.default_rng(315144)
    results = []
    interior = prob[h.diagonal() == h[1, 1]].mean(axis=0)
    for tag, v in [("first_interior", prob[1]), ("uniform_interior", interior)]:
        weights = np.outer(v, v)
        results.append(
            dict(
                sharp_source=tag,
                records=[
                    read_filter(weights, pair, 2 * h[1, 1], s, rng)
                    for s in (0, 13, 19, 23)
                ],
            )
        )
    q, D = np.array([[r[key] for r in data["rows"]] for key in ("q", "energy")])
    m, c = data["rest"], data["contact_speed"]
    idx = int(np.argmin(abs(q - m / c)))
    return dict(
        protocol_sha256=PROTOCOL,
        single_sharp=single_records,
        paired_sharp=pairs,
        filters=results,
        fixed_law=dict(
            gap=3.5,
            upper_bandwidth=1.0,
            bulk_compton_length_cells=2**0.5 / 3.5,
            bandwidth_over_rest=1 / 3.5,
            cross_regime_index=idx,
            measured_q=q[idx],
            pair_energy=D[idx],
            continuum_energy=(m * m + (q[idx] * c) ** 2) ** 0.5,
            continuum_squared_error=D[idx] ** 2 - m**2 - (q[idx] * c) ** 2,
            measured_lattice_energy=(m * m + 4 * c * c * np.sin(q[idx] / 2) ** 2)
            ** 0.5,
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    output = (
        json.dumps(run(), sort_keys=True, indent=2, default=lambda v: v.item()) + "\n"
    )
    target = HERE / "copy_source.json"
    if parser.parse_args().write:
        target.write_text(output)
    elif target.read_text() != output:
        raise AssertionError("Gateway selector replay mismatch")
    print("Verified sharp-gateway selection and fixed-coupling scale test.")
