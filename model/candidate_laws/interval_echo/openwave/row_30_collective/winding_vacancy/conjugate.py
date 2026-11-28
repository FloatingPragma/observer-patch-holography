"""Same-inventory addition and removal interferometry with an independently formed gateway."""

import argparse
from hashlib import sha256
import io
import json
from itertools import product
from pathlib import Path
import networkx as nx
import numpy as np

from model.experiments.interval_echo_lepton_register.pole_meter import measure as pole
from . import density, geometry, measure, transfer

HERE = measure.HERE
PROTOCOL = "7db0688f2b1edbcfb9b56066df43e5ab172c77522c1f026e29038fc2d66d4e52"


def local_symmetry():
    law = measure.source.measure.LAW
    image = []
    for slot in range(law.DIM):
        head = slot // 9
        a, b = law.SRC[slot % 9]
        mapped = 0
        if head:
            c, d = law.LABELS[head - 1]
            mapped = 1 + law.LABELS.index((1 - c, -d))
        image.append(9 * mapped + law.SRC.index((2 - a, 2 - b)))
    ew, ec = law.local_term_parts()
    for op in (law.echo(), ew, ec):
        np.testing.assert_array_equal(op[np.ix_(image, image)], op)
    blocks = measure.source.measure.blocks()
    assert not np.any(blocks[0, 1])
    return image


def excess_gateway():
    words, ix, terms = transfer.receiver.build(True)
    tables = measure.source.measure.blocks()
    loads = measure.source.measure.LAW.SRC
    positive = np.zeros_like(terms)
    for k, (left, right) in enumerate(transfer.receiver.EDGES):
        marked = k < 6
        local = tables[1, -1 if k == 0 else 1] if marked else tables["blank"]
        for j, (site, word) in enumerate(words):
            coin = (1 ^ ((word >> k) & 1)) if marked else 0
            column = 9 * coin + loads.index(
                (1 + int(site == left), 1 + int(site == right))
            )
            filled = 9 * coin + loads.index((1, 1))
            positive[k, j, j] -= local[filled, filled]
            for end in np.flatnonzero(local[:, column]):
                out = loads[end % 9]
                new = left if out[0] == 2 else right if out[1] == 2 else site
                newword = word ^ (((end // 9) ^ coin) << k) if marked else word
                positive[k, ix[new, newword], j] += local[end, column]
    np.testing.assert_array_equal(positive, terms)
    return positive.sum(axis=0) / 8, [
        ix[3, a + 8 * b] for a, b in product((0, 7), repeat=2)
    ]


def clock(chart):
    board = measure.source.material.federation(2)
    g = nx.Graph(board.seams)
    used = set(
        chart["signal_ports"]
        + chart["reference_ports"]
        + chart["interferometer_route"]
        + [chart["compensator"]]
    )
    length = nx.multi_source_dijkstra_path_length(g, used)
    a, b = next((a, b) for a, b in board.seams if min(length[a], length[b]) >= 3)
    other = nx.multi_source_dijkstra_path_length(g, used | {a, b})
    donor = max(g, key=lambda v: (other[v], -v))
    assert other[donor] >= 3
    # A blank pair has material-relative onsite 8.5 and exchange 1/2, giving unit gap.
    hb, _ = measure.source.contacts()["blank"]
    full = measure.source.measure.blocks()["blank"]
    inds = [measure.source.measure.LAW.SRC.index(w) for w in ((0, 1), (1, 0))]
    assert full[inds[0], inds[1]] == -4 and full[inds[0], inds[0]] == hb
    p = np.sin(16 / (64 * 2)) ** 2
    n = 2**24
    hit = int(np.random.default_rng(312144).binomial(n, p))
    original = json.loads((Path(pole.__file__).parent / "results.json").read_text())
    return dict(
        dimer=[a, b],
        compensator=donor,
        quadrature="vacancy swap",
        shots=n,
        swaps=hit,
        stays=n - hit,
        native_probability=p,
        measured_step=2 * np.arcsin(np.sqrt(hit / n)) / 16,
        source_meter_step=original["short_wait_clock_estimate"],
        source_meter_counts=original["short_wait_clock_swaps_and_stays"],
    )


def run():
    assert sha256((HERE / "CONJUGATE.md").read_bytes()).hexdigest() == PROTOCOL
    image = local_symmetry()
    chart = transfer.embedded()
    positive, ports = geometry.microscopic(chart["signal_ports"], charge=1)
    negative, _ = geometry.microscopic(chart["signal_ports"])
    np.testing.assert_array_equal(positive, negative)
    e, v = np.linalg.eigh(positive / 8)
    native = (ports, e, v)
    ref, sub = excess_gateway()
    amp, p_ref = transfer.ref_evolution(ref, sub, 1 / 1024)
    t = transfer.TIMES
    rng = np.random.default_rng(311144)
    raw = []
    results = []
    for label, energy, wave in transfer.probes(native):
        # Source rotations and their conjugates act on the corresponding sharp charge.
        np.testing.assert_allclose(positive @ wave / 8, energy * wave, atol=3e-13)
        cross = np.exp(-1j * energy * t[:, None]) * amp.conj()
        counts = []
        for comp in (cross.real, cross.imag):
            p = np.stack(
                [
                    (1 + abs(amp) ** 2 + 2 * comp) / 4,
                    (1 + abs(amp) ** 2 - 2 * comp) / 4,
                    (p_ref - abs(amp) ** 2) / 2,
                    (1 - p_ref) / 2,
                ],
                axis=-1,
            )
            assert p.min() > -1e-12
            counts.append(density.samples(np.maximum(p, 0), 2**20, rng))
        r, i = counts
        observed = ((r[:, :, 0] - r[:, :, 1]) + 1j * (i[:, :, 0] - i[:, :, 1])).mean(
            axis=1
        ) / 2**20
        err = max(abs(np.angle(observed * np.exp(1j * (energy - 4.5) * t))))
        raw.append(counts)
        results.append(
            dict(
                label=label,
                formation=energy,
                worst_phase_residual=err,
                short_measured_phase=-np.angle(observed[0]),
            )
        )
    timer = clock(chart)
    stream = io.BytesIO()
    np.savez_compressed(stream, positive_charge_counts=np.array(raw))
    blob = stream.getvalue()
    neg = json.loads((HERE / "transfer.json").read_text())
    free = neg["gateway"]["modes"][0]["phase_rows"]
    out = dict(
        protocol_sha256=PROTOCOL,
        counts_sha256=sha256(blob).hexdigest(),
        conjugated_local_columns=len(image),
        guide_dimension=len(e),
        gateway_dimension=len(ref),
        native_generator_max_integer_error=0,
        timer=timer,
        passive_negative_short_slope_max_sampled_error=max(
            abs(
                r["short_read_slope"] / 64 / timer["measured_step"]
                - r["single_branch_expected"]
            )
            for r in free
        ),
        positive_modes=results,
        max_positive_phase_error=max(r["worst_phase_residual"] for r in results),
    )
    return json.dumps(out, indent=2, sort_keys=True) + "\n", blob


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    write = p.parse_args().write
    text, blob = run()
    for name, data in [
        ("conjugate.json", text.encode()),
        ("conjugate_counts.npz", blob),
    ]:
        path = HERE / name
        if write:
            path.write_bytes(data)
        else:
            assert path.read_bytes() == data, f"Conjugate transfer mismatch: {name}"
    print("Conjugate transfers written." if write else "Conjugate transfers verified.")
