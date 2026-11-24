"""Test local dipole polarization in the registered blank-seam material sector."""

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier

HERE = Path(__file__).resolve().parent
SHOTS = 2**26
SEED = 290611
TIMES = np.array([0, 1, 2, 4, 8, 16, 32, 64])
LABELS = ("x", "y", "z", "circular_plus", "circular_minus")


def geometry():
    points, seams = carrier.icosahedron_from_coordinates()
    mapping = carrier.isomorphism(carrier.CELL_SEAMS, seams)
    points = points[mapping]
    return points / np.linalg.norm(points[0])


def frame(port):
    pts = geometry()
    z = pts[port]
    x = next(row - z * (row @ z) for row in pts if abs(row @ z) < 0.99)
    x /= np.linalg.norm(x)
    axes = np.stack([x, np.cross(z, x), z], axis=1)
    orbitals = pts @ axes / 2
    np.testing.assert_allclose(orbitals.T @ orbitals, np.eye(3), atol=1e-12)
    return orbitals


def path_and_gates(board):
    graph = nx.Graph()
    ports = {}
    for p, q in board.external:
        a, b = p // 12, q // 12
        graph.add_edge(a, b)
        ports[a, b], ports[b, a] = p % 12, q % 12
    paths = [
        path
        for node, length in nx.single_source_shortest_path_length(graph, 0).items()
        if length == 3
        for path in nx.all_shortest_paths(graph, 0, node)
    ]
    path = min(paths)
    return path, [ports[0, path[1]], *(ports[b, a] for a, b in zip(path, path[1:]))]


def one_bridge(ports, enabled=True):
    seams = list(carrier.CELL_SEAMS) + [(12 + p, 12 + q) for p, q in carrier.CELL_SEAMS]
    if enabled:
        seams.append((ports[0], 12 + ports[1]))
    return carrier.Carrier("two spherical cells", 24, seams)


def generator(board, material, k):
    kinetic = board.laplacian().toarray() / 2
    if material == "vacancy":
        kinetic -= np.diag(2 * k * (6 - board.degree))
    return kinetic  # After subtraction of the single-defect common 6*k.


def spectral_groups(energy):
    return np.split(np.arange(len(energy)), np.flatnonzero(np.diff(energy) > 1e-9) + 1)


def response(matrix):
    s = np.linalg.svd(matrix[:2, :2], compute_uv=False)
    longitudinal = max(np.linalg.norm(matrix[:, 2]), np.linalg.norm(matrix[2, :]))
    return [s[0], s[1], longitudinal]


def qualifies(data):
    a, b, longitudinal = np.moveaxis(np.asarray(data), -1, 0)
    return (a >= 0.001) & (b >= 0.5 * a) & (longitudinal <= 0.1 * b)


def draw(rng, probs):
    flat = np.maximum(np.asarray(probs), 0).round(13)
    np.testing.assert_allclose(flat.sum(), 1, atol=2e-10)
    return rng.multinomial(SHOTS, flat / flat.sum())


def detectors(at_cell, local_frame, admitted, rng):
    # Columns are x, y, z, circular+ and circular- preparations. In addition
    # to the local outcomes keep remote quiet and the clock's rejected port.
    circular = np.array([[1, 1, 0], [1j, -1j, 0], [0, 0, np.sqrt(2)]]) / np.sqrt(2)
    inputs = np.column_stack([np.eye(3), circular[:, :2]])
    visible = local_frame.T @ at_cell
    detected = at_cell @ inputs
    accepted = np.real(np.diag(inputs.conj().T @ admitted @ inputs))
    cell_prob = np.sum(abs(detected) ** 2, axis=0)
    modes = [visible @ inputs, circular.conj().T @ visible @ inputs]
    counts = []
    for amp in modes:
        channels = abs(amp) ** 2
        events = np.vstack(
            [
                channels,
                cell_prob - channels.sum(axis=0),
                accepted - cell_prob,
                1 - accepted,
            ]
        )
        counts.append(np.stack([draw(rng, p) for p in events.T]))
    quadratures = np.stack([visible.real, visible.imag])
    bright = rng.binomial(SHOTS, np.clip((1 + quadratures) / 2, 0, 1).round(13))
    return np.array(counts), bright, response(visible)


def board_run(board, owners, ports, material, k, rng, records, prefix):
    energy, v = eigh(generator(board, material, k), driver="evd")
    groups = spectral_groups(energy)
    local = frame(ports[0])
    source = v[:12].T @ local
    weights = [np.sum(source[g, :2] ** 2) / 2 for g in groups]
    # Break numerical ties using the lower frequency.
    selected = next(i for i, w in enumerate(weights) if w >= max(weights) - 1e-12)
    center = energy[groups[selected][0]]
    filter_factor = np.exp(-0.25j * np.outer(np.arange(256), energy - center)).mean(
        axis=0
    )
    spectral = np.stack(
        (energy[[g[0] for g in groups]], np.array([len(g) for g in groups]), weights),
        axis=1,
    )
    records[prefix + "_lines"] = spectral.round(10)
    result = []
    for steps, (owner, port) in enumerate(zip(owners[1:], ports[1:], strict=True), 1):
        det = v[12 * owner : 12 * owner + 12]
        out = frame(port)
        projected = out.T @ det
        residues = np.array([response(projected[:, g] @ source[g]) for g in groups])
        key = f"{prefix}_d{steps}"
        records[key + "_residues"] = residues.round(11)
        counts, pointers, diagnostics = [], [], []
        for propagation in [*(np.exp(-1j * t * energy) for t in TIMES), filter_factor]:
            waves = propagation[:, None] * source
            c, q, d = detectors(det @ waves, out, waves.conj().T @ waves, rng)
            counts.append(c)
            pointers.append(q)
            diagnostics.append(d)
        counts, pointers = np.array(counts), np.array(pointers)
        records[key + "_detections"] = counts
        records[key + "_quadratures"] = pointers
        decoded = 2 * pointers / SHOTS - 1
        observed = np.array([response(a + 1j * b) for a, b in decoded])
        accepted = np.flatnonzero(qualifies(residues))
        result.append(
            dict(
                steps=steps,
                owner=owner,
                facing_port=port,
                key=key,
                line_pass=accepted.tolist(),
                selected_residue=residues[selected].tolist(),
                free_and_filtered=np.array(diagnostics).round(9).tolist(),
                from_interferometers=observed.round(9).tolist(),
            )
        )
    return dict(
        material=material,
        contrast=k,
        prefix=prefix,
        source_port=ports[0],
        selected_line=int(selected),
        selected_frequency=float(center),
        source_weight=float(weights[selected]),
        number_lines=len(groups),
        paths=result,
    )


def run():
    rng = np.random.default_rng(SEED)
    records, boards = {}, []
    for level in (0, 1):
        full = carrier.federation(level)
        owners, ports = path_and_gates(full)
        for name, board, path, axes in [
            ("full", full, owners, ports),
            ("single", one_bridge(ports), [0, 1], ports[:2]),
            ("shut", one_bridge(ports, False), [0, 1], ports[:2]),
        ]:
            for material, k in (
                ("unit", candidate.KAPPA),
                ("unit", 0),
                ("vacancy", candidate.KAPPA),
                ("vacancy", 0),
            ):
                prefix = f"level{level}_{name}_{material}_{k}"
                row = board_run(board, path, axes, material, k, rng, records, prefix)
                row.update(cells=full.cells, bridge_setting=name, path=owners)
                boards.append(row)
    measurements = sum(
        a.size for key, a in records.items() if key.endswith("_quadratures")
    )
    # Each +/- pointer mean is in [-1,+1]. Union across the real and imaginary pointers.
    eps = np.sqrt(2 * np.log(2 * measurements / 0.01) / SHOTS)
    # With the addressed defect emitter idle the uniform host is an eigenword
    # of every native seam. A degree-reference term is diagonal on that word.
    exchange, contrast = candidate.native_law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    dark = []
    for host in (0, 1):
        idx = candidate.native_law.SRC.index((host, host))
        departure = float(np.linalg.norm(np.delete(h[:, idx], idx)))
        if departure != 0:
            raise ValueError("The homogeneous dark source can generate defects")
        records[f"dark_host_{host}"] = draw(rng, [0.0, 1.0])
        dark.append(
            dict(
                host=host,
                local_departure_norm=departure,
                receiver_clicks=int(records[f"dark_host_{host}"][0]),
            )
        )
    summary = dict(
        shots=SHOTS,
        seed=SEED,
        times=TIMES.tolist(),
        filter_times=256,
        filter_step=0.25,
        labels=LABELS,
        dark_sources=dark,
        quadratures=measurements,
        simultaneous_quadrature_error=eps,
        boards=boards,
    )
    return json.loads(json.dumps(summary)), records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report, arrays = run()
    path, bank = HERE / "results.json", HERE / "counts.npz"
    if args.write:
        path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        np.savez_compressed(bank, **arrays)
    else:
        stored = json.loads(path.read_text())
        assert report.keys() == stored.keys()
        for key in report.keys() - {"boards"}:
            assert report[key] == stored[key]
        assert len(stored["boards"]) == len(report["boards"])
        for old, fresh in zip(stored["boards"], report["boards"], strict=True):
            assert old.keys() == fresh.keys()
            for k in old:
                if k not in ("paths", "selected_frequency", "source_weight"):
                    assert fresh[k] == old[k]
            np.testing.assert_allclose(
                [old["selected_frequency"], old["source_weight"]],
                [fresh["selected_frequency"], fresh["source_weight"]],
                atol=1e-9,
            )
            for a, b in zip(old["paths"], fresh["paths"], strict=True):
                for key in a:
                    if key in (
                        "free_and_filtered",
                        "from_interferometers",
                        "selected_residue",
                    ):
                        np.testing.assert_allclose(a[key], b[key], atol=1e-8)
                    else:
                        assert a[key] == b[key]
        with np.load(bank) as saved:
            assert set(arrays) == set(saved)
            for k, a in arrays.items():
                if a.dtype.kind in "iu":
                    np.testing.assert_array_equal(saved[k], a)
                else:
                    np.testing.assert_allclose(saved[k], a, atol=2e-9)
    print(f"{len(arrays)} spectral and detector arrays checked.")


if __name__ == "__main__":
    main()
