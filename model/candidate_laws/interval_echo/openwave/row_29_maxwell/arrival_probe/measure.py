"""Local face current and deflection of an independently launched material double."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

from model.experiments.material_inertia.neutral_defect import generators
from .. import experiment as dipole
from ..polarization_channel.measure import INPUTS

HERE = Path(__file__).resolve().parent
SHOTS = 2**38
SEED = 294164
TIMES = (0.125, 0.5)
NAMES = ("x", "y", "diag+", "diag-", "circle+", "circle-", "z", "blank")
PROBES = ("plus", "minus", "real")
# 2 laws, 2 times, 3 probes, 2 controls, 8 sources, 4 projections.
EPS = float(np.sqrt(2 * np.log(2 * 2 * 2 * 3 * 2 * 8 * 4 / 0.01) / SHOTS))


def draw(rng, rates):
    p = np.asarray(rates, dtype=float)
    if np.min(p) < -2e-12 or abs(np.sum(p) - 1) > 3e-11:
        raise ValueError("invalid complete outcome law")
    p = np.maximum(0, p).round(12)
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


def oriented_faces(cell, points):
    faces = []
    for a, b, c in cell.faces():
        normal = np.cross(points[b] - points[a], points[c] - points[a])
        if normal @ points[a] < 0:
            b, c = c, b
            normal *= -1
        faces.append(((a, b, c), normal / np.linalg.norm(normal)))
    return faces


def read_faces(psi, faces, rng, blocked=None):
    arrived = float(np.vdot(psi, psi).real)
    totals = []
    for face, _ in faces:
        amplitudes = psi[list(face)]
        # +i/2 hops forward in the oriented current; e^(+2pi i j/3)
        # has its positive eigenvalue sqrt(3)/2.
        modes = np.fft.fft(amplitudes) / np.sqrt(3)
        plus, minus = abs(modes[1:]) ** 2
        other = [1 - arrived] if blocked is None else [1 - arrived - blocked, blocked]
        totals.append(draw(rng, [plus, minus, arrived - plus - minus, *other]))
    return totals


def inputs():
    board = dipole.carrier.federation(1)
    route, gates = dipole.path_and_gates(board)
    h = dipole.generator(board, "unit", dipole.candidate.KAPPA)
    energy, v = eigh(h, driver="evd")
    source = v[:12].T @ dipole.frame(gates[0])
    d = 12 * route[-1]
    cell = v[d : d + 12] @ (np.exp(-64j * energy)[:, None] * source)
    transverse = cell[:, :2] @ INPUTS.T
    words = np.column_stack((transverse, cell[:, 2], np.zeros(12)))
    return route, gates, words.T


def probe_banks(psi, probes, a, b, x, y, one_gate, two_gate, rng):
    psi = psi.copy()
    remote = 1 - np.vdot(psi, psi).real
    blocked = sum(abs(psi[[a, b]]) ** 2)
    psi[[a, b]] = 0
    admit = np.vdot(psi, psi).real
    bank = {}
    for label, launch in zip(PROBES, probes, strict=True):
        alone = abs(one_gate @ launch) ** 2
        both = two_gate @ (launch[x] * psi[y])
        interacting = np.bincount(x, weights=abs(both) ** 2, minlength=12)
        bank[label] = {
            tag: draw(rng, [blocked, *(remote * alone), *local])
            for tag, local in (("contact", interacting), ("stored", admit * alone))
        }
    return bank


def contrast(bank, points, u):
    diffs = {}
    for mode in PROBES:
        case = bank[mode]
        local = np.subtract(case["contact"][13:], case["stored"][13:])
        diffs[mode] = local @ points / SHOTS
    odd = (diffs["plus"] - diffs["minus"]) / 2
    even = (diffs["plus"] + diffs["minus"]) / 2
    return {
        "odd": odd.round(10).tolist(),
        "even": even.round(10).tolist(),
        "real": diffs["real"].round(10).tolist(),
        "odd_parallel": round(float(odd @ u), 10),
    }


def run():
    route, gates, words = inputs()
    cell = dipole.carrier.cell()
    points = dipole.geometry()
    a = gates[-1]
    b = min(cell.neighbours[a])
    u = points[b] - points[a]
    u /= np.linalg.norm(u)
    probes = np.zeros((3, 12), complex)
    probes[:, a] = 1 / np.sqrt(2)
    probes[:, b] = np.array([1j, -1j, 1]) / np.sqrt(2)
    faces = oriented_faces(cell, points)
    x, y, adj, kinetic, _, _ = generators(cell)
    rng = np.random.default_rng(SEED)
    out = dict(
        shots=SHOTS,
        seed=SEED,
        route=route,
        gateways=gates,
        probe_ports=[a, b],
        direction=u.tolist(),
        face_order=[list(f) for f, _ in faces],
        cohort_projection_half_width=EPS,
        contrast_half_width=2 * EPS,
        face_counts={},
        axial_reading={},
        records={},
    )
    normals = np.array([n for _, n in faces])
    for name, psi in zip(NAMES, words, strict=True):
        faces_read = read_faces(psi, faces, rng)
        out["face_counts"][name] = faces_read
        counts = np.asarray(faces_read)
        axial = (counts[:, 0] - counts[:, 1]) @ normals * np.sqrt(3) / (40 * SHOTS)
        out["axial_reading"][name] = axial.round(10).tolist()
    matched_rng = np.random.default_rng(SEED + 1)
    out["matched_seed"] = SEED + 1
    out["matched_face_counts"] = {}
    out["matched_axial_reading"] = {}
    for name, psi in zip(NAMES, words, strict=True):
        admitted = psi.copy()
        admitted[[a, b]] = 0
        blocked = np.sum(abs(psi[[a, b]]) ** 2)
        counts = np.array(read_faces(admitted, faces, matched_rng, blocked))
        out["matched_face_counts"][name] = counts.tolist()
        axial = (counts[:, 0] - counts[:, 1]) @ normals * np.sqrt(3) / (40 * SHOTS)
        out["matched_axial_reading"][name] = axial.round(10).tolist()
    # Subtract the reference-vacancy sector constant. With shut bridges
    # every single receives the same diagonal correction.
    for law, kappa in (("native", dipole.candidate.KAPPA), ("zero_contrast", 0)):
        out["records"][law] = {}
        pair = kinetic.toarray() - np.diag(4 * kappa * adj)
        for t in TIMES:
            one_gate = expm((-0.5j * t) * cell.laplacian().toarray())
            two_gate = expm((-1j * t) * pair)
            panel = {}
            for name, psi in zip(NAMES, words, strict=True):
                bank = probe_banks(psi, probes, a, b, x, y, one_gate, two_gate, rng)
                panel[name] = dict(counts=bank, **contrast(bank, points, u))
            out["records"][law][str(t)] = panel
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    dst = HERE / "records.json"
    if args.write:
        dst.write_text(json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n")
    elif data != json.loads(dst.read_text()):
        raise ValueError("archived counts differ")
    print("Full counts verified. Simultaneous contrast half-width:", 2 * EPS)
    for law, panels in data["records"].items():
        for t, panel in panels.items():
            print(law, t, {n: v["odd_parallel"] for n, v in panel.items()})


if __name__ == "__main__":
    main()
