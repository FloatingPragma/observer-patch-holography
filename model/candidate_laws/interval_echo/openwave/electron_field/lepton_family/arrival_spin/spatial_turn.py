"""Steer the parked pair around its physical hub and compare a single mass anchor to data."""

import argparse
from collections import Counter, defaultdict
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
# The archive replay launches this script by filename from a temporary directory.
# ruff: noqa: E402
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
    rest_partner as rest,
)


def orbit_faces(fixture):
    _, hub, first, second = fixture
    origin = hub // 12 * 12
    turns = rest.wind.circ.tree.rotations()
    for rotate in turns:

        def turn(face):
            return tuple(sorted(origin + rotate[p - origin] for p in face))

        if turn((hub,)) != (hub,) or turn(turn(first)) != second:
            continue
        faces = [first]
        for _ in range(4):
            faces.append(turn(faces[-1]))
        if (
            len(set(first) & set(faces[1])) != 2
            or len(set(faces)) != 5
            or turn(faces[-1]) != first
        ):
            continue
        return hub, faces, tuple(rotate)
    raise ArithmeticError("No adjacent proper 72-degree turn")


class Wheel:
    def __init__(self, fixture):
        self.hub, self.faces, self.turn = orbit_faces(fixture)
        self.edges = sorted({e for f in self.faces for e in rest.wind.face_edges(f)})
        self.ports = sorted({v for e in self.edges for v in e})
        self.masks = []
        self.cw = []
        # Transport the orientation convention with the physical turn.
        face = self.faces[0]
        for _ in range(5):
            ring = list(zip(face, (*face[1:], face[0])))
            self.masks.append(
                sum(1 << self.edges.index(rest.wind.edge(x, y)) for x, y in ring)
            )
            self.cw.append(
                sum(
                    1 << self.edges.index(rest.wind.edge(x, y))
                    for x, y in ring
                    if y < x
                )
            )
            base = self.hub // 12 * 12
            face = tuple(base + self.turn[v - base] for v in face)
        self.groups = defaultdict(list)
        for word in range(1 << len(self.edges)):
            heads = [
                e[0] if word & (1 << j) else e[1] for j, e in enumerate(self.edges)
            ]
            self.groups[tuple(heads.count(p) for p in self.ports)].append(word)

    def directed(self, w, j):
        return (w & self.masks[j]) in (self.cw[j], self.cw[j] ^ self.masks[j])

    def pose_words(self, words, j):
        return [
            w for w in words if self.directed(w, j) and self.directed(w, (j + 2) % 5)
        ]

    def frame(self, words, mass):
        n = len(words)
        rows = []
        tags = []
        for i, w in enumerate(words):
            for j, mask in enumerate(self.masks):
                other = w ^ mask
                if self.directed(w, j) and w < other:
                    row = np.zeros(n)
                    row[i], row[words.index(other)] = -np.sqrt(0.5), np.sqrt(0.5)
                    rows.append(row)
                    tags.append(j)
        rows.extend(np.sqrt(8) * np.eye(n))
        tags.extend([5] * n)
        b = np.array(rows)
        tags = np.array(tags)

        def generator(j):
            return rest.guide.generator(
                b * np.isin(tags, [j, (j + 2) % 5, 5])[:, None], np.array([[mass]])
            )

        def pair(j):
            selected = self.pose_words(words, j)
            if len(selected) != 4:
                raise ArithmeticError("Orientation does not have four source words")
            faces = (j, (j + 2) % 5)
            chi = np.zeros(n)
            for w in selected:
                chi[words.index(w)] = (-1) ** sum(
                    (w & self.masks[f]) != self.cw[f] for f in faces
                ) / 2
            columns = []
            for tag in (*faces, 5):
                column = (b @ chi) * (tags == tag)
                columns.append(column / np.linalg.norm(column))
            normal = np.array([1, 1, np.sqrt(8)]) / np.sqrt(10)
            axis = np.eye(3)[0] - normal[0] * normal
            axis /= np.linalg.norm(axis)
            target = np.column_stack(columns) @ np.column_stack(
                (axis, np.cross(normal, axis))
            )
            state = np.vstack((np.zeros((n, 2)), target))
            np.testing.assert_allclose(generator(j) @ state, -mass * state, atol=2e-14)
            return state

        return generator, pair


def steer(eigens, source, duration, ticks):
    wave = source.copy()

    def hold(wave, setting, time):
        e, v = eigens[setting]
        return v @ (np.exp(-1j * time * e)[:, None] * (v.T @ wave))

    for i in range(ticks):
        new = (1 - np.cos(np.pi * (i + 0.5) / ticks)) / 2
        dt = duration / ticks
        wave = hold(wave, 0, dt * (1 - new) / 2)
        wave = hold(wave, 1, dt * new)
        wave = hold(wave, 0, dt * (1 - new) / 2)
    np.testing.assert_allclose(wave.conj().T @ wave, np.eye(2), atol=2e-11)
    return wave


def wheel_experiment(fixture, mass):
    wheel = Wheel(fixture)
    records = []
    coverage = Counter()
    for heads, words in sorted(wheel.groups.items()):
        poses = [j for j in range(5) if wheel.pose_words(words, j)]
        if not poses:
            continue
        coverage[len(poses)] += 1
        if 0 not in poses:
            continue
        row = {"inward_counts": heads, "reachable_poses_including_initial": poses}
        records.append(row)
        if 1 not in poses:
            continue
        g, p = wheel.frame(words, mass)
        source, target = p(0), p(1)
        eigens = [eigh(g(j)) for j in (0, 1)]
        panels = []
        for duration in (4, 16, 64):
            out = steer(eigens, source, duration, 128)
            k = target.T @ out
            singular = np.linalg.svd(k, compute_uv=False)
            e, v = eigens[0]
            sham = (v * np.exp(-1j * duration * e)) @ (v.T @ source)
            np.testing.assert_allclose(
                sham, np.exp(1j * mass * duration) * source, atol=3e-12
            )
            panels.append(
                {
                    "exposure": duration,
                    "transfer_real": k.real.tolist(),
                    "transfer_imag": k.imag.tolist(),
                    "minimum_and_maximum_target_probability": (
                        singular[::-1] ** 2
                    ).tolist(),
                    "failure_for_mixed_input": float(1 - np.sum(singular**2) / 2),
                    "initial_to_target_sham_mixed": float(
                        np.linalg.norm(target.T @ sham) ** 2 / 2
                    ),
                }
            )
        row["steering"] = panels
    return {
        "hub": wheel.hub,
        "actual_faces": wheel.faces,
        "proper_port_turn": wheel.turn,
        "wheel_edges": wheel.edges,
        "head_ports": wheel.ports,
        "all_words": sum(map(len, wheel.groups.values())),
        "pose_counts_across_sectors": dict(sorted(coverage.items())),
        "complete_five_step_sectors": coverage[5],
        "source_sectors": records,
    }


def rigid_half_turns():
    ident = tuple(range(12))
    half = [
        a
        for a in rest.wind.circ.tree.rotations()
        if a != ident and tuple(a[i] for i in a) == ident
    ]
    seams = rest.wind.circ.tree.cell().seams
    records = []
    for a, b in combinations(half, 2):
        ab = tuple(a[i] for i in b)
        if ab != tuple(b[i] for i in a):
            continue
        group = (ident, a, b, ab)
        port_orbits = {tuple(sorted({g[p] for g in group})) for p in ident}
        edge_orbits = {
            tuple(sorted({rest.wind.edge(g[u], g[v]) for g in group})) for u, v in seams
        }
        sizes = list(map(len, edge_orbits))
        total = [
            sum(s for j, s in enumerate(sizes) if mask & (1 << j))
            for mask in range(1 << len(sizes))
        ]
        # One vacancy gives sum(internal heads minus signed load) = number of written edges + 1.
        allowed = sum((n + 1) % 4 == 0 for n in total)
        if {len(o) for o in port_orbits} != {4} or allowed:
            raise ArithmeticError("Unexpected conserved-count orbit")
        records.append(
            {
                "turns": [a, b],
                "port_orbit_lengths": sorted(map(len, port_orbits)),
                "edge_orbit_lengths": sorted(sizes),
                "invariant_masks": len(total),
                "compatible_masks": allowed,
            }
        )
    if len(records) != 15:
        raise ArithmeticError("Incomplete perpendicular-pair census")
    return records


def anchor_test(masses):
    # 2026 PDG lepton summary, pages 1 and 2. All entries are mass energies in MeV.
    observed = np.array([0.51099895069, 105.6583755, 1776.93])
    errors = np.array([0.00000000016, 0.0000023, 0.09])
    ratios = masses / masses[0]
    predicted = observed[0] * ratios
    return {
        "physical_reference": "https://pdg.lbl.gov/2026/tables/rpp2026-sum-leptons.pdf",
        "unit": "MeV",
        "observed": observed.tolist(),
        "reported_uncertainties": errors.tolist(),
        "native_signed_levels": (-masses).tolist(),
        "anchored_magnitudes": predicted.tolist(),
        "literal_511_keV_anchor_magnitudes": (0.511 * ratios).tolist(),
        "prediction_minus_measurement": (predicted - observed).tolist(),
        "reverse_muon_anchor_electron_keV": float(observed[1] / ratios[1] * 1e3),
    }


def run():
    board, route = rest.guide.longest_unpinned_route()
    fixture = min(rest.geometries(board, route))
    _, masses = rest.guide.register_effect()
    return {
        "wheel": wheel_experiment(fixture, float(masses[0])),
        "rigid_perpendicular_pairs": rigid_half_turns(),
        "energy_scale": anchor_test(masses),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")
