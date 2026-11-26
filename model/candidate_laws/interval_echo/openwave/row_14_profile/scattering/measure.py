"""Coherent address and balanced material probes on the guarded, filtered source."""

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np

from model.candidate_laws.interval_echo.openwave.row_14_profile import measure as source
from model.executable_model import federation

from .readout import GRID, STARTS, decode

HERE = Path(__file__).resolve().parent
SEED = 140326
DURATION = np.pi / 32
APERTURES = {"all": (0, 1, 2, 3), "center": (0,), "outer": (2, 3)}
KAPPA = source.candidate.KAPPA


def routes(level):
    """Reversible address transport on blank internal seams to one external probe."""
    board = federation(level)
    graph = nx.Graph(board.seams)
    gateway = 1
    (external,) = [i for i in graph[gateway] if i >= 12]
    spectator = min(i for i in graph[external] if i // 12 == external // 12)
    cell = graph.subgraph(range(12))
    paths = [nx.shortest_path(cell, i, gateway) for i in range(12)]
    for i, path in enumerate(paths):
        # Reverse every exchange after the contact, including in the spill branch.
        word = list(range(12))
        for x, y in zip(path, path[1:]):
            word[x], word[y] = word[y], word[x]
        if word[gateway] != i:
            raise ValueError("Failed address transport")
        for x, y in zip(path[-1:0:-1], path[-2::-1]):
            word[x], word[y] = word[y], word[x]
        if word != list(range(12)):
            raise ValueError("Failed inverse transport")
    gateways = list(map(int, np.flatnonzero(board.degree[:12] == 6)))
    if gateways != [1, 3, 10]:
        raise ValueError("Different private-reference geometry")
    return dict(paths=paths, probe_contact=[gateway, external], spectator=spectator)


def gain(kappa):
    # Matched order stays, while the opposite order has mismatch two.
    no_exchange = np.exp(-1j * (4 * kappa + 0.5) * DURATION) * np.cos(DURATION / 2)
    return (1 - no_exchange) / 2


def detector(mean, power, g, bypass):
    """Probabilities summed over the selected material ensemble."""
    elastic = abs(g) ** 2 * power
    if bypass:
        sign = g.real * mean
        return np.array(
            [
                (1 + elastic + 2 * sign) / 4,
                (1 + elastic - 2 * sign) / 4,
                (1 - elastic) / 2,
            ]
        )
    return np.array([elastic, 1 - elastic])


def acquire(rng, blocks, acceptance):
    """Round all outcomes on a fixed grid, assigning its residue to the largest cell."""
    original = np.append(acceptance * np.asarray(blocks).ravel(), 1 - acceptance)
    if np.min(original) < 0 or abs(original.sum() - 1) > 2e-11:
        raise ArithmeticError("Incomplete scattering instrument")
    largest = int(np.argmax(original))
    weights = np.rint(original * GRID).astype(np.int64)
    weights[largest] += GRID - weights.sum()
    if np.min(weights) < 0:
        raise ArithmeticError("Negative sampling weight")
    sampled = weights / GRID
    if np.abs(sampled - original).sum() > len(weights) / GRID + 2e-11:
        raise ArithmeticError("Excessive probability quantization")
    counts = rng.multinomial(STARTS, sampled)
    return dict(
        probability=original.tolist(),
        grid_weights=weights.tolist(),
        counts=list(map(int, counts)),
    )


def radial_moments(weights, centers, shells, connected, aperture):
    # A sharp support read makes the center and the topology available locally.
    occupancy = centers * weights[:, None]
    size = np.isin(source.DIST[0], aperture).sum()
    total = shells[:, :, aperture].sum(axis=2) / size
    return [
        (
            occupancy[subset].sum(),
            (occupancy[subset] * total[subset]).sum(),
            (occupancy[subset] * total[subset] ** 2).sum(),
        )
        for subset in (~connected, connected)
    ]


def probe_bank(rng, moments, g, bypass, acceptance):
    blocks = [p * detector(m / p, s / p, g, bypass) for p, m, s in moments]
    return acquire(rng, blocks, acceptance)


def calibration(rng, k):
    result = {}
    for n in range(3):
        result[str(n)] = acquire(rng, [detector(n - 1, (n - 1) ** 2, gain(k), True)], 1)
    return result


def panel(rng, name, private, standards):
    words, swaps = source.sector(*source.NAMES[name])
    mask, _ = source.native_geometry()
    kinetic, contrast, reference = source.operators(words, swaps, mask)
    rho, spec = source.filtered(kinetic, contrast, reference, KAPPA, private)
    weights = np.diag(rho).copy()
    del rho, kinetic, swaps
    if np.min(weights) < 0 or abs(weights.sum() - 1) > 2e-12:
        raise ArithmeticError("Invalid filtered source")
    med, shells, flags = source.graph_read(words)
    connected = np.asarray(flags)
    centers = {"intrinsic": med}
    for i in np.flatnonzero(~mask):
        centers[f"gateway_{i}"] = np.broadcast_to(np.eye(12)[i], med.shape)
    acceptance = spec["acceptance"]
    result = {}
    for center, selection in centers.items():
        for label, aperture in APERTURES.items():
            if label == "all" and center != "intrinsic":
                continue
            moments = radial_moments(weights, selection, shells, connected, aperture)
            banks = {}
            for law, k in [("native", KAPPA), ("exchange", 0)]:
                for mode in ("fringe", "elastic"):
                    bypass = mode == "fringe"
                    bank = probe_bank(rng, moments, gain(k), bypass, acceptance)
                    if bypass:
                        # Detector report receives only terminal port counts,
                        # admission and the independently stocked sharp standard.
                        size = int(np.isin(source.DIST[0], aperture).sum())
                        standard = standards[law]["2"]
                        bank["signed_shell_interval"] = decode(bank, standard, size)
                    banks[f"{law}_{mode}"] = bank
            if label == "all":
                # Recording the selected address kills cross terms. Its match
                # to the uniform return mode contributes another factor of 1/12.
                q2 = np.sum(words[0] ** 2) / 12**2
                incoherent = [(p, 0, p * q2) for p, _, _ in moments]
                banks["recorded_address"] = probe_bank(
                    rng, incoherent, gain(KAPPA), False, acceptance
                )
            result[f"{center}/{label}"] = banks
    return dict(acceptance=acceptance, settings=result)


def run():
    rng = np.random.default_rng(SEED)
    standards = {
        law: calibration(rng, k) for law, k in [("native", KAPPA), ("exchange", 0)]
    }
    panels = {
        label: panel(rng, name, private, standards)
        for label, name, private in (
            ("neutral", "neutral", True),
            ("positive", "positive", True),
            ("neutral_reference_shut", "neutral", False),
        )
    }
    return dict(
        starts=STARTS,
        grid=GRID,
        seed=SEED,
        routes={str(20 * 4**level): routes(level) for level in range(3)},
        outcomes={
            "fringe": ["plus", "minus", "other"],
            "elastic": ["match", "other"],
            "order": ["disconnected", "connected", "source_refusal"],
            "calibration": ["plus", "minus", "other", "source_refusal"],
        },
        standards=standards,
        panels=panels,
    )


def verify(got, saved):
    if isinstance(saved, dict):
        if got.keys() != saved.keys():
            raise ValueError("Receipt keys differ")
        for key, val in saved.items():
            verify(got[key], val)
    elif isinstance(saved, list):
        for a, b in zip(got, saved, strict=True):
            verify(a, b)
    elif isinstance(saved, float):
        if not np.isclose(got, saved, rtol=2e-9, atol=2e-11):
            raise ValueError("Computed probability or interval differs")
    elif got != saved:
        raise ValueError("Recorded count or setting differs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    path = parser.parse_args().output
    record = run()
    if path:
        with path.open("x") as out:
            json.dump(record, out, indent=1, allow_nan=False)
            out.write("\n")
    else:
        verify(record, json.loads((HERE / "results.json").read_text()))
        print("Balanced scattering receipt verified")
