"""Apply a position-dependent scalar impulse to the same parked vacancy and a mobile control."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (  # noqa: E402
    anchor_spectrum,
    rest_partner as rest,
)


def run():
    board, route = rest.guide.longest_unpinned_route()
    patch = rest.Patch(board, route, min(rest.geometries(board, route)))
    columns, pair = patch.local_space()
    source = columns[:, 0]
    parked = columns @ pair
    distance = shortest_path(
        board.adjacency(), indices=patch.hub, directed=False, unweighted=True
    )[patch.ports]
    outside = patch.ports != patch.hub
    outgoing = patch.bridge_flux != 0
    if np.any(outgoing & ~outside):
        raise ValueError("Outgoing bridge head lies at the initial hub")
    effect, masses = rest.guide.register_effect()
    # The light selector acts only on the family record. Its accepted spatial mode is the same for all characters.
    # Recover its setting from archived discovery counts, then recompute its yield from the source.
    archive = json.loads((HERE / "anchor_spectrum.json").read_text())["discovery"]
    fits, _, _ = anchor_spectrum.infer(
        anchor_spectrum.unpack(archive["survey_counts"]),
        anchor_spectrum.unpack(archive["refinement_counts"]),
    )
    omegas, weights, _, _ = anchor_spectrum.instrument(effect)
    for reject in fits[1:]:
        t = np.pi / (reject - fits[0])
        weights = weights * np.cos(t * (omegas - fits[0]) / 2) ** 2
    success = weights.sum()
    rng = np.random.default_rng(178032)
    shots = 2**19
    records = []
    selected = {}
    for i, m in enumerate(masses):
        h = rest.guide.generator(patch.b, np.array([[m]]))
        np.testing.assert_allclose(h @ parked, -m * parked, atol=2e-13)
        np.testing.assert_allclose(parked[outside], 0, atol=2e-13)
        e, v = eigh(h)

        def evolve(state, time):
            return v @ (np.exp(-1j * time * e)[:, None] * (v.T @ state))

        initial = evolve(np.column_stack([parked, source]), 1.0)
        for strength in (0, np.pi / 4, -np.pi / 4):
            impulse = np.exp(-1j * strength * distance)
            np.testing.assert_allclose(impulse[:, None] * parked, parked, atol=2e-13)
            for t in (0.5, 2, 4):
                wave = evolve(impulse[:, None] * initial, t)
                probs = np.array(
                    [
                        np.sum(abs(wave[k]) ** 2, axis=0)
                        for k in (~outside, outside & ~outgoing, outgoing)
                    ]
                )
                np.testing.assert_allclose(probs.sum(axis=0), 1, atol=2e-13)
                np.testing.assert_allclose(
                    probs[:, :2], [[1, 1], [0, 0], [0, 0]], atol=3e-13
                )
                row = {
                    "character": i,
                    "impulse": float(strength),
                    "read_wait": t,
                    "vertex_control_probabilities": probs[:, 2].tolist(),
                    "vertex_control_counts": rng.multinomial(
                        shots, np.round(probs[:, 2], 13)
                    ).tolist(),
                    "max_parked_departure_probability": float(
                        probs[1:, :2].sum(axis=0).max()
                    ),
                }
                selected.setdefault((strength, t), np.zeros(3))[:] += (
                    weights[i] * probs[:, 0]
                )
                records.append(row)
    mixed = []
    for (strength, t), probs in selected.items():
        outcome = np.r_[probs, 1 - success]
        mixed.append(
            {
                "impulse": float(strength),
                "read_wait": t,
                "counts_hub_elsewhere_outgoing_failed": rng.multinomial(
                    shots, outcome / outcome.sum()
                ).tolist(),
            }
        )
    return {
        "shots": shots,
        "seed": 178032,
        "hub": patch.hub,
        "bridge_port": patch.port,
        "light_filter_yield": float(success),
        "acquisitions": records,
        "unlabeled_selected_source": mixed,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    with args.output.open("x") as f:
        json.dump(run(), f, indent=2, allow_nan=False)
        f.write("\n")
