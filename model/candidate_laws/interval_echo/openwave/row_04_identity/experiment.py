"""Compare the free row-4 preparation under a guarded SU(3) centre word and sham."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate, currents
from model.candidate_laws.interval_echo.openwave.row_04.__main__ import (
    evolve,
    verify_local_reader,
)
from model.executable_model import cell, federation
from model.experiments.interval_echo_gauge.__main__ import (
    coordinates,
    ordered_word,
    portable,
)


HERE = Path(__file__).resolve().parent
EXPOSURES = (0, 32, 150)
ARMS = ("sham", "centre", "centre_and_inverse")
SHOTS = 262144
STEPS = 2048
WORD_TOLERANCE = 0.0001


def calibrated_words():
    local, sources, _, _, bands, _, _, _, addresses, _ = currents._service(
        cell(), True, 1.0
    )
    embedding = np.column_stack(bands)
    axis = bands[0].T @ sources[:, 0]
    axis /= np.linalg.norm(axis)
    generator = block_diag(np.eye(3) / 3 - np.outer(axis, axis), np.zeros((3, 3)))
    coefficients = np.linalg.lstsq(
        coordinates(addresses).T, coordinates([generator])[0], rcond=1e-12
    )[0]
    if (
        np.linalg.norm(np.einsum("i,ijk->jk", coefficients, addresses) - generator)
        > 1e-10
    ):
        raise ArithmeticError(
            "The native address expansion failed to reproduce the traceless axis"
        )

    turn = ordered_word(addresses, coefficients, 2 * np.pi, STEPS)
    inverse = ordered_word(addresses, coefficients, -2 * np.pi, STEPS)
    phase = np.exp(2j * np.pi / 3)
    ideal = block_diag(phase * np.eye(3), np.eye(3))
    word_error = np.linalg.norm(turn - ideal, 2)
    undo_error = np.linalg.norm(inverse @ turn - np.eye(6), 2)
    if max(word_error, undo_error) > WORD_TOLERANCE:
        raise ArithmeticError("The guarded pulse does not resolve the proposed centre")

    # Reconstruct both full operators independently from the colour turn and central current.
    centre_projector = embedding @ ((ideal - np.eye(6)) / (phase - 1)) @ embedding.T
    primitive_abelian = embedding @ (np.sum(addresses, axis=0) / 2) @ embedding.T
    support_error = np.linalg.norm(centre_projector - primitive_abelian, 2)
    spectrum = np.linalg.eigvalsh(primitive_abelian)
    if support_error > 1e-10 or not np.allclose(
        spectrum, [0.0] * 63 + [1.0] * 3, atol=1e-10
    ):
        raise ArithmeticError(
            "The proposed colour/abelian support relation fails on this cell"
        )
    audit = dict(
        port_coefficients=coefficients,
        symmetric_steps=STEPS,
        centre_operator_distance=word_error,
        undo_operator_distance=undo_error,
        centre_projector_minus_primitive_abelian_norm=support_error,
        primitive_abelian_rank=int((spectrum > 0.5).sum()),
        local_pair_dimension=local.P,
    )
    return local, embedding, (np.eye(6), turn, inverse @ turn), audit


def run():
    verify_local_reader()
    local, embedding, words, audit = calibrated_words()
    screen = federation(0)
    sector = PairSector(screen)
    home = np.array(
        [sector.index(i, j) for i, j in zip(local.aa, local.bb, strict=True)]
    )
    outside = np.ones(sector.P, dtype=bool)
    outside[home] = False
    source = np.zeros(sector.P)
    source[sector.index(0, 1)] = 1.0
    states = evolve(
        sector.hamiltonian(candidate.KAPPA), candidate.KAPPA, source, EXPOSURES
    )
    baseline = json.loads((HERE.parent / "row_04/results.json").read_text())["boards"][
        0
    ]
    expected = {r["time"]: r["home_internal_seam"] for r in baseline["arms"]["sharp"]}
    if baseline["original_seam"] != [0, 1] or baseline["home_internal_addresses"] != [
        list(e) for e in cell().seams
    ]:
        raise ValueError(
            "Released addresses differ from the archived sharp source and reader"
        )
    n_outcomes = len(local.seam_states) + 1
    count_tolerance = np.sqrt((n_outcomes * np.log(2) + np.log(6e6)) / SHOTS)
    cutoff = count_tolerance + WORD_TOLERANCE
    trials = []
    for t, state in zip(EXPOSURES, states, strict=True):
        local_pair = state[home]
        unselected = np.linalg.norm(state[outside]) ** 2
        kept = embedding.T @ local_pair
        arms = {}
        for arm_index, (arm, word) in enumerate(zip(ARMS, words, strict=True)):
            result = local_pair + embedding @ ((word - np.eye(6)) @ kept)
            occupation = abs(result) ** 2
            probabilities = np.r_[
                occupation[local.seam_states],
                unselected + occupation[~local.adjacent].sum(),
            ]
            if abs(probabilities.sum() - 1) > 1e-10 or min(probabilities) < 0:
                raise ArithmeticError("The patch instrument lost a preparation")
            # Remove floating point norm error for the multinomial sampler only.
            seed = [41937, t, arm_index]
            counts = np.random.default_rng(seed).multinomial(
                SHOTS, probabilities / probabilities.sum()
            )
            arms[arm] = dict(
                probabilities=probabilities, counts=counts.tolist(), seed=seed
            )
        sham = arms["sham"]
        if abs(sham["probabilities"][:-1].sum() - expected[t]) > 1e-9:
            raise ArithmeticError(
                "The sham changed the row-4 event on its original sharp preparation"
            )
        comparisons = {}
        for arm in ARMS[1:]:
            exact = np.abs(arms[arm]["probabilities"] - sham["probabilities"]).sum() / 2
            observed = np.abs(
                np.asarray(arms[arm]["counts"]) - sham["counts"]
            ).sum() / (2 * SHOTS)
            comparisons[arm] = dict(
                probability_tv=exact,
                count_tv=observed,
                resolves_difference=bool(observed > cutoff),
            )
        trials.append(
            dict(
                free_exposure=t,
                both_loads_in_home_probability=float(
                    np.vdot(local_pair, local_pair).real
                ),
                model_weight_in_first_local_triplet=float(
                    np.linalg.norm(kept[:3]) ** 2
                ),
                arms=arms,
                comparisons_to_sham=comparisons,
            )
        )
    primary = trials[-1]["comparisons_to_sham"]
    resolved = (
        primary["centre"]["resolves_difference"]
        and not primary["centre_and_inverse"]["resolves_difference"]
    )
    return portable(
        dict(
            candidate=candidate.NAME,
            kappa=candidate.KAPPA,
            source=[0, 1],
            cells=screen.cells,
            active_seams_in_flight=len(screen.seams),
            clock="supplied native exposure",
            outcomes=[list(e) for e in cell().seams] + ["quiet"],
            shots_per_arm=SHOTS,
            per_family_type_i_bound=1e-6,
            decision_cutoff_tv=cutoff,
            pulse_audit=audit,
            trials=trials,
            centre_invariance_rejected_at_150=resolved,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, help="Write a fresh receipt to a new file"
    )
    output = parser.parse_args().output
    data = run()
    if output:
        with output.open("x") as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write("\n")
    elif data != json.loads((HERE / "results.json").read_text()):
        raise SystemExit("Guarded colour response differs from the archived receipt")
    for trial in data["trials"]:
        print(trial["free_exposure"], trial["comparisons_to_sham"])
    print(
        "Primary centre invariance rejected:", data["centre_invariance_rejected_at_150"]
    )
