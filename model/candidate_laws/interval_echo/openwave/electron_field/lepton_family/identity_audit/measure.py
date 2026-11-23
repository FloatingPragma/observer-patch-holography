"""Replay the complete five-pointer charge record on the seven-cell Dirac comparison."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.dirac_read import measure as source

TIMES = (1, 2, 4, 8, 12)
SEED = 396121
SHOTS = 64


def distal_readers(route, sign_sector):
    """Build the receiving pointer independently from material address and seam-head records."""
    ends = list(route) + list(route[:-1])
    receiver = route[-1] // 12
    material = np.array([address // 12 == receiver for address in ends])
    region = {p for p in route if p // 12 == receiver}
    reference = sum(v for p, v in zip(route, sign_sector["reference_head_counts_minus_charge"]) if p in region)
    edges = [sorted(e) for e in zip(route, route[1:])]
    words = sign_sector["written_words_in_port_order"]
    net = []
    for word in words + words[:-1]:
        inward = sum((a if word & (1 << i) else b) in region for i, (a, b) in enumerate(edges))
        net.append(inward - reference)
    flux = np.array(net)
    np.testing.assert_array_equal(flux, sign_sector["sign"] * material)
    return material, flux


def propagators(matrix):
    levels, states = np.linalg.eigh(matrix)
    intervals = np.diff((0, *TIMES))
    return [(states * np.exp(-1j * dt * levels)) @ states.conj().T for dt in intervals]


def record_law(pulses, inside, initial):
    """Each column collects an unnormalized conditional state for one pointer word."""
    amplitudes = initial.astype(complex)[:, None]
    for advance in pulses:
        amplitudes = advance @ amplitudes
        outside = (~inside)[:, None] * amplitudes
        occupied = inside[:, None] * amplitudes
        amplitudes = np.stack([outside, occupied], axis=2).reshape(len(initial), -1)
    records = (abs(amplitudes)**2).sum(axis=0)
    np.testing.assert_allclose(records.sum(), 1, atol=3e-13)
    return records


def score(laws):
    pairwise = []
    for a, b in combinations(range(len(laws)), 2):
        tv = float(abs(laws[a] - laws[b]).sum() / 2)
        pairwise.append({"sources": [a, b], "total_variation": tv,
                         "equal_prior_correct_from_one_word": (1 + tv)/2})
    return {"pairwise_records": pairwise,
            "three_way_equal_prior_correct_from_one_word": float(laws.max(axis=0).sum()/len(laws))}


def sampled(laws, rng):
    items = []
    for law in laws:
        codes = rng.choice(len(law), size=SHOTS, p=law/law.sum())
        words = [format(code, f"0{len(TIMES)}b") for code in codes.tolist()]
        items.append(words)
    return items


def measure():
    board, route = source.longest_unpinned_route()
    bonds = source.path_incidence(board, route)
    effect, _ = source.register_effect()
    masses, characters = np.linalg.eigh(effect)
    full = source.generator(bonds, effect)
    n = len(route)
    space = n + len(bonds)
    recv_minus, flux_minus = distal_readers(route, source.head_paths(route)[0])
    recv_plus, flux_plus = distal_readers(route, source.head_paths(route)[1])
    np.testing.assert_array_equal(recv_plus, recv_minus)
    np.testing.assert_array_equal(flux_plus, -flux_minus)
    inside = np.repeat(recv_minus, 3)
    assert recv_minus.sum() == 1 and route[-1] // 12 != route[0] // 12
    pulses = propagators(full)
    scalar_laws = {}
    gram_errors = []
    eigenspectrum = np.linalg.eigvalsh(full)
    kinetics = np.linalg.eigvalsh(bonds.T @ bonds)
    np.testing.assert_array_equal(bonds @ np.ones(n), 0)
    predicted_first = 2*np.sqrt(source.HOP)*np.sin(np.pi/(2*n))
    np.testing.assert_allclose(kinetics[1], predicted_first**2, atol=1e-14)
    assert 0 < masses[0] < masses[1] < masses[2] < predicted_first
    multiplicities = [int(np.sum(abs(eigenspectrum-f) < 1e-9)) for f in masses]
    assert multiplicities == [1, 1, 1]
    # Complex conjugation in the real scalar blocks pulls back to U K.
    # Its choice of phase in each nondegenerate character is conventional.
    real_structure = np.kron(np.eye(space), characters @ characters.T)
    np.testing.assert_allclose(real_structure @ full.conj() @ real_structure.conj().T, full, atol=1e-13)
    np.testing.assert_allclose(real_structure @ real_structure.conj(), np.eye(len(full)), atol=1e-13)
    spatial_sources = {"uniform_rest": np.r_[np.ones(n)/np.sqrt(n), np.zeros(len(bonds))],
                       "sharp_at_first_port": np.eye(space)[0]}
    for tag, spatial in spatial_sources.items():
        rows = []
        for j, mass in enumerate(masses):
            # Intertwiner covers every state of this orbital, not just the given source.
            lift = np.kron(np.eye(space), characters[:,j:j+1])
            orbital = source.generator(bonds, np.array([[mass]]))
            gram_errors.append(float(np.linalg.norm(full @ lift - lift @ orbital, 2)))
            np.testing.assert_allclose(inside[:,None]*lift, lift*recv_minus[None,:], atol=1e-14)
            complete = record_law(pulses, inside, lift @ spatial)
            alternative = record_law(propagators(orbital), recv_minus, spatial)
            np.testing.assert_allclose(complete, alternative, atol=3e-13)
            rows.append(complete)
        scalar_laws[tag] = np.array(rows)
    assert max(gram_errors) < 1e-12
    # A constant formation cost in the one-charge sector leaves a phase on
    # every trajectory. Its neutral sector can independently retain zero energy.
    formation_cost = 1.
    shift_pulses = propagators(full + formation_cost*np.eye(len(full)))
    offset_error = 0.
    for tag, spatial in spatial_sources.items():
        for j, v in enumerate(characters.T):
            shifted = record_law(shift_pulses, inside, np.kron(spatial, v))
            offset_error = max(offset_error, float(abs(shifted-scalar_laws[tag][j]).max()))
    assert offset_error < 3e-13
    # At rest all characters share the same stationary real-space distribution.
    for row in scalar_laws["uniform_rest"]:
        np.testing.assert_allclose(row[len(row)//2:].sum(), 1/n, atol=2e-13)
    mean = np.trace(effect).real/len(effect)
    flat_pulses = propagators(source.generator(bonds, mean*np.eye(3)))
    flat = np.array([record_law(flat_pulses, inside, np.kron(spatial_sources["uniform_rest"], v))
                     for v in characters.T])
    np.testing.assert_allclose(flat, np.repeat(flat[:1], len(masses), axis=0), atol=2e-13)

    projectors = [np.kron(np.eye(space), np.outer(v, v.conj())) for v in characters.T]
    zero = max(np.linalg.norm(full @ projector - projector @ full, 2) for projector in projectors)
    assert zero < 1e-13
    for projector in projectors:
        np.testing.assert_allclose(inside[:,None]*projector, projector*inside[None,:], atol=1e-15)
    # The same detector registers a supplied local conversion on its occupied branch.
    z = np.diag(np.exp(2j*np.pi*np.arange(3)/3))
    addr_kick = [z, z.conj().T]
    kick = max(addr_kick, key=lambda op: abs(np.vdot(characters[:,0], op @ characters[:,2])))
    np.testing.assert_allclose(abs(np.vdot(characters[:,0], kick @ characters[:,2])), 1., atol=1e-13)
    parent = np.kron(spatial_sources["uniform_rest"], characters[:,2])
    kicked = parent.reshape(space, 3).copy()
    kicked[recv_minus] = (kick @ kicked[recv_minus].T).T
    kicked = kicked.ravel()
    daughter = float(np.linalg.norm(inside * (projectors[0] @ kicked))**2)
    work = float(np.vdot(kicked, full @ kicked).real - np.vdot(parent, full @ parent).real)
    np.testing.assert_allclose(daughter, 1/n, atol=2e-14)
    np.testing.assert_allclose(work, (masses[0] - masses[2])/n, atol=2e-14)
    # A spin spectator couples identically in these instruments for every spin density matrix.
    test_spin = np.array([1., 1j]) / np.sqrt(2)
    spin_lift = np.kron(np.eye(space), test_spin[:,None])
    scalar = source.generator(bonds, np.array([[masses[0]]]))
    np.testing.assert_allclose(np.kron(scalar,np.eye(2)) @ spin_lift, spin_lift @ scalar, atol=1e-14)
    rng = np.random.default_rng(SEED)
    return {"fixture": "incidence Dirac comparison with the historical cyclic read as mass operator",
            "physical_path_ports": list(route), "receiver_cell": route[-1]//12,
            "signs": [-1, 1], "native_energy_levels": masses.tolist(),
            "rest_rank_by_family_and_sign": multiplicities, "first_nonuniform_spatial_singular_value": predicted_first,
            "complete_one_sign_dimension": len(eigenspectrum),
            "real_block_conjugation_squared": 1,
            "rest_doublet_under_antiunitary_squaring_to_minus_one_preserving_H_and_sign": False,
            "largest_scalar_intertwining_residual": max(gram_errors),
            "max_family_commutator_norm": float(zero), "time_in_native_units": list(TIMES),
            "pointer_words_in_order": [format(i,"05b") for i in range(2**len(TIMES))],
            "pointer_seed": SEED, "fresh_trials_per_input": SHOTS,
            "sources": {tag: {"complete_word_probabilities": laws.tolist(), "scores": score(laws),
                              "simulated_pointer_words": sampled(laws,rng)} for tag,laws in scalar_laws.items()},
            "equal_mass_control": {"input_independent_uniform_rest_word_probabilities": flat[0].tolist(),
                                   "three_way_equal_prior_correct": score(flat)["three_way_equal_prior_correct_from_one_word"]},
            "constant_charged_formation_cost_control": {
                "neutral_energy_change": 0, "one_charge_energy_change": formation_cost,
                "modified_ratios": ((masses + formation_cost)/(masses[0]+formation_cost)).tolist(),
                "max_observer_word_probability_error": offset_error},
            "rest_one_time_distal_probability_at_every_wait": 1/n,
            "driven_distal_heavy_to_light_probability": daughter, "driven_mean_work": work,
            "undriven_heavy_to_light_probability_all_times_and_charge_reads": 0,
            "source_sha256": sha256(Path(source.__file__).read_bytes()).hexdigest(),
            "code_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--output", type=Path, required=True)
    path = args.parse_args().output
    report = measure()
    with path.open("x") as f:
        json.dump(report, f, indent=2, allow_nan=False)
        f.write("\n")


if __name__ == "__main__":
    main()
