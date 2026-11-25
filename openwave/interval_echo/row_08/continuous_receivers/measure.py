"""Terminal coincidences in one uninterrupted ring, wing and detector evolution."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.row_08_annihilation import walls
from model.candidate_laws.interval_echo.openwave.row_08_annihilation.passive_receiver import measure as meter

HERE = Path(__file__).resolve().parent
TIMES = (0, 8, 16, 32, 48, 64)
METER_LOAD = np.array([0, 1, 1, 2])
BASES = ("ordered_sham", "opposite_walls")


def generator(cut, omega, coupling):
    host = sparse.csr_matrix(walls.operators(8, cut, 3))
    size = host.shape[0] * 4
    h = sparse.kron(host, sparse.eye(4), format="csr")
    h += sparse.diags(omega * np.tile(METER_LOAD, host.shape[0]))
    # Ten dimer bits followed by the left and right receiver bits.
    # Flip the far dimer of one wing and its own receiver in the same term.
    idx = np.arange(size)
    for mask in ((1 << 5) | 2, (1 << 2) | 1):
        h += sparse.csr_matrix((np.full(size, coupling), (idx, idx ^ mask)), shape=h.shape)
    return h


def wall_charges(word):
    bits = [int(x) for x in f"{word:04b}"]
    order = np.array([(-1)**j * (2*b - 1) for j, b in enumerate(bits)])
    return ((np.roll(order, -1) - order) // 2).tolist()


def panel(state, ordered, omega):
    joint = (abs(state.reshape(16, 64, 4))**2).sum(axis=1)
    outcomes = joint.sum(axis=0)
    return dict(ring_word_by_receiver=joint.tolist(), receiver=outcomes.tolist(),
                both=float(outcomes[3]), ordered_both=float(joint[ordered, 3].sum()),
                ordered=float(joint[ordered].sum()),
                stored_detector_energy=float(omega * (outcomes @ METER_LOAD)))


def run():
    boundary, host, ordered = meter.parent.reference()
    levels = eigh(meter.ex.wing(8, 3, int(boundary[0])), eigvals_only=True)
    omega = float(levels[2] - levels[0])
    size = host.shape[0] * 4
    initial = np.zeros((size, 2))
    initial[::4] = host
    source_words = np.argmax(np.sum(host.reshape(16, 64, 2)**2, axis=1), axis=0)
    incoming = [wall_charges(word) for word in source_words]
    assert sorted(incoming[0]) == [0, 0, 0, 0] and sorted(incoming[1]) == [-1, 0, 0, 1]
    results = {}
    for label, cut, coupling in (("connected", False, meter.COUPLING),
                                 ("cut_source_links", True, meter.COUPLING),
                                 ("uncoupled_receivers", False, 0)):
        full = generator(cut, omega, coupling)
        if coupling == 0:
            # Prove the initially occupied detector block is invariant, then evolve it.
            assert full[np.arange(size) % 4 != 0][:, ::4].nnz == 0
            h, seed = full[::4, ::4], host
        else:
            h, seed = full, initial
        energy, vec = eigh(h.toarray(order="F"), overwrite_a=True, driver="evd")
        overlap = vec.T @ seed
        e0 = np.sum(seed * (h @ seed), axis=0)
        check_time = 1/32
        def evolve(time):
            return vec @ (np.exp(-1j * time * energy)[:, None] * overlap)
        independent = expm_multiply(-1j * check_time * h, seed)
        np.testing.assert_allclose(evolve(check_time), independent, rtol=0, atol=2e-11)
        if label == "connected":
            # Cross-check the primary endpoint using sparse propagation directly on
            # the two launch vectors, with no energy eigenbasis or intermediate read.
            center = h.diagonal().mean()
            shifted = h - sparse.eye(size, format="csr") * center
            direct = expm_multiply(-64j * shifted, seed, traceA=-64j * shifted.diagonal().sum())
            np.testing.assert_allclose(abs(direct)**2, abs(evolve(64))**2, rtol=1e-8, atol=2e-11)
        readings = {}
        for time in TIMES:
            flight = evolve(time)
            expectation = np.sum(flight.conj() * (h @ flight), axis=0).real
            np.testing.assert_allclose(expectation, e0, atol=2e-8, rtol=0)
            np.testing.assert_allclose(np.sum(abs(flight)**2, axis=0), 1, atol=2e-11, rtol=0)
            if coupling == 0:
                expanded = np.zeros_like(initial, dtype=complex)
                expanded[::4] = flight
                flight = expanded
            readings[str(time)] = {name: panel(flight[:, i], ordered, omega)
                                   for i, name in enumerate(BASES)}
        results[label] = dict(initial_energy=e0.tolist(), observations=readings)
    for time in TIMES:
        for name in BASES:
            assert results["uncoupled_receivers"]["observations"][str(time)][name]["stored_detector_energy"] == 0
        cut = results["cut_source_links"]["observations"][str(time)]
        np.testing.assert_allclose(cut[BASES[0]]["receiver"], cut[BASES[1]]["receiver"],
                                   atol=1e-10, rtol=0)
    primary = results["connected"]["observations"]["64"]
    diffs = {k: primary[BASES[1]][k] - primary[BASES[0]][k]
             for k in ("both", "ordered_both", "stored_detector_energy")}
    excess = float(np.diff(results["connected"]["initial_energy"]).item())
    diffs["fraction_of_initial_excess_energy_in_receivers"] = diffs["stored_detector_energy"] / excess
    return dict(exposures=TIMES, primary_exposure=64, meter_order=["00", "01", "10", "11"],
                ring_words=[f"{w:04b}" for w in range(16)], ordered_words=ordered,
                incoming_wall_charges=dict(zip(BASES, incoming)), receiver_gap=omega,
                coupling=meter.COUPLING, initial_excess_energy=excess,
                primary_wall_minus_sham=diffs,
                passes_uninterrupted_transfer=bool(min(diffs["both"], diffs["ordered_both"]) > 1e-10),
                arms=results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        # JSON normalizes the exposure tuple to a list before comparing.
        saved = json.loads(path.read_text())
        meter.parent.compare(json.loads(json.dumps(result)), saved)
        for arm, data in result["arms"].items():
            for t, by_source in data["observations"].items():
                for name, obs in by_source.items():
                    other = saved["arms"][arm]["observations"][t][name]
                    np.testing.assert_allclose(obs["ring_word_by_receiver"], other["ring_word_by_receiver"],
                                               rtol=1e-8, atol=2e-13)
    print(json.dumps(result["primary_wall_minus_sham"], indent=2))
