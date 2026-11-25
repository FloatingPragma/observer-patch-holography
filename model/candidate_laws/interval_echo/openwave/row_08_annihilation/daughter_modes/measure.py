"""Identify the guide normal modes and read joint coherence on the collision source."""

import argparse
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from .. import experiment, walls
from ..daughter_spectrum.measure import check, compare, reference

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "results.json"
SPECTRA = HERE.parent / "daughter_spectrum" / "results.json"
SHOTS, SEED = 2**23, 81169


def constants(kappa):
    matrix = experiment.term(kappa)
    law = experiment.LAW
    opposite, equal, turned = [law.SRC.index(x) for x in ((0, 2), (0, 0), (2, 0))]
    coupling = (matrix[opposite, opposite] - matrix[equal, equal])/2
    hopping = -matrix[opposite, turned]
    return float(coupling), float(hopping)


def frequencies(n, kappa, cut):
    coupling, hopping = constants(kappa)
    if cut:
        bonds = hopping*np.eye(n) + coupling*np.eye(n, k=1)
        return np.sort(2*np.linalg.svd(bonds, compute_uv=False))
    momenta = np.pi*np.arange(1, n+1)/(n+1)
    return 2*np.sqrt((coupling-hopping)**2 + 4*coupling*hopping*np.sin(momenta/2)**2)


def energy_table(n, kappa, cut):
    modes = frequencies(n, kappa, cut)
    occupations = np.array(list(product((0, 1), repeat=n)), dtype=int)
    return occupations @ modes, occupations


def decode_spectra():
    receipt = json.loads(SPECTRA.read_text())
    decoded = []
    for arm in receipt["arms"]:
        k, cut = arm["contrast"], arm["attachments_cut"]
        levels, words = energy_table(3, k, cut)
        check(np.min(np.diff(np.sort(levels))) > 1e-5, "spectral read is ambiguous")
        for panel in arm["readings"]:
            for bank in panel["spectra"]:
                by_number = np.zeros((4,4), dtype=np.int64)
                by_mode = np.zeros((3,3), dtype=np.int64)
                payload = bank["joint_spectrum"]["counts"]
                for event, count in zip(bank["spectrum_outcome_labels"], payload[:-1], strict=True):
                    pair = []
                    for key in ("left_energy", "right_energy"):
                        difference = abs(levels-event[key])
                        selected = int(np.argmin(difference))
                        check(difference[selected] < 3e-8, "frequency has no unique mode label")
                        pair.append(words[selected])
                    sizes = tuple(word.sum() for word in pair)
                    by_number[sizes] += count
                    if sizes == (1,1):
                        by_mode[tuple(word.argmax() for word in pair)] += count
                check(by_number.sum()+payload[-1] == receipt["shots"], "energy counts lost")
                decoded.append(dict(contrast=k, cut=cut, wall=bank["wall_source"],
                                    exposure=panel["exposure"], ring_miss=payload[-1],
                                    mode_number_counts=by_number.tolist(),
                                    one_each_mode_counts=by_mode.tolist()))
    return decoded


@lru_cache(maxsize=2)
def calibrated_modes(boundary):
    _, frame = eigh(experiment.wing(8, 3, boundary))
    outer_flip = np.kron(np.eye(4), [[0.,1.],[1.,0.]])
    selected = frame[:, 1:3]
    overlap = frame[:, 0] @ outer_flip @ selected
    check(np.min(abs(overlap)) > .1, "local phase calibration unresolved")
    return selected * np.sign(overlap)


def outcomes(state, left, right, axes, dephase):
    """Pooled nine terminal labels, including local rejects, for one fixed ring phase."""
    frames = [left, right]
    if axes == "XX" and not dephase:
        hadamard = np.array([[1.,1.],[1.,-1.]])/np.sqrt(2)
        frames = [f @ hadamard for f in frames]
    left, right = frames
    p = np.zeros((3,3))
    p[:2,:2] = abs(left.T @ state @ right)**2
    p[:2,2] = np.sum(abs(left.T @ state)**2, axis=1) - p[:2,:2].sum(axis=1)
    p[2,:2] = np.sum(abs(state @ right)**2, axis=0) - p[:2,:2].sum(axis=0)
    p[2,2] = np.linalg.norm(state)**2 - p.sum()
    if axes == "XX" and dephase:
        # Dephasing in each two-mode energy chart makes both X outputs fair.
        # This preserves the observed two-mode, single-mode and rejection weights.
        p[:2,:2] = p[:2,:2].sum()/4
        p[:2,2] = p[:2,2].sum()/2
        p[2,:2] = p[2,:2].sum()/2
    check(np.min(p) > -2e-13, "negative terminal branch")
    return p.clip(min=0)


def acquisition(wave, boundary, rng, dephase):
    reads = {}
    for axes in ("XX","ZZ"):
        weights = []
        for phase, state in enumerate(wave):
            frames = [calibrated_modes(int(b if not phase else 2-b)) for b in boundary]
            weights.extend(outcomes(state, *frames, axes, dephase).ravel().tolist())
        weights.append(1-float(np.linalg.norm(wave)**2))
        probs = np.asarray(weights).clip(min=0)
        check(abs(probs.sum()-1) < 2e-10, "incomplete local or central outcome")
        probs /= probs.sum()
        counts = rng.multinomial(SHOTS, probs)
        joint = counts[:-1].reshape(2,3,3).sum(axis=0)[:2,:2]
        hits = int(joint.sum())
        correlation = (float((joint.trace()-(hits-joint.trace()))/hits) if hits else None)
        reads[axes] = dict(all_19_counts=counts.tolist(), admitted=hits, correlation=correlation)
    # Two settings, twelve arms, two-sided Hoeffding bounds for variables in [-1,1].
    resolved = [r for r in reads.values() if r["admitted"]]
    lower = sum(max(0., abs(r["correlation"])-np.sqrt(2*np.log(48/.01)/r["admitted"]))
                for r in resolved)
    return dict(reads=reads, simultaneous_99_witness_lower=float(lower))


def run():
    boundary, starts, ordered = reference()
    rng = np.random.default_rng(SEED)
    panels = []
    for kappa, cut in ((8,False),(8,True),(0,False)):
        h = sparse.csr_matrix(walls.operators(kappa, cut, 3))
        step = -32j*(h - sparse.eye(len(starts))*np.mean(h.diagonal()))
        wave = expm_multiply(step, starts).reshape(16,8,8,2)
        for source, name in enumerate(("sham","wall")):
            for dephase in (False, True):
                reading = acquisition(wave[ordered,:,:,source], boundary, rng, dephase)
                panels.append(dict(contrast=kappa, cut=cut, source=name,
                                   local_energy_dephasing=dephase, **reading))
    coupling, hopping = constants(8)
    return dict(geometry="four-dimer source ring and two three-dimer receivers",
                seed=SEED, shots=SHOTS, exposure=32, standing_wavenumbers_in_pi=[.25,.5,.75],
                native_one_mode_frequencies=frequencies(3, 8, False).tolist(),
                native_uniform_bulk_gap=2*abs(coupling-hopping),
                native_uniform_bulk_max_group_speed=2*min(coupling,hopping),
                energy_counts=decode_spectra(), coherence=panels)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run()
    if args.output:
        encoded = json.dumps(report, indent=2, allow_nan=False)+"\n"
        with args.output.open("x") as stream:
            stream.write(encoded)
    else:
        compare(report, json.loads(ARCHIVE.read_text()), "daughter modes")
        print("Replayed daughter mode assignment and both joint-coherence settings.")
