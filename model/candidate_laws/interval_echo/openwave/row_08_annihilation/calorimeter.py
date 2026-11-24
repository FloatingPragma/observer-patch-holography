"""Far-end material spectroscopy after monitored and unmonitored conversion."""

import argparse
import json
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply
from scipy import sparse
from . import experiment as ex, converter

SHOTS = 2**30
BINS = 2**20


def sample(prob, rng):
    p = np.maximum(0, np.asarray(prob))
    order = np.argsort(p.round(10), kind="stable")
    q = p[order].round(12)
    q[-1] = 1 - q[:-1].sum()
    counts = np.empty(len(q), np.int64)
    counts[order] = rng.multinomial(SHOTS, q)
    return dict(probabilities=p.round(12).tolist(), counts=counts.tolist())


def terminal_sensors(wing):
    e, v = eigh(wing)
    e = e - e[0]
    d = len(e)
    x = np.kron(np.eye(d // 2), np.array([[0, 1], [1, 0]]))
    z = np.diag([1 - 2 * (i % 2) for i in range(d)])
    choices = [("soft", e[1], z)]
    hard = sorted(
        {round(e[i] - e[j], 10) for i in range(d) if 16 < e[i] < 48 for j in (0, 1)}
    )
    choices += [("hard", delta, x) for delta in hard]
    effects = []
    settings = []
    for kind, target, op in choices:
        # Each delay conjugates a local unitary; the selected clock amplitude averages them.
        theta = 2 * np.pi / (BINS * e[1]) * (e[:, None] - e[None, :] + target)
        weight = (
            np.exp(0.5j * (BINS - 1) * theta)
            * np.sinc(0.5 * BINS * theta / np.pi)
            / np.sinc(0.5 * theta / np.pi)
        )
        kraus = weight * (v.T @ op @ v)
        m = kraus.conj().T @ kraus
        assert np.linalg.eigvalsh(m)[-1] < 1 + 1e-10
        effects.append(v @ m @ v.T)
        settings.append(
            dict(
                type=kind,
                frequency=float(target),
                responses_by_level=m.diagonal().real.round(12).tolist(),
            )
        )
    return effects, settings


def acquire(psi, wing, rng):
    d = len(wing)
    e, v = eigh(wing)
    psi = psi.reshape(-1, d, d)
    joint = np.sum(abs(v.T @ psi @ v) ** 2, axis=0)
    success = joint.sum()
    effects, settings = terminal_sensors(wing)
    paired = []
    for i, left in enumerate(effects):
        for j, right in enumerate(effects):
            p = ex.read_joint(psi, left, right)
            paired.append(dict(left=i, right=j, **sample(np.r_[p, 1 - success], rng)))
    return dict(
        wing_levels=(e - e[0]).round(12).tolist(),
        calibrations=settings,
        full_spectrum=sample(np.r_[joint.ravel(), 1 - success], rng),
        paired_detectors=paired,
    )


def run():
    rng = np.random.default_rng(81126)
    h, hs, _, _, g = ex.basis(8, 3)
    incoming, _ = ex.launch(3, 1, h, g)
    h, wing, _ = converter.system(8, 3, 0.5)
    d = len(g)
    psi = np.zeros(len(h.diagonal()))
    i = ex.PAIRS.index((0, 3))
    psi[i * d * d : (i + 1) * d * d] = np.kron(g, g)
    h = h - sparse.eye(h.shape[0]) * h.diagonal().mean()
    reaction = expm_multiply(-64j * h, psi).reshape(14, d, d)[12:]
    return dict(
        seed=81126,
        shots=SHOTS,
        clock_bins=BINS,
        registered=acquire(incoming, hs[1], rng),
        constant_converter=acquire(reaction, wing, rng),
        neutral_ground=acquire(np.outer(g, g), wing, rng),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = json.loads(json.dumps(run()))
    path = ex.HERE / "calorimeter.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from .verify import compare

        compare(data, json.loads(path.read_text()), "calorimeter")
        print(
            "Verified the complete two-wing spectra and locally controlled end-detector coincidences."
        )


if __name__ == "__main__":
    main()
