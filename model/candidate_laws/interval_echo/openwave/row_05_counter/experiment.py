"""Joint angular readouts on two native five-port rings with guarded compositions."""

import argparse
import json
from pathlib import Path
import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import belt
from model.executable_model import carrier

HERE = Path(__file__).parent
D = 5
F = np.exp(2j * np.pi * np.outer(np.arange(D), np.arange(D)) / D) / np.sqrt(D)
TIMES = np.arange(129) / 8
SHOTS = 2**22
GEOMETRIES = ("all", "one", "cut")


def board_h(probe, strength, geometry):
    car = carrier.cell()
    order = next(iter(belt.belts().values()))
    left, right = order[::2], order[1::2]
    words = []
    for a in left:
        for b in right:
            occ = np.ones(12, dtype=int)
            occ[a] = 0
            occ[b] = 1 + probe
            words.append(tuple(occ))
    loc = {w: i for i, w in enumerate(words)}
    law = candidate.native_law
    h_local = law.seam_term(strength)
    rims = [e for e in car.seams if set(e) <= set(left) or set(e) <= set(right)]
    links = [
        e for e in car.seams if len(set(e) & set(left)) == len(set(e) & set(right)) == 1
    ]
    assert len(links) == len(rims) == 10
    active = rims + {"all": links, "one": links[:1], "cut": []}[geometry]
    h = np.zeros((D**2, D**2))
    for i, w in enumerate(words):
        h[i, i] = strength * np.dot(6 - car.degree, np.square(w))
        for a, b in active:
            p, q = w[a], w[b]
            k = law.SRC.index((p, q))
            h[i, i] += h_local[k, k]
            if p == q:
                continue
            z = list(w)
            z[a], z[b] = q, p
            j = loc.get(tuple(z))
            if j is not None:
                h[j, i] += h_local[law.SRC.index((q, p)), k]
    return h, left, right, links


def eigenspaces(h):
    e, v = np.linalg.eigh(h)
    indices = np.split(np.arange(len(e)), np.flatnonzero(np.diff(e) > 1.0e-9) + 1)
    return [(float(e[g].mean()), v[:, g] @ v[:, g].conj().T) for g in indices]


def rotate(x):
    return F.conj().T @ x.reshape(D, D) @ F.conj()


def joint(a):
    left, right = a[:, 1], a[:, -1]
    miss = (abs(a[:, [0, 2, 3]]) ** 2).sum(axis=1)
    return np.array(
        [
            np.stack(
                (abs(left + c * right) ** 2 / 2, abs(left - c * right) ** 2 / 2, miss),
                axis=-1,
            )
            for c in (1, 1j)
        ]
    )


def draw(probs, rng):
    p = np.clip(probs.ravel().round(12), 0, 1)
    major = p.argmax()
    order = [i for i in range(len(p)) if i != major] + [major]
    return rng.multinomial(SHOTS, p[order])[np.argsort(order)].reshape(probs.shape)


def census(h, levels):
    # The two-ring fifth-turn commutes with the all-contact Hamiltonian.
    form = np.kron(F, F)
    hf = form.conj().T @ h @ form
    sectors = (np.arange(25) // 5 + np.arange(25) % 5) % 5
    assert np.max(abs(hf) * (sectors[:, None] != sectors[None, :])) < 3e-10
    rows = []
    for m in range(D):
        sel = np.flatnonzero(sectors == m)
        for energy, project in eigenspaces(hf[np.ix_(sel, sel)]):
            weights = {
                f"{x},{y}": float(
                    project[
                        np.where(sel == (x % 5 * 5 + y % 5))[0][0],
                        np.where(sel == (x % 5 * 5 + y % 5))[0][0],
                    ].real
                )
                for x in (1, -1)
                for y in (1, -1)
                if (x + y) % 5 == m
            }
            rows.append(
                dict(
                    total_m=m,
                    energy=energy,
                    multiplicity=int(round(np.trace(project).real)),
                    weights=weights,
                )
            )
    # The coherent probe difference between the two source signs is a sum over energy gaps.
    terms = {}
    standing = (F[:, 1] + F[:, -1]) / np.sqrt(2)
    for sign in (1, -1):
        input_state = np.kron(F[:, sign], standing)
        projected = [rotate(p @ input_state) for _, p in levels]
        for (ep, _), a in zip(levels, projected):
            for (eq, _), b in zip(levels, projected):
                gap = round(ep - eq, 9)
                terms[gap] = terms.get(gap, 0j) + 2 * sign * np.vdot(b[:, -1], a[:, 1])
    poles = [
        dict(gap=g, coefficient=[c.real, c.imag])
        for g, c in sorted(terms.items())
        if abs(c) > 1e-9
    ]
    return dict(levels=rows, odd_poles=poles)


def run():
    bank = {}
    cases = []
    rng = np.random.default_rng(511241)
    for probe in (-1, 1):
        for k in (candidate.KAPPA, 0):
            for geom in GEOMETRIES:
                h, left, right, links = board_h(probe, k, geom)
                h -= h[0, 0] * np.eye(len(h))
                levels = eigenspaces(h)
                shared = (F[:, 1] + F[:, -1]) / np.sqrt(2)
                params = []
                retention = []
                coherences = []
                records = []
                for sign in (1, -1):
                    isometry = np.kron(F[:, sign, None], np.eye(D))
                    params.append(isometry.conj().T @ h @ isometry)
                    src = isometry @ shared
                    waves = np.exp(
                        -1j * np.outer(TIMES, [v for v, _ in levels])
                    ) @ np.array([p @ src for _, p in levels])
                    modes = np.array([rotate(x) for x in waves])
                    retention.append(np.sum(abs(modes[:, sign, :]) ** 2, axis=-1))
                    c = 2 * np.sum(modes[:, :, 1] * modes[:, :, -1].conj(), axis=-1)
                    coherences.append(c)
                    reads = np.empty((len(TIMES), 2, D, 3), np.int64)
                    for itime, a in enumerate(modes):
                        for axis, p in enumerate(joint(a)):
                            reads[itime, axis] = draw(p, rng)
                    records.append(reads)
                key = f"{probe}_{int(k)}_{geom}"
                counts = np.array(records)
                bank[key] = counts
                tot = counts.sum(axis=3)
                obs = (tot[:, :, :, 0] - tot[:, :, :, 1]) / SHOTS
                delta = (np.array(coherences)[0] - np.array(coherences)[1]).imag
                retained = np.min(retention, axis=0)
                observed = obs[0, :, 1] - obs[1, :, 1]

                def peak(eligible):
                    # Symmetry-related peaks and exact nulls can differ at roundoff.
                    magnitude = np.abs(delta[eligible])
                    i = eligible[np.flatnonzero(magnitude >= magnitude.max() - 1.e-10)[0]]
                    return dict(
                        time=float(TIMES[i]),
                        exact=float(delta[i]),
                        observed=float(observed[i]),
                        both_source_retention=float(retained[i]),
                        visibility=[float(abs(c[i])) for c in coherences],
                    )

                cases.append(
                    dict(
                        probe=probe,
                        contrast=k,
                        geometry=geom,
                        key=key,
                        peak=peak(np.arange(len(TIMES))),
                        retained_peak=peak(np.flatnonzero(retained >= 0.95 - 1e-12)),
                        fixed_source_difference=float(
                            np.linalg.norm(params[0] - params[1])
                        ),
                        spectral=census(h, levels) if geom == "all" else None,
                    )
                )
    return dict(
        shots=SHOTS,
        times=TIMES.tolist(),
        rings=[left, right],
        contacts=links,
        simultaneous_error=float(
            np.sqrt(4 * np.log(2 * len(cases) * len(TIMES) / 0.01) / SHOTS)
        ),
        cases=cases,
    ), bank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data, bank = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
        np.savez_compressed(HERE / "counts.npz", **bank)
    else:
        compare(json.loads(json.dumps(data)), json.loads(path.read_text()), "counter")
        with np.load(HERE / "counts.npz") as stored:
            assert set(stored) == set(bank)
            for key, n in bank.items():
                np.testing.assert_array_equal(n, stored[key])
    print(
        "Verified two-ring stationary source controls, energy sectors and all joint readings."
    )


if __name__ == "__main__":
    main()
