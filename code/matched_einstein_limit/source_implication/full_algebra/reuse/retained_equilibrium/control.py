"""Construct the retained-head controller and test its spectral conditions."""

from itertools import permutations

import numpy as np
from scipy.optimize import linprog

from .probe import density_center
from . import certificate, gates

TWIRL = 1 / 8
I4 = np.eye(4)
PERMS = np.array(list(permutations(range(4))))


def eigen(a):
    w, v = np.linalg.eigh((a + a.conj().T) / 2)
    return w[::-1], v[:, ::-1]


def spectrum_checked(a):
    """Weyl enclosure from a dyadic eigenframe, including its polar correction."""
    q = 2**34

    def grid(v):
        return np.rint(v * q).astype(np.int64).astype(object)

    w, v = eigen(a)
    ar, ai = grid(a.real), grid(a.imag)
    vr, vi = grid(v.real), grid(v.imag)
    d = grid(w)
    np.testing.assert_array_equal(ar, ar.T)
    np.testing.assert_array_equal(ai, -ai.T)
    rr = ar @ vr - ai @ vi - vr * d
    ri = ar @ vi + ai @ vr - vi * d
    gr = vr.T @ vr + vi.T @ vi - np.eye(4, dtype=object) * q**2
    gi = vr.T @ vi - vi.T @ vr
    for x, y in ((rr, ri), (gr, gi)):
        assert (np.sum(x * x) + np.sum(y * y)) * 10**16 < q**4
    assert np.sum(ar * ar) + np.sum(ai * ai) <= q * q and max(abs(d)) <= q
    # ||V-polar(V)||<=||V*V-I||. Both matrix norms are <=1, giving
    # an error below 4e-8 per eigenvalue including the matrix grid.
    return w


def state_outputs(blocks, x):
    return np.einsum("hij,j->hi", blocks, x.ravel()).reshape(-1, 4, 4)


def slack(blocks, state):
    eigenvalues = sum(spectrum_checked(z) for z in blocks)
    return (eigenvalues - spectrum_checked(state)).cumsum()[:3]


def calibrate(blocks, state):
    """Choose one physical program among the admissible eigenvector phases and lotteries."""
    pieces = state_outputs(blocks, state)
    mux, vx = eigen(state)
    eigenpairs = [eigen(z) for z in pieces]
    lam = sum(w for w, _ in eigenpairs)
    flat = state.trace().real / 4
    desired = (mux - TWIRL * flat) / (1 - TWIRL)
    # Permuting already aligned spectra implements majorization of their sum.
    system = np.vstack((np.ones(24), lam[PERMS].T[:3]))
    fit = linprog(
        np.zeros(24),
        A_eq=system,
        b_eq=np.r_[1, desired[:3]],
        bounds=(0, None),
        method="highs",
    )
    assert fit.success
    return np.array([vx] + [v for _, v in eigenpairs]), np.asarray(fit.x)


def build(blocks, state, frames, trust):
    """Replay the deposited program, including its eigenvector phases."""
    vx = frames[0]
    pieces = state_outputs(blocks, state)
    for frame, y in zip(frames, [state, *pieces], strict=True):
        diag = frame.conj().T @ y @ frame
        np.testing.assert_allclose(diag, np.diag(eigen(y)[0]), atol=2e-13, rtol=0)
        np.testing.assert_allclose(frame.conj().T @ frame, I4, atol=2e-13, rtol=0)
    cz = gates.native_cz()
    packed = np.zeros_like(blocks)
    compiler_error = 0.0
    for k, (basis, phi) in enumerate(zip(frames[1:], blocks, strict=True)):
        restored = np.zeros_like(phi)
        for perm, wt in zip(PERMS, trust, strict=True):
            if wt < 1e-12:
                continue
            unitary = vx @ I4[perm] @ basis.conj().T
            physical = gates.compile_gate(unitary, cz)
            error = np.linalg.norm(physical - unitary, 2)
            compiler_error = max(compiler_error, error)
            restored += wt * (np.kron(physical, physical.conj()) @ phi)
        flat_channel = np.outer(I4.ravel(), I4.ravel() @ phi) / 4
        packed[k] = TWIRL * flat_channel + (1 - TWIRL) * restored
    return packed, dict(
        permutations=PERMS[trust > 1e-12].tolist(),
        probabilities=trust[trust > 1e-12].tolist(),
        maximum_gate_error=compiler_error,
    )


def fixed(ch, weights):
    # The rank-one term fixes the trace and vanishes on all traceless tangents.
    return np.array(
        [
            (
                np.linalg.solve(
                    np.eye(16) - c + np.outer(I4.ravel() / 4, I4.ravel()),
                    w * I4.ravel() / 4,
                )
            ).reshape(4, 4)
            for c, w in zip(ch, weights, strict=True)
        ]
    )


def compare():
    x, *_ = density_center()
    bank = certificate.channels()
    weights = x.trace(axis1=1, axis2=2).real
    with np.load(certificate.HERE / "program.npz") as program:
        stored_frames = program["frames"]
        stored_weights = program["weights"]
    schedules = []
    report = []
    for b in range(2):
        # Branch states are unnormalized. Dividing by 24 assigns route/charge tickets.
        frames = bank[:, b].reshape(24, 2, 16, 16) / 24
        spl = frames.reshape(-1, 16, 16)
        source = x[b]
        pieces = state_outputs(spl, source).reshape(24, 2, 4, 4)
        fine = slack(pieces.reshape(-1, 4, 4), source)
        blind = slack(pieces.sum(axis=1), source)
        bit_only = slack(pieces.sum(axis=0), source)
        flats = np.arange(1, 4) * weights[b] / 4
        mu = eigen(source)[0].cumsum()[:3]
        threshold = fine / (fine + mu - flats)
        pc, receipt = build(spl, source, stored_frames[b], stored_weights[b])
        schedules.append(pc)
        report.append(
            dict(
                with_bit=fine.tolist(),
                without_bit=blind.tolist(),
                bit_only=bit_only.tolist(),
                maximum_uniform_twirl=float(min(threshold)),
                **receipt,
            )
        )
    return x, np.array(schedules), report


if __name__ == "__main__":
    # Optional commissioning. Replacing eigenvector phases gives the same
    # reference density but can change future correlations and held slopes.
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", help="destination for a newly commissioned gate program"
    )
    args = parser.parse_args()
    x, *_ = density_center()
    bank = certificate.channels()
    panel = [calibrate(bank[:, b].reshape(-1, 16, 16) / 24, x[b]) for b in range(2)]
    np.savez_compressed(
        args.output,
        frames=np.array([t[0] for t in panel]),
        weights=np.array([t[1] for t in panel]),
    )
