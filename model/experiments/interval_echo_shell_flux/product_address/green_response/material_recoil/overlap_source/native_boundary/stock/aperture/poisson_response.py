"""Test the Poisson response and screening of the native ordinary medium."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh, expm_multiply, spsolve

from . import isolated_medium as lab

HERE = lab.HERE
K = lab.lab.K
T = 65536


def weak_density(lap, m):
    """Derivative of N times site probability at the blank rest eigenmode."""
    rhs = 4 * K * (m - m.mean())
    v = np.r_[0.0, spsolve(lap[1:, 1:], rhs[1:])]
    return v - v.mean()


def current_sums(board, shells, field):
    a, b = np.array(board.seams).T
    du = (field[a] - field[b]) * np.sign(shells[b] - shells[a])
    cuts = np.minimum(shells[a], shells[b])
    return np.bincount(cuts, weights=du, minlength=shells.max() + 1)[:-1]


def lowest(h, count=1):
    # A fixed nonsymmetric start also sees odd modes on reflection-symmetric pairs.
    e, v = eigsh(
        h,
        k=count,
        which="SA",
        tol=2e-13,
        v0=np.random.default_rng(7307).normal(size=h.shape[0]),
    )
    ix = np.argsort(e)
    return e[ix], v[:, ix]


def dephased(rho, e, v):
    # Distinct clusters in this finite numerical comparison are well separated.
    gap = np.diff(e)
    assert np.min(gap[gap > 1e-9]) > 1e-6
    transformed = v.T @ rho @ v
    return v @ (transformed * (abs(e[:, None] - e) < 1e-9)) @ v.T


def pulse(h, hw, start, duty, steps, time):
    result = start.astype(complex)
    for _ in range(steps):
        result = expm_multiply((-1j * time * (1 - duty) / steps) * h, result)
        result = expm_multiply((-1j * time * duty / steps) * hw, result)
    return result


def panel(level, tape):
    board = lab.lab.cold.federation(level)
    if level == 0:
        marks, _ = lab.lab.cold.read_map(board, tape["spectral_tape"])
    else:
        edge = tuple(sorted(min(board.external)))
        reading = lab.lab.cold.source_tape(
            board, list(edge), 2, K, np.random.default_rng(7410 + level)
        )
        marks, _ = lab.lab.cold.read_map(board, reading)
    shells = lab.lab.cold.distances(board, marks[0])
    h, lap, _ = lab.lab.cold.motion(board, [], 1, K)
    hw, live, m = lab.lab.cold.motion(board, marks, 1, K)
    n = board.n
    chi = weak_density(lap, m)
    flux = current_sums(board, shells, chi)
    volumes = np.bincount(shells).cumsum()[:-1]
    np.testing.assert_allclose(flux, 4 * K * (m.sum() - m.mean() * volumes), atol=5e-11)
    cold_zero, _, _ = lab.lab.cold.motion(board, marks, 1, 0)
    np.testing.assert_array_equal(cold_zero @ np.ones(n), np.zeros(n))
    e, v = lowest(hw, 2)
    degree = np.array(live.diagonal())
    assert e[1] <= -5.5 and e[0] >= -8
    adj = sparse.diags(degree) - live
    # For each exact eigenmode the positive tail bound uses |E|>=5.5 and
    # the row maximum principle. Also evaluate the two bound modes as a group.
    p = (v * v).sum(axis=1) / 2
    f = 4 * K * (p + adj @ p)
    pflux = current_sums(board, shells, p)
    fflux = current_sums(board, shells, f) / (n - 1)
    blankforce = current_sums(board, shells, (1 + board.degree) * 4 * K / n) / (n - 1)
    viable = (0.01 * (6 / 17) ** (2 * np.arange(len(flux)))) > 2e-15
    assert np.all(pflux[viable] > 0) and np.all(
        fflux[viable & (np.arange(len(flux)) >= 1)] > 0
    )
    seed_fraction = float((v[list(marks[0]), :] ** 2).sum() / 2)
    regular = current_sums(board, shells, 4 * K * (chi + board.adjacency() @ chi))
    spectra = {}
    for duty in (1e-5, 1e-3, 1 / 32, 1):
        en, mode = lowest(h + duty * (hw - h))
        density = n * mode[:, 0] ** 2 - 1
        spectra[str(duty)] = dict(
            energy=float(en[0]),
            relative_linear_error=float(
                np.linalg.norm(density - duty * chi) / np.linalg.norm(duty * chi)
            ),
            endpoint_probability=float((mode[m > 0, 0] ** 2).sum()),
        )
    result = dict(
        ports=n,
        degree_range=[int(degree.min()), int(degree.max())],
        written=marks[0],
        energy=e.tolist(),
        endpoint_seed_bound_fraction=seed_fraction,
        sources_per_port=float(m.mean()),
        sites=np.bincount(shells).tolist(),
        cut_edges=np.bincount(
            np.minimum(
                shells[np.array(board.seams)[:, 0]], shells[np.array(board.seams)[:, 1]]
            )[
                shells[np.array(board.seams)[:, 0]]
                != shells[np.array(board.seams)[:, 1]]
            ]
        ).tolist(),
        relative_density_flux=flux.tolist(),
        exterior_release_radii=list(range(1, len(flux))),
        zero_stiffness_rest_response=0,
        weak_probe_flux_times_N_times_Nminus1=regular[1:].tolist(),
        bound_mass_per_shell=np.bincount(shells, weights=p).tolist(),
        bound_probability_flux=pflux.tolist(),
        bound_probe_acceleration=fflux[1:].tolist(),
        bound_minus_blank_acceleration=(fflux - blankforce)[1:].tolist(),
        finite_exposure=spectra,
    )
    if level > 0:
        result["source_tape"] = reading
    if level == 0:
        psi = np.ones(n) / np.sqrt(n)
        duty = 1 / 32
        reference = expm_multiply(-1j * 4 * (h + duty * (hw - h)), psi)
        result["native_switch_error_at_time_4"] = {
            str(steps): float(
                np.linalg.norm(pulse(h, hw, psi, duty, steps, 4) - reference)
            )
            for steps in (8, 32, 128)
        }
    return result


def endpoint_experiment(tape):
    board = lab.lab.cold.federation(0)
    marks, _ = lab.lab.cold.read_map(board, tape["spectral_tape"])
    shells = lab.lab.cold.distances(board, marks[0])
    rng = np.random.default_rng(549128)
    start = np.diag(np.isin(np.arange(board.n), marks[0]).astype(float) / 2)
    raw, stationary, fields = {}, [], []
    for tag, mask in (("blank", []), ("written", marks)):
        h, meters = lab.instruments(board, shells, mask)
        e, v = eigh(h.toarray())
        rho = lab.averaged(start, (e, v), T)
        stationary.append(
            np.array([np.sum(dephased(start, e, v) * o) for _, _, o in meters])
            / (board.n - 1)
        )
        p = rho.diagonal().real
        raw[tag] = dict(
            sites=rng.multinomial(lab.READS, p / p.sum()).tolist(),
            fringes=[lab.fringe_counts(c, rho, rng) for _, c, _ in meters],
        )
        fields.append(
            np.array([np.sum(rho.real * o) for _, _, o in meters]) / (board.n - 1)
        )
    record = dict(observations={"endpoint": {str(T): raw}})
    record["decoded"] = lab.decode(record, tape)
    lo, hi = np.array(record["decoded"]["endpoint"][str(T)]).T
    assert np.all((lo <= fields[1] - fields[0]) & (fields[1] - fields[0] <= hi))
    record["stationary_change"] = (stationary[1] - stationary[0]).tolist()
    return record


def ordinary_geometry(tape):
    """Shared calibrated geometry with a clamped ordinary side port."""
    b = lab.lab.cold.federation(0)
    marks, _ = lab.lab.cold.read_map(b, tape["spectral_tape"])
    shells = lab.lab.cold.distances(b, marks[0])
    y = int(np.argmax(shells))
    keep = np.arange(b.n) != y
    h, lap, _ = lab.lab.cold.motion(b, [], 1, K)
    h0 = h[keep][:, keep].tocsr()
    hw = lab.lab.cold.motion(b, marks, 1, K)[0][keep][:, keep]
    near = np.asarray(b.adjacency()[y, :].toarray()).ravel()[keep]
    ordinary = lab.lab.ordinary
    joint, i, j = ordinary.two_loads(b, K)
    pinned = np.where(i == y)[0]
    np.testing.assert_array_equal(j[pinned], np.where(keep)[0])
    # Joint native energies are measured relative to both isolated costs.
    ho = joint[pinned][:, pinned] - sparse.eye(keep.sum()) * b.degree[y] / 2
    np.testing.assert_allclose(
        (ho - h0).toarray(), np.diag(-(4 * K + 0.5) * near), atol=0
    )
    e, v = lowest(h0)
    phi = np.abs(v[:, 0])
    conductance = sparse.diags(phi) @ (-h0) @ sparse.diags(phi)
    weighted = sparse.diags(np.asarray(conductance.sum(axis=1)).ravel()) - conductance
    gain = lab.ordinary_gain(tape).mean()
    m = np.isin(np.arange(b.n), marks[0]).astype(float)[keep]
    diff = lap - lab.lab.cold.motion(b, marks, 1, K)[1]
    cut = diff[keep][:, keep]

    def source(name, psi, g):
        return (
            (g + 0.5) * near * phi * psi
            if name == "ordinary"
            else phi * (g * m * psi / 4 + (cut @ psi) / 2)
        )

    rows = {}
    for name, changed in (("ordinary", ho), ("written", hw)):
        src = source(name, phi, gain)
        rhs = src - src.sum() * phi**2
        u = np.r_[0.0, spsolve(weighted[1:, 1:], rhs[1:])]
        u -= phi**2 @ u
        derivative = 2 * phi**2 * u
        per = {}
        for duty in (1e-5, 1e-3, 1):
            en, mode = lowest(h0 + duty * (changed - h0))
            psi = abs(mode[:, 0])
            exact_source = duty * source(name, psi, 4 * K)
            residual = weighted @ (psi / phi) - (
                exact_source + (en[0] - e[0]) * phi * psi
            )
            np.testing.assert_allclose(residual, 0, atol=3e-13)
            per[str(duty)] = dict(
                energy=float(en[0]),
                relative_linear_error=float(
                    np.linalg.norm(psi**2 - phi**2 - duty * derivative)
                    / np.linalg.norm(duty * derivative)
                ),
                exact_equation_max_residual=float(max(abs(residual))),
            )
        rows[name] = dict(
            linear_source=src.tolist(),
            occupation_derivative=derivative.tolist(),
            exposure=per,
        )
    return dict(
        clamp=y,
        ports=int(keep.sum()),
        blank_gap_reference=float(e[0]),
        ordinary_gain=gain,
        rows=rows,
    )


def run():
    origin = HERE / "mediated_release.json"
    tape = json.loads(origin.read_bytes())
    lab.lab.ordinary.check_joint_columns(K)
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (origin, HERE / "POISSON_RESPONSE.md")
        },
        ordinary=ordinary_geometry(tape),
        boards=[panel(level, tape) for level in range(3)],
        endpoint=endpoint_experiment(tape),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    modes = cli.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "poisson_response.json"
    if args.counts_only:
        out = json.loads(path.read_bytes())
        src = json.loads((HERE / "mediated_release.json").read_bytes())
        lab.lab.cold.apparatus.check(
            out["endpoint"]["decoded"], lab.decode(out["endpoint"], src)
        )
    else:
        out = run()
        if args.write:
            path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
        else:
            # The finite equation residuals and deep tails approach machine precision.
            saved = json.loads(path.read_bytes())
            lab.lab.cold.apparatus.check(saved, json.loads(json.dumps(out)))
    print(
        [
            (b["ports"], b["energy"], b["endpoint_seed_bound_fraction"])
            for b in out["boards"]
        ]
    )
    print(out["endpoint"]["decoded"])
