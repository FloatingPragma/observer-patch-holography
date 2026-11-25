"""Conserved density of a resonant source on the addressed cubic Interval Echo guide."""

import json
from math import pi

import numpy as np
from scipy import sparse
from scipy.linalg import eigh, expm
from scipy.optimize import brentq

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.experiments.seam_program import cubic_resonance as cubic
from model.experiments.seam_program.range_contact import native_check


def axial_kernel(side, mass, offsets):
    """Invert the periodic massive difference operator along one axis for each transverse mode."""
    if (
        side < 4
        or side % 2
        or mass <= 0
        or any(r < 0 or r > side // 2 for r in offsets)
    ):
        raise ValueError(
            "need positive screening mass, an even side and offsets in the half box"
        )
    angles = 2 * pi * np.arange(side) / side
    split = 4 * np.sin(angles / 2) ** 2
    transverse = mass**2 + split[:, None] + split[None, :]
    decay = 2 * np.arcsinh(np.sqrt(transverse) / 2)
    weights = 1 / (2 * np.sinh(decay) * (-np.expm1(-side * decay)))
    return np.array(
        [
            np.mean(weights * (np.exp(-r * decay) + np.exp(-(side - r) * decay)))
            for r in offsets
        ]
    )


def binding(side, kappa):
    if kappa <= 0:
        raise ValueError("a strictly attractive contact is required")

    def mismatch(mass):
        center, neighbor = axial_kernel(side, mass, (0, 1))
        return center / neighbor - 1 - 2 * kappa

    return brentq(mismatch, 1e-7 / side, 2 * np.sqrt(3 + 2 * kappa), xtol=5e-15)


def normalized_mode(side, mass):
    frequencies = 4 * np.sin(pi * np.arange(side) / side) ** 2
    symbol = (
        frequencies[:, None, None]
        + frequencies[None, :, None]
        + frequencies[None, None, :]
    )
    green = np.fft.ifftn(1 / (mass**2 + symbol)).real
    exterior = green.ravel()[1:]
    return exterior / np.linalg.norm(exterior), float(np.vdot(exterior, exterior))


def local_energy(kinetic, shell, psi, kappa, hopping, offset):
    """Compare two site energy meters and compute exact number and energy continuity residuals."""
    diagonal = offset - 2 * hopping * kappa * shell
    h = hopping * kinetic + sparse.diags(diagonal)
    phi = h @ psi
    p = abs(psi) ** 2
    sym = (psi.conj() * phi).real
    bonds = sparse.triu(kinetic, k=1, format="coo")
    x, y = bonds.row, bonds.col
    charge = hopping * abs(psi[x] - psi[y]) ** 2 / 2
    bond = (
        diagonal * p
        + np.bincount(x, charge, minlength=len(p))
        + np.bincount(y, charge, minlength=len(p))
    )
    np.testing.assert_allclose(bond - sym, -hopping * (kinetic @ p) / 2, atol=3e-14)
    current = -2 * hopping * bonds.data * (psi[x].conj() * psi[y]).imag
    energy_current = (
        -hopping * bonds.data * (psi[x].conj() * phi[y] + phi[x].conj() * psi[y]).imag
    )

    def divergence(j):
        return np.bincount(x, j, minlength=len(p)) - np.bincount(y, j, minlength=len(p))

    pdot = 2 * (psi.conj() * phi).imag
    edot = (1j * phi.conj() * phi - 1j * psi.conj() * (h @ phi)).real
    return dict(
        sym=sym,
        bond=bond,
        number_error=float(np.max(abs(pdot + divergence(current)))),
        energy_error=float(np.max(abs(edot + divergence(energy_current)))),
        max_bond_current=float(np.max(abs(current))),
        max_energy_current=float(np.max(abs(energy_current))),
    )


def check_finite(side, kappa, check_spectrum=False):
    mass = binding(side, kappa)
    psi, normalization = normalized_mode(side, mass)
    kinetic, shell, _ = cubic.punctured_cube(side)
    h = kinetic - sparse.diags(2 * kappa * shell)
    residual = np.linalg.norm(h @ psi + mass**2 * psi)
    np.testing.assert_allclose(residual, 0, atol=3e-13)
    if check_spectrum:
        eigs, vecs = eigh(h.toarray(), subset_by_index=[0, 1])
        np.testing.assert_allclose(eigs[0], -(mass**2), atol=5e-14)
        np.testing.assert_allclose(abs(vecs[:, 0] @ psi), 1.0, atol=1e-13)
        assert eigs[1] > eigs[0]
    hopping, offset = 0.5, 6 * KAPPA
    energy = local_energy(kinetic, shell, psi, kappa, hopping, offset)
    np.testing.assert_allclose(
        energy["sym"], (offset - hopping * mass**2) * psi**2, atol=2e-13
    )
    r = max(1, int(np.sqrt(side) / 2))
    i = r * side**2 - 1
    return dict(
        side=side,
        contact=kappa,
        screening_mass=mass,
        green_norm_squared=normalization,
        unitary_eigen_residual=float(residual),
        probability_at_axial_probe=float(psi[i] ** 2),
        bond_energy_per_load_at_probe=float(energy["bond"][i] / psi[i] ** 2),
        symmetrized_energy_per_load=offset - hopping * mass**2,
        probe_offset=r,
    )


def profile(side, kappa):
    m = binding(side, kappa)
    r = int(np.sqrt(side) / 2)
    near, far = axial_kernel(side, m, (r, 2 * r))
    return dict(
        side=side,
        separation=r,
        screening_mass=m,
        mass_times_side=side * m,
        density_ratio=float((near / far) ** 2),
    )


def hardware(contact, duration=1 / 64):
    """Blank single-load hopping, degree-phase restoration, and occupied-anchor attraction."""
    h = native_law.seam_term(KAPPA)
    ids = [native_law.SRC.index(w) for w in ((1, 0), (0, 1))]
    np.testing.assert_array_equal(
        h[np.ix_(ids, ids)], (KAPPA + 0.5) * np.eye(2) - 0.5 * np.fliplr(np.eye(2))
    )
    # In the parked reduction, the desired well on six neighbors is -contact times anchor occupancy.
    # Its conditional phase uses the verified native four-column operation; the sham uses the same helpers.
    tau = contact * duration / (2 * KAPPA + 1)
    t = KAPPA * duration / (KAPPA + 0.5)
    z = np.diag([-1, 1])
    half = expm(-0.5j * t * h[np.ix_(ids, ids)])
    restored = (z @ half) @ (z @ half)
    restoration_error = np.linalg.norm(
        restored - np.exp(-1j * KAPPA * duration) * np.eye(2)
    )
    np.testing.assert_allclose(restoration_error, 0, atol=5e-14)
    return dict(
        occupied=native_check(tau),
        sham=native_check(tau, coupled=False),
        degree_restoration_error=float(restoration_error),
    )


def run():
    g0 = cubic.green(0.0)
    critical = 1 / (12 * (g0 - 1 / 6))
    boards = [check_finite(n, critical, n == 8) for n in (8, 16, 32)]
    return dict(
        critical_contact=critical,
        hardware=hardware(critical),
        boards=boards,
        critical_profile=[profile(n, critical) for n in (64, 128, 256, 512)],
        unadjusted=check_finite(8, KAPPA, True),
        unadjusted_profile=profile(64, KAPPA),
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
