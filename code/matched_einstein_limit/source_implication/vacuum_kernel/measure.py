"""Ball-averaged and canonical Gaussian reads for the supplied routed scalar frequencies."""

from math import pi, sqrt

import numpy as np
from scipy.linalg import eigh

from model.experiments.seam_program.clocked_field import eigenphase
from model.experiments.seam_program.scalar_vacuum import SPEED


def frequencies(sites, mesh):
    """The real-transform half-cube, with the routing program's private angle fixed to zero."""
    full = 2 * np.pi * np.fft.fftfreq(sites)
    half = 2 * np.pi * np.fft.rfftfreq(sites)
    angles = np.zeros((sites, len(half), 4))
    angles[:, :, 1] = full[:, None]
    angles[:, :, 2] = half[None, :]
    result = np.empty((sites, sites, len(half)))
    for idx, angle in enumerate(full):
        angles[:, :, 0] = angle
        result[idx] = eigenphase(angles) / mesh
    return result


def ball_sites(radius_mesh):
    limit = int(np.floor(radius_mesh))
    cube = np.array(np.meshgrid(*([np.arange(-limit, limit + 1)] * 3), indexing="ij"))
    points = cube.reshape(3, -1).T
    return points[(points**2).sum(axis=1) <= radius_mesh**2]


def covariances(sites, mesh):
    """Unscaled site canonical variables q=h^(3/2) phi, p=h^(3/2) pi with the constant removed."""
    omega = frequencies(sites, mesh)
    inv = np.divide(0.5, omega, out=np.zeros_like(omega), where=omega > 0)
    dims = (sites,) * 3
    return (np.fft.irfftn(inv, s=dims, axes=(0, 1, 2)),
            np.fft.irfftn(omega / 2, s=dims, axes=(0, 1, 2)))


def averaged_read(sites, mesh, radius=1):
    points = ball_sites(radius / mesh)
    if 2 * radius >= sites * mesh:
        raise ValueError("the ball must fit in the periodic box")
    omega = frequencies(sites, mesh)
    weights = np.zeros((sites,) * 3)
    weights[tuple(points.T)] = 1 / len(points)
    spectrum = np.abs(np.fft.rfftn(weights))**2
    multiplicity = np.full(omega.shape[-1], 2)
    multiplicity[0] = 1
    if sites % 2 == 0:
        multiplicity[-1] = 1
    variance = float(np.sum(multiplicity * np.divide(spectrum, 2 * omega,
                            out=np.zeros_like(omega), where=omega > 0)) / (sites * mesh)**3)
    ideal = 9 / (16 * pi**2 * SPEED * radius**2)
    return variance, 1 / (2 * variance), variance / ideal


def modular_read(sites, mesh, radius=1):
    """Uniform *local* field or momentum displacement. A zero-average continuation exists outside."""
    points = ball_sites(radius / mesh)
    if 2 * radius >= sites * mesh:
        raise ValueError("the ball must fit in the periodic box")
    coord = tuple((points[:, None] - points).transpose(2, 0, 1))
    full_x, full_p = covariances(sites, mesh)
    xmat, pmat = full_x[coord], full_p[coord]
    count = len(points)
    # [q_i,p_j]=i(delta_ij-1/N). Use the canonical regional momentum D_R^{-1}p_R.
    inv_comm = np.eye(count) + np.ones((count, count)) / (sites**3 - count)
    pcan = inv_comm @ pmat @ inv_comm
    xval, xb = eigh(xmat)
    root = (xb * np.sqrt(xval)) @ xb.T
    inv_root = (xb / np.sqrt(xval)) @ xb.T
    val, rot = eigh(root @ pcan @ root)
    nu = np.sqrt(val)
    if np.min(nu) - 0.5 < 1e-12:
        raise ValueError("near-pure mode requires more precision for the logarithm")
    eps = np.log1p(1 / (nu - 0.5))
    delta = np.full(count, mesh**1.5)
    q = rot.T @ inv_root @ delta
    p = rot.T @ root @ inv_comm @ delta
    field = float((nu * eps * q**2).sum() / 2)
    momentum = float((eps / nu * p**2).sum() / 2)
    return dict(points=count, symplectic_gap=nu.min()-0.5, field=field, momentum=momentum,
                field_ratio=field / (4*pi**2*SPEED*radius**2/3),
                momentum_ratio=momentum / (4*pi**2*radius**4/(15*SPEED)))


if __name__ == "__main__":
    for radius_mesh in (1, 1.5, 2, 2.5, 3):
        try:
            print(f"R/h={radius_mesh:g}, modular read {modular_read(48, 1/radius_mesh)}")
        except ValueError as error:
            print(f"R/h={radius_mesh:g}: {error}")
    for sites, mesh in ((32,1/4), (64,1/8), (96,1/8), (128,1/8)):
        print(f"N={sites}, h={mesh}, ball average {averaged_read(sites, mesh)}")
