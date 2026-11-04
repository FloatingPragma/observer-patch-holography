"""Finite memory response evaluated on the unchanged source law.

Floating trajectories are diagnostic evaluations. The dimension, signs,
spectral window and convergence statements are rational certificates.
"""
import numpy as np

from ..native_equilibration import reciprocal_vacuum as source


def projected_trajectory(populations, steps, *, model=None):
    if type(steps) is not int or not 0 <= steps <= 10000:
        raise ValueError('integer 0<=steps<=10000 required')
    model = source.finite_model() if model is None else model
    q = np.asarray(populations)
    if (q.shape != (16,) or q.dtype.kind not in 'fi' or not np.isfinite(q).all()
            or np.min(q) < 0 or abs(q.sum()-1) > 1e-12
            or model.field_dimension != 8 or len(model.active) != 16):
        raise ValueError('all sixteen normalized source populations are required')
    pi, P = model.probabilities, model.transition.T
    p = float(pi[8:].sum())
    f = np.repeat([0., 1.], 8)-p
    variance = p*(1-p)
    Q = np.eye(16)-np.outer(np.ones(16), pi)-np.outer(f, pi*f)/variance
    C = Q@P@Q
    b = Q@P@f
    a = float((pi*f)@(P@f)/variance)
    u = q/pi-1
    z0 = Q@u
    kernels, forces, powers_b, powers_z = [], [], [], []
    v, z = b.copy(), z0.copy()
    for _ in range(steps):
        kernels.append(float((pi*b)@v/variance))
        forces.append(float((pi*b)@z))
        powers_b.append(v.copy())
        powers_z.append(z.copy())
        v, z = C@v, C@z
    powers_z.append(z.copy())
    x = [float((pi*f)@u)]
    for n in range(steps):
        x.append(a*x[n]+sum(kernels[j]*x[n-1-j] for j in range(n))+forces[n])
    centered_energy = model.energies-pi@model.energies
    energy_slope = float((pi*centered_energy)@f/variance)
    energy = []
    for n in range(steps+1):
        eliminated = powers_z[n].copy()
        for j in range(n):
            eliminated += powers_b[j]*x[n-1-j]/variance
        energy.append(float(pi@model.energies+energy_slope*x[n]
                            +(pi*centered_energy)@eliminated))
    geometry = p+np.asarray(x)
    area = 20/9+(20/11-20/9)*geometry
    energy = np.asarray(energy)
    return {'geometry1_probability': geometry, 'area': area, 'energy': energy,
            'system_heat_per_update': np.diff(energy),
            'memory_kernel': np.asarray(kernels), 'matter_forcing': np.asarray(forces),
            'instantaneous_coefficient': a}
