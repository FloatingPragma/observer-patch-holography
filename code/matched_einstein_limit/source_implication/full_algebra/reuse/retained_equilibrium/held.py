"""Introduce signed perturbations to the persistent bath and match its measured volume."""

import numpy as np
from scipy.optimize import brentq

from .control import TWIRL, I4, fixed
from .probe import AXES
from ..stationary_baths.forecast import apparatus as app
from ... import station


ORDINARY_GAIN = 23.104236196474325


def panel(channels, x):
    """Return a held stationary-state preparation, directed POVMs and linear responses."""
    _, clocks = app.instruments(8)
    eff = np.einsum("brtij,t->rbij", clocks, app.SCORES)
    wt = x.trace(axis1=1, axis2=2).real
    cap = channels.sum(axis=1)
    xstar = fixed(cap, wt)

    def survey(z):
        return np.einsum("rbij,bji->r", eff, z).real

    means = survey(xstar)

    def geo_delta(z):
        return np.imag(station.geometry(means + 1e-17j * survey(z))) / 1e-17

    resolvent = [np.eye(16) - v + np.outer(I4.ravel() / 4, I4.ravel()) for v in cap]

    def tangent(dz):
        return np.array(
            [
                np.linalg.solve(k, d.ravel()).reshape(4, 4)
                for k, d in zip(resolvent, dz, strict=True)
            ]
        )

    heat = tangent((wt[:, None, None] * I4 / 4 - xstar) / (1 - TWIRL))
    candidates = []
    for i, j in ((0, 1), (0, 2), (0, 3), (1, 0), (2, 0), (3, 0)):
        h = np.kron(AXES[i], AXES[j])
        rotation = tangent(-0.5j * (h @ xstar - xstar @ h))
        rot_v = geo_delta(rotation)[1]
        candidates.append((abs(rot_v), i, j, rotation, rot_v))
    _, i, j, rot, rv = max(candidates, key=lambda e: e[0])
    h = np.kron(AXES[i], AXES[j])
    dtheta = -geo_delta(heat)[1] / rv
    drift = heat + dtheta * rot
    area, volume = geo_delta(drift)
    entropy = -np.einsum("bij,bji", drift, station.log_blocks(xstar)).real
    flat = np.outer(I4.ravel() / 4, I4.ravel())
    restore = (cap - TWIRL * flat) / (1 - TWIRL)

    def endpoint(t):
        # The fit uses only the reference volume and clock readings.
        def trial(z):
            u = np.cos(z / 2) * I4 - 1j * np.sin(z / 2) * h
            c = np.kron(u, u.conj()) @ ((TWIRL + t) * flat + (1 - TWIRL - t) * restore)
            return fixed(c, wt)

        vol_ref = station.geometry(means)[1]
        angle = brentq(
            lambda a: station.geometry(survey(trial(a)))[1] - vol_ref,
            -4 * abs(t * dtheta),
            4 * abs(t * dtheta),
            xtol=2e-13,
        )
        return trial(angle)

    return (
        endpoint,
        clocks,
        dict(
            fixed_state_trace_error=float(np.abs(np.linalg.eigvalsh(xstar - x)).sum()),
            phase_axes=[i, j],
            phase_per_twirl=dtheta,
            entropy=entropy,
            area=area,
            volume=volume,
            linear_residual=entropy + ORDINARY_GAIN * area,
        ),
    )


def read(channels, x):
    state_at, clocks, report = panel(channels, x)

    def records(t):
        y = state_at(t)
        means = app.census(y, clocks) @ app.SCORES
        return np.array([station.entropy(y), *station.geometry(means)])

    finite = []
    for step in (1 / 100000, 1 / 200000):
        ds, da, dv = (records(step) - records(-step)) / (2 * step)
        finite.append(
            dict(
                step=step,
                entropy=ds,
                area=da,
                volume=dv,
                residual=ds + ORDINARY_GAIN * da,
            )
        )
    return dict(**report, finite=finite)
