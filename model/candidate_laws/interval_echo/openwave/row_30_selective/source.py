"""Compile and replay native reference-phase controls on the free arrow component."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30_arrows.experiment import (
    local_check, sector,
)
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
STEP = 0.173
LAYERS = 96


class Source:
    def __init__(self):
        local_check()
        car = federation(0)
        self.car = car
        h, self.sites, inventory = sector(car, list(range(car.n)))
        if inventory["omitted_states"] or h.shape != (32, 32):
            raise ValueError("Source requires the complete ascending-arrow sector")
        self.h = h.toarray()
        self.energies, self.modes = np.linalg.eigh(self.h)
        self.modes *= np.where(self.modes[0] >= 0, 1., -1.)
        self.step = expm(-1j * STEP * self.h)
        self.initial = np.eye(len(self.sites), dtype=complex)[:, 0]
        self.ports = sorted(int(p) for p in set(self.sites) if car.degree[p] == 5)
        self.masks = np.array([self.sites == p for p in self.ports], dtype=float)
        self.indices = np.argsort(-self.modes[0]**2, kind="stable")[:6].tolist()
        self.reference = self.indices[0]
        self.targets = {}
        for j in self.indices:
            self.targets[f"mode_{j}"] = self.modes[:, j]
            if j != self.reference:
                self.targets[f"beat_{j}"] = (self.modes[:, self.reference] + self.modes[:, j]) / np.sqrt(2)
        # Only the contrast service is applied to the blank private degree reference.
        law = candidate.native_law
        local = np.eye(law.DIM)[:, [law.SRC.index((0, 0)), law.SRC.index((1, 0))]]
        _, phase_operator = law.local_term_parts()
        if np.linalg.norm(phase_operator @ local - local @ np.diag([0., 1.])) > 1e-12:
            raise ArithmeticError("Blank reference does not apply the required local phase")

    def trajectory(self, angles):
        phases = np.exp(-1j * (angles @ self.masks))
        path = np.empty((len(angles) + 1, len(self.sites)), dtype=complex)
        path[0] = self.initial
        for i, phase in enumerate(phases):
            path[i + 1] = phase * (self.step @ path[i])
        return path, phases

    def compile(self, target, seed):
        def loss_gradient(flat):
            angles = flat.reshape(LAYERS, len(self.ports))
            path, phases = self.trajectory(angles)
            amplitude = np.vdot(target, path[-1])
            backward = target * amplitude
            gradient = np.empty_like(angles)
            for j in reversed(range(len(phases))):
                gradient[j] = -2 * (self.masks @ (backward.conj() * path[j + 1]).imag)
                backward = self.step.conj().T @ (phases[j].conj() * backward)
            return float(1. - abs(amplitude)**2), gradient.ravel()

        initial = np.random.default_rng(seed).normal(0, .2, (LAYERS, len(self.ports)))
        result = minimize(loss_gradient, initial.ravel(), method="L-BFGS-B", jac=True,
                          options=dict(maxiter=450, maxcor=8, ftol=1e-13, gtol=2e-7))
        return dict(angles=(result.x.reshape(initial.shape) % (2*np.pi)).round(11).tolist(),
                    iterations=int(result.nit), infidelity=float(result.fun))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Choose a fresh path for the new pulse settings")
    source = Source()
    controls = {}
    for index, (name, target) in enumerate(source.targets.items()):
        controls[name] = source.compile(target, 93181 + index)
        print(name, controls[name]["infidelity"], flush=True)
    with args.output.open("x") as output:
        json.dump(dict(ports=source.ports, step=STEP, layers=LAYERS, pulses=controls), output,
                  separators=(",", ":"), allow_nan=False)
        output.write("\n")
