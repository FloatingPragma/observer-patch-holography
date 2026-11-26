"""Intrinsic areas on the tilted exponential rest slice of a reconstructed station."""
import argparse
import json
from itertools import combinations
import math
from pathlib import Path

import numpy as np
from mpmath import iv
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from ..geometry import ALL, observed_squared, raw_bounds

HERE = Path(__file__).resolve().parent
ROUND_RADIUS = 0.1 / (math.pi / 2 - math.pi / 20)
CUTOFF_PLATEAU = 0.02
TARGETS = tuple(1 / n for n in (128, 256, 512, 1024))
B = 2
S = 1


def angular_rule(order):
    z, weights = leggauss(order)
    phi = 2 * np.pi * np.arange(2*order)/(2*order)
    z, phi = np.broadcast_arrays(z[:, None], phi[None, :])
    z, phi = z.ravel(), phi.ravel()
    s = np.sqrt(1-z*z)
    n = np.column_stack([z, s*np.cos(phi), s*np.sin(phi)])
    first = np.column_stack([-s, z*np.cos(phi), z*np.sin(phi)])
    second = np.column_stack([np.zeros_like(z), -np.sin(phi), np.cos(phi)])
    return n, np.stack([first, second], axis=1), np.repeat(weights, 2*order)*np.pi/order


class Slice:
    def __init__(self, gamma, anisotropy, order):
        self.n, self.tangents, self.angular = angular_rule(order)
        self.lift = np.diag([gamma, 1, 1])
        self.dt_squared = np.diag([gamma*gamma-1, 0, 0])
        self.B = np.diag([B+anisotropy, B-anisotropy, B])
        self.radial, self.weights = leggauss(order)

    def metric(self, r):
        """Induced three-metric in the orthonormal normal coordinates y=r*n."""
        x = r*self.n @ self.lift
        length = np.linalg.norm(x, axis=1)
        if np.max(length) >= CUTOFF_PLATEAU:
            raise ValueError("quadrature leaves the round, unit-cutoff neighborhood")
        direction = x/length[:, None]
        radial = direction[:, :, None]*direction[:, None, :]
        sinc2 = np.sinc(length/(np.pi*ROUND_RADIUS))**2
        round_metric = radial + sinc2[:, None, None]*(np.eye(3)-radial)
        cross = np.stack([np.cross(x, e) for e in np.eye(3)], axis=2)
        h = round_metric + cross.transpose(0, 2, 1) @ self.B @ cross
        return self.lift @ h @ self.lift - self.dt_squared

    def volume(self, r):
        total = 0
        for u, weight in zip(self.radial, self.weights, strict=True):
            radius = r*(u+1)/2
            density = np.sqrt(np.linalg.det(self.metric(radius)))
            total += weight*radius**2 * (self.angular @ density)
        return float(r*total/2)

    def area(self, r):
        surface = self.tangents @ self.metric(r) @ self.tangents.transpose(0, 2, 1)
        return float(r*r*(self.angular @ np.sqrt(np.linalg.det(surface))))

    def ball(self, ell):
        target = 4*math.pi*ell**3/3
        r = brentq(lambda x: self.volume(x)-target, 0.98*ell, 1.04*ell,
                   xtol=2e-18, rtol=1e-15)
        return dict(proper_radius=r, volume=self.volume(r), boundary_area=self.area(r))


def outward(interval):
    return [math.nextafter(float(interval.a), -math.inf), math.nextafter(float(interval.b), math.inf)]


def measure(order=10):
    square = observed_squared(raw_bounds())
    pairs = [square[i, j] for i, j in combinations(ALL, 2)]
    assert len(pairs) == 15 and all(v.a > (13/5)**2 and v.b < (31/10)**2 for v in pairs)
    station_ports = (0, 1, 38, 39, 5, 6)
    cylinders = []
    for i, j in combinations(ALL, 2):
        length = iv.sqrt(square[i, j])
        tube = length-iv.mpf("0.2")
        cylinders.append(dict(ports=[station_ports[i], station_ports[j]],
                              station_distance=outward(length), cylinder_length=outward(tube)))
    tables = {}
    for name, gamma in (("rest", 1), ("axial", 5/3)):
        plus, minus = (Slice(gamma, s*S, order) for s in (1, -1))
        rows = []
        for ell in TARGETS:
            a, b = plus.ball(ell), minus.ball(ell)
            rows.append(dict(ell=ell, target_volume=4*math.pi*ell**3/3, plus=a, minus=b,
                             signed_half_difference_per_s_ell4=(a["boundary_area"]-b["boundary_area"])/(2*S*ell**4)))
        tables[name] = rows
    return dict(chamber_radius=ROUND_RADIUS, cutoff_plateau=CUTOFF_PLATEAU, b=B, s=S,
                signed_axial_limit=-64*math.pi/45, native_station_pairs=len(pairs), cylinders=cylinders, tables=tables)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = measure()
    if parser.parse_args().write:
        (HERE / "results.json").write_text(json.dumps(result, indent=2)+"\n")
    else:
        prior = json.loads((HERE / "results.json").read_text())
        assert result["cylinders"] == prior["cylinders"]
        for name, rows in result["tables"].items():
            for a, b in zip(rows, prior["tables"][name], strict=True):
                for sign in ("plus", "minus"):
                    np.testing.assert_allclose(list(a[sign].values()), list(b[sign].values()), rtol=1e-12, atol=0)
                np.testing.assert_allclose(a["signed_half_difference_per_s_ell4"], b["signed_half_difference_per_s_ell4"], atol=1e-7, rtol=0)
        print("Verified proper areas of sixteen volume-matched balls.")
