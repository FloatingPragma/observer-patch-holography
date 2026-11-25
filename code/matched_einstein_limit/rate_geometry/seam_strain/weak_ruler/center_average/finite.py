"""Numerical positive-source replay of the individually volume-matched area and count."""

import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

from .. import aperture, flight, hinge
from ... import live_cadence, live_radar, live_source, live_survey
from ...source import measure as source


HERE = Path(__file__).parent


def range_rebates():
    """Timeout-subtracted ten fair orientation- and charge-pooled range effects."""
    divisor = 16
    t = np.arange(flight.DEADLINE * divisor + 1) / divisor
    effects = np.zeros((10, 8, 8), complex)
    for q in (-1, 1):
        h, rows = flight.hardware(q)
        kernels = live_cadence.kernel_grid(h, rows[:5], t)
        for edge, (a, b) in enumerate(hinge.EDGES):
            returns = live_cadence.effects(kernels, a, b, -np.expm1(-1/divisor))
            effects[edge] -= np.einsum("t,tij->ij", flight.DEADLINE - t, returns) / 4
    return effects


def area_measure(node_count):
    # u=sin(latitude) gives a smooth integrand through both endpoints.
    nodes, weights = np.polynomial.legendre.leggauss(node_count)
    grid = []
    for z, w in zip(nodes, weights, strict=True):
        t = (mp.mpf(float(z)) + 1) * mp.pi / 4
        u, h = mp.sin(t), mp.cos(t)
        grid.append((u, h, mp.mpf(float(w))*mp.pi/4*u*h))

    def individual(alpha):
        eta = mp.pi - alpha/2
        sn, cs = mp.sin(eta), mp.cos(eta)
        total = 0
        for u, h, weight in grid:
            y, v = u*sn, u*cs
            lam, psi = mp.atan2(sn*h, cs), mp.atan2(h, v)
            p = 1 - (lam - y*psi)/mp.pi
            b = 1 - (lam - y*(3-y*y)*psi/2 + y*v*h/2)/mp.pi
            total += weight*p*p*b**(-mp.mpf(7)/3)
        return 2*mp.pi*alpha*(1 - 2*total)

    def measure(base, shift):
        lengths = mp.matrix(5)
        for i, j in hinge.EDGES:
            lengths[i, j] = lengths[j, i] = mp.mpf(float(base[i,j])) + float(shift[i,j])

        def gram(i, j):
            return (lengths[4,i]**2 + lengths[4,j]**2 - lengths[i,j]**2)/2

        angles = [0]*4
        for face in hinge.BASES:
            for i in face:
                j, k = [a for a in face if a != i]
                g, x, y = gram(i,i), gram(i,j), gram(i,k)
                cosine = (gram(j,k) - x*y/g) / mp.sqrt(
                    (gram(j,j)-x*x/g) * (gram(k,k)-y*y/g)
                )
                angles[i] += mp.acos(cosine)
        return sum(lengths[4,i] * individual(angle) for i, angle in enumerate(angles))/3

    return measure


def run(repair=1e-7):
    if not 0 < repair < 1:
        raise ValueError("the re-pairing probability must be interior")
    board = source.formation.federation(1)
    route = source.formation.region(board)
    formed, product = source.form(board, route, -1)
    g = 2 * live_survey.rotations(formed)["right"]
    w0 = np.asarray(source.RAILS["reference"])
    dw = np.asarray(source.RAILS["ordinary"]) - w0
    start = np.kron(np.diag(w0), formed)
    kraus = live_source.head_channel(True)
    wave = kraus.reshape(-1, 2, 4, 8).transpose(1, 0, 2, 3).reshape(-1, 8)
    ports = aperture.enroll(board, 1)[0]
    inside = live_radar.station_rows(board, ports).ravel()
    outside = np.delete(kraus, ports, axis=0)
    number = wave[inside].conj().T @ wave[inside]
    hw = source.system(board, route, -1)[0] @ wave / 40
    ph = wave[inside].conj().T @ hw[inside]
    energy = (ph + ph.conj().T) / 2
    eigen, basis = np.linalg.eigh(number)
    root_number = (basis * np.sqrt(eigen)) @ basis.conj().T

    def read(operator, state):
        return float(np.trace(operator @ state).real)

    def region_entropy(state):
        blocks = (root_number @ state @ root_number,
                  np.einsum("pai,ij,pbj->ab", outside, state, outside.conj(), optimize=True))
        spec = np.concatenate([np.linalg.eigvalsh(a) for a in blocks])
        if min(spec) <= 0 or abs(sum(spec) - 1) > 1e-11:
            raise ArithmeticError("the regional state must be normalized and faithful")
        return float(-spec @ np.log(spec))

    # Pull back the actual radar to incoming rail/heads through the four-tick
    # pump. Exact state contrasts are traceless, including the complete timeout.
    radar = np.einsum("pai,eab,pbj->eij", kraus.conj(), range_rebates(), kraus, optimize=True)
    base = live_survey.square(3 + live_survey.read(radar, start), 5)
    zero = np.zeros((5, 5))
    paired = formed + repair * (product - formed)

    def preparation(theta):
        # Unit rotation has generator g/2, g²=I. Compute tiny changes directly.
        dh = (repair*(product-formed) - np.sin(theta/2)**2*(paired-g@paired@g)
              - 0.5j*np.sin(theta)*(g@paired-paired@g))
        change_head = np.kron(np.diag(w0), dh)
        change_rail = np.kron(np.diag(dw), formed + dh)
        rail = -read(number, change_head) / read(number, change_rail)
        return change_head + rail*change_rail, rail

    def delta_lengths(change):
        return live_survey.square(live_survey.read(radar, change), 5)

    with mp.workdps(40):
        area = area_measure(32)
        initial_area = area(base, zero)

        def residual(rate):
            change, _ = preparation(rate*repair)
            return float(area(base, delta_lengths(change)) - initial_area) / repair

        # Solve only from the measured area and count. The dyadic interval
        # transfer and the independent energy readout play no role in this root.
        at_zero = residual(0)
        upper = 1
        while upper < 8192 and residual(upper) * at_zero > 0:
            upper *= 2
        angle_rate = brentq(residual, 0, upper, xtol=1e-10)
        theta = angle_rate * repair
        change, rail = preparation(theta)
        shifted = delta_lengths(change)
        check = area_measure(64)
        difference = float(check(base, shifted) - check(base, zero))
        # Assemble the two positive lotteries and the unitary directly to
        # check the cancellation-stable variation used by the constraint solve.
        unitary = np.cos(theta/2)*np.eye(4) - 1j*np.sin(theta/2)*g
        weights = w0 + rail*dw
        assembled = np.kron(np.diag(weights), unitary @ paired @ unitary.conj().T)
        s = region_entropy(assembled) - region_entropy(start)
        e = read(energy, change)
        n = read(number, change)
        unchanged_controls, _ = preparation(0)
        return dict(
            repair_probability=repair, right_angle=theta, rail_parameter=rail,
            rail_weights=weights.tolist(), enrolled_ports=ports,
            minimum_input_eigenvalue=float(np.linalg.eigvalsh(assembled)[0]),
            input_trace=float(np.trace(assembled).real),
            independently_assembled_input_discrepancy=float(
                np.max(abs(assembled - (start + change)))
            ),
            reference_individual_area_coefficient=float(initial_area),
            reference_number=read(number, start), reference_energy=read(energy, start),
            constrained_area_change=difference, constrained_number_change=n,
            area_change_at_count_match_without_rotation=float(
                check(base, delta_lengths(unchanged_controls)) - check(base, zero)
            ),
            regional_energy_change=e, regional_entropy_change=s,
            directly_subtracted_count=read(number, assembled) - read(number, start),
            directly_subtracted_energy=read(energy, assembled) - read(energy, start),
            energy_secant=e/repair, entropy_secant=s/repair,
            rail_secant=rail/repair, angle_secant=angle_rate,
            quadrature_nodes_for_solve=32, quadrature_nodes_for_check=64,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "finite.json").write_text(data)
    else:
        print(data, end="")
