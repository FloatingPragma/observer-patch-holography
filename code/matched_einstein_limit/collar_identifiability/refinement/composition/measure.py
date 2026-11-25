"""Local seam arrays, conserved transport controls, and phase-sensitive comparison records."""

import json
from functools import cache
from math import pi, sqrt

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_einstein_station.native_tidal_response.measure import local_calibration
from model.experiments.seam_program import cubic_resonance as cubic
from matched_einstein_limit.coherent_halo import measure as halo
from .. import admission, transport


def edge_probabilities(r1, r2, pairing, fresh, formed):
    """Both-hit and neither-hit inputs get fresh, balanced comparisons."""
    positive = admission.compare_returns(fresh + pairing*(formed-fresh), fresh)
    negative = 1 - positive
    only_left = r1*(1-r2)
    only_right = r2*(1-r1)
    return only_left*positive + only_right*negative + (1-only_left-only_right)/2


def edge_array(level, fresh, formed, u=1/4096, reference=.5):
    board = federation(level)
    owners = np.asarray(board.external)//12
    sites = 4**level
    rates = np.zeros(board.cells)
    rates[:sites] = u/3
    left, right = rates[owners[:, 0]], rates[owners[:, 1]]
    p = edge_probabilities(left, right, reference, fresh, formed)
    y = (formed-fresh)*(left-right)
    np.testing.assert_allclose(p, (1+reference*y)/2, atol=3e-15, rtol=0)
    charge = float(y @ np.arctanh(reference*y))
    boundary = int(np.count_nonzero(left-right))
    return dict(cells=board.cells, enrolled=sites, return_channels=len(owners), boundary=boundary,
                charge_per_u_squared=charge/u**2, leading=reference*(formed-fresh)**2*boundary/9)


def ball_panel():
    # Integrate the three-dimensional ball by transverse disks at fixed x.
    nodes, weights = np.polynomial.legendre.leggauss(6)
    rows = []
    for radius in (1., .5, .25):
        x = radius*nodes
        disk_squared_radius = radius**2 - x**2
        uniform = weights*disk_squared_radius
        modular = weights*disk_squared_radius**2
        amplitude = 1 + x/2
        mean = float(uniform @ amplitude / uniform.sum())
        local_square = float(modular @ amplitude**2 / modular.sum())
        rows.append(dict(radius=radius, pooled_square=mean**2, normalized_kernel_integral=local_square,
                         scaled_difference=(local_square-mean**2)/radius**2))
    return rows


def periodic_step(values):
    """The same eight fair-bit directions as the absorbing transport, on a closed chart."""
    total = values/4
    for axis in range(3):
        total = total + (np.roll(values, 1, axis=axis) + np.roll(values, -1, axis=axis))/8
    return total


def predecessor_variance(values):
    after = periodic_step(values)
    variance = (values-after)**2/4
    for axis in range(3):
        for delta in (-1, 1):
            variance += (np.roll(values, delta, axis=axis)-after)**2/8
    return variance


def transport_panel():
    transport.native_fair_bit()
    mean = np.zeros((9, 9, 9))
    mean[4, 4, 4] = 1
    second = mean**2
    rows = []
    for step in range(1, 5):
        loss = predecessor_variance(mean)
        transported = periodic_step(mean**2)
        after = periodic_step(mean)
        variance = second - mean**2
        next_second = periodic_step(second)
        np.testing.assert_allclose(after**2, transported-loss, atol=1e-15)
        np.testing.assert_allclose(next_second-after**2, periodic_step(variance)+loss, atol=1e-15)
        rows.append(dict(step=step, squared_sum=float(np.sum(after**2)),
                         new_unresolved_variance=float(np.sum(loss)), full_second_moment=float(next_second.sum())))
        mean, second = after, next_second
    return rows


@cache
def phase_meter():
    law = candidate.native_law
    h0, h1 = law.local_term_parts()
    indices = [law.SRC.index(pair) for pair in ((1, 0), (0, 1))]
    h = h0 + candidate.KAPPA*h1
    u = expm(-1j*(pi/2)*h)[np.ix_(indices, indices)]
    np.testing.assert_allclose(u.conj().T @ u, np.eye(2), atol=2e-14)
    _, contacts = local_calibration()
    shift = contacts["coin_even"]["native_potential"][0]
    np.testing.assert_allclose(shift, -32.)
    phases = np.exp(-1j*shift*np.array([pi/64, pi/32]))
    return u, phases


def measured_quadratures(signal, signal_weight, coherence, reference_weight=None):
    """One-unit gate on the parked reference and query address, with an unbiased empty outcome."""
    u, phases = phase_meter()
    if reference_weight is None:
        reference_weight = 1-signal_weight
    if not 0 <= reference_weight <= 1-signal_weight + 1e-15 or not 0 <= coherence <= 1:
        raise ValueError("incompatible signal, reference, and coherence weights")
    records = []
    for phase in phases:
        # The rest of the signal lies outside this pair and produces an unbiased tie.
        a = phase*sqrt(reference_weight)
        b = sqrt(signal_weight)*signal
        local = np.array([[abs(a)**2, coherence*a*b.conjugate()],
                          [coherence*b*a.conjugate(), abs(b)**2]])
        after = u @ local @ u.conj().T
        records.append(float((after[0, 0]-after[1, 1]).real))
    return np.asarray(records)


def interferometry(contrast):
    sample = .3 + .4j
    rows = []
    for weight in (.25, .5, .75):
        full = measured_quadratures(sample, weight, 1.)
        faded = measured_quadratures(sample, weight, .5)
        target = 2*sqrt(weight*(1-weight))*np.array([sample.real, sample.imag])
        np.testing.assert_allclose(full, target, atol=2e-14)
        np.testing.assert_allclose(faded, full/2, atol=2e-14)
        occupancy = weight*abs(sample)**2
        sparse = 1/4096
        local_delta = measured_quadratures(sample, weight+1e-6, 1)-measured_quadratures(sample, weight-1e-6, 1)
        observed = -(contrast*sparse*local_delta/2e-6) @ np.arctanh(contrast*sparse*full)/sparse**2
        rows.append(dict(signal_weight=weight, occupation=occupancy, quadratures=full.tolist(),
                         half_visibility_quadratures=faded.tolist(), response_per_population=float(observed),
                         population_slope_limit=(4*weight-2)*contrast**2*abs(sample)**2))
    flat_gain = []
    for weight in (.25, .5, .75):
        def record(t, visibility):
            return measured_quadratures(sample, t, visibility, reference_weight=(1/16-t/32)/t)

        full = record(weight, 1.)
        delta = (record(weight+1e-6, 1.)-record(weight-1e-6, 1.))/2e-6
        sparse = 1/4096
        slope = -(contrast*sparse*delta) @ np.arctanh(contrast*sparse*full)/sparse**2
        flat_gain.append(dict(signal_weight=weight, reference_weight=(1/16-weight/32)/weight,
                              response_per_population=float(slope),
                              population_slope_limit=contrast**2*abs(sample)**2/16,
                              zero_visibility_max_bias=float(np.max(abs(record(weight, 0))))))
    # A phase-gradient two-address signal supplies a competing density-increasing tangent.
    n = abs(sample)**2
    packet = np.array([sample, -1j*sample*sqrt((1-n)/n)])
    h0, h1 = candidate.native_law.local_term_parts()
    indices = [candidate.native_law.SRC.index(word) for word in ((1, 0), (0, 1))]
    h = (h0 + candidate.KAPPA*h1)[np.ix_(indices, indices)]
    t = .75
    speed = 2*(packet[0].conjugate()*(h @ packet)[0]).imag
    clock_factor = n/(t*speed)
    assert speed > 0
    v = (1/16-t/32)/t
    epsilon = 1e-6
    def shape_read(offset):
        moved = expm(-1j*h*clock_factor*offset) @ packet
        return measured_quadratures(moved[0], t, 1, v), t*abs(moved[0])**2
    plus, occ_plus = shape_read(epsilon)
    minus, occ_minus = shape_read(-epsilon)
    derivative = (plus-minus)/(2*epsilon)
    shape_slope = -(contrast*sparse*derivative) @ np.arctanh(contrast*sparse*shape_read(0.)[0])/sparse**2
    comparison = dict(same_local_occupation_derivative=(occ_plus-occ_minus)/(2*epsilon),
                      population_slope=flat_gain[-1]["response_per_population"], native_shape_slope=float(shape_slope),
                      shape_slope_limit=-2*v*n*contrast**2)
    return dict(two_rail=rows, three_rail=flat_gain, matched_tangents=comparison)


def occupancy_meter(q, p):
    delta = q-p
    def gain(hit):
        return float(-delta @ np.log(p+hit*delta))
    rows = []
    for hit in (0., .25, .5, .75, 1.):
        rows.append(dict(hit=hit, gain=gain(hit), slope=float(-(delta**2/(p+hit*delta)).sum())))
    # Constant density on a unit-volume comparison torus. All radii below ensure a physical lottery.
    gamma = 1/gain(0.)
    radii = []
    for radius in (1/64, 1/128, 1/256):
        kernel_volume = 8*pi**2*radius**4/15
        assert gamma*pi*radius < 1
        background = gamma*kernel_volume
        radii.append(dict(radius=radius, peak_lottery_probability=gamma*pi*radius,
                          normalized_charge=gain(background)*gamma))
    # One common position and two material return registers. The second control overlaps the first ball.
    probability = np.array([.25, .2, .3, .25])
    tangent = np.array([.2, -.4, .1, .1])
    joint = []
    for name, lotteries in (("disjoint", [[.03, .05, 0, 0], [0, 0, .01, .07]]),
                            ("overlapping", [[.03, .05, 0, 0], [0, .04, .01, .07]])):
        weights = np.array(lotteries)
        chances = p+weights[:, :, None]*delta
        conditional = chances[0, :, :, None]*chances[1, :, None, :]
        law = np.einsum("x,xij->ij", probability, conditional)
        direction = np.einsum("x,xij->ij", tangent, conditional)
        hits, dh = weights @ probability, weights @ tangent
        np.testing.assert_allclose(law.sum(axis=1), p + delta * hits[0])
        np.testing.assert_allclose(law.sum(axis=0), p + delta * hits[1])
        score = float(-np.sum(direction*np.log(law)))
        linear = gain(0.)*sum(dh)
        ratio = delta/p
        envelope = np.max(abs(ratio))**2 * hits.sum() * np.sum(abs(dh)) / np.min(np.r_[1., q/p])
        if name == "disjoint":
            np.testing.assert_allclose(law/(p[:, None]*p[None, :]), 1+hits[0]*ratio[:, None]+hits[1]*ratio[None, :])
            assert abs(score-linear) <= envelope
        joint.append(dict(regions=name, raw_joint_slope=score, linear_marginal_sum=linear,
                          disjoint_error_bound=float(envelope) if name == "disjoint" else None))
    return dict(formed=p.tolist(), independent=q.tolist(), gains=rows, inverse_gain=gamma, flat=radii,
                same_state_regions=joint)


def energy_panel(q, p):
    """One randomized matching read followed by a native, positive three-outcome collar read."""
    side = 6
    radius = 1/131072
    shape = (side,)*3
    coord = np.stack(np.unravel_index(np.arange(1, side**3), shape), axis=1)
    critical = 1/(12*(cubic.green(0.)-1/6))
    mass = halo.binding(side, critical)
    stationary, _ = halo.normalized_mode(side, mass)
    kinetic, shell, _ = cubic.punctured_cube(side)
    hopping, offset = .5, 6*candidate.KAPPA
    hamiltonian = hopping*kinetic.toarray()+np.diag(offset-2*hopping*critical*shell)
    matchings = []
    for axis in range(3):
        for parity in (0, 1):
            a = np.flatnonzero(coord[:, axis] % 2 == parity)
            other = coord[a].copy()
            other[:, axis] = (other[:, axis]+1) % side
            b = np.ravel_multi_index(other.T, shape) - 1
            matchings.append((a[b>=0], b[b>=0]))
    displacement = (coord - np.array([1, 0, 0]) + side//2) % side - side//2
    kernel = pi*radius*np.maximum(0, 1-np.sum(displacement**2, axis=1)/4)
    regional = hamiltonian*(kernel[:, None]+kernel[None, :])/2
    u, phases = phase_meter()
    ratio = q-p
    gain = -ratio @ np.log((p+q)/2)
    gamma = 1/gain
    diagonal_scores = 7*kernel*hamiltonian.diagonal()
    rows = []
    for name, phase in (("stationary", np.ones_like(stationary)),
                        ("staggered", (-1.)**coord.sum(axis=1)),
                        ("phase wave", np.exp(2j*pi*coord.sum(axis=1)/side))):
        packet = stationary*phase
        pos = abs(packet)**2
        score = pos @ diagonal_scores / 7
        peak = float(max(abs(diagonal_scores)))
        for a, b in matchings:
            pair = np.stack([packet[a]*phases[0], packet[b]])
            detected = abs(u @ pair)**2
            outcomes = -7*hopping*(kernel[a]+kernel[b])/2
            peak = max(peak, float(max(abs(outcomes))))
            score += (detected[0]-detected[1]) @ outcomes / 7
        expected = (packet.conj() @ regional @ packet).real
        np.testing.assert_allclose(score, expected, atol=2e-16)
        assert gamma*peak < .5
        hit = .5+gamma*score
        rows.append(dict(source=name, weighted_occupation=float(kernel @ pos), weighted_energy=float(expected),
                         observed_signed_estimator=float(score), terminal_hit=float(hit),
                         finite_gain=float(-gamma*ratio @ np.log(p+hit*ratio))))
    return dict(local_radius=radius, side=side, settings=7, gamma=gamma,
                eigenenergy=offset-hopping*mass**2, largest_lottery_offset=gamma*peak, rows=rows)


def run():
    material = admission.material_panel(1, 1, 8)["records"][-1]
    fresh, z = material["q"][-1], material["d"][-1]
    q = np.asarray(material["q"])
    p = q + np.asarray(material["d"])
    return dict(material_contrast=z, local_arrays=[edge_array(i, fresh, fresh+z) for i in range(3)],
                ball_smearing=ball_panel(), transport=transport_panel(), phase_readout=interferometry(z),
                weighted_occupancy_read=occupancy_meter(q, p), weighted_energy_read=energy_panel(q, p))


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
