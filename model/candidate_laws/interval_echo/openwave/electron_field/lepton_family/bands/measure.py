"""Bloch spectrum and independent signed pulse releases on a verified Gauss ring."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
from model.executable_model import federation

HOP = 4
DEGREE_FIVE_COST = 8
WIDTH = 0.23


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_route():
    route = json.loads((HERE / 'route.json').read_text())
    board = federation(route['federation_level'])
    cells, ports = np.array(route['cells']), np.array(route['ports'])
    n = route['cells_on_ring']
    assert n == len(cells) == len(set(cells)) and len(set(ports)) == 3*n
    np.testing.assert_array_equal(ports // 12, np.repeat(cells, 3))
    np.testing.assert_array_equal(board.degree[ports], np.tile([6, 5, 6], n))
    active = {tuple(sorted((int(a), int(b)))) for a, b in zip(ports, np.roll(ports, -1))}
    assert len(active) == len(ports) and active.issubset(board.seams)
    # The only cycle in a simple ring uses every seam. Neither cell triangles
    # nor elementary intercell faces (5/6 cells) can fit this 48-cell cycle.
    assert n > 6
    return n, {'sphere_cells': board.cells, 'ring_cells': n,
               'active_seams': len(active), 'active_elementary_field_loops': 0,
               'motif_degrees': [6, 5, 6], 'ordered_pair_dimension': 3*n*(3*n-1)}


def bloch(k, cost=DEGREE_FIVE_COST):
    h = np.array([[2*HOP, -HOP, -HOP*np.exp(-1j*k)],
                  [-HOP, 2*HOP+cost, -HOP],
                  [-HOP*np.exp(1j*k), -HOP, 2*HOP]], dtype=complex)
    d = np.zeros_like(h)
    d[0, 2] = -1j*h[0, 2]
    d[2, 0] = d[0, 2].conjugate()
    d2 = np.zeros_like(h)
    d2[0, 2], d2[2, 0] = -h[0, 2], -h[2, 0]
    energies, u = np.linalg.eigh(h)
    u *= np.exp(-1j*np.angle(u[0, :]))
    slope = u.conj().T @ d @ u
    curvature = np.diag(u.conj().T @ d2 @ u).real.copy()
    for b in range(3):
        for j in range(3):
            if j != b and abs(energies[b] - energies[j]) > 1e-12:
                curvature[b] += 2*abs(slope[j, b])**2/(energies[b]-energies[j])
    return h, energies, u, np.diag(slope).real, curvature


def report_spectrum(n, k, vectors):
    zero, edge = bloch(0), bloch(np.pi)
    nz, ne = zero[1], edge[1]
    # det(H-EI) has a single momentum term, -2*HOP**3*cos(k). The open
    # gaps and lack of degeneracy place every extremum at 0 or pi.
    for x in np.linspace(0, np.pi, 65):
        h = bloch(x)[0]
        for e in (1., 8.5, 17.):
            np.testing.assert_allclose(np.linalg.det(h - e*np.eye(3)),
                                       np.linalg.det(zero[0]-e*np.eye(3)) + 2*HOP**3*(1-np.cos(x)),
                                       atol=2e-12)
    bands = []
    floors = np.minimum(nz, ne)
    for b in range(3):
        rest = zero if nz[b] < ne[b] else edge
        bands.append({'band': b+1, 'minimum': float(floors[b]), 'minimum_momentum_over_pi': int(rest is edge),
                      'zero_cell_momentum_energy': float(nz[b]), 'maximum': float(max(nz[b], ne[b])),
                      'minimum_ratio_to_lowest': float(floors[b]/floors[0]),
                      'zero_momentum_ratio_to_lowest': float(nz[b]/nz[0]),
                      'curvature_at_minimum': float(rest[-1][b])})
    comm = 0.
    for idx in range(n):
        h = bloch(k[idx])[0]
        for b in range(3):
            p = np.outer(vectors[idx, :, b], vectors[idx, :, b].conj())
            comm = max(comm, float(np.linalg.norm(h@p-p@h, ord=2)))
    homogeneous_error = max(float(np.max(abs(bloch(q, 0)[1]-np.sort(
        [2*HOP*(1-np.cos((q+2*np.pi*j)/3)) for j in range(3)])))) for q in k)
    assert comm < 1e-12 and homogeneous_error < 1e-12
    return {'complete_sorted_bands': bands, 'neutral_vacuum_energy': 0,
            'koide_Q_of_band_minima': float(floors.sum()/np.sqrt(floors).sum()**2),
            'koide_Q_at_zero_cell_momentum': float(nz.sum()/np.sqrt(nz).sum()**2),
            'max_norm_isolated_band_commutator': comm,
            'homogeneous_folded_chain_max_energy_error': homogeneous_error,
            'homogeneous_folded_energies_at_0_and_pi': [bloch(0, 0)[1].tolist(), bloch(np.pi, 0)[1].tolist()]}


def pair_space(n, cut=False):
    length = 3*n
    a, b = np.where(~np.eye(length, dtype=bool))
    ids = np.full((length, length), -1, int)
    ids[a, b] = np.arange(len(a))
    rows, cols = [], []
    for moving_a in (True, False):
        p, other = (a, b) if moving_a else (b, a)
        for direction in (-1, 1):
            target = (p+direction) % length
            allowed = target != other
            if cut:
                allowed &= target//3 == p//3
            dst = ids[target, other] if moving_a else ids[other, target]
            rows.append(dst[allowed])
            cols.append(np.flatnonzero(allowed))
    i, j = np.concatenate(rows), np.concatenate(cols)
    off = sparse.csr_matrix((np.full(len(i), -HOP, float), (i, j)), shape=(len(a), len(a)))
    kinetic = off + sparse.diags(-np.asarray(off.sum(axis=1)).ravel(), dtype=float, format='csr')
    np.testing.assert_array_equal(kinetic @ np.ones(len(a)), 0)
    cost = DEGREE_FIVE_COST*((a%3 == 1).astype(int)+(b%3 == 1))
    h = kinetic + sparse.diags(cost, dtype=float, format='csr')
    exchange = ids[b, a]
    assert (h-h.T).nnz == 0 and (h[exchange][:, exchange]-h).nnz == 0
    if not cut:
        # A consecutive arc of reversed heads has charge +1 and -1 only at
        # its ends. Vectorize the check for each plus address and every minus.
        link = np.arange(length)
        for start in range(length):
            ends = np.delete(link, start)
            span = (ends-start) % length
            reversed_heads = ((link-start) % length < span[:, None]).astype(int)
            net = reversed_heads - np.roll(reversed_heads, 1, axis=1)
            loads = np.zeros_like(net)
            loads[:, start] = 1
            loads[np.arange(len(ends)), ends] = -1
            np.testing.assert_array_equal(net, loads)
    return a, b, h, exchange


def packet(k, vectors, velocities, band, momentum, center, exact_band):
    dk = np.angle(np.exp(1j*(k-momentum)))
    weights = np.exp(-dk**2/(4*WIDTH**2)-1j*k*center)
    weights /= np.linalg.norm(weights)
    spinor = vectors[:, :, band] if exact_band else bloch(momentum)[2][:, band]
    pulse = np.fft.ifft(weights[:, None]*spinor, axis=0, norm='ortho').ravel()
    return pulse, float(np.sum(abs(weights)**2*velocities[:, band]))


def paired(a, b, plus, minus):
    state = plus[a]*minus[b]
    return state/np.linalg.norm(state)


def band_weights(a, b, state, vectors):
    n = len(vectors)
    full = np.zeros((3*n, 3*n), complex)
    full[a, b] = state
    wave = np.fft.fft(full.reshape(n, 3, 3*n), axis=0, norm='ortho')
    amplitudes = np.einsum('kab,kac->kbc', vectors.conj(), wave)
    return np.sum(abs(amplitudes)**2, axis=(0, 2))


def evolve(h, state, stop, steps):
    return expm_multiply(-1j*h, state, start=0., stop=stop, num=steps, traceA=-1j*h.diagonal().sum())


def direction_releases(n, k, vectors, velocities, a, b, h, exchange):
    source, far = n//4, 3*n//4
    spectator = packet(k, vectors, velocities, 0, 0, far, True)[0]
    offsets = np.stack([((pos//3-center+n//2) % n)-n//2 for pos, center in [(a, source), (b, far)]], axis=1)
    gap = abs(a//3-b//3)
    near = np.minimum(gap, n-gap) <= 1
    setups = []
    for exact in (False, True):
        for sign in (1, -1):
            for band in range(3):
                pulse, speed = packet(k, vectors, velocities, band, sign*np.pi/2, source, exact)
                setups.append((exact, sign, band, speed, paired(a, b, pulse, spectator)))
    # The three local port profiles are independent at each launch momentum.
    overlap = max(float(np.max(abs(np.column_stack([row[-1] for row in setups[i:i+3]]).conj().T
                                   @ np.column_stack([row[-1] for row in setups[i:i+3]]) - np.eye(3))))
                  for i in (0, 3))
    assert overlap < 1e-10
    cut = pair_space(n, cut=True)[2]
    releases = []
    for exact, sign, band, speed, state in setups:
        waves = evolve(h, state, 4, 3)
        probabilities = abs(waves)**2
        np.testing.assert_allclose(probabilities.sum(axis=1), 1, atol=2e-12, rtol=0)
        measured = (probabilities[-1]-probabilities[0]) @ offsets / 4
        p0, p4 = (band_weights(a, b, w, vectors) for w in (waves[0], waves[-1]))
        assert p0[band] > .90 and p4[band] > .90
        assert measured[0]*speed > 0 and abs(measured[0]) > .05
        if exact:
            np.testing.assert_allclose(measured, [speed, 0], atol=1e-6, rtol=0)
        # By the verified charge-exchange symmetry these exactly swapped states
        # solve the other sign's release. Read their individual physical loads.
        swapped = abs(waves[:, exchange])**2
        swapped_offsets = np.stack([((pos//3-center+n//2) % n)-n//2 for pos, center in [(a, far), (b, source)]], axis=1)
        swapped_speed = (swapped[-1]-swapped[0]) @ swapped_offsets / 4
        np.testing.assert_allclose(swapped_speed, measured[::-1], atol=2e-12, rtol=0)
        cut_speed = None
        if not exact:
            cw = evolve(cut, state, 4, 2)
            cut_speed = float(((abs(cw[1])**2-abs(cw[0])**2) @ offsets[:, 0])/4)
            assert abs(cut_speed) < 1e-12
        releases.append({'source': 'band_wavepacket' if exact else 'local_three_port_profile',
                         'band': band+1, 'central_momentum_over_pi': .5*sign,
                         'local_load_speeds_plus_projectile': measured.tolist(),
                         'local_load_speeds_minus_projectile': swapped_speed.tolist(),
                         'isolated_pure_band_prediction': speed,
                         'projectile_band_weights_at_0_and_4': [p0.tolist(), p4.tolist()],
                         'encounter_probability_at_0_2_4': (probabilities @ near).tolist(),
                         'cut_guide_local_cell_speed': cut_speed})
    return {'time_samples': [0, 2, 4], 'momentum_standard_deviation': WIDTH,
            'max_local_source_gram_error': overlap, 'records': releases}


def collisions(n, k, vectors, velocities, a, b, h):
    records = []
    gap = abs(a//3-b//3)
    distance = np.minimum(gap, n-gap)
    for band in (1, 2):
        for approach in (True, False):
            oriented = 1 if approach else -1
            kp = np.sign(bloch(np.pi/2)[3][band])*np.pi/2*oriented
            km = -np.pi/2*oriented
            plus = packet(k, vectors, velocities, band, kp, n//4, True)[0]
            minus = packet(k, vectors, velocities, 0, km, n//4+10, True)[0]
            initial = paired(a, b, plus, minus)
            waves = evolve(h, initial, 12, 7)
            probabilities = abs(waves)**2
            error_norm = float(np.max(abs(probabilities.sum(axis=1)-1)))
            energy = [float(np.vdot(w, h@w).real) for w in waves]
            error_energy = float(np.ptp(energy))
            assert error_norm < 3e-12 and error_energy < 2e-10
            records.append({'projectile_band': band+1, 'direction': 'toward' if approach else 'apart',
                            'momentum_over_pi_plus_minus': [float(kp/np.pi), float(km/np.pi)],
                            'inclusive_projectile_band_weights': [band_weights(a, b, w, vectors).tolist() for w in waves],
                            'mean_minimal_cell_distance': (probabilities @ distance).tolist(),
                            'encounter_probability_within_one_cell': (probabilities @ (distance <= 1)).tolist(),
                            'total_energy': energy[0], 'max_energy_drift': error_energy, 'max_norm_error': error_norm})
    return {'time_samples': [0, 2, 4, 6, 8, 10, 12], 'records': records}


def contact_pulse(wave, port, theta):
    left, right = wave[port], wave[port+1]
    c, s = np.cos(theta), 1j*np.sin(theta)
    phase = np.exp(-1j*theta)
    wave[port] = phase*(c*left+s*right)
    wave[port+1] = phase*(s*left+c*right)


def synthesize(desired, center):
    # The same finite window, independent of source band or direction.
    start, stop = 3*(center-6), 3*(center+7)
    sites = np.arange(start, stop)
    target = np.zeros(len(desired), complex)
    target[sites] = desired[sites]
    retained = float(np.linalg.norm(target)**2)
    target /= np.sqrt(retained)
    wave = np.zeros_like(target)
    wave[start] = 1.
    durations, swaps, holds = [], 0, 0
    for left in sites[:-1]:
        fraction = min(1., abs(target[left])/abs(wave[left]))
        theta = np.arccos(fraction)
        contact_pulse(wave, left, theta)
        durations.append(theta/HOP)
    np.testing.assert_allclose(abs(wave), abs(target), atol=2e-12, rtol=0)
    # Two empty sinks just outside the packet have opposite degrees. Pure
    # contact swap trains transpose source and sink and restore all in between.
    cost = DEGREE_FIVE_COST*(np.arange(len(wave))%3 == 1).astype(int)
    period = 2*np.pi/DEGREE_FIVE_COST
    for port in sites:
        if abs(target[port]) < 1e-15:
            continue
        sink = stop if cost[port] else stop+1
        delta = np.angle(target[port]/wave[port])
        wait = float((delta/(cost[sink]-cost[port])) % period)
        path = [*range(port, sink), *range(sink-2, port-1, -1)]
        # P exp(+i V t) P, followed by exp(-i V t). The inverse hold
        # is a positive-duration hold because every local cost is 0 or 8.
        for hold_time in (period-wait, wait):
            for left in path:
                contact_pulse(wave, left, np.pi/2)
                swaps += 1
            wave *= np.exp(-1j*cost*hold_time)
            holds += 1
            durations.append(hold_time)
    err = float(np.linalg.norm(wave-target))
    assert err < 2e-12
    return wave, {'sharp_seed_port_on_route': start, 'output_ports_first_last': [start, stop-1],
                  'empty_phase_sink_ports': [stop, stop+1], 'retained_probability_before_normalization': retained,
                  'native_fractional_contacts': len(sites)-1, 'native_full_swaps': swaps,
                  'global_degree_holds': holds, 'total_addressed_duration': float(sum(durations) + swaps*np.pi/8),
                  'prepared_vector_l2_error': err}


def addressed_sources(n, k, vectors, velocities, a, b, h):
    center, far = n//4, 3*n//4
    spectator, compensation = synthesize(packet(k, vectors, velocities, 0, 0, far, True)[0], far)
    offsets = np.stack([((p//3-c+n//2)%n)-n//2 for p, c in ((a, center), (b, far))], axis=1)
    out = []
    for sign in (1, -1):
        for band in range(3):
            desired = packet(k, vectors, velocities, band, sign*np.pi/2, center, False)[0]
            pulse, circuit = synthesize(desired, center)
            # The two pulse trains address disjoint intervals. At every contact
            # the spectator is clear of the addressed seam. For every phase
            # gadget the two global holds sum to one full degree period,
            # giving spectator phase exactly one at every address.
            assert circuit['empty_phase_sink_ports'][1] < compensation['sharp_seed_port_on_route']
            init = paired(a, b, pulse, spectator)
            np.testing.assert_allclose(np.linalg.norm(pulse[a]*spectator[b]), 1, atol=2e-12, rtol=0)
            waves = evolve(h, init, 4, 3)
            prob = abs(waves)**2
            np.testing.assert_allclose(prob.sum(axis=1), 1, atol=2e-12, rtol=0)
            weights = [band_weights(a, b, w, vectors) for w in (waves[0], waves[-1])]
            speed = (prob[-1]-prob[0])@offsets/4
            assert min(row[band] for row in weights) > .95
            assert speed[0]*bloch(sign*np.pi/2)[3][band] > 0 and abs(speed[0]) > .05
            out.append({'band': band+1, 'central_momentum_over_pi': .5*sign, 'source_circuit': circuit,
                        'target_band_weight_at_times_0_and_4': [float(row[band]) for row in weights],
                        'local_load_speeds_plus_projectile': speed.tolist()})
    return {'initial_condition': 'one sharp compensated pair with its open arc of oriented heads',
            'addressing': 'one native contact at a time, or the entire native degree term alone',
            'compensator_circuit': compensation, 'releases': out}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output', '--out', dest='out', type=Path, default=HERE/'results.json')
    args = cli.parse_args()
    n, geometry = verify_route()
    k = 2*np.pi*np.fft.fftfreq(n)
    grid = [bloch(q) for q in k]
    vectors = np.stack([row[2] for row in grid])
    velocities = np.stack([row[3] for row in grid])
    a, b, h, exchange = pair_space(n)
    result = {'design_sha256': digest(HERE/'design.md'), 'route_sha256': digest(HERE/'route.json'),
              'executable_sha256': digest(Path(__file__)), 'law': 'written positive contact Laplacian, load-one reference',
              'geometry': geometry, 'energy_unit': 'native generator energy', 'speed_unit': 'route cells per native time',
              'spectrum': report_spectrum(n, k, vectors),
              'releases': direction_releases(n, k, vectors, velocities, a, b, h, exchange),
              'encounters': collisions(n, k, vectors, velocities, a, b, h),
              'addressed_sources': addressed_sources(n, k, vectors, velocities, a, b, h)}
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
