"""Ordinary unguarded native energy of a daughter enrolled on a decaying source."""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model.carrier import federation

from ...row_20_electric import experiment as prior
from .apparatus import Guide
from .emitter import launch, pure_overlap, reflected, selected_parent
from .mode_reader import analyzer

SHOTS = 2**42
LEVELS = np.array([0, 8, 9, 32, 33])


def calibrate():
    h = law.seam_term(8)
    e, v = np.linalg.eigh(h)
    times = 2*np.pi*np.arange(64)/64
    bins = np.arange(64)
    amp = np.exp(1j*np.outer(bins, times)) @ np.exp(-1j*np.outer(times, e))/64
    np.testing.assert_allclose(abs(amp)**2, (bins[:,None] == np.rint(e)), atol=2e-13)
    pairs = {}
    for ab in [(0,0),(0,1),(1,2),(0,2)]:
        start = law.SRC.index(ab)
        pairs[str(ab)] = np.round(abs(amp[LEVELS])**2 @ abs(v[start])**2, 13).tolist()
    for sign in (1, -1):
        head = (v[law.SRC.index((0,2))]+sign*v[law.SRC.index((2,0))])/np.sqrt(2)
        pairs[f'balanced {sign:+}'] = np.round(abs(amp[LEVELS])**2 @ abs(head)**2, 13).tolist()
    return pairs


def pfaffian(matrix):
    """Pivoted elimination of a real skew-symmetric matrix in its original order."""
    a = matrix.copy()
    result = 1.
    for i in range(0, len(a), 2):
        pivot = i + 1 + np.argmax(abs(a[i, i+1:]))
        if a[i, pivot] == 0:
            return 0.
        if pivot != i + 1:
            a[[i+1,pivot]] = a[[pivot,i+1]]
            a[:,[i+1,pivot]] = a[:,[pivot,i+1]]
            result = -result
        entry = a[i,i+1]
        result *= entry
        u, v = a[i,i+2:], a[i+1,i+2:]
        a[i+2:,i+2:] += (np.outer(v,u)-np.outer(u,v))/entry
    return result


def local_bins(cov, u, v):
    """u,v are the positions of both physical ports within the selected path."""
    lo, hi = sorted((u, v))
    j, k = lo//2, hi//2
    if j == k:
        x = cov[2*j,2*j+1]
        return np.array([0, 0, 0, (1+x)/2, (1-x)/2])
    block = cov[2*j+1:2*k+1,2*j+1:2*k+1]
    zz = pfaffian(block)
    unlike = (1 - (-1)**(lo%2+hi%2)*zz)/2
    return np.array([1-unlike, 0, 0, unlike/2, unlike/2])


def enroll(g):
    r = np.linalg.matrix_power(g.r, 4*g.n)
    start = launch(Guide(g.n,'native',0),'parent')
    pg, post = selected_parent(r @ start @ r.T, 1)
    one = reflected(g.c, g.modes[:,0].real)
    pm = pure_overlap(post, one)
    return one, np.array([1-pg, pg*(1-pm), pg*pm])


def measure_contact(p, admission, rng):
    refuse_ground, refuse_mode, admitted = admission
    probs = np.concatenate(([refuse_ground, refuse_mode], admitted*p))
    if min(probs) < -2e-10:
        raise ValueError('invalid local spectral weight')
    probs = probs.clip(0.)
    counts = rng.multinomial(SHOTS, probs/probs.sum())
    m = int(sum(counts[2:]))
    if m == 0:
        raise ArithmeticError('no admitted specimens')
    return counts.tolist(), float(LEVELS @ (counts[2:]/m)), m


def acquire_mount(board, ports, mount, daughter, host, admission, rng):
    slot = {port: i for i, port in enumerate(ports)}
    path_edges = list(zip(ports, ports[1:]))
    live = path_edges if mount == 'path' else board.seams
    rows = []
    measured = exact = error_sum = escape = 0.
    mixed = outside = 0
    for a, b in live:
        if a not in slot and b not in slot:
            outside += 1
            continue
        if a not in slot or b not in slot:
            mixed += 1
            continue
        u, v = slot[a], slot[b]
        signal = local_bins(daughter,u,v)
        blank = local_bins(host,u,v)
        cs, es, ns = measure_contact(signal,admission,rng)
        cb, eb, nb = measure_contact(blank,(0.,0.,1.),rng)
        measured += es-eb
        exact += LEVELS @ (signal-blank)
        width = 1 if u//2 == v//2 else 33
        error_sum += width**2*(1/ns+1/nb)
        if mount == 'path' and u//2 != v//2:
            escape += (signal[-1]+signal[-2])/4
        rows.append(dict(ports=[int(a),int(b)], path_positions=[u,v], daughter=cs, host=cb))
    r = float(np.sqrt(np.log(2*4/.01)/2*error_sum))
    weights = 6-board.degree
    outside_refs = 0 if mount == 'path' else int(8*sum(weights[v] for v in range(board.n) if v not in slot))
    for a,b in zip(ports[::2],ports[1::2],strict=True):
        assert weights[a] == weights[b]
    return dict(mount=mount, measured_increment=measured, exact_increment=exact,
                joint_99_radius=r, interval=[measured-r,measured+r], contacts=rows,
                unchanged_data_host_contacts=mixed, unchanged_host_host_contacts=outside,
                common_external_contact_mean=8.5*mixed, unchanged_outside_private_energy=outside_refs,
                unchanged_dimer_private_energy=int(32*sum(weights[a] for a in ports[::2])),
                path_departure_quadratic_coefficient=escape if mount == 'path' else None)


def low_weight_certificate(n, measured_mount):
    """Block-separation constants depend only on the physical graph and fixed trial cut."""
    cuts = 0
    links = Counter()
    for term in measured_mount['contacts']:
        u, v = sorted(term['path_positions'])
        cuts += int(u < n <= v)
        if u//2 != v//2:
            links[u//2, v//2] += 1
    mixed = measured_mount['unchanged_data_host_contacts']
    numerator = (mixed + sum(m*m for m in links.values()))/4
    separation = 32*n - 32.5*cuts - .1
    if separation <= 0:
        raise ValueError('fixed block trial fails to separate the spectral intervals')
    return dict(trial_cut_edges=cuts, inter_dimer_multiplicities=dict(Counter(links.values())),
                mixed_edges=mixed, coupling_norm_squared_upper=numerator,
                separated_gap_lower=separation, low_weight_upper=numerator/separation**2)


def acquisition():
    board = federation(2)
    _, path, _ = prior.guide(board)
    rng = np.random.default_rng(90129)
    panels = []
    for n in (23,45):
        ports = path[1:1+2*n]
        g = Guide(n,'native',.5)
        one, admission = enroll(g)
        reader, vacuum_bits, holds = analyzer(g)
        single_bits = vacuum_bits.copy()
        single_bits[0] *= -1
        observed = reader @ one @ reader.T
        np.testing.assert_allclose(observed[::2,1::2],np.diag(single_bits),atol=3e-12)
        e = g.phi/g.period
        panels.append(dict(dimers=n, native_mode_reader=dict(hold_count=len(holds),
                           vacuum_word=vacuum_bits.tolist(), daughter_word=single_bits.tolist(),
                           inverse_restoration_error=float(np.abs(reader.T @ observed @ reader-one).max())),
                           lowest_two_driven_lines=e[:2].tolist(),
                           driven_cutoff=.1/(1+2/65), one_mode_admission=admission.tolist(),
                           path_low_weight_upper=(n-1)/(4*(32*n-32.5-.1)**2),
                           mounts=[acquire_mount(board,ports,m,one,g.c,admission,rng)
                                   for m in ('path','full')]))
    for panel in panels:
        for mount in panel['mounts']:
            mount['static_ground_spectral_bound'] = low_weight_certificate(panel['dimers'], mount)
    return dict(seed=90129, original_starts_per_contact_and_arm=SHOTS, instrument=calibrate(),
                clock_levels=LEVELS.tolist(), clock_bins=64,
                count_order=['parent refused','one-mode refused']+LEVELS.tolist(), panels=panels)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    parser.parse_args().output.write_text(json.dumps(acquisition(),indent=1,allow_nan=False)+'\n')
