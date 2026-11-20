"""Independent finite cross-check of the primary A-contrast-memory-spherical-v1 contrast law.

This law is separate from the scalar CAR source. Its complete 99-state seam
law and reducing blank-memory sector are checked before a two-load cell run.
"""
from fractions import Fraction
from itertools import combinations_with_replacement
import numpy as np
from scipy.linalg import eigh, norm
from .scalar_channel.carrier import load_source


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def seam(kappa=8.0):
    require(np.isfinite(kappa) and kappa >= 0, 'Contrast coefficient must be nonnegative')
    W = np.zeros((90, 9))
    for a in range(3):
        for b in range(3):
            s, x = a+b, a-b
            for c in (0, 1):
                y = (s//2+c*(s % 2), s//2+(1-c)*(s % 2))
                record = c*5 + x+2
                W[record*9 + 3*y[0]+y[1], 3*a+b] = 1/np.sqrt(2)
    U = np.block([[np.zeros((9, 9)), W.T], [W, np.eye(90)-W@W.T]])
    R = np.zeros((10, 10))
    for c in (0, 1):
        for x in range(-2, 3):
            R[5*c-x+2, 5*c+x+2] = 1
    G = np.zeros((99, 99))
    G[9:, 9:] = np.kron((np.eye(10)-R)/2 +
                         kappa*np.diag([x*x for c in (0, 1) for x in range(-2, 3)]), np.eye(9))
    H = U @ G @ U
    swap = np.zeros((9, 9))
    for a, b in ((a,b) for a in range(3) for b in range(3)):
        swap[3*b+a, 3*a+b] = 1
    expected = (np.eye(9)-swap)/2+kappa*np.diag([(a-b)**2 for a in range(3) for b in range(3)])
    require(norm(W.T@W-np.eye(9), 2)<1e-12, 'Write loses a coherent input')
    require(norm(U@U-np.eye(99), 2)<1e-12, 'Complete seam extension is not involutive')
    require(norm(H[:9,:9]-expected, 2)<1e-12, 'Blank-memory action differs')
    require(norm(H[9:,:9], 2)<1e-12, 'Blank memory leaks under free evolution')
    return H, W, U


def cell(kappa):
    require(np.isfinite(kappa) and kappa >= 0, 'Contrast coefficient must be nonnegative')
    graph = load_source('cell')
    # The whole qutrit total-two sector includes twelve double occupancies.
    # They reduce under site swaps, but remain outcomes of the same detector.
    states = list(combinations_with_replacement(range(12), 2))
    index = {s:i for i,s in enumerate(states)}
    edges = {tuple(sorted(map(int, edge))) for edge in graph['edges']}
    H = np.zeros((len(states), len(states)))
    for col, state in enumerate(states):
        for a, b in edges:
            moved = tuple(sorted(b if v==a else a if v==b else v for v in state))
            H[col,col] += .5
            H[index[moved],col] -= .5
        occupation = np.bincount(state, minlength=12)
        contrast = 6 * int(occupation @ occupation) - 2 * sum(
            int(occupation[a] * occupation[b]) for a, b in edges)
        H[col,col] += kappa * contrast
    return H, states, np.array([s in edges for s in states])


def family_pair_bound(kappa):
    """Exact all-time theorem for the declared compensated graph family.

    Complete icosahedral cells joined by an external port partial matching,
    with d0=6, have adjacent/nonadjacent blocks A <= (10*kappa+6) I,
    C >= 12*kappa I and ||B||**2 <= 10. The last bound follows from
    cross-block absolute column sums <=5 and row sums <=2; the first uses
    the internal-pair line graph's minimum adjacency eigenvalue -2.
    For kappa>3 the ordered gap gamma=2*kappa-6 is positive. The spectral
    subspace Sylvester bound gives ||Q0 Pi||, ||Q P0|| <= sqrt(10)/gamma.
    Splitting exp(-itH) at Pi bounds every initially adjacent state's escape
    probability by 40/gamma**2, uniformly in t and with any reference.
    The twelve double occupancies belong to the complementary block.

    This is the imported contrast-law theorem, not a scalar-CAR transfer or
    a conclusion inferred from the finite time grid in run().
    """
    value = Fraction(str(kappa))
    require(value > 3, 'The ordered pair-band bound requires kappa > 3')
    gamma = 2 * value - 6
    return min(Fraction(1), Fraction(40) / gamma**2)


def run():
    H, W, U = seam()
    values = eigh(H, eigvals_only=True)
    checks = dict(write_isometry_error=float(norm(W.T@W-np.eye(9), 2)),
                  full_involution_error=float(norm(U@U-np.eye(99), 2)),
                  blank_memory_leakage=float(norm(H[9:,:9], 2)),
                  full_seam_spectrum=values.tolist())
    records = []
    for kappa in (8.0, 0.0):
        h, states, adjacent = cell(kappa)
        vals, vecs = eigh(h)
        initial = states.index((0, 1))
        outputs = []
        for t in (0., .125, .5, 1., 2., 4., 8.):
            psi = (vecs*np.exp(-1j*t*vals))@vecs[initial,:].conj()
            weights = abs(psi)**2
            total = float(weights.sum())
            require(np.isfinite(weights).all() and abs(total-1)<1e-12,
                    'Incomplete two-load instrument')
            # Normalize numerical state norm; retain the unmodified weights
            # and their error so this cannot hide a failed evolution.
            law = weights / total
            outputs.append(dict(action=t, probabilities=law.tolist(),
                                unnormalized_weights=weights.tolist(),
                                state_norm_error=abs(total-1),
                                numerical_normalization=total,
                                no_adjacent_pair_probability=float(law[~adjacent].sum())))
        if kappa == 8:
            require(max(r['no_adjacent_pair_probability'] for r in outputs) <= .4,
                    'Pair-persistence screen failed')
        records.append(dict(kappa=kappa, eigenvalues=vals.tolist(), outcomes=states, records=outputs))
    # Meaningful alternative law, same preparation and complete detector.
    require(max(r['no_adjacent_pair_probability'] for r in records[1]['records']) > .4,
            'Same-apparatus zero-contrast control was not distinguished')
    bound = family_pair_bound(8)
    return dict(schema='oph.contrast_memory_replay.v1', construction='A-contrast-memory-spherical-v1',
                seam=checks, cell=records,
                scope='one complete seam; two-unit-load invariant sector of one cell',
                all_time_bound=dict(kappa=8, numerator=bound.numerator,
                                    denominator=bound.denominator,
                                    domain='initially adjacent total-two qutrit states; d0=6; '
                                           'complete icosahedral cells with external port partial matching',
                                    hypothesis='kappa > 3',
                                    source='model/contrast_memory_spherical_v1/PROOF.md',
                                    grid_role='finite numerical screen, not the all-time proof'),
                source_join='blank-memory reducing operator equality; no scalar-CAR transfer')


def seam_orientation_chart(reverse=False):
    """Permutation from old to relabelled ordered seam coordinates.

    Reversing endpoints swaps both loads and sends a nonblank record
    (coin, contrast) to (1-coin, -contrast). Blank memory stays blank.
    """
    require(type(reverse) is bool, 'Seam orientation must be boolean')
    indices = np.arange(99)
    if reverse:
        for a in range(3):
            for b in range(3):
                indices[3*a+b] = 3*b+a
                for c in (0, 1):
                    for d in range(-2, 3):
                        old = 9 + 9*(5*c+d+2) + 3*a+b
                        new = 9 + 9*(5*(1-c)-d+2) + 3*b+a
                        indices[old] = new
    require(np.array_equal(np.sort(indices), np.arange(99)), 'Incomplete seam chart')
    return indices


def primary_join():
    """Full seam and total-two-sector transfer to the primary contrast code.

    This is a relabelling and reference shift within the contrast law only.
    It neither identifies the scalar CAR dynamics nor transfers its clocks.
    """
    from .contrast_memory_spherical_v1 import carrier, law, sectors
    port_map = (0, 1, 5, 7, 11, 9, 10, 8, 4, 6, 2, 3)
    original = load_source('cell')
    mapped = {tuple(sorted((port_map[a], port_map[b]))) for a,b in original['edges']}
    require(mapped == set(carrier.CELL_SEAMS), 'Carrier labels do not intertwine')
    charts = {name: seam_orientation_chart(reverse) for name, reverse in
              (('preserved', False), ('reversed', True))}
    seam_charts = []
    for a, b in original['edges']:
        aa, bb = port_map[a], port_map[b]
        seam_charts.append(dict(source_edge=[int(a), int(b)],
                               primary_edge=sorted((aa, bb)),
                               orientation='reversed' if aa > bb else 'preserved'))
    require(any(item['orientation'] == 'reversed' for item in seam_charts),
            'Reversed seam control has no applicable edge')
    reversed_chart = np.eye(99)[:, charts['reversed']]
    load_swap = np.zeros((9, 9))
    for a in range(3):
        for b in range(3):
            load_swap[3*b+a, 3*a+b] = 1
    incomplete_chart = np.kron(np.eye(11), load_swap)
    p = sectors.PairSector(carrier.cell())
    checks = []
    for kappa in (0, 8):
        h, states, _ = cell(kappa)
        primary = np.zeros_like(h)
        primary[:66,:66] = p.hamiltonian(kappa).toarray()+12*kappa*np.eye(66)
        primary[66:,66:] = sectors.doubled_hamiltonian(carrier.cell(), kappa).toarray()+12*kappa*np.eye(12)
        indices = []
        for a,b in states:
            aa,bb = sorted((port_map[a],port_map[b]))
            indices.append(p.index(aa,bb) if aa!=bb else 66+aa)
        error = float(norm(h-primary[np.ix_(indices,indices)],2))
        require(error == 0, 'Complete contrast sector transfer failed')
        native, _, _ = seam(kappa)
        h0,h1 = law.local_term_parts()
        seam_error = float(norm(native-h0-kappa*h1,2))
        require(seam_error < 1e-12, 'Full memory complement law differs')
        reversed_error = float(norm(native - reversed_chart.T @ (h0+kappa*h1)
                                    @ reversed_chart, 2))
        require(reversed_error < 1e-12, 'Reversed raw-memory seam does not intertwine')
        omitted_memory_error = float(norm(native - incomplete_chart.T @ (h0+kappa*h1)
                                          @ incomplete_chart, 2))
        require(omitted_memory_error > .5,
                'Load-only relabelling control failed to distinguish raw memory')
        checks.append(dict(kappa=kappa, full_78_state_operator_error=error,
                           full_99_state_seam_error=seam_error,
                           full_99_reversed_seam_error=reversed_error,
                           load_only_orientation_control_error=omitted_memory_error,
                           primary_energy_shift=12*kappa))
    return dict(port_map=list(port_map), all_78_outcomes_mapped=True,
                all_99_memory_states_retained=True, checks=checks,
                seam_charts=seam_charts,
                seam_coordinate_permutations={name: indices.tolist()
                                             for name, indices in charts.items()},
                reversed_record_map='(coin,contrast) -> (1-coin,-contrast)',
                scalar_law_equivalence_claimed=False)


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2))
