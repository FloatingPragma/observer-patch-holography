"""Certified local readout and a static Coulomb refinement of native records.

The local arithmetic is integer arithmetic. Floating matrices are independent
diagnostics; rational bounds certify the ideal graph comparison, including
fixed-point rounding. No reserved acceptance configuration is executed.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as F
from itertools import product
from pathlib import Path

import numpy as np
import scipy.linalg as la

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BRIDGE = ROOT / 'openwave-tuning/auditor/native_field_bridge'
N20 = ROOT / 'openwave-push/issues/N20_native_green_field/n20lib.py'
FROZEN_MANIFEST = '4597dc9c226d203d890bc3d3a929bbe1be80e2dc57ec03898c77da4fe5b954e1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def dependencies():
    paths = [BRIDGE / 'FINAL_MANIFEST.json', N20,
             ROOT / 'reverse-engineering-reality/model/observer/native.py',
             ROOT / 'reverse-engineering-reality/model/native_spatial_refinement/source.py',
             ROOT / 'reverse-engineering-reality/model/native_spatial_refinement/PROOF.md',
             ROOT / 'reverse-engineering-reality/Lean/Screen/DiscreteCoulombGreen.lean']
    manifest = json.loads((BRIDGE / 'FINAL_MANIFEST.json').read_text())
    assert sha(BRIDGE / 'FINAL_MANIFEST.json') == FROZEN_MANIFEST
    for relative, expected in manifest['files'].items():
        path = BRIDGE / relative
        assert sha(path) == expected, relative
        paths.append(path)
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(paths))}


def sources():
    fixed = load_module('field_feedback', BRIDGE / 'field_feedback.py')
    n20 = load_module('controlled_limit_n20', N20)
    return fixed, n20


def round_ratio(num, den=12):
    """Nearest integer, ties toward positive infinity; error at most 1/2."""
    return (2 * num + den) // (2 * den)


def local_fixed_point(state, charge, neighbours, rounds, bits=40):
    scale = 1 << bits
    out = list(state)
    for _ in range(rounds):
        out = [round_ratio(6 * out[x] + sum(out[int(y)] for y in neighbours[x])
                           + scale * int(charge[x])) for x in range(len(out))]
    return out


def residual_squared_certificate(state, charge, neighbours, bits=40):
    """||B*(phi-G rho)||² <= 2 ||rho-L phi||²; lambda2 >= 1/2."""
    scale = 1 << bits
    residual = [scale * int(charge[x]) - 6 * state[x]
                + sum(state[int(y)] for y in neighbours[x]) for x in range(len(state))]
    assert sum(residual) == 0
    return F(2 * sum(v * v for v in residual), scale * scale)


def rounded_upper(value, places=15):
    scale = 10 ** places
    return str(F(-((-value.numerator * scale) // value.denominator), scale))


def track_path(path, fixed, edges, neighbours, green, lap, rounds=32, bits=40):
    """Only known initial positions and observed currents enter the field."""
    assert path['terminal']['kind'] in ('no_further_event', 'cap_exhausted')
    n = len(neighbours)
    counters = [[int(i == x) for i in range(n)] for x in path['initial']]
    charge = np.array(counters[0]) - np.array(counters[1])
    state = [0] * n
    contraction = F(23, 24) ** rounds
    floor = F(408, 1 << bits)
    bound = F(2)  # sqrt(rho*G rho) <= ||rho||/sqrt(1/2) <= 2.
    maximum_ratio = maximum_error = maximum_gauss = 0.
    source = [tuple(c) for c in counters]
    native = fixed.literal_native()
    snapshots = []
    updates = [None, *path['records']]
    for index, record in enumerate(updates):
        moved = 0
        if record is not None:
            side, seam = record['copy'], record['seam']
            a, b = edges[seam]
            source[side], (_, current) = native(source[side], (a, b), record['coin'])
            assert current == record['current'] == record['observer_current']
            counters[side][a] -= current
            counters[side][b] += current
            assert source[side] == tuple(counters[side])
            moved = abs(current)
            charge = np.array(counters[0]) - np.array(counters[1])
        state = local_fixed_point(state, charge, neighbours, rounds, bits)
        bound = contraction * (bound + moved) + floor * (1 - contraction)
        # Round outward after every block to keep exact integer sizes bounded.
        bound = F(rounded_upper(bound, 18))
        phi = np.array(state, dtype=float) / (1 << bits)
        error_vector = phi - green @ charge
        error = math.sqrt(max(0., float(error_vector @ lap @ error_vector)))
        gauss = float(la.norm(lap @ phi - charge))
        certificate = residual_squared_certificate(state, charge, neighbours, bits)
        assert error <= float(bound) + 2e-13
        assert error * error <= float(certificate) + 1e-24
        assert gauss <= math.sqrt(12) * float(bound) + 2e-13
        maximum_error = max(maximum_error, error)
        maximum_gauss = max(maximum_gauss, gauss)
        maximum_ratio = max(maximum_ratio, error / float(bound))
        if moved or index in (0, len(updates)-1):
            snapshots.append({'after_record_count': index,
                              'actual_field_error_diagnostic': error,
                              'tracking_bound_exact_upper': str(bound),
                              'residual_squared_bound_exact_upper': rounded_upper(certificate, 24)})
    assert counters == path['final_counters']
    return {'initial_positions': path['initial'],
            'source_kappa': path.get('kappa'), 'source_seed': path.get('seed'),
            'native_records': len(path['records']),
            'native_holds': sum(r['current'] == 0 for r in path['records']),
            'native_moves_by_copy': [sum(r['copy']==side and r['current']!=0
                                         for r in path['records']) for side in (0,1)],
            'terminal': path['terminal'], 'local_rounds_per_record': rounds,
            'fixed_point_bits': bits, 'local_rounds': rounds * len(updates),
            'max_field_error_diagnostic': maximum_error,
            'max_Gauss_residual_diagnostic': maximum_gauss,
            'max_error_to_proved_bound_ratio': maximum_ratio,
            'snapshot_certificates': snapshots,
            'final_charge': charge.tolist()}


def transient_checks(edges, neighbours, lap, green):
    """Uncoupled kappa=0 branch only, in the same joint Poisson clock."""
    rho = np.zeros(len(neighbours)); rho[1] = 1; rho[35] = -1
    target = green @ rho
    rows = []
    for time in (0.5, 8., 24.):
        mean_rho = la.expm(-time * lap / 2) @ rho
        finite = green @ (rho - mean_rho)
        difference = finite - target
        error = math.sqrt(float(difference @ lap @ difference))
        bound = 2 * math.exp(-time / 4)
        assert error <= bound + 1e-13
        rows.append({'Poisson_time': time, 'mean_record_field_error': error,
                     'analytic_upper': bound,
                     'formula': '2 exp(-t/4)', 'kappa': 0,
                     'law': 'uncapped expectation; no capped output substituted'})
    return rows


def pi_interval():
    def atan(q, terms):
        value = sum((F((-1)**j, (2*j+1)*q**(2*j+1)) for j in range(terms)), F())
        error = F(1, (2*terms+1)*q**(2*terms+1))
        return value - error, value + error
    a, b = atan(5, 36), atan(239, 10)
    return 16*a[0]-4*b[1], 16*a[1]-4*b[0]


def mode_bounds(size, mode):
    h = F(1, size)
    pi_lo, pi_hi = pi_interval()
    r = 2*pi_hi*h*max(abs(v) for v in mode)
    assert 0 < r <= F(1, 2)
    q = (2*pi_lo*h)**2 * sum(v*v for v in mode) / 8
    remainder = F(3, 8) * r**4
    lower = q - remainder
    assert lower > 0
    optical = 4-r
    defect = r*r/(4-r-r*r/2)**2
    projected = remainder/(q*lower) + defect/lower + defect/optical
    cross = r/(optical*lower)
    # For any signed pair: |u*rho_hat| <= 2/sqrt(12), ||Q rho_hat|| <= 2.
    # 2/sqrt(12) = 1/sqrt(3) < 3/5.
    point_error = h*h*(projected/6 + F(3, 5)*cross)
    uniform_error = h*h*projected/6
    return point_error, uniform_error, projected


def symbol(n20, k, wrong_dimension=False):
    h = 6*np.eye(12, dtype=complex)
    for a, b in n20.frame_seams(): h[a,b] -= 1; h[b,a] -= 1
    for axis, links in enumerate(n20.cubic_links()):
        value = 0. if wrong_dimension and axis == 2 else k[axis]
        phase = np.exp(1j*value)
        for a, b in links: h[a,b] -= phase; h[b,a] -= phase.conjugate()
    return h


def continuum_checks(n20, final_charge):
    occupied = [(i, int(v)) for i, v in enumerate(final_charge) if v]
    assert sum(v for _,v in occupied) == 0 and sum(abs(v) for _,v in occupied) <= 2
    # The finite Fejer product is a declared observation filter, not a preparation.
    modes = [m for m in product((-1,0,1), repeat=3) if any(m)]
    points = tuple(product((0., 1/8, 3/8, 5/8), repeat=3))
    unit = np.ones(12)/math.sqrt(12)
    rows = []
    for size in (16, 32, 64):
        discrete = []; continuum = []; uniform = []
        error_bound = uniform_bound = field_bound = F()
        for mode in modes:
            weight = F(1, 2**sum(v != 0 for v in mode))
            xi = 2*np.pi*np.array(mode); k = xi/size
            rho = np.zeros(12, dtype=complex)
            for port, sign in occupied:
                cell, internal = divmod(port, 12)
                location = np.array(np.unravel_index(cell, (2,2,2)))/2
                rho[internal] += sign*np.exp(-1j*float(xi@location))
            total = rho.sum()
            h = symbol(n20, k)
            exact = float(weight) * (unit@la.solve(h, rho))/(math.sqrt(12)*size**2)
            uniform_source = unit*(unit@rho)
            averaged = float(weight) * (unit@la.solve(h, uniform_source))/(math.sqrt(12)*size**2)
            limit = float(weight) * F(2,3) * total / float(xi@xi)
            discrete.append(complex(exact)); uniform.append(complex(averaged)); continuum.append(complex(limit))
            point, flat, _ = mode_bounds(size, mode)
            error_bound += weight*point; uniform_bound += weight*flat
            # |xi|_2 <= |xi|_1 supplies a rational Euclidean vector bound.
            field_bound += weight*point*2*pi_interval()[1]*sum(abs(v)for v in mode)
        phase_matrix = np.exp(2j*np.pi*np.array(points)@np.array(modes).T)
        observed = phase_matrix@np.array(discrete)
        expected = phase_matrix@np.array(continuum)
        projected = phase_matrix@np.array(uniform)
        error = float(max(abs(observed-expected)))
        flat_error = float(max(abs(projected-expected)))
        derivatives = -2j*np.pi*np.array(modes)
        electric = phase_matrix@(derivatives*np.array(discrete)[:,None])
        electric_limit = phase_matrix@(derivatives*np.array(continuum)[:,None])
        field_error = float(max(la.norm(electric-electric_limit,axis=1)))
        assert error < float(error_bound) and flat_error < float(uniform_bound)
        assert field_error < float(field_bound)
        # A separate dimension intervention uses a genuine z-source mode.
        k = np.array((1.,0.,1.))*2*np.pi/size
        correct = float((unit@la.solve(symbol(n20,k),unit)).real)/size**2
        dim2 = float((unit@la.solve(symbol(n20,k,True),unit)).real)/size**2
        target = 8/(8*np.pi*np.pi)
        _,_,resolvent = mode_bounds(size,(1,0,1))
        assert abs(correct-target) <= float(resolvent)/size**2
        dimension_rejects = bool(abs(dim2-target) > float(resolvent)/size**2)
        # Replacing acoustic response by a gapped optical mode loses its pole.
        optical = np.linalg.eigvalsh(symbol(n20,k))[1]
        wrong_mode = 1/(size**2*optical)
        mode_rejects = bool(abs(wrong_mode-target) > float(resolvent)/size**2)
        assert dimension_rejects and mode_rejects
        rows.append({'cells_per_axis': size, 'matrix_order': 12,
                     'Fourier_modes': len(modes), 'diagnostic_points': len(points),
                     'actual_record_filter_max_potential_error': error,
                     'all_position_supnorm_bound_exact_upper': rounded_upper(error_bound),
                     'actual_record_filter_max_electric_vector_error': field_error,
                     'all_position_electric_bound_exact_upper': rounded_upper(field_bound),
                     'uniform_internal_comparison_error': flat_error,
                     'uniform_internal_supnorm_bound_exact_upper': rounded_upper(uniform_bound),
                     'wrong_dimension_rejected': dimension_rejects,
                     'wrong_optical_mode_rejected': mode_rejects})
    assert rows[-1]['actual_record_filter_max_potential_error'] < rows[0]['actual_record_filter_max_potential_error']
    return rows


def cap_checks():
    receipt = json.loads((BRIDGE/'CAPPED_LAW.json').read_text())
    rows=[]
    for row in receipt['count_controls']:
        if row['copies'] != 2 or row['kappa'] != 0: continue
        mu=576*row['horizon']; cap=row['cap']
        conclusive=math.exp(-mu)*sum(mu**k/math.factorial(k) for k in range(cap))
        assert abs(conclusive-row['success_probability']) < 2e-13
        assert abs(1-conclusive-row['cap_probability']) < 2e-13
        assert row['cap_probability'] > .01
        rows.append({'cap':cap,'native_event_rate':576,'time':row['horizon'],
                     'conclusive':conclusive,'cap_probability':1-conclusive,
                     'omitted_cap_rejected':abs(conclusive-1)>.01})
    assert len(rows)==2
    return rows


def run():
    before=dependencies(); fixed,n20=sources(); edges,nb=fixed.graph()
    assert set(edges)=={tuple(sorted(e)) for e in n20.cubic_graph(2).seams}
    coefficient,_,eig8=n20.cubic_quadratic_coefficient()
    assert eig8 and coefficient==[[F(int(i==j),8)for j in range(3)]for i in range(3)]
    # The exact rational Fourier certificate is reused on its stated size-two domain.
    gap=fixed.exact_gap_certificate()
    lap=6*np.eye(len(nb))
    for a,b in edges: lap[a,b]-=1; lap[b,a]-=1
    green=la.inv(lap+np.ones_like(lap)/len(nb))-np.ones_like(lap)/len(nb)
    frozen=json.loads((BRIDGE/'RECIPROCAL_FEEDBACK.json').read_text())
    selected=[p for p in frozen['complete_paths']
              if p['initial']==[1,35] and p['kappa'] in (0.,1.) and p['seed']==804]
    selected.append(frozen['cap_failure'])
    assert len(selected)==3
    tracking=[track_path(p,fixed,edges,nb,green,lap)for p in selected]
    # Use a nonzero actual present-charge record, including its internal port labels.
    observed=next(t['final_charge']for t in tracking if any(t['final_charge']))
    continuum=continuum_checks(n20,observed)
    result={'schema':'oph.controlled_native_field_limit.v1',
      'checks':{'frozen_reciprocal_source_verified':True,
                'same_literal_native_records_replayed':True,
                'fixed_point_rounding_has_rational_error_bound':True,
                'local_tracking_and_Gauss_bounds':True,
                'uniform_native_transient_bound':True,
                'cubic_static_Fourier_refinement':True,
                'dimension_mode_cap_mutations_rejected':True,
                'source_bytes_unchanged':dependencies()==before},
      'source_pins':before,'side_two_gap_certificate':gap,
      'local_tracking':tracking,'uncoupled_transient':transient_checks(edges,nb,lap,green),
      'static_refinement':continuum,'complete_cap_controls':cap_checks(),
      'scope':{'source':'two literal classical one-load copies, observed signed charge',
        'clock':'frozen joint Poisson source time; local arithmetic rounds are separately counted readout work',
        'preparation':'exact frozen known initial positions; no uniform creation prior is assumed',
        'continuum_domain':'unit periodic three-torus; Fejer-filtered nonzero modes |m|_infinity<=1',
        'normalization':'-Delta Phi=(2/3)rho; infinite-space fundamental coefficient 1/(6 pi)',
        'refinement_transfer':'static charge lifting preserves physical cell coordinates and port labels; no trajectory/clock equivalence between different carrier sizes is asserted',
        'C10_transfer':False,'Maxwell_waves':False,'fine_structure_constant_derived':False,
        'closed_Hamiltonian_force_dynamics':False,'official_rows_changed':[]}}
    assert all(result['checks'].values())
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'tracking_paths':len(tracking),
                      'native_records':sum(t['native_records']for t in tracking),
                      'refinement':continuum,'cap':result['complete_cap_controls']},indent=2))


if __name__=='__main__': run()
