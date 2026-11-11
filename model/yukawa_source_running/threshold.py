"""Pauli-blocking thresholds, with complete source-carrier controls.

The no-hopping reduction is exact for the declared finite Higgs join. Its
addition frequencies are not Yukawa coefficients or physical particle masses.
"""
from __future__ import annotations

from functools import lru_cache
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]


def exact_source_commutators():
    """Every column of the rational retained operator packet, independently."""
    base = ROOT/'model/interacting_extension/higgs_join'
    packet = json.loads((base/'joined_exact_operators.json').read_text())['operators']
    basis = json.loads((base/'joined_basis.json').read_text())
    norms = {q:list(map(Fraction,basis[str(q)]['raw_norms'])) for q in (0,1)}

    def matrix(name, target):
        return {(i,j):Fraction(value)/norms[target][i] for i,j,value in packet[name]}

    def multiply(a,b):
        rows={}
        for (i,j),value in b.items():
            rows.setdefault(i,[]).append((j,value))
        result={}
        for (i,k),value in a.items():
            for j,other in rows.get(k,()):
                result[i,j]=result.get((i,j),Fraction(0))+value*other
        return {key:value for key,value in result.items() if value}

    def equal(a,b):
        if a != b:
            raise ArithmeticError('exact complete source intertwiner changed')

    loops={q:matrix(f'q{q}_loop',q) for q in (0,1)}
    loops_adjoint={q:{(j,i):value*norms[q][i]/norms[q][j]
                      for (i,j),value in loops[q].items()} for q in (0,1)}
    electric={q:{(i,i):Fraction(value) for i,value in enumerate(
                  packet[f'q{q}_diagonal']['electric_native_trace']) if Fraction(value)} for q in (0,1)}
    identities=0
    for q in (0,1):
        nh=list(map(int,packet[f'q{q}_diagonal']['N_H']))
        for left in range(3):
            for right in range(3):
                w=matrix(f'q{q}_W_{left}_{right}',q)
                if any(nh[i]-nh[j] != 1 for i,j in w):
                    raise ArithmeticError('literal scalar-raising Higgs operator')
                for field in (electric[q],loops[q],loops_adjoint[q]):
                    equal(multiply(field,w),multiply(w,field));identities+=1
    for family in range(3):
        for spin in range(2):
            b=matrix(f'addition_{family}_{spin}',1)
            for fields in (electric,loops,loops_adjoint):
                equal(multiply(fields[1],b),multiply(b,fields[0]));identities+=1
    return {'complete_exact_operator_identities':identities,
            'Gauss_dimensions':{str(q):len(norms[q]) for q in (0,1)},
            'all_Higgs_family_entries_and_addition_probes_retained':True,
            'field_commutes_with_Higgs_exact':True,'addition_intertwines_field_exact':True,
            'Higgs_raises_scalar_occupation_exact':True}


def source_spectral_bounds(y):
    """Exact Sylvester inertias bracket each source eigenvalue."""
    brackets=((sp.Rational(264,10**6),sp.Rational(265,10**6)),
              (sp.Rational(516,10**4),sp.Rational(517,10**4)),
              (sp.Rational(948,1000),sp.Rational(9481,10000)))
    counts=[]
    for j,pair in enumerate(brackets):
        row=[]
        for shift in pair:
            previous=sp.Integer(1);negative=0
            for k in (1,2,3):
                minor=sp.simplify((y-shift*sp.eye(3))[:k,:k].det())
                pivot=sp.simplify(minor/previous)
                if pivot.is_positive is True:
                    pass
                elif pivot.is_negative is True:
                    negative+=1
                else:
                    raise ArithmeticError('unresolved exact inertia sign')
                previous=minor
            row.append(negative)
        if row != [j,j+1]:
            raise ArithmeticError('source eigenvalue bracket changed')
        counts.append(row)
    lower=[sp.factor((brackets[j+1][0]/brackets[j][1])**2) for j in (0,1)]
    return {'eigenvalue_brackets':[[str(x) for x in b] for b in brackets],
            'negative_inertias_at_endpoints':counts,
            'no_hopping_adjacent_threshold_ratio_lower_exact':list(map(str,lower)),
            'no_hopping_adjacent_threshold_ratio_lower_numeric':list(map(float,lower)),
            'domain':'g != 0, scalar > 0, reference offset zero, hopping zero'}


def _create(mask, mode):
    if mask & (1 << mode):
        return None
    return mask | (1 << mode), (-1) ** ((mask & ((1 << mode)-1)).bit_count())


@lru_cache(maxsize=1)
def exact_local_car():
    """Independent literal CAR calculation on vacuum and all six e inputs.

    Scalar labels are orthogonal, both weak components and both spin terms
    are retained. The source tensor has one unit charged-lepton entry per
    Higgs component. The signs are recovered from the actual tensor packet.
    """
    path = ROOT / 'model/observer_particles/results/unified_yukawa_operators.json'
    vertices = json.loads(path.read_text())['vertices']
    vertex = next(v for v in vertices if tuple(v['constituent_q6Y']) == ('-3', '6'))
    entries = [v['entries'] for v in vertex['coefficient_matrices']]
    supports = [{int(i) for i,j,v in part} for part in entries]
    right, = supports[0] & supports[1]
    gamma = []
    for part in entries:
        value, = [sp.sympify(v) for i,j,v in part if j == right]
        if abs(value) != 1:
            raise ArithmeticError('literal unit lepton Higgs tensor changed')
        gamma.append(value)
    yy = sp.Matrix(3, 3, sp.symbols('y0:9', complex=True))

    def acted(mask):
        out = {}
        for h in range(2):
            for left in range(3):
                for right_family in range(3):
                    for left_spin, right_spin, epsilon in ((0,1,1),(1,0,-1)):
                        one = _create(mask, 2*right_family+right_spin)
                        if one is None:
                            continue
                        two = _create(one[0], 6+6*h+2*left+left_spin)
                        if two is not None:
                            key = h, two[0]
                            out[key] = out.get(key, 0)+gamma[h]*epsilon*one[1]*two[1]*yy[left,right_family]
        return out

    def inner(a, b):
        return sp.expand(sum(sp.conjugate(value)*b.get(key, 0) for key,value in a.items()))

    vacuum = acted(0)
    inputs = [acted(1 << mode) for mode in range(6)]
    gram = sp.Matrix(6, 6, lambda i,j: inner(inputs[i], inputs[j]))
    total = sp.trace(yy.H*yy)
    expected = 4*total*sp.eye(6)-2*sp.kronecker_product((yy.H*yy).T, sp.eye(2))
    if sp.expand(inner(vacuum, vacuum)-4*total) != 0:
        raise ArithmeticError('vacuum CAR normalization')
    if any(sp.expand(x) != 0 for x in gram-expected):
        raise ArithmeticError('complete one-particle blocking Gram')
    return {'arbitrary_complex_coefficients': 9, 'one_particle_gram_entries_checked': 36,
            'weak_components': 2, 'spin_terms_per_component': 2,
            'source_tensor_signs': list(map(str, gamma)),
            'vacuum_Gram': '4 Tr(Ydagger Y)',
            'one_particle_Gram': '4 Tr(Ydagger Y) I6 - 2 (Ydagger Y)^T tensor I2',
            'identity_exact': True}


def analytic_records(yukawa, coupling=1., scalar=1., reference=0.):
    """Two bright lines per family; every fixed-coordinate probe is retained."""
    y = np.asarray(yukawa, complex)
    controls=(coupling,scalar,reference)
    if any(isinstance(x,(bool,np.bool_)) or not isinstance(x,(int,float,np.integer,np.floating))
           or not np.isfinite(x) for x in controls):
        raise ValueError('finite real coupling, scalar energy and reference required')
    if y.shape != (3,3) or not np.isfinite(y).all() or not scalar > 0:
        raise ValueError('finite 3x3 coefficient and positive scalar energy required')
    k = (y.conj().T@y).T
    squared, vectors = np.linalg.eigh(k)
    squared = np.maximum(squared, 0)
    total = float(np.trace(k).real)
    d0 = np.sqrt(scalar*scalar+16*coupling*coupling*total)
    if not np.isfinite(d0):
        raise ValueError('finite resolvable threshold discriminant required')
    a0 = np.sqrt((1+scalar/d0)/2)
    b0 = np.sqrt((1-scalar/d0)/2)
    arrays = {'squared_singular_values': squared}
    rows, residues, frequencies = [], [], []
    for j, value in enumerate(squared):
        ci_squared = 4*total-2*value
        di = np.sqrt(scalar*scalar+4*coupling*coupling*ci_squared)
        ai, bi = np.sqrt((1+scalar/di)/2), np.sqrt((1-scalar/di)/2)
        relative = np.sqrt(ci_squared/(4*total)) if total else 1.
        # Rationalized subtraction resolves small positive frequencies.
        low = 4*coupling*coupling*value/(d0+di)+reference
        high = (d0+di)/2+reference
        low_weight = (a0*ai+b0*bi*relative)**2
        high_weight = (-a0*bi+b0*ai*relative)**2
        blocked = b0*b0*value/(2*total) if total else 0.
        projector = np.kron(np.outer(vectors[:,j], vectors[:,j].conj()), np.eye(2))
        for frequency, weight in ((low,low_weight),(high,high_weight)):
            frequencies.append(frequency)
            residues.append(weight*projector)
        rows.append({'family': j, 'squared_singular_value': float(value),
                     'low_frequency': float(low), 'upper_frequency': float(high),
                     'low_weight_per_spin': float(low_weight),
                     'upper_weight_per_spin': float(high_weight),
                     'blocked_probability_per_spin': float(blocked)})
    arrays.update(frequencies=np.asarray(frequencies), residues=np.asarray(residues),
                  blocked_effect=np.eye(6)-np.sum(residues, axis=0))
    if np.max(abs(np.diag(arrays['blocked_effect']).imag)) > 1e-13:
        raise ArithmeticError('complex blocked probability')
    return {'source_norm_squared': total, 'vacuum_discriminant': float(d0),
            'rows': rows, 'complete_success_plus_blocked_max': float(max(
                abs(r['low_weight_per_spin']+r['upper_weight_per_spin']+
                    r['blocked_probability_per_spin']-1) for r in rows))}, arrays


def correlation(frequencies, residues, time):
    return np.einsum('n,nab->ab', np.exp(-1j*time*frequencies), residues)


def reconstruct_coupling(paired_frequencies, blocked_effect):
    """Infer |g| singular(Y) and a common reference offset from full records.

    The input pairs are resolved lower/upper lines with the same family
    residue projector. No source matrix, g, scalar energy or reference shift
    is an input. Pair resolution and the stated finite model are premises.
    """
    pairs=np.asarray(paired_frequencies,float).reshape(3,2)
    blocked=np.asarray(blocked_effect,complex)
    if (blocked.shape!=(6,6) or not np.isfinite(pairs).all()
            or not np.isfinite(blocked).all()):
        raise ValueError('finite paired lines and six complete probe failures required')
    if np.max(abs(blocked-blocked.conj().T))>2e-11:
        raise ValueError('a Hermitian blocked effect is required')
    eigenvalues=np.linalg.eigvalsh(blocked)
    if eigenvalues.min()<-2e-11 or eigenvalues.max()>1+2e-11:
        raise ValueError('blocked effect must lie between zero and identity')
    trace=float(np.trace(blocked).real)
    kappa=1-2*trace
    differences=pairs[:,1]-pairs[:,0]
    if np.any(differences<=0) or not -1e-12<=kappa<=1+1e-12:
        raise ValueError('outside positive-scalar resolved two-line domain')
    d0_squared=2*np.dot(differences,differences)/(5+kappa*kappa)
    d0=np.sqrt(d0_squared)
    squares=(d0_squared-differences*differences)/8
    if squares.min() < -1e-11:raise ValueError('inconsistent spectral and failure records')
    offsets=(pairs.sum(axis=1)-d0)/2
    if np.ptp(offsets)>1e-10:raise ValueError('reference shift is not common')
    return {'coupling_times_singular_values':np.sqrt(np.maximum(squares,0)).tolist(),
            'scalar_energy':float(kappa*d0),'reference_energy':float(offsets.mean()),
            'input_source_matrix_or_coupling_or_scalar_energy':False,
            'weak_upper_line_and_all_blocked_outcomes_required':True,
            'finite_shot_error_bound_provided':False}


def complete_control(y, sigma, *, coupling=1., hopping=0., reference=0.):
    """Recompute both complete carriers; compare fixed-coordinate record laws."""
    from ..interacting_extension.higgs_join.interacting_carrier import build_join
    from ..interacting_extension.higgs_join.run_joined_response import block_eigh, addition_data, maxabs
    joined = build_join(family0=y, family1=sigma, yukawa=y, higgs_coupling=coupling,
                        hopping=hopping, reference_energy=reference)
    eigen = {q: block_eigh(joined['sectors'][q]) for q in (0,1)}
    arrays = {}
    record = addition_data(joined, eigen, 'complete', arrays)
    for q in (0,1):
        arrays[f'q{q}_all_energies'] = eigen[q][0]
        arrays[f'q{q}_all_spin_labels'] = eigen[q][2]
    gram = arrays['complete_addition_probe_gram']
    blocked = np.eye(6)-gram
    arrays['complete_blocked_effect'] = blocked
    if np.linalg.eigvalsh(blocked).min() < -2e-11:
        raise ArithmeticError('creation completion is not positive')
    field_errors, parity_errors = [], []
    for q,x in joined['sectors'].items():
        w, field = x['W'], x['H_field']
        field_errors.append(maxabs(field@w-w@field))
        sign = np.asarray(1-2*x['observables']['N_H'].diagonal().real)
        conjugate = x['H'].multiply(sign[:,None]).multiply(sign[None,:])
        parity_errors.append(maxabs(conjugate-(x['H']-2*x['H_Higgs'])))
    result = {'coupling': coupling, 'hopping': hopping, 'reference_energy': reference,
              'sectors': {str(q): eigen[q][3] for q in (0,1)},
              'complete_addition_record': record,
              'field_Higgs_commutator_max': max(field_errors),
              'Higgs_sign_parity_conjugation_max': max(parity_errors),
              'lowest_six_addition_frequencies': arrays['complete_addition_frequencies'][:6].tolist(),
              'blocked_probabilities_fixed_probes': np.diag(blocked).real.tolist()}
    if max(field_errors+parity_errors) > 2e-12:
        raise ArithmeticError('source field or parity identity')
    if hopping == 0:
        analytic, aa = analytic_records(y, coupling=coupling, reference=reference)
        times = (0., .25, 1., 3., 11.)
        errors = [float(np.max(abs(correlation(arrays['complete_addition_frequencies'],
                    arrays['complete_addition_residues'], t)-correlation(aa['frequencies'],aa['residues'],t))))
                  for t in times]
        blocked_error = float(np.max(abs(blocked-aa['blocked_effect'])))
        if max([*errors,blocked_error]) > 2e-10:
            raise ArithmeticError('analytic complete source records differ')
        result['analytic'] = analytic
        reconstructed=reconstruct_coupling(aa['frequencies'],aa['blocked_effect'])
        expected=abs(coupling)*np.linalg.svd(y,compute_uv=False)[::-1]
        if not np.allclose(reconstructed['coupling_times_singular_values'],expected,rtol=2e-6,atol=1e-7):
            raise ArithmeticError('complete-record coupling inverse failed')
        result['complete_record_coupling_inverse']=reconstructed
        result['analytic_fixed_probe_correlation_checks'] = dict(zip(map(str,times), errors))
        result['analytic_blocked_effect_error'] = blocked_error
        arrays.update({'analytic_'+k:v for k,v in aa.items()})
    return result, arrays
