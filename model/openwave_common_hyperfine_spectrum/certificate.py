"""Complete exact spectra and common native-current calibration records."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb
import json
import numpy as np
import scipy.sparse as sp
import sympy as sy
from scipy.linalg import eigh
from model.openwave_six_quark import source as native6
from model.openwave_colour_string.fixed_witness import spectral_enclosure, integers
from model.openwave_axial_state.certificate import exact_action, sqrt_upper, reduced_line
from model.openwave_spin_gradient_nuclei.certificate import reduce_exact
from model.openwave_spin_gradient_nuclei.fixed_nuclear_basis import channel as fixed_nuclear_channel
from model.openwave_calibrated_hadrons.source import angular4
from model.openwave_nuclear_distance import nuclear
from model.openwave_string_positive_cone.source import sqrt_interval
from . import source as s

UNIT = F(s.INPUTS['electron_unit_GeV'])
WITNESSES = json.loads((s.HERE / 'spectral_witnesses.json').read_text())


@lru_cache(None)
def baryon_exact(S2, I2):
    d = s.components('baryon')
    eye = sp.eye(len(d['basis']), dtype=np.int64, format='csr')
    S4 = angular4(d['modes'], d['basis'], (0, 1), True)
    I4 = angular4(d['modes'], d['basis'], (0, 1), False)
    ids = [i for i, row in enumerate(d['labels'])
           if sum(f == 0 for x, z, f in row) == 1 and sum(z == 0 for x, z, f in row) == 2]
    a = 15 * eye - S4 if S2 == 1 else S4 - 3 * eye
    b = 15 * eye - I4 if I2 == 1 else I4 - 3 * eye
    B = native6.exact_product(a, native6.exact_product(b, d['Binteger'][:, ids]))
    G = sy.Matrix(native6.exact_product(B.T, B).toarray().tolist())
    _, cols = G.rref()
    B = B[:, list(cols)].tocsr()
    G = sy.Matrix(native6.exact_product(B.T, B).toarray().tolist())
    assert (S4 @ B - S2 * (S2 + 2) * B).nnz == 0
    assert (I4 @ B - I2 * (I2 + 2) * B).nnz == 0
    H, T = reduce_exact(B, G, d['parts'], s.parameters())
    return H, T, B


@lru_cache(None)
def baryon_spectra():
    rows = []
    for S, I in product((1, 3), repeat=2):
        H, _, _ = baryon_exact(S, I)
        spec = spectral_enclosure(H, WITNESSES[f'baryon_{S}_{I}'])
        center = [[str(F(a) - s.parameters()['delta'] / 2), str(F(b) - s.parameters()['delta'] / 2)]
                  for a, b in spec['all_eigenvalue_intervals']]
        rows.append({'twice_spin': S, 'twice_isospin': I, 'fixed_Ms_half_Mi_minus_half_dimension': H.rows,
                     'full_magnetic_dimension': H.rows * (S + 1) * (I + 1),
                     'complete_spectrum': spec, 'center_intervals': center})
    assert sum(x['full_magnetic_dimension'] for x in rows) == 120
    return rows


def threshold(S, I):
    choices = []
    for a, b in product(baryon_spectra(), repeat=2):
        sa, sb = F(a['twice_spin'], 2), F(b['twice_spin'], 2)
        ia, ib = F(a['twice_isospin'], 2), F(b['twice_isospin'], 2)
        if abs(sa - sb) <= S <= sa + sb and abs(ia - ib) <= I <= ia + ib:
            x, y = a['center_intervals'][0], b['center_intervals'][0]
            choices.append((F(x[0]) + F(y[0]), F(x[1]) + F(y[1])))
    return min(x for x, y in choices), min(y for x, y in choices)


@lru_cache(None)
def nuclear_exact(S, I, parity):
    d, c = nuclear.carrier(), fixed_nuclear_channel(S, I, parity)
    if not c['rank']:
        return None
    parts = {**d['parts'], 'G4': 2 * d['parts']['S4'] - 4 * d['Sg']}
    H, T = reduce_exact(c['Binteger'], c['Gram'], parts, s.parameters())
    return H, T, c['Binteger'], parts, d


@lru_cache(None)
def nuclear_spectra():
    rows = []
    for S, I, parity in product(range(4), range(4), (-1, 1)):
        data = nuclear_exact(S, I, parity)
        if data is None:
            continue
        spec = spectral_enclosure(data[0], WITNESSES[f'nuclear_{S}_{I}_{parity}'])
        lo, hi = map(F, spec['all_eigenvalue_intervals'][0])
        tl, th = threshold(S, I)
        rows.append({'spin': S, 'isospin': I, 'cell_parity': parity,
                     'fixed_Ms0_Mi0_dimension': data[0].rows,
                     'full_magnetic_dimension': data[0].rows * (2 * S + 1) * (2 * I + 1),
                     'complete_spectrum': spec, 'ground_GeV_interval': list(map(str, (lo * UNIT, hi * UNIT))),
                     'binding_GeV_interval': list(map(str, ((tl - hi) * UNIT, (th - lo) * UNIT))),
                     'strictly_bound': tl > hi, 'strictly_above_compatible_breakup': lo > th,
                     'all_Iz_shifts': [str(-s.parameters()['delta'] * m) for m in range(-I, I + 1)]})
    assert sum(x['fixed_Ms0_Mi0_dimension'] for x in rows) == 240
    assert sum(x['full_magnetic_dimension'] for x in rows) == 2520
    return rows


@lru_cache(None)
def axial_certificate():
    old, d = reduced_line(), s.components('baryon')
    H, T = reduce_exact(old['B'], old['Gram'], d['parts'], s.parameters())
    v = exact_action(old['B'], integers(WITNESSES['axial']['trial_basis_integer_coefficients'], old['B'].shape[1], 2**45))
    norm = sum(x * x for x in v)
    hv = [F(0) for _ in v]
    for k, coefficient in s.weights(s.parameters()).items():
        for i, x in enumerate(exact_action(d['parts'][k], v)):
            hv[i] += coefficient * x
    mean = sum(x * y for x, y in zip(v, hv)) / norm
    variance = sum((x - mean * y)**2 for x, y in zip(hv, v)) / norm
    spec = spectral_enclosure(H, WITNESSES['axial'])
    gap = F(spec['all_eigenvalue_intervals'][1][0]) - mean
    assert gap > 0
    sine = sqrt_upper(variance) / gap
    gA = F(sum(x * y for x, y in zip(v, exact_action(old['X'], v))), norm)
    error = 12 * sine
    assert abs(gA - F(s.INPUTS['axial_calibration_target'])) + error < F(1, 10**6)
    # Compare actual certified vectors, not only the same axial label.
    parent = __import__('json').loads((s.HERE.parent / 'openwave_spin_gradient_nuclei/results.json').read_text())['native_axial_ground']
    w = list(map(int, parent['integer_native_trial']))
    wn = sum(x * x for x in w)
    overlap2 = F(sum(x * y for x, y in zip(v, w))**2, norm * wn)
    distance = sqrt_upper(1 - overlap2) + sine + F(parent['ground_vector_sine_bound'])
    return {'gV_exact': '1', 'gA_interval': [str(gA - error), str(gA + error)],
            'gA_interval_decimal': [float(gA - error), float(gA + error)],
            'integer_native_trial': list(map(str, v)), 'trial_norm': str(norm), 'trial_mean': str(mean), 'trial_variance': str(variance),
            'ground_vector_sine_bound': str(sine), 'complete_selected_spectrum': spec,
            'five_input_ground_projector_norm_distance_upper': str(distance),
            'arbitrary_reference_complete_instrument_TV_transfer_upper': str(2 * distance),
            'transfer_scope': 'The same full native CAR carrier, currents and reference frame; only the stated ground preparation is compared. Other strong eigenstates and all nuclear states are recalculated.'}


@lru_cache(None)
def nuclear_trial():
    H, T, B, parts, d = nuclear_exact(1, 0, 1)
    v = exact_action(B, integers(WITNESSES['nuclear_1_0_1']['trial_basis_integer_coefficients'], B.shape[1], 2**45))
    norm = sum(x * x for x in v)
    hv = [F(0) for _ in v]
    for k, c in s.native.weights(s.parameters()).items():
        for i, x in enumerate(exact_action(parts[k], v)):
            hv[i] += c * x
    mean = sum(x * y for x, y in zip(v, hv)) / norm
    variance = sum((x - mean * y)**2 for x, y in zip(hv, v)) / norm
    spec = next(x for x in nuclear_spectra() if (x['spin'], x['isospin'], x['cell_parity']) == (1, 0, 1))
    gap = F(spec['complete_spectrum']['all_eigenvalue_intervals'][1][0]) - mean
    sine = sqrt_upper(variance) / gap
    Q = F(sum(x * x * int(n * (6 - n)) for x, n in zip(v, d['n0'])), 18 * norm)
    return {'quadrupole_graph_interval': [str(max(F(0), Q - sine)), str(Q + sine)],
            'trial_energy_GeV': str(mean * UNIT), 'ground_vector_sine_bound': str(sine),
            'native_integer_trial': list(map(str, v)), 'trial_norm': str(norm)}


@lru_cache(None)
def full_fock_local_bound(kappa):
    """All169 local quark/antiquark counts, with no Fock-state truncation."""
    fundamental = [(1, 0), (-1, 1), (0, -1)] * 4
    character = {(0, 0, 0, 0): 1}
    for species, weights in enumerate((fundamental, [(-a, -b) for a, b in fundamental])):
        for u, v in weights:
            out = dict(character)
            for (nq, na, a, b), m in character.items():
                key = (nq + int(species == 0), na + int(species == 1), a + u, b + v)
                out[key] = out.get(key, 0) + m
            character = out
    roots = ((2, -1), (-1, 2), (1, 1))
    denominator = [(sum(x for x, y in z), sum(y for x, y in z), (-1)**n)
                   for n in range(4) for z in combinations(roots, n)]
    p, rows, dimension = s.parameters(), [], 0
    c = max(F(0), F(kappa) / 4 - p['U'])
    for nq, na in product(range(13), repeat=2):
        irreps, count, maximum = [], 0, F(0)
        for a, b in product(range(9), repeat=2):
            m = sum(sign * character.get((nq, na, a + x, b + y), 0) for x, y, sign in denominator)
            assert m >= 0
            if m:
                casimir = F(a*a + b*b + a*b + 3*a + 3*b, 3)
                dim = (a+1)*(b+1)*(a+b+2)//2
                irreps.append([a, b, m])
                count += m * dim
                maximum = max(maximum, casimir)
        assert count == comb(12, nq) * comb(12, na)
        dimension += count
        pair_floor = nq * na - min(nq, na)
        value = ((p['mbar'] - abs(p['delta']) / 2) * (nq + na)
                 + p['A'] * F((nq - na)**2, 9) + p['pair_triplet'] * pair_floor - c * maximum)
        assert value > 0 if nq + na else value == 0
        rows.append({'n_quark': nq, 'n_antiquark': na, 'dimension': count,
                     'complete_irreps': irreps, 'Casimir_maximum': str(maximum),
                     'pair_triplet_lower': pair_floor, 'local_energy_lower': str(value)})
    assert dimension == 2**24
    floor = min(F(x['local_energy_lower']) / (x['n_quark'] + x['n_antiquark'])
                for x in rows if x['n_quark'] + x['n_antiquark'])
    return {'full_local_Fock_dimension': dimension, 'complete_count_sectors': rows,
            'uniform_number_floor': str(floor), 'kappa': str(kappa),
            'all_chain_lengths_and_all_local_occupations': True,
            'scope': 'One-generation zero-boundary-flux Gauss carrier, with global positive pair-triplet interaction; every positive cross-cell pair is retained.'}


@lru_cache(None)
def shared_length_attachment():
    qlo, qhi = map(F, nuclear_trial()['quadrupole_graph_interval'])
    target = F('.285699')
    alo, ahi = sqrt_interval(target / qhi)[0], sqrt_interval(target / qlo)[1]
    klo, khi = alo * 3 / (4 * UNIT), ahi * 3 / (4 * UNIT)
    # Select a rational point, carrying the interval instead of claiming exact empirical data.
    k = (klo + khi) / 2
    a = 4 * k * UNIT / 3
    return {'kappa': str(k), 'kappa_required_interval': [str(klo), str(khi)],
            'a_fm': str(a), 'a_fm_diagnostic': float(a),
            'quadrupole_fm2_interval': [str(a*a*qlo), str(a*a*qhi)],
            'illustrative_tension_GeV_per_fm': '1',
            'full_Fock_positivity': full_fock_local_bound(khi),
            'ancestry': 'Quadrupole comparison0.285699fm^2 and illustrative tension1GeV/fm select this conditional length and kappa on the seven-input source; neither is an independent prediction.'}


def complete_records():
    rows = []
    for kind, singlet, full in (('baryon', 120, 2024), ('meson', 64, 576), ('nuclear', 2520, comb(24, 6))):
        rows.append({'carrier': kind, 'original_occupation_dimension': full, 'complete_colour_singlet_dimension': singlet,
                     'singlet_preparation_success': str(F(singlet, full)), 'coloured_preparation_failure': str(1-F(singlet, full)),
                     'each_resolved_rank_one_spectral_record_probability': str(F(1, full)),
                     'conditional_probability_per_rank': str(F(1, singlet)), 'outside_number_sector_probability': '0',
                     'outside_number_sector_effect': 'I minus the declared full native number-sector projector',
                     'all_magnetic_and_multiplicity_records_retained': True,
                     'reference': 'Untouched dual native occupation reference of the maximally entangled input; every outcome retains it.',
                     'nonselective_mean_work_exact': '0 on the full uniform native occupation input',
                     'instrument': 'Actual joint colour, number, commuting spin/isospin/charge and complete Hamiltonian spectral projectors; degenerate multiplicity spaces retained without imposing a species label.'})
    return rows


def calculate():
    p = s.parameters()
    assert min(p.values()) > 0 and p['mbar'] > p['delta']/2
    baryons, nuclei, axial = baryon_spectra(), nuclear_spectra(), axial_certificate()
    H, d = s.exact_matrix('meson'), s.components('meson')
    mesons = []
    for charge in sorted(set(np.rint(np.diag(d['q']).real).astype(int))):
        ids = np.flatnonzero(np.rint(np.diag(d['q']).real).astype(int) == charge).tolist()
        mesons.append({'charge': int(charge), 'complete_spectrum': spectral_enclosure(H.extract(ids, ids), WITNESSES[f'meson_{charge}'])})
    from model.openwave_hadron_spin_constraints.constraints import block
    fitted = []
    for name, kind, charge, spin in (('proton','baryon',1,1), ('neutron','baryon',0,1), ('charged_pion','meson',1,0), ('neutral_vector','meson',0,2)):
        P = block(kind, charge, spin)
        value = float(eigh(P.conj().T @ s.matrix(kind) @ P, eigvals_only=True)[0]) * float(UNIT)
        target = float(s.INPUTS['mass_calibration_targets_GeV'][name])
        assert abs(value-target) < 2e-10
        fitted.append({'candidate': name, 'energy_GeV': value, 'target_GeV': target, 'role': 'empirical calibration input'})
    delta_spec = next(x for x in baryons if (x['twice_spin'],x['twice_isospin']) == (3,3))
    delta_interval = [(F(x)-p['delta'])*UNIT for x in delta_spec['complete_spectrum']['all_eigenvalue_intervals'][0]]
    target = F(s.INPUTS['mass_calibration_targets_GeV']['Delta_plus'])
    assert delta_interval[0] <= target <= delta_interval[1]
    fitted.append({'candidate':'Delta_plus', 'energy_GeV_interval': list(map(str,delta_interval)), 'target_GeV': str(target), 'role':'empirical calibration input; native S=I=3/2,Q=1'})
    even = next(x for x in nuclei if (x['spin'],x['isospin'],x['cell_parity']) == (1,0,1))
    lo,hi = map(F,even['binding_GeV_interval'])
    assert lo <= F(s.INPUTS['deuteron_binding_calibration_target_GeV']) <= hi
    for S,I in ((0,0),(1,1)):
        assert next(x for x in nuclei if (x['spin'],x['isospin'],x['cell_parity']) == (S,I,-1))['strictly_above_compatible_breakup']
    neutral = block('meson',0,0)
    neutral_mass = float(eigh(neutral.conj().T@s.matrix('meson')@neutral,eigvals_only=True)[0])*float(UNIT)
    attachment = shared_length_attachment()
    from .gauss import fusion_certificate
    fusion = fusion_certificate()
    fusion['selected_local_colour_lower_coefficient'] = str(p['U'] - F(2,11)*F(attachment['kappa_required_interval'][1]))
    fusion['strict_positive_coefficient_shortcut'] = F(fusion['selected_local_colour_lower_coefficient']) >= 0
    return {'schema':'oph.common_native_hyperfine_spectrum.v1',
            'verdict':'PASS: seven calibrated native-current outputs and complete retained finite spectra',
            'parameters_rational':{k:str(v) for k,v in p.items()}, 'input_ancestry':s.INPUTS,
            'native_pair_identity':{'full_baryon_occupation_dimension':len(s.components('baryon')['basis']), 'full_meson_occupation_dimension':len(s.components('meson')['basis']), 'Pair4_zero_on_all_pure_quark_sectors':True, 'Pair4_equals_total_spin4_on_one_pair':True, 'positive_full_Fock_pair_sum':True},
            'five_mass_calibrations':fitted, 'native_axial_ground':axial,
            'complete_baryon_multiplets':baryons, 'complete_meson_charge_spectra':mesons,
            'complete_nuclear_channels':nuclei, 'deuteron_binding_calibration':even['binding_GeV_interval'],
            'all_unassigned_bound_nuclear_channels':[x for x in nuclei if x['strictly_bound'] and (x['spin'],x['isospin'],x['cell_parity']) != (1,0,1)],
            'nuclear_ground_record':nuclear_trial(), 'complete_preparation_and_spectral_records':complete_records(),
            'common_length_attachment':attachment, 'sharp_Gauss_fusion_bound':fusion,
            'held_out_neutral_pion_control':{'candidate_GeV':neutral_mass, 'comparison_GeV':s.INPUTS['held_out_neutral_pion_GeV'], 'difference_GeV':neutral_mass-float(s.INPUTS['held_out_neutral_pion_GeV']), 'matches':abs(neutral_mass-float(s.INPUTS['held_out_neutral_pion_GeV']))<.00001},
            'boundaries':['Seven target agreements are calibrated, not held-out predictions.','Finite Hermitian energies calibrated to Breit-Wigner conventions do not supply rho/Delta resonance poles or widths.','Pair-triplet action includes all q-antiquark pairs across the two-cell source; no four-dimensional local field completion is inferred.','The common physical length is reattached on this source, with every local Fock occupation retained in its positivity proof.']}
