"""Charged poles of a specified interacting native-exchange architecture.

The open-chain law is H=J sum_b(I-S_b)+m sum_x(-1)^x(n_x-1/2).
J and each m are the existing native correlation coefficients. Alternating
the sign of the onsite phase is an explicit coherent control choice. It
uses the inverse onsite word, without introducing empirical mass values.
The density interaction in ordinary carrier exchange is retained.

Compatibility is the spectrum and CAR Lehmann law of this fixed action,
independently of a state-preparation or repair algorithm. In two cells the
global ground state has charge one. Both charged poles are
sqrt(J^2+m^2)-J, with nonzero local CAR residues. They differ from the bare
coefficient m; the squared hierarchy still fails the lepton comparison.

For an even open path of L cells, the charge-N unperturbed ground vector is
the uniform Dicke vector. H0/J is minus the symmetric-exclusion generator.
For w_x=(-1)^x, solve the path Poisson equation Lap f=w with sum f=0.
Then w.f=L/2 and H0(sum f_x n_x)|D_N>=J(sum w_x n_x)|D_N>.
The second-order energy is -m^2 N(L-N)/(2J(L-1)); reflection makes the
energy even in m. This coefficient shrinks with L. Its perturbative domain
is finite-volume and is not uniform at fixed m, so it proves no interchange
of the weak-field and thermodynamic limits or positive bulk particle mass.

Nearest reused theorem: Caputo, Liggett and Richthammer, arXiv:0906.1238,
Theorem1.1 and Section4.1.1 Eq4.1. The actual fixed-N exchange matrix is the
finite symmetric-exclusion Laplacian with edge rate J. Its gap is the path
random-walk gap, 2J(1-cos(pi/L)). Added here are its native parameter mapping,
charge-resolved response, explicit Poisson coefficient and bounded remainder.
The matrix is also used as a coherent Hamiltonian; Markov elapsed time is
not identified with its declared coherent evolution parameter.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
import sympy as s

from tools import strict_json
from . import repair_field_attachment as field
from ..calibrated_measurement.bounds import root_interval

ROOT = Path(__file__).resolve().parents[2]
SIZES = (2, 4, 6, 8, 10)
EIGENPAIR_ADMISSION = 1e-11


def _length(cells):
    if type(cells) is not int or cells not in SIZES:
        raise ValueError('an admitted even open-chain size is required')
    return cells


def _positive(value):
    if type(value) not in (int, F) or value <= 0:
        raise ValueError('an original positive exact rational is required')
    return F(value)


def native_parameters():
    p = field.logical_parameters()
    return F(p['exchange_strength']), tuple(F(x) for x in p['positive_occupation_gaps'])


def sector_matrix(cells, charge, mass, exchange):
    """Full adjacent-carrier exchange in a fixed occupation basis, exactly."""
    cells = _length(cells)
    if type(charge) is not int or not 0 <= charge <= cells:
        raise ValueError('a charge in the full finite CAR space is required')
    mass, exchange = _positive(mass), _positive(exchange)
    basis = tuple(mask for mask in range(1 << cells) if mask.bit_count() == charge)
    index = {mask: j for j, mask in enumerate(basis)}
    h = s.zeros(len(basis))
    for col, mask in enumerate(basis):
        h[col, col] += s.Rational(mass)*sum((-1)**x*((mask >> x)&1) for x in range(cells))
        for x in range(cells-1):
            if ((mask >> x) ^ (mask >> (x+1))) & 1:
                other = mask ^ (3 << x)
                h[col, col] += s.Rational(exchange)
                h[index[other], col] -= s.Rational(exchange)
    return basis, s.ImmutableMatrix(h)


def dimer_certificate(mass, exchange):
    mass, exchange = _positive(mass), _positive(exchange)
    # Normalize before the fixed-resolution rational square-root enclosure.
    # Its radicand lies in[1,2], so even arbitrarily small exact source
    # coefficients keep a positive lower root for the residue divisions.
    scale = max(mass, exchange)
    lo, hi = root_interval((exchange/scale)**2+(mass/scale)**2)
    lo, hi = scale*lo, scale*hi
    # Rationalize the difference so small charged gaps are never obtained
    # by subtracting two nearly equal floating eigenvalues.
    gap = (mass*mass/(hi+exchange), mass*mass/(lo+exchange))
    return {'mass_coefficient': str(mass), 'exchange': str(exchange),
        'charged_gap_exact': 'sqrt(J^2+m^2)-J',
        'charged_gap_interval': list(map(str, gap)),
        'charge_one_ground_energy_interval': [str(-gap[1]), str(-gap[0])],
        'charge_zero_and_two_energies': ['0', '0'],
        'charge_one_excited_energy_interval': [str(exchange+lo), str(exchange+hi)],
        'site_zero_removal_residue_interval': [str((1-mass/lo)/2), str((1-mass/hi)/2)],
        'site_zero_addition_residue_interval': [str((1+mass/hi)/2), str((1+mass/lo)/2)],
        'ground_occupation': 1,
        'both_charged_poles_are_positive': gap[0] > 0,
        'pole_conversion_is_a_relativistic_mass_identification': False}


def dimer_control(mass, exchange):
    """Independent full-Fock CAR spectral residues, including the interaction."""
    mass, exchange = _positive(mass), _positive(exchange)
    a = [np.asarray(field.annihilator(j, 2), complex) for j in (0, 1)]
    n = [op.conj().T@op for op in a]
    h = float(exchange)*(np.eye(4)-np.asarray(field.carrier_swap(0, 1, 2), complex))
    h += float(mass)*(n[0]-n[1])
    values, vectors = eigh(h)
    omega = vectors[:, 0]
    number = n[0]+n[1]
    if abs(np.vdot(omega, number@omega).real-1) > EIGENPAIR_ADMISSION:
        raise ArithmeticError('the global dimer ground state lost its charge-one domain')
    gap = float(F(dimer_certificate(mass, exchange)['charged_gap_interval'][0]))
    rows = []
    for site, op in enumerate(a):
        removed, added = op@omega, op.conj().T@omega
        for label, state in (('removal', removed), ('addition', added)):
            residue = float(np.vdot(state, state).real)
            residual = float(np.linalg.norm((h-values[0]*np.eye(4))@state-gap*state))
            if not residue > 0 or residual > EIGENPAIR_ADMISSION:
                raise ArithmeticError('a claimed charged pole lacks its full-Fock residue')
            rows.append({'site': site, 'channel': label, 'residue': residue})
    # A free second quantization of the one-particle block misses the
    # -2J n0*n1 term: the charge-two energy would be2J instead of0.
    free = h+2*float(exchange)*n[0]@n[1]
    return {'local_CAR_poles': rows,
        'free_one_particle_lift_charge_two_error_exact': str(2*exchange),
        'full_interaction_retained': True}, {'dimer_hamiltonian': h,
        'dimer_eigenvalues': values, 'dimer_ground_projector': np.outer(omega, omega.conj()),
        'dimer_free_lift': free}


def poisson_coefficient(cells, charge):
    cells = _length(cells)
    if type(charge) is not int or not 0 <= charge <= cells:
        raise ValueError('a finite-chain charge is required')
    lap = s.zeros(cells)
    for i in range(cells-1):
        lap[i, i] += 1; lap[i+1, i+1] += 1
        lap[i, i+1] -= 1; lap[i+1, i] -= 1
    w = s.Matrix([(-1)**i for i in range(cells)])
    f = (lap+s.ones(cells)/cells).inv()*w
    if lap*f != w or sum(f) != 0 or (w.T*f)[0] != s.Rational(cells, 2):
        raise ArithmeticError('the exact open-path Poisson coefficient failed')
    return s.Rational(charge*(cells-charge), 2*(cells-1))


def perturbative_bound(cells, charge, mass, exchange):
    """Explicit finite-volume even-series remainder; no uniform bulk claim.

    For the path,2(1-cos(pi/L))>=4/L^2. The perturbation norm is at most
    |m|*min(N,L-N). On the complex circle radius g/(4v), the isolated
    ground eigenvalue has modulus at most g/4. Reflection removes odd
    Taylor coefficients. Cauchy's estimate and a geometric sum bound all
    terms after second order. This is conservative and can be inconclusive.
    """
    cells = _length(cells)
    mass, exchange = _positive(mass), _positive(exchange)
    coefficient = F(poisson_coefficient(cells, charge))/exchange
    v = min(charge, cells-charge)
    if v == 0:
        return {'admitted': True, 'second_order_energy': '0', 'remainder_upper': '0',
                'isolation_ratio': '0', 'gap_lower': str(4*exchange/(cells*cells))}
    g = 4*exchange/(cells*cells)
    ratio = 4*mass*v/g
    admitted = ratio < 1
    remainder = g*ratio**4/(4*(1-ratio*ratio)) if admitted else None
    return {'admitted': admitted, 'second_order_energy': str(-mass*mass*coefficient),
        'remainder_upper': str(remainder) if admitted else None,
        'isolation_ratio': str(ratio), 'gap_lower': str(g)}


def finite_spectrum(cells, mass, exchange):
    """Numerical finite-volume diagnostics, separate from exact certificates."""
    cells = _length(cells)
    rows, arrays = [], {}
    for charge in range(cells+1):
        basis, exact = sector_matrix(cells, charge, mass, exchange)
        h = np.asarray(exact, float)
        values, vectors = eigh(h)
        residual = np.linalg.norm(h@vectors-vectors*values[None, :], 2)
        defect = np.linalg.norm(vectors.T@vectors-np.eye(len(basis)), 2)
        if residual > EIGENPAIR_ADMISSION or defect > EIGENPAIR_ADMISSION:
            raise ArithmeticError('finite charge spectrum exceeds numerical admission')
        bound = perturbative_bound(cells, charge, mass, exchange)
        if bound['admitted'] and abs(values[0]-float(F(bound['second_order_energy']))) > float(F(bound['remainder_upper']))+EIGENPAIR_ADMISSION:
            raise ArithmeticError('finite spectrum contradicts the analytic remainder')
        rows.append({'charge': charge, 'dimension': len(basis), 'ground_energy': float(values[0]),
                     'perturbative_certificate': bound})
        arrays[f'N{charge}_energies'] = values
    ground = min(row['ground_energy'] for row in rows)
    half = rows[cells//2]['ground_energy']
    if half-ground > EIGENPAIR_ADMISSION:
        raise ArithmeticError('half filling is not a global ground sector in this finite case')
    addition = rows[cells//2+1]['ground_energy']-half
    removal = rows[cells//2-1]['ground_energy']-half
    return {'cells': cells, 'charge_sectors': rows, 'half_filled_ground_energy': half,
        'addition_gap_diagnostic': addition, 'removal_gap_diagnostic': removal,
        'second_order_charged_gap_exact': str(mass*mass/(2*exchange*(cells-1))),
        'a_positive_thermodynamic_mass_is_certified': False}, arrays


def build():
    exchange, masses = native_parameters()
    certificates = [dimer_certificate(m, exchange) for m in masses]
    ordered = sorted((tuple(map(F, row['charged_gap_interval'])) for row in certificates), key=lambda x: x[0])
    ratios = [(lo/ordered[0][1], hi/ordered[0][0]) for lo, hi in ordered]
    arrays, dimers, finite = {}, [], []
    for sector, mass in enumerate(masses):
        result, values = dimer_control(mass, exchange)
        dimers.append(result)
        arrays.update({f's{sector}_{key}': value for key, value in values.items()})
        for size in SIZES:
            result, values = finite_spectrum(size, mass, exchange)
            finite.append({'native_sector': sector, **result})
            arrays.update({f's{sector}_L{size}_{key}': value for key, value in values.items()})
    prediction = {'schema': 'oph.native_charged_spectrum.v1',
        'source_law': 'J sum_open_bonds(I-S)+Delta_s sum_x(-1)^x(n_x-1/2)',
        'native_exchange_exact': str(exchange), 'native_mass_coefficients_exact': list(map(str, masses)),
        'empirical_inputs_to_prediction': [], 'dimer_certificates': certificates,
        'ascending_dimer_gap_ratio_intervals': [list(map(str, interval)) for interval in ratios],
        'full_Fock_dimer_controls': dimers, 'finite_charge_spectra': finite,
        'numerical_eigenpair_admission': EIGENPAIR_ADMISSION,
        'finite_volume_theorem': {'second_order_energy': '-m^2*N*(L-N)/(2J*(L-1))',
            'second_order_half_filled_charged_gap': 'm^2/(2J*(L-1))',
            'evenness': 'open-path reflection maps m to -m at fixed charge',
            'finite_volume_remainder': '(g/4)*r^4/(1-r^2), r=4|m|min(N,L-N)/g<1, g=4J/L^2',
            'limits_are_interchanged': False,
            'nearest_theorem': 'Caputo-Liggett-Richthammer, arXiv0906.1238, Theorem1.1 and Eq4.1',
            'primary_source': 'https://arxiv.org/abs/0906.1238'},
        'scope': {'onsite_sign_schedule_is_a_declared_control': True,
            'native_local_terms_inherited_from_executed_attachment': True,
            'new_complete_native_many_body_history_executed_here': False,
            'finite_matrix_eigensolvers_are_outward_certificates': False,
            'physical_mass_clock_dictionary_or_3D_interacting_limit_derived': False,
            'observed_particle_spectrum_derived': False}}
    # Empirical comparison is performed only after all action predictions.
    windows = ((F(206), F(207)), (F(3470), F(3485)))
    prediction['external_comparison'] = {
        'source': 'https://pdg.lbl.gov/2025/download/db2024.pdf',
        'inherited_broad_charged_lepton_ratio_windows': [list(map(str, pair)) for pair in windows],
        'ratio_overlaps': [not (hi < low or lo > high) for (lo, hi), (low, high) in zip(ratios[1:], windows, strict=True)],
        'windows_are_not_quoted_standard_uncertainties': True,
        'no_parameter_or_geometry_selected_by_these_windows': True}
    return prediction, arrays


def source_bindings():
    paths = set(field.source_bindings()) | {'model/common_source/repair_field_attachment.py',
        'model/common_source/native_charged_spectrum.py',
        'model/common_source/test_native_charged_spectrum.py',
        'model/calibrated_measurement/bounds.py'}
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def _compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual):
        raise ValueError(path+' original type differs')
    if type(expected) is dict:
        if set(expected) != set(actual):
            raise ValueError(path+' keys differ')
        for key in expected:
            _compare(expected[key], actual[key], path+'.'+key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError(path+' length differs')
        for i, (a, b) in enumerate(zip(expected, actual, strict=True)):
            _compare(a, b, path+f'[{i}]')
    elif type(expected) is float:
        exact = path.endswith('numerical_eigenpair_admission')
        if not math.isfinite(actual) or (expected != actual if exact else not math.isclose(expected, actual, rel_tol=5e-12, abs_tol=5e-12)):
            raise ValueError(path+' numerical comparison failed')
    elif expected != actual:
        raise ValueError(path+' differs')


def check(path):
    path = Path(path)
    report = strict_json.load(path)
    if type(report) is not dict or report.get('source_bindings') != source_bindings():
        raise ValueError('charged spectrum source bindings differ')
    artifact = path.with_suffix('.npz')
    if report.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('charged spectrum binary custody differs')
    fresh, outputs = build()
    _compare(fresh, {k: v for k, v in report.items() if k not in ('source_bindings', 'artifact')})
    with np.load(artifact, allow_pickle=False) as stored:
        if len(stored.files) != len(set(stored.files)) or set(stored.files) != set(outputs):
            raise ValueError('charged spectrum complete array keys differ')
        for key, value in outputs.items():
            old = stored[key]
            if old.shape != value.shape or old.dtype != value.dtype or not np.isfinite(old).all():
                raise ValueError('charged spectrum array schema differs: '+key)
            distance = np.linalg.norm(old-value, 2) if value.ndim == 2 else np.max(np.abs(old-value))
            if distance > 5e-12:
                raise ValueError('charged spectrum arrays differ: '+key)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check', type=Path)
    mode.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('Native interacting charged spectrum and finite-volume controls: PASS')
        return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists():
        parser.error('choose fresh receipt and artifact paths')
    report, arrays = build()
    report['source_bindings'] = source_bindings()
    with artifact.open('xb') as output:
        np.savez_compressed(output, **arrays)
    report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open('x') as output:
        json.dump(report, output, indent=2, sort_keys=True, allow_nan=False)
        output.write('\n')


if __name__ == '__main__':
    main()
