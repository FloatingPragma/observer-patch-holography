"""Reuse native current forms and the existing product field for energy tests.

The native quadratic form is derived from complete C/E records. Rest gap and
radial symbol curvature belong to the existing declared protected product
mode, not a localized defect or a physical momentum measurement. Finite pulse
subdivision has a controlled mass-shell discrepancy; no units are fitted.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
import sympy as sp

from . import clock_wave
from .joins_contract import binding, source_bindings, check_receipt
from ..observer import native, record_recovery as records
from ..observer.source_scalar import compiler
from ..observer_time_fields import continuum_dynamics as continuum
from ..observer_time_fields import spatial_dynamics as fields
from ..common_source import slow_geometry

HERE = Path(__file__).resolve().parent
Q = 630


def registered_incidence_map():
    """Map actual source labels and current orientations to the Lean instance."""
    checked = slow_geometry.registered_binding()
    permutation = checked['source_to_registered_lean_port']
    path = HERE.parents[1]/'Lean/Screen/SeamCurrentCarrierQuotient.lean'
    text = path.read_text()
    arrays = []
    for name in ('seamLeft', 'seamRight'):
        found = re.search(r'def '+name+r'[^=]*:=\s*!\[([^]]+)\]', text)
        if found is None:
            raise ValueError('registered Lean seam table missing')
        arrays.append(tuple(map(int, found.group(1).split(','))))
    edges = tuple(zip(*arrays, strict=True))
    if len(edges) != 30 or len(set(edges)) != 30:
        raise ValueError('registered seam table is not complete')
    maps = []
    for a, b in records.seams():
        pa, pb = permutation[a], permutation[b]
        edge = tuple(sorted((pa, pb)))
        sign = 1 if pa < pb else -1
        index = edges.index(edge)
        # One incidence column maps both charge and current, including sign.
        actual = [int(i == pb)-int(i == pa) for i in range(12)]
        target = [sign*(int(i == edge[1])-int(i == edge[0])) for i in range(12)]
        if actual != target:
            raise ArithmeticError('native-to-Lean current orientation failed')
        maps.append({'source_seam': [a, b], 'registered_seam_index': index,
                     'current_orientation': sign})
    return {**checked, 'source_seam_maps': maps,
            'charge_and_current_incidence_intertwine': True}


def encoded(value):
    if isinstance(value, sp.MatrixBase):
        return [[str(x) for x in row] for row in value.tolist()]
    if isinstance(value, (F, sp.Basic)):
        return str(value)
    if isinstance(value, dict):
        return {key: encoded(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encoded(item) for item in value]
    return value


def current_forms():
    """Stationary source covariance and zero-inclusive current quadratic form."""
    basis = native.states(binary=True)
    centered = [sp.Matrix([sp.Rational(x)-sp.Rational(1, 6) for x in state]) for state in basis]
    covariance = sum((x*x.T for x in centered), sp.zeros(12))/66
    laplacian = sp.Matrix(fields.native_density()['laplacian'])
    increments = {}
    for law in ('C', 'E'):
        total = sp.zeros(12)
        columns = 0
        for state in basis:
            for edge in records.seams():
                for coin in (0, 1):
                    output, record = records.recorded_step(law, state, edge, coin)
                    current = records.record_current(law, record)
                    delta = sp.Matrix([b-a for a, b in zip(state, output, strict=True)])
                    expected = sp.zeros(12, 1)
                    expected[edge[0]], expected[edge[1]] = -current, current
                    if delta != expected or sum(delta):
                        raise ArithmeticError('native current continuity failed')
                    total += delta*delta.T
                    columns += 1
        increments[law] = total/columns
    if (increments['C'] != increments['E'] or increments['C'] != laplacian/198
            or covariance != sp.Rational(5, 33)*(sp.eye(12)-sp.ones(12)/12)):
        raise ArithmeticError('native record form differs from source covariance or Laplacian')
    profile = sp.Matrix([1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    if laplacian*profile != 6*profile:
        raise ArithmeticError('the retained elementary source mode changed')
    variance = (profile.T*covariance*profile)[0]
    increment = (profile.T*increments['C']*profile)[0]
    if increment/(2*variance) != sp.Rational(1, 10):
        raise ArithmeticError('record-derived elementary frequency changed')
    return {'source_dimension': 66, 'complete_columns_per_law': columns,
            'stationary_reference': 'uniform binary66 under uniform seam/coin selection',
            'parameter': 'native attempted events, including zero current',
            'covariance': covariance, 'increment_second_moment': increments['C'],
            'C_and_E_forms_equal': True, 'current_continuity_all_columns': True,
            'registered_finite_continuity_map': registered_incidence_map(),
            'dirichlet_form': 'D(v)=E[(F_v(after)-F_v(before))^2]/2=v^T L v/396',
            'mode_profile': list(profile), 'mode_eigenvalue': 6,
            'mode_variance': variance, 'mode_increment_second_moment': increment,
            'normalized_dirichlet_rate': increment/(2*variance),
            'elementary_wave_phase_relation': '1-cos(theta)=D(v)/Var(F_v)=1/10',
            'bare_excess_movement_redefined': False,
            'physical_energy_identified': False}


def native_mode():
    """Verify the existing full-port fixture; do not identify distinct carriers."""
    witness = continuum.product_native_mode_control()
    if (witness['ports'] != 240 or witness['seams'] != Q
            or witness['native_eigenvalue'] != 6):
        raise ArithmeticError('the existing full-port product mode changed')
    return witness


def ratio_bounds(subdivisions):
    """Enclose chi*cot(chi), chi=acos(209/210)/(2*n), n a power of two.

    The half-angle recurrence is exact: cot(x/2)=cot(x)+sqrt(cot(x)^2+1).
    The defect bound uses sin(x)-x*cos(x)=integral_0^x t*sin(t)dt.
    """
    if type(subdivisions) is not int or subdivisions not in (1, 2, 4, 8):
        raise ValueError('retained pulse subdivisions are 1,2,4,8')
    theta = clock_wave._wave_frequency_bounds((F(2, 7), F(2, 7)))
    cot = clock_wave._sqrt_bounds(F(419))
    count = 1
    while count < subdivisions:
        cot = (cot[0]+clock_wave._sqrt_bounds(cot[0]**2+1)[0],
               cot[1]+clock_wave._sqrt_bounds(cot[1]**2+1)[1])
        count *= 2
    chi = theta[0]/(2*subdivisions), theta[1]/(2*subdivisions)
    ratio = chi[0]*cot[0], chi[1]*cot[1]
    upper_defect = chi[1]**2/(3*(1-chi[1]**2/6))
    if not 0 < 1-ratio[1] <= 1-ratio[0] <= upper_defect:
        raise ArithmeticError('finite split inertia bound failed')
    return {'subdivisions': subdivisions,
            'rest_to_independent_inertia_ratio': clock_wave._outward(ratio, 14),
            'mass_shell_defect_upper': clock_wave._outward((upper_defect, upper_defect), 14)[1],
            'pulses_per_original_trigger': 2*subdivisions,
            'source_load_export_calls_per_trigger': 4*subdivisions}


def source_word(radial_frequency, subdivisions=1):
    """Existing product pulses, optionally subdivided at fixed elapsed parameter.

    Radial frequency is a declared symbol parameter, not an inferred native
    momentum. Both original primitive requests, including a zero pulse, remain.
    """
    ratio_bounds(subdivisions)
    p = clock_wave._fraction(radial_frequency)
    if not 0 <= p <= 1:
        raise ValueError('retained radial symbol parameter lies in [0,1]')
    theta = sp.acos(sp.Rational(209, 210))
    angles = (theta/(2*subdivisions), sp.Rational(p)/(2*subdivisions*sp.sqrt(Q)))
    word = [compiler.pulse(axis, float(continuum._source_dyadic_scalar(angle)), [0])
            for axis, angle in zip(('Z', 'X'), angles, strict=True)]
    return word*subdivisions


def pulse_checks():
    """Original source-array execution, including an untouched Bell reference."""
    rows = []
    x = np.array([[0, 1], [1, 0]], complex)
    z = np.diag([1., -1.])
    theta = float(sp.N(sp.acos(sp.Rational(209, 210)), 60))
    a = np.sqrt(Q)*theta
    bell = np.array([1., 0., 0., 1.], complex)/np.sqrt(2)
    for n in (1, 2, 4, 8):
        for p in (F(0), F(2, 5), F(4, 5)):
            word = source_word(p, n)
            actual = compiler.execute(word, 1)
            alpha, beta = theta/(2*n), float(p)/(2*n*np.sqrt(Q))
            # Independent closed Pauli multiplication, no compiler helpers.
            elementary = ((np.cos(beta)*np.eye(2)+1j*np.sin(beta)*x)
                          @(np.cos(alpha)*np.eye(2)+1j*np.sin(alpha)*z))
            exact = np.linalg.matrix_power(elementary, n)
            residual = float(np.linalg.norm(actual-exact, 2))
            omega = np.sqrt(a*a+float(p)**2)
            normal = (a*z+float(p)*x)/omega
            ideal = (np.cos(omega/(2*np.sqrt(Q)))*np.eye(2)
                     +1j*np.sin(omega/(2*np.sqrt(Q)))*normal)
            error = float(np.linalg.norm(actual-ideal, 2))
            split_bound = a*float(p)/(4*Q*n)
            numeric_budget = 1e-10
            if residual > numeric_budget or error > split_bound+numeric_budget:
                raise ArithmeticError('source product word exceeded matrix or split bound')
            output = compiler.execute(word, 2, bell)
            expected = np.kron(np.eye(2), actual)@bell
            if np.linalg.norm(output-expected) > numeric_budget:
                raise ArithmeticError('product word changed the untouched reference')
            probabilities = abs(output)**2
            if abs(sum(probabilities)-1) > numeric_budget:
                raise ArithmeticError('source output normalization changed')
            # E is the gap between the two eigenphases, not either half-angle.
            gap = 2*np.sqrt(Q)*np.arccos(np.clip(np.trace(actual).real/2, -1, 1))
            def diagnostic_interval(value):
                # Coarse outward display of a floating diagnostic, not an
                # interval-arithmetic certificate. Independent analytic bounds
                # are recorded separately. The grid also avoids BLAS last bits.
                grid = 1 << 26
                center = round(value*grid)
                # Exact dyadic probabilities (notably1/2) must not cross a
                # display-bin boundary when BLAS changes the last few bits.
                return [str(F(max(0, center-1), grid)), str(F(center+1, grid))]
            rows.append({'subdivisions': n, 'radial_symbol_parameter': str(p),
                         'word_sha256': hashlib.sha256(json.dumps(word, sort_keys=True).encode()).hexdigest(),
                         'source_pulses': len(word),
                         'matrix_residual_diagnostic_interval': diagnostic_interval(residual),
                         'target_generator_unitary_error_diagnostic_interval': diagnostic_interval(error),
                         'analytic_split_bound_formula': 'a*p/(4*Q*n)',
                         'split_bound_display_interval': diagnostic_interval(split_bound),
                         'floating_execution_budget': '1/10000000000',
                         'spectral_gap_diagnostic_interval': diagnostic_interval(float(gap)),
                         'complete_system_reference_probability_intervals': [diagnostic_interval(p) for p in probabilities]})
    return {'rows': rows,
            'source_compiler_proper_code_and_failure_contract_reused': True,
            'floating_display_intervals_are_rigorous_interval_arithmetic': False,
            'ideal_proper_code_failure_probability': 0,
            'new_hardware_fidelity_certificate': False}


def symbolic_response():
    """Differentiate the executed split symbol before extracting inertia."""
    p = sp.Symbol('p', real=True)
    a, h = sp.symbols('a h', positive=True)
    # Its local branch is valid for 0<a*h/2<pi/2. Differentiate before p=0.
    gap = 2*sp.acos(sp.cos(a*h/2)*sp.cos(p*h/2))/h
    curvature = sp.diff(gap, p, 2).subs(p, 0)
    expected = h*sp.cos(a*h/2)/(2*sp.sqrt(1-sp.cos(a*h/2)**2))
    if sp.simplify(curvature-expected) != 0:
        raise ArithmeticError('independent split-symbol differentiation failed')
    return {'executed_gap': 'E_n(p)=2*n*sqrt(Q)*acos(cos(theta/(2*n))*cos(p/(2*n*sqrt(Q))))',
            'rest_gap': 'E_n(0)=sqrt(Q)*theta=a',
            'independent_radial_curvature': 'E_n_second(0)=cot(theta/(2*n))/(2*n*sqrt(Q))',
            'kinetic_mass_from_curvature': 'm_n=1/E_n_second(0)=2*n*sqrt(Q)*tan(theta/(2*n))',
            'massless_branch_signal_speed': '1 in the declared radial symbol units, below phase aliasing',
            'rest_to_inertia_ratio': '(theta/(2*n))*cot(theta/(2*n))',
            'ideal_product_gap': 'sqrt(a^2+p^2)',
            'ideal_independent_curvature': '1/a',
            'ideal_mass_shell_relation': 'E0=m_kin*c^2 with c=1',
            'radial_symbol_is_native_physical_momentum': False,
            'pulse_subdivision_is_carrier_refinement': False,
            'gravitational_stress_constructed': False}


def theorem_maps():
    return [
        {'path': 'Lean/Geometry/SourceActionMeasure.lean',
         'theorem': 'normalized_symmetrizer_unique',
         'model_map': 'The connected native twelve-port Laplacian is symmetric with negative reciprocal edge support; uniform positive capacities normalized to1 are its unique diagonal symmetrizer.',
         'verified_inputs': 'Exact complete-record covariance, current form and native Laplacian; rank11 and edge-square factorization checked independently.',
         'remaining': 'A relative quadratic action does not select physical kinetic normalization, volume or an action unit.'},
        {'path': 'Lean/Screen/RegionalContinuity.lean',
         'theorem': 'regional_continuity',
         'model_map': 'q=x, qNext=native_output and source=0, transported by the checked sourceToRERPort permutation; C/E current is mapped to the registered seam with its required orientation sign.',
         'verified_inputs': 'Pointwise continuity on every supported binary source/seam/coin column implies every regional sum.',
         'remaining': 'Conserved occupation is not identified with stress-energy. WardLimitManifest still requires a stress realization and controlled refinement.'},
        {'path': 'Lean/Geometry/MassShellKinematics.lean',
         'theorem': 'fourMomentum_time_sq_eq_mass_sq_add_spatial',
         'model_map': 'For the ideal existing product symbol, E=sqrt(a^2+p^2), frame=(E/a,p/a,0,0), declared rest parametera>0 and c=1.',
         'verified_inputs': 'Actual source eigenmode fixes theta anda; anticommuting product axes give the ideal symbol. Independent curvature1/a agrees. Finite split words have the retained controlled deviation.',
         'remaining': 'Radial apparatus and momentum dictionary are declared. No localized defect or physical Lorentz frame is inferred.'},
        {'path': 'Lean/Geometry/InternalEnergyInertia.lean',
         'theorem': 'compositeStationary_iff_eom',
         'model_map': 'The theorem locates the action-shape premise behind equal rest and inertial energy; it is not imported from bare record activity.',
         'verified_inputs': 'The product gap and its radial-symbol curvature are computed independently, so their equality or finite splitting mismatch is tested rather than defined.',
         'remaining': 'The theorem worldline/composite action is not yet derived for the native history defect.'},
        {'path': 'Lean/Geometry/ProperTimeInternalAction.lean',
         'theorem': 'slope_selection',
         'model_map': 'Reuse its refinement-invariance and proper-time criteria as the action-selection contract for any defect action.',
         'verified_inputs': 'The finite product response tends to its own ideal mass-shell relation as pulse subdivision error vanishes.',
         'remaining': 'Pulse subdivision is not the theorem midpoint worldline refinement and does not discharge its proper-time premise.'},
        {'path': 'Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean',
         'theorem': 'modularEnergy_gibbs_sub',
         'model_map': 'The existing product-mode faithful Gibbs target uses energy levels0,omega and common beta; the modular ledger difference equals beta times energy difference.',
         'verified_inputs': 'The existing source product preparation/readout test is reused; the finite pulse word is separately compared to the ideal generator.',
         'remaining': 'An arbitrary transported three-charge target or finite prepared state cannot silently replace this Gibbs target.'},
    ]


def report():
    paths = set(continuum.SOURCE_PATHS) | {
        'model/clock_defect_energy/joins_contract.py', 'model/clock_defect_energy/contract.py',
        'model/clock_defect_energy/mode_energy.py', 'model/clock_defect_energy/test_mode_energy.py',
        'model/clock_defect_energy/MODE_ENERGY.md',
        'model/clock_defect_energy/clock_wave.py', 'model/observer/record_recovery.py',
        'Lean/Geometry/ProperTimeInternalAction.lean', 'Lean/Geometry/InternalEnergyInertia.lean',
        'Lean/Geometry/MassShellKinematics.lean', 'Lean/Geometry/SourceActionMeasure.lean',
        'Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean',
        'Lean/Screen/RegionalContinuity.lean', 'Lean/Screen/SeamChargeContinuity.lean',
        'Lean/Screen/SeamCurrentCarrierQuotient.lean',
        'Lean/Screen/PrimitivePortFrameQuotient.lean', 'Lean/Screen/PortFrameGram.lean',
        'model/common_source/slow_geometry.py',
        'Lean/Dynamics/WardLimitManifest.lean', 'tools/strict_json.py'}
    return encoded({'schema': 'oph.clock_defect_energy.mode_energy.v1', 'campaign': binding(),
                    'source_sha256': source_bindings(paths), 'native_record_forms': current_forms(),
                    'existing_full_source_mode': native_mode(), 'response': symbolic_response(),
                    'theorem_maps': theorem_maps(),
                    'finite_inertia_bounds': [ratio_bounds(n) for n in (1, 2, 4, 8)],
                    'source_execution': pulse_checks(),
                    'scope': {'physical_mass_derived': False, 'localized_defect_identified': False,
                              'physical_comparison_count_added': 0, 'integrated_gate_accepted': False}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = report()
    output = HERE/'mode_energy_results.json'
    if args.check:
        check_receipt(output, value)
    if args.write:
        output.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    print(json.dumps(value['finite_inertia_bounds'], indent=2))


if __name__ == '__main__':
    main()
