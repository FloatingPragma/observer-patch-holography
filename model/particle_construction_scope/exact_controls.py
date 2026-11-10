"""Small exact controls supplementing, without altering, the Fable producers."""
import ast
from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as s


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def controls(rer, work):
    sys.path.insert(0, str(rer))
    from model.coherent_seam_repair.program import Observer, State
    from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
    from model.native_weak_current_interface.model import scalar_controls

    # Execute the actual optional controller on a vector supported on every
    # basis coordinate. Its diagonal action verifies all coordinates at once.
    observer = Observer()
    state = State((1,) * observer.dimension, observer.dimension)
    backward = observer.addresses.index((3, 2))
    targets = [observer.index(port, 1, backward) for port in range(12)]
    record = []
    for target in targets:
        for _ in range(2):
            state, event = observer.quarter_turn(state, target)
            record.append(event)
    wanted = tuple(-1 if i in targets else 1 for i in range(observer.dimension))
    require(state.amplitudes == wanted, 'literal native V controller differs')
    require(state.norm_squared() == 1, 'native V lost norm')

    # General second-order kept-clock self-energy: E_f=E_f0+d and
    # d=w/(d-omega) for each flag. All three residues are equal.
    w, omega = s.symbols('w omega', positive=True)
    shift = (omega-s.sqrt(omega**2+4*w))/2
    require(s.simplify(shift*(shift-omega)-w) == 0, 'kept-clock equation')
    residue = s.simplify(1/(1+w/(shift-omega)**2))

    # HT's actual exact register checker, extracted by AST without running its
    # separate calibrated D11 flow. The function body is executed unchanged.
    ht = work/'HT/part1_c1.py.raw'
    parsed = ast.parse(ht.read_text())
    function = next(node for node in parsed.body if isinstance(node, ast.FunctionDef)
                    and node.name == 'common_eigvecs')
    namespace = {'sp': s}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(ht), 'exec'), namespace)
    check = namespace['common_eigvecs']
    memory = s.Matrix([[0,1],[0,0]])
    cycle = s.Matrix([[0,0,1],[1,0,0],[0,1,0]])
    require(check(memory) == [], 'finite lowering common eigenstate changed')
    require(set(check(cycle)) == {s.Integer(1)}, 'finite cyclic modulus changed')
    # No finite exact CCR; a finite oscillator nevertheless has exact CCR on
    # its retained interior. The edge defect retains the zero trace.
    ccr = []
    for dimension in (2, 3, 5, 8):
        a = s.zeros(dimension)
        for k in range(1, dimension):
            a[k-1,k] = s.sqrt(k)
        comm = a*a.H-a.H*a
        expected = s.eye(dimension)
        expected[-1,-1] = 1-dimension
        require(comm == expected and s.trace(comm) == 0, 'finite CCR edge control')
        ccr.append({'dimension': dimension, 'full_commutator_trace': '0',
                    'interior_commutator_equals_identity': True,
                    'retained_edge_value': str(1-dimension)})

    # Execute the public same-source scalar countercontrol and independently
    # check literal n-fold states through three copies.
    public = json.loads((rer/'model/native_weak_current_interface/results.json').read_text())
    # The function's bound parameter does not enter the expectation/variance
    # control. A rational positive test value avoids identifying v with a mass.
    scalar = scalar_controls(Fraction(1,2))
    phi = higgs_memory()[0]
    eta = s.Matrix([s.sqrt(2),1,0,1,0])/2
    samples = []
    for n in (1,2,3):
        psi = eta
        for _ in range(n-1):
            psi = s.kronecker_product(psi, eta)
        averaged = s.zeros(5**n)
        for i in range(n):
            factors = [phi[0] if j == i else s.eye(5) for j in range(n)]
            op = factors[0]
            for factor in factors[1:]:
                op = s.kronecker_product(op, factor)
            averaged += op/n
        residual = averaged*psi-s.sqrt(2)*psi/2
        squared = s.simplify((residual.H*residual)[0])
        require(squared == s.Rational(1,4*n), 'collective scalar variance differs')
        samples.append({'copies': n, 'literal_dimension': 5**n, 'squared_residual': str(squared)})
    require(public['fixed_scalar_carrier_obstruction']['explicit_expectation_control']
            == scalar['explicit_expectation_control'], 'public collective control differs')

    paths = [ht,
             rer/'model/native_weak_current_interface/results.json',
             rer/'model/native_higgs_spectral_join/results.json',
             rer/'model/interacting_block_quadrature/PROOF.md',
             rer/'model/interacting_block_quadrature/results.json']
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path and Path(path).suffix == '.py' and Path(path).resolve().is_relative_to(rer.resolve()) and not Path(path).resolve().is_relative_to(rer/'.venv'):
            paths.append(Path(path).resolve())
    return {
        'input_paths': sorted(set(map(str, paths))),
        'native_V': {'all_full_observer_coordinates_checked': observer.dimension,
                     'selected_diagonal_coordinates': len(targets), 'exact_equal': True,
                     'retained_controller_events': record,
                     'optional_apparatus_intervention_not_idle_evolution': True},
        'kept_clock_identity': {'common_shift': str(shift), 'common_residue': str(residue),
                                'identity_exact_for_the_declared_second_order_self_energy': True,
                                'all_orders_source_theorem_claimed': False},
        'HT_exact_function': {'source': 'retained/HT/part1_c1.py.raw:common_eigvecs',
                              'unchanged_function_body_executed': True,
                              'two_level_common_eigenvectors': [], 'Z3_modulus': ['1'],
                              'D11_flow_replayed_in_this_audit': False},
        'finite_CCR_controls': ccr,
        'collective_scalar': {'literal_checks': samples,
                             'general_identity': 'For eta^tensor n and centered X_i, cross expectations vanish; ||n^-1 sum X_i eta^tensor n||^2=Var_eta(phi)/n=1/(4n).',
                             'all_complementary_scalar_states_retained': True,
                             'actual_interacting_ground_or_stationary_background': False,
                             'local_field_silently_replaced_by_average': False,
                             'proves_a_Higgs_or_SM_limit': False},
    }


def summarize(fresh, exact):
    cm = fresh['CM/r2.json']
    comparison = fresh['CM/r2_compare.json']
    summary = comparison['summary']
    require(summary['native_all_four_pass'] == summary['conditional_all_four_pass'] == 0,
            'CM match verdict changed')
    fits = [r for group in comparison['fit'].values() for r in group['smallest_four']]
    require(all(not r['verdict']['all_four_pass'] for r in fits), 'a fitted CM reading passes')
    w3 = cm['native']['W3|A|N']
    require(w3['J'] != 0 and cm['native']['W3|A|N0']['J'] == 0, 'CM clock control changed')
    dr = fresh['DR/register.json']
    require(dr['looks']['passes'] == [], 'DR registered dressing passes')
    kept = dr['rows']['P1_link_kept_(S)']
    require([r[0] for r in kept['weights']['H']] == ['1/16'], 'hold transport absent')
    require(all([r[0] for r in kept['weights'][flag]] == ['1/16'] for flag in ('R+','R-')),
            'transport flag weights unequal')
    for row in dr['rows'].values():
        require(all(x == '0' for x in row['flag_mixing_second_order'].values()), 'flag mixing changed')
    l1 = fresh['SIM/derive/lnat3/lnat3_t1.json']
    l2 = fresh['SIM/derive/lnat3/lnat3_t2.json']
    l4 = fresh['SIM/derive/lnat3/lnat3_t4.json']
    require(l1['B_matching_code']['every_logical_state_every_seam_both_colorings_has_one_pair_per_endpoint'],
            'source-balanced preparation changed')
    require(l2['separability_theorem']['product_analyzer_gives_scalar_times_C_operator'], 'separability identity')
    require(l4['Phi_Schmidt_rank_C|E'] == 2, 'FR analyzer rank changed')
    require(l4['multi_tick_unread_K_n']['n=2']['K_n']
            == {'a (I)':'sqrt(3)/9','b (R)':'0','c (R^-1)':'0'}, 'literal repeated tick differs')
    return {
        'CM2': {'native_W3_A': {k:w3[k] for k in ('V12','V23','V13','J','unitarity_defect')},
                'native_readings_including_theta_zero_controls': summary['native_readings'],
                'conditional_readings': summary['conditional_readings'],
                'fitted_assignments': len(comparison['fit']),
                'fitted_roots_retained': len(fits), 'all_four_passes': 0,
                'comparison_rule': 'All four relative deviations <=10 percent against the report brief; not a confidence interval.',
                'X8_scope': 'First-order site perturbation and 300 numerical controls; no arbitrary-strength theorem.'},
        'DR2': {'literal_second_order_flag_weights': kept['weights'],
                'kept_result': {k:kept[k] for k in ('mu_over_e','theta_eff','phi_even','balance_factors_Z+/ZH_Z-/ZH')},
                'absent_result': dr['rows']['P2_link_absent_(S)'],
                'needed_target_informed_weight': dr['needed_NOT_BLIND_absent_reading'],
                'registered_passes': [], 'clock_choice_gap_coupling_and_dressed_compression_declared': True},
        'LNAT3': {'balanced_preparation_cases': l1['B_matching_code']['cases'],
                  'joint_analyzer_Schmidt_rank': l4['Phi_Schmidt_rank_C|E'],
                  'optional_native_V_verified_exactly': exact['native_V']['exact_equal'],
                  'literal_second_tick': l4['multi_tick_unread_K_n']['n=2']['K_n'],
                  'product_readout_obstruction_only_for_listed_local_grammar': True,
                  'joint_heralded_filter_exception_retained': l2['A5_heralded_filter'],
                  'entangling_extension_impossible_claimed': False},
        'HT': {'finite_exact_CCR_excluded': True,
               'finite_local_lowering_scalar_background_excluded': True,
               'growing_register_free_scalar_limit_excluded': False,
               'fixed_modulus_to_quartic_identification_is_imported': True,
               'new_Higgs_mass_prediction_claimed': False},
        'source_derived_SM_GeV_mass': False,
    }
