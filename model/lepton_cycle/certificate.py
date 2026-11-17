"""Recompute native source, symbolic channel identities and finite-cycle consequences."""
from __future__ import annotations
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from . import composition as c, instrument as i, sensitivity

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def symbolic_checks():
    z, w = sp.symbols('z w', nonzero=True)
    p = sp.Rational(2, 9)
    weak = (1-p+p*z)*(1-p+p/w)+p*(1-p)*(z-1)*(1/w-1)
    e, f = sp.symbols('e f', real=True)
    u, ub = z-1, 1/z-1
    complete_diagonal = (1+u*e)*(1+ub*e)+u*ub*(e-e**2)
    complete_offdiagonal = (1+u*e)*(1+ub*f)-u*ub*e*f
    omega = -sp.Rational(1, 2)+sp.sqrt(3)*sp.I/2
    return {
        'weak_complete_channel_including_nonreturns': sp.simplify(weak-(1-p+p*z/w)) == 0,
        'mass_diagonal_normalization': sp.simplify(complete_diagonal-1) == 0,
        'mass_offdiagonal_multiplier': sp.simplify(complete_offdiagonal-(1+u*e+ub*f)) == 0,
        'equal_analyzer_flags_make_pulled_vectors_orthogonal': all(
            sp.simplify(sum(omega**(r*d) for r in range(3))/3) == int(d == 0) for d in (0,1,2)),
        'positive_nonreturn_gram_weight': i.TRACE_EFFECT < 1,
        'native_pair_to_double_load_rate_factor': i.TRACE_EFFECT / F(2,9) == F(5,4),
    }


def report():
    records, weights = i.source_data()
    checks = symbolic_checks()
    checks['native_address_flags_50_5_5'] = [records.count(r) for r in (0,2,3)] == [50,5,5]
    h, theta = 1/100, 2/9
    weak = i.weak_instrument(h)
    mass = i.mass_instrument(theta,h)
    source, e, _ = i.analyzer_data(theta)
    source_matrix = source.reshape(-1,3)
    checks['actual_source_gram_equals_diagonal_effect'] = bool(np.max(abs(source_matrix.conj().T@source_matrix-np.diag(e))) < 2e-12)
    checks['actual_effect_trace_5_18'] = bool(abs(sum(e)-float(i.TRACE_EFFECT)) < 2e-12)
    H = i.address_preparation_unitary()
    checks['explicit_preparation_is_unitary_and_uniform'] = bool(np.max(abs(H.T@H-np.eye(60))) < 2e-13 and
        np.max(abs(H[:,0]-1/np.sqrt(60))) < 2e-13)
    for name,parts,expected in (('weak',weak,i.unconditional_weak_multiplier(h)),
                               ('mass',mass,i.unconditional_mass_multiplier(theta,h))):
        grams = i.gram_parts(parts)
        checks[name+'_complete_actual_channel_matches_formula'] = bool(np.max(abs(sum(grams.values())-expected)) < 2e-12)
        checks[name+'_all_branch_grams_positive'] = all(np.linalg.eigvalsh(g).min() > -2e-12 for g in grams.values())
        checks[name+'_four_outcome_history_sums_to_channel'] = bool(np.max(abs(sum(c.history_multipliers(grams,3).values())-expected**3)) < 2e-12)
        checks[name+'_stationary_diagonal_register'] = bool(np.max(abs(np.diag(sum(grams.values()))-1)) < 2e-12)
        checks[name+'_nontrivial_finite_cycle_not_unitary'] = bool(np.max(abs(expected-np.ones((3,3)))) > 1e-8 and abs(expected[0,1]) < 1-1e-10)
        checks[name+'_reset_preserves_exported_source_and_reference_map'] = all(
            np.array_equal(i.reset_branch(label, block).archived_source_register, block)
            and np.array_equal(i.reset_branch(label, block).fresh_source, i.preparation())
            for label,block in parts.items())
    checks['weak_no_reset_control_changes_count'] = rounds_without_reset()
    if not all(checks.values()):
        raise ArithmeticError('failed lepton-cycle checks: '+str([k for k,v in checks.items() if not v]))
    with mp.workdps(60):
        lam=sensitivity.spectrum(mp.mpf(2)/9,mp.sqrt(2),mp.sqrt(mp.mpf(5)/108))
        rates=[mp.nstr(x*x,25) for x in lam]
    # Bind this producer's inputs explicitly; independent comparison extensions
    # have their own receipts and must not introduce a circular receipt graph.
    files = [HERE/p for p in ('PROOF.md','README.md','__init__.py','__main__.py',
        'certificate.py','composition.py','instrument.py','sensitivity.py','test_cycle.py')]
    files += [ROOT/p for p in ('model/lepton_register_round/rounds.py','model/lepton_register_round/walk.py',
        'model/lepton_register_round/carrier.py','model/coherent_seam_repair/common.py',
        'model/coherent_seam_repair/program.py','Lean/EventAlgebra/SlotLocalitySurface.lean')]
    return {'schema':'oph.complete_lepton_cycle.v1', 'checks':checks,
            'source_bindings':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            'instrument':{'source_address_dimension':3600,'register_dimension':3,'outcomes':list(i.LABELS),
                          'weak_exclusive_probability':'2/9',
                          'weak_projected_X_norms':{k:str(weights[k]) for k in i.LABELS},
                          'mass_effect_trace':'5/18','mass_effect_eigenvalues_k0_k1_k2':rates,
                          'mass_nonreturn_gram':'diag(e)-e e^T',
                          'discarded_or_renormalized_outcomes':False,
                          'reset':'SWAP with a freshly prepared pair; old source and its correlations exported'},
            'reference':{'fixed_register_states':'Every density diagonal in the indicated three-character/current basis.',
                         'unique_stationary_reference':False,
                         'archive_stationary':False,'physical_Gibbs_or_vacuum_identification':False},
            'finite_budgets':[c.error_budget(n,cycle=kind) for kind in ('weak','mass') for n in (10,100,1000)],
            'rounds_for_ideal_trace_distance_1e_minus_6':{kind:c.required_rounds(F(1,10**6),cycle=kind) for kind in ('weak','mass')},
            'implementation_designs_1e_minus_6':{kind:c.implementation_design(F(1,10**6),cycle=kind) for kind in ('weak','mass')},
            'resources_for_1000_complete_cycles':c.resources(1000),
            'sensitivity':sensitivity.report(),
            'viewer_projection':{
                'weak':{'probability':'2/9','current_eigenvalues':[-1,0,1],
                    'return_probability':'1 - (14/81) * 4*sin(h*j/2)^2',
                    'nonreturn_probabilities':{'x_no':'(11/324) * 4*sin(h*j/2)^2',
                        'y_no':'(125/1296) * 4*sin(h*j/2)^2','both_no':'(55/1296) * 4*sin(h*j/2)^2'},
                    'limiting_phase_per_unit':'2/9','complete_channel_half_diamond_bound':'14*t^2/(81*N)'},
                'mass':{'effect_trace':'5/18','native_pair_effect_trace':'2/9',
                    'double_load_to_native_pair_common_rate_factor':'5/4',
                    'rates_k0_k1_k2':rates,'complete_channel_half_diamond_bound':'349*t^2/(1296*N)',
                    'nonreturn_environment_rank':3},
                'reset':{'fresh_pair_required_each_cycle':True,'all_nonreturn_states_exported':True,
                         'archive_stationary':False},
                'sensitivity':sensitivity.report(),
                'physical_pole_mass_identification':False},
            'scope':{'derived':'Complete finite instruments, declared-clock weak limits, reference families and input sensitivities.',
                     'supplied':'Fresh double loads, uniform-address preparation gate, analyzer, controlled kick, inverse ticks, return apparatus, SWAP/archive services and declared clock.',
                     'not_claimed':'Autonomous preparation/reset, a physical mass pole, GeV scale, interacting relativistic limit or new Lean formalization.'}}


def rounds_without_reset():
    from model.lepton_register_round import carrier, rounds
    return rounds.without_reset(carrier.carrier_walk())['counts'] == (F(2,9),F(7,9),F(7,9),F(2,9))
