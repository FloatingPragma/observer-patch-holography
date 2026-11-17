"""Replay complete taste records and the finite same-action thermal map."""
from pathlib import Path
import hashlib
from fractions import Fraction as F
from . import tastes,thermal

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


def source_bindings():
    paths=list(HERE.glob('*.py'))+[HERE/'README.md']
    paths.extend(ROOT/p for p in ('model/bound_pair/fiber.py','model/bound_pair/spectrum.py',
        'model/interacting_field/lattice.py','model/observer/scalar_response/experiment.py'))
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def run():
    taste=tastes.evolution_report()
    positive=thermal.report()
    probability=taste['complete_evolutions'][0]['unconditioned_outside_probability']
    bound=taste['short_time_bound']
    checks={
        'all_eight_nodes_retained':len(taste['all_nodes'])==8 and taste['net_Weyl_node_chirality']==0,
        'single_pocket_has_nonzero_exact_obstruction':F(taste['commutator']['commutator_operator_norm_squared_exact'])>0,
        'complete_space_reduces_action':taste['all_pockets_control']['interacting_pocket_reduces_full_action'],
        'complete_probabilities_normalized':all(r['normalization_residual']<2e-12 for r in taste['complete_evolutions']),
        'short_time_positive_leakage_enclosed':float(F(bound['unconditioned_outside_probability_lower']))<probability<float(F(bound['unconditioned_outside_probability_upper'])),
        'all_seven_removed_pockets_observable':all(p>1e-5 for p in taste['complete_evolutions'][-1]['all_eight_pocket_probabilities'][1:]),
        'physical_fermion_exchange_preserved':all(r['physical_exchange_residual']<2e-12 for r in taste['complete_evolutions']),
        'same_action_faithful_Gibbs':positive['smallest_Gibbs_eigenvalue']>0,
        'reflection_positivity':positive['reflection_Gram_minimum_eigenvalue']>-2e-12 and positive['independent_trace_formula_residual']<2e-12,
        'KMS_and_wrong_reference_control':positive['KMS_residual']<2e-12 and positive['wrong_tracial_reference_KMS_residual']>1e-6}
    return {'schema':'oph.same_field_tastes.v1','source_bindings':source_bindings(),
            'owner_issues':[1142,1137],'taste_interaction':taste,'finite_thermal_positivity':positive,
            'checks':checks,'all_pass':all(checks.values()),
            'integrated_gate_accepted':False,'physical_comparison_count_added':0,
            'exclusion_scope':'Proper exchange-symmetric momentum-pocket projection with nonzero and nonunit contact weight at the fixed supplied contact action.',
            'not_excluded':['Other regulators','Dressed invariant subspaces','Different interactions','Controlled decoupling or continuum limits with an actual error map']}
