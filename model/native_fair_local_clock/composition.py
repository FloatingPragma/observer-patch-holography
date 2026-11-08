"""One source programme for a prepared reference and complete local readout."""
from fractions import Fraction as F

from model.native_fair_local_preparation import preparation as prep
from model.native_fair_local_readout import readout
from model.weighted_source_family import family as native_graph
from . import clock as c


def unitality_certificate():
    """All actual fair E raw slots on the full endpoint occupation factor."""
    graph = native_graph.Family((1, 1), ((0, 1),), (1,))
    identity = {(x, 0, x, 0): F(1) for x in range(4)}
    rows = prep.readout.cone.fixed_instrument(graph, identity, (0, 1), law='E')
    output = {}
    for _, block in rows:
        for key, value in block.items():
            output[key] = output.get(key, F(0))+value
        if any(x.bit_count() != y.bit_count() for x, _, y, _ in block):
            raise ArithmeticError('an endpoint output left its number block')
    if output != identity or len(rows) != 25:
        raise ArithmeticError('complete fair E is not exactly unital')
    # Each nonzero Kraus map separately preserves number, stronger than the
    # diagonal identity control. Enumerate its four input columns directly.
    checks = 0
    for x in range(4):
        column = prep.readout.cone.fixed_instrument(graph, {(x, 0, x, 0): F(1)}, (0, 1), law='E')
        for _, block in column:
            if any(y.bit_count() != x.bit_count() for y, _, _, _ in block):
                raise ArithmeticError('a complete endpoint branch changed particle count')
            checks += 1
    return {'raw_slots': 25, 'endpoint_input_columns': 4,
            'number_preserving_column_checks': checks,
            'sum_K_K_adjoint_exact_identity': True,
            'fixed_global_number_identity_preserved_by_tensor_extension': True,
            'external_reference_independence_required': False,
            'canonical_source_marginal_required_for_canonical_conclusion': True}


def combined_error(preparation_error, coverage_error, ambiguous_weight, *, global_constant=False):
    """D J_b to O J_* after marginalizing only mixing-record values.

    Ambiguous weight is the actual source marginal at readout start (equal
    to its initial value when preparation is internal to the detector).
    It is not inferred from an arbitrary input or an archive factorization.
    """
    for value in (preparation_error, coverage_error, ambiguous_weight):
        c.limit.rational(value, upper=F(1))
    if type(global_constant) is not bool:
        raise ValueError('an explicit global singleton-sector flag required')
    oracle = readout.error_bound(coverage_error, ambiguous_weight, global_constant=global_constant)
    return {'preparation_reference_halfdiamond_upper': preparation_error,
            'actual_preparation_to_its_own_oracle_upper': oracle,
            'actual_protocol_to_reference_preparation_oracle_upper': min(F(1), preparation_error+oracle)}


def canonical_protocol(level=6, blocks=64):
    """Bound one admitted long protocol, including the finite routed detector."""
    c.native._integer(level); prep.integer(blocks, 0, prep.MAX_BLOCKS)
    V, M = 240*4**level, 630*4**level
    word = tuple((edge, 'E') for edge in prep.edges(blocks)); L = len(word)
    partition = (F(0), F(L, M), F(L, M)+64)
    plan = c.compile_plan(level, V//3, partition, (word, ()), ('E', 'E'))
    if plan['program']['commands'][:-1] != prep.commands(blocks):
        raise ArithmeticError('the actual preparation prefix changed')
    if plan['epochs'][0]['uniform_attempts'] or plan['epochs'][1]['fixed_attempts']:
        raise ArithmeticError('the protocol inserted an uncounted background command')
    if plan['program']['commands'][-1] != readout.uniform_program(level, V//3, 64*M)['commands'][0]:
        raise ArithmeticError('the complete readout epoch differs from its owning source')
    outer = c.routing.routed_window(level, 5)
    inner = c.routing.routed_window(level, 4)
    keep = set(inner['source_ports'])
    indices = tuple(i for i, source_port in enumerate(outer['source_ports']) if source_port in keep)
    detector = readout.rooted_plan(outer, indices)
    if not {v for edge in prep.SEAMS for v in edge} <= set(detector['ports']):
        raise ArithmeticError('preparation is not internal to the complete detector')
    if len(detector['ports']) != 348 or len(detector['tree']) != 347:
        raise ArithmeticError('the authenticated finite detector changed')
    coverage = readout.coverage(M, len(detector['tree']), 64*M)
    eta = readout.canonical_erasure(V, V//3, len(detector['ports']), use_global_total=True)
    error = combined_error(F(1, 2**blocks), coverage['upper'], eta['product_upper'])
    if blocks == 64 and (error['actual_protocol_to_reference_preparation_oracle_upper'] >= F(1, 10**16)
                          or eta['product_upper'] >= F(1, 10**60)):
        raise ArithmeticError('the composed complete-instrument error exceeds its stated gate')
    return {'kind': 'canonical_source_marginal_with_arbitrary_reference_extension',
            'plan': plan, 'population': V, 'global_total': V//3,
            'detector_source_plan': detector,
            'detector_plan_digest': c.native.digest(detector),
            'detector_ports': len(detector['ports']), 'detector_tree_edges': len(detector['tree']),
            'detector_incident_seams': len(detector['incident_seams']),
            'prep_tree_records_retained': len(prep.TREE),
            'prep_mixing_record_values_marginalized': 60*blocks,
            'source_deletes_mixing_records': False,
            'readout_refers_to_postpreparation_occupations': True,
            'initial_source_marginal_is_declared_canonical': True,
            'external_reference_may_be_correlated': True,
            'prep_archive_independence_asserted': False,
            'unread_source_marginal_preserved_exactly': True,
            'coverage_upper': coverage['upper'],
            'initial_and_postprep_ambiguous_weight_upper': eta['product_upper'],
            'error': error, 'short_time_Poisson_budget_used': False,
            'expanded_host_completion_asserted': False}


def sparse_protocol(level, total=2, blocks=64):
    """Exact fresh root-tree readout when the global sector fixes constants.

    This supports arbitrary source/reference inputs in that global sector.
    Its number of attempts is independent of capacity; per-seam duration is
    consequently capacity dependent.
    """
    c.native._integer(level); c.native._integer(total)
    V, M = 240*4**level, 630*4**level
    c.native._integer(total, 0, V)
    feasible = readout.feasibility(V, total, 12)
    if feasible['ambiguous_constant']:
        raise ValueError('this global sector does not remove local constant ambiguity')
    word = tuple((edge, 'E') for edge in prep.edges(blocks))
    tree = tuple((edge, 'E') for edge in prep.TREE)
    L, J = len(word), len(tree)
    plan = c.compile_plan(level, total, (F(0), F(L, M), F(L+J, M)), (word, tree), ('E', 'E'))
    if (plan['program']['commands'][:L] != prep.commands(blocks) or
            plan['uniform_attempts'] != 0 or plan['total_attempts'] != 22+60*blocks):
        raise ArithmeticError('the composed internal E word changed')
    return {'kind': 'known_global_sector_with_unambiguous_local_constants',
            'plan': plan, 'population': V, 'global_total': total,
            'constant_feasibility': feasible,
            'feasible_local_numbers': tuple(range(max(0, total-(V-12)), min(12, total)+1)),
            'prep_tree_records_retained': 11,
            'prep_mixing_record_values_marginalized': 60*blocks,
            'fresh_readout_attempts': 11, 'all_fresh_readout_records_retained': True,
            'readout_oracle_error': F(0),
            'actual_protocol_to_reference_preparation_oracle_upper': F(1, 2**blocks),
            'arbitrary_source_reference_input_in_declared_global_sector': True,
            'same_positive_per_seam_duration_across_capacity_asserted': False,
            'canonical_input_or_selected_density_required': False,
            'expanded_host_completion_asserted': False}
