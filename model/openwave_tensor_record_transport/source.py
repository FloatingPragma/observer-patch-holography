"""Exact interpolation of classical tensor records, with its retained error."""
from fractions import Fraction as F
from itertools import product
import sympy as s
from model.openwave_staggered_ward import source as staggered
from model.openwave_einstein_error import source as decoder


def move(p, axis, distance):
    q = list(p); q[axis] += distance
    return tuple(q)


def sample_sites(x, a, b):
    if a == b:
        return (x,)
    return tuple(move(move(x, a, i), b, j) for i, j in product((-1, 1), repeat=2))


def centers(tensor):
    """Every nonzero interpolation output, including the complete outside halo."""
    out = set()
    for (a, b), field in tensor.items():
        for p in field:
            out.update(sample_sites(p, a, b))
    assert all(all(v % 2 == 0 for v in p) for p in out)
    return sorted(out)


def colocate(tensor, x):
    if any(v % 2 for v in x):
        raise ValueError('a scalar-lattice record location is required')
    return s.Matrix(4, 4, lambda a, b: s.Rational(component(tensor, x, a, b)))


def component(tensor, x, a, b):
    locations = sample_sites(x, a, b)
    return sum((tensor[a,b].get(y, F()) for y in locations), F()) / len(locations)


def variation(field, axis):
    domain = set(field) | {move(p, axis, -2) for p in field}
    return max((abs(field.get(move(p, axis, 2), F()) - field.get(p, F()))
                for p in domain), default=F())


def transport_error(tensor):
    """Actual finite variation, not an assumed continuum smoothness constant."""
    rows = []
    for a in range(4):
        for b in range(a + 1, 4):
            da, db = variation(tensor[a, b], a), variation(tensor[a, b], b)
            rows.append({'component': [a, b], 'axis_variations_exact': [str(da), str(db)],
                         'transport_error_upper_exact': str((da + db) / 2)})
    return max(F(row['transport_error_upper_exact']) for row in rows), rows


def exact_transport_check(tensor):
    staggered.validate_tensor(tensor)
    delta, rows = transport_error(tensor)
    maximum = F()
    for a in range(4):
        for b in range(a + 1, 4):
            # Include zero source values that border a nonzero interpolation.
            points = set(tensor[a, b])
            for x in centers(tensor):
                points.add(move(move(x, a, 1), b, 1))
            bound = (variation(tensor[a, b], a) + variation(tensor[a, b], b)) / 2
            for p in points:
                x = move(move(p, a, -1), b, -1)
                difference = abs(tensor[a, b].get(p, F()) - component(tensor, x, a, b))
                assert difference <= bound
                maximum = max(maximum, F(difference))
    return delta, maximum, rows


def kernel_control():
    """The interpolation is lossy; the variation budget prevents false inversion."""
    n = 4
    values = {(i, j): F((-1)**(i+j)) for i, j in product(range(n), repeat=2)}
    averages = [sum(values[i % n, j % n] for i, j in
                    ((x, y), (x-1, y), (x, y-1), (x-1, y-1))) / 4
                for x, y in product(range(n), repeat=2)]
    assert set(averages) == {0}
    assert max(abs(v) for v in values.values()) == 1
    return {'periodic_cells_per_axis': n, 'all_interpolated_values_exact': ['0'] * (n*n),
            'original_maximum_absolute_value_exact': '1',
            'two_axis_variation_budget_exact': '2',
            'zero_null_records_without_transport_budget_do_not_force_zero_tensor': True}


def run():
    examples = []
    for seed in (1, 2):
        tensor = staggered.einstein(staggered.fixture(seed))
        assert all(not staggered.divergence(tensor, b) for b in range(4))
        delta, error, rows = exact_transport_check(tensor)
        sites = centers(tensor)
        # The existing exact inverse identity proves every tensor case.
        # Reuse it instead of repeating symbolic inversions at every site.
        assert decoder.design()[2]*decoder.design()[1] == s.eye(9)
        for x in sites[:3] + sites[-3:]:
            matrix = colocate(tensor, x)
            recovered = decoder.tracefree_matrix(decoder.design()[2] * decoder.charges(matrix))
            assert s.simplify(recovered - matrix + decoder.metric_coefficient(matrix)*decoder.ETA) == s.zeros(4)
        examples.append({'fixture': seed, 'complete_colocated_support_size': len(sites),
                         'maximum_actual_transport_error_exact': str(error),
                         'uniform_transport_budget_exact': str(delta), 'components': rows,
                         'all_four_staggered_divergences_exact_zero': True,
                         'six_direct_decoder_samples_checked': True,
                         'every_colocated_tensor_decodes_by_the_existing_exact_inverse_identity': True})
    return {'schema': 'oph.finite_tensor_record_transport.v1',
        'verdict': 'PASS_COMPLETE_RECORD_TRANSPORT_WITH_EXPLICIT_VARIATION_BUDGET',
        'construction': {'input': 'jointly available classical tensor records on the declared staggered lattice',
            'diagonal_rule': 'C_aa(x)=B_aa(x)',
            'off_diagonal_rule': 'C_ab(x)=sum_{i,j=+-1} B_ab(x+i e_a+j e_b)/4',
            'complete_halos_and_zero_records_retained': True,
            'record_channel': 'deterministic finite arithmetic on orthogonal classical record labels; preserve the inputs and write the result into initialized output registers',
            'arithmetic': 'rational interpolation and exact Q(sqrt(3)) null charges; finite bounded input words give finite output words',
            'memory_resource': 'initialized output registers and retained inputs are supplied; no erasure or zero-work physical apparatus theorem is claimed',
            'causal_scope': 'postprocessing at a common successor of the input records, not simultaneous acquisition of incompatible quantum fields'},
        'theorem': {'charge_error_bound': 'epsilon bounds each of the nine null contractions of C, including measurement error',
            'component_variation': 'd_i(ab)=sup_y |B_ab(y+2e_i)-B_ab(y)| on its actual staggered sublattice',
            'transport_error': 'tau=max_{a!=b}(d_a(ab)+d_b(ab))/2',
            'staggered_remainder_bound': 'delta=(2+sqrt(3))*epsilon+tau',
            'metric_term': 'lambda(x)=tr(eta*C(x))/4 on the scalar lattice',
            'connected_edge_bound': '|lambda(x+2e_j)-lambda(x)| <= a*w+8*delta',
            'path_bound': 'sum the actual edge budgets along the connected scalar graph',
            'existing_Lean_reused': ['Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean',
                                    'Lean/ObserverPatchHolography/EinsteinResidualBound.lean'],
            'averaging_inverse_assumed': False},
        'complete_exact_examples': examples, 'lost_high_frequency_control': kernel_control(),
        'scope': {'finite_staggered_to_null_decoder_error_premise_closed': True,
            'quantum_curvature_acquisition_from_native_source_proved': False,
            'physical_stress_or_metric_identification_proved': False,
            'nonlinear_or_continuum_Einstein_equation_proved': False,
            'cosmological_metric_term_determined_by_null_records': False,
            'new_empirical_fit': False, 'new_Nature_validation': False}}
