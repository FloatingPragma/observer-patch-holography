"""Exact finite transfer, reflection forms and source-clock gap bounds.

These statements concern the stationary population marginal of the complete
native instrument. They do not replace the instrument by its population law.
"""
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
import json

from model.openwave_native_process import source as native
from model.openwave_native_response.source import population_column
from model.openwave_native_event_clock.source import population_and_record_check


def distances(neighbours, start=0):
    result = [-1] * len(neighbours)
    result[start] = 0
    queue = deque([start])
    while queue:
        i = queue.popleft()
        for j in neighbours[i]:
            if result[j] < 0:
                result[j] = result[i] + 1
                queue.append(j)
    if min(result) < 0:
        raise ValueError('connected finite source required')
    return result


def minus_log_bounds(value, terms=48):
    """Rational bounds for -log(value), with the complete positive tail."""
    value = F(value)
    if not 0 < value <= 1 or type(terms) is not int or terms < 1:
        raise ValueError('0 < value <= 1 and positive integer terms required')
    x = 1 - value
    partial = sum((x**k / k for k in range(1, terms + 1)), F())
    return partial, partial + x**(terms+1) / ((terms+1) * (1-x))


def apply(rows, values):
    return [sum((weight*values[j] for j, weight in row.items()), F())
            for row in rows]


def reflection_forms(rows, diagonal_functions):
    """Conditional contraction of a cylinder, retaining every Markov outcome."""
    n = len(rows)
    g = [F(1)] * n
    for values in reversed(diagonal_functions):
        g = apply(rows, g)
        g = [a*b for a, b in zip(values, g, strict=True)]
    Tg = apply(rows, g)
    site = sum((x*x for x in g), F()) / n
    bond = sum((a*b for a, b in zip(g, Tg, strict=True)), F()) / n
    return site, bond


def finite_certificate(neighbours, rows, denominator):
    """P=I-L/denominator, with exact stochastic/symmetry checks."""
    n = len(rows)
    d = distances(neighbours)
    diameter_bound = min(n-1, 2*max(d))
    if n < 2 or denominator <= 2*max(map(len, neighbours)):
        raise ValueError('strictly positive lazy connected transfer required')
    edges = [(i, j) for i, ns in enumerate(neighbours) for j in ns if i < j]
    for i, row in enumerate(rows):
        expected = {j: F(1, denominator) for j in neighbours[i]}
        expected[i] = 1-F(len(neighbours[i]), denominator)
        assert dict(row) == expected and sum(row.values()) == 1
        assert all(rows[j].get(i, F()) == p for j, p in row.items())
    floor = 1-F(2*max(map(len, neighbours)), denominator)
    lap_lower = F(4, n*diameter_bound)
    mean = F(sum(d), n)
    centered = [x-mean for x in d]
    variance_sum = sum(x*x for x in centered)
    energy = sum((centered[i]-centered[j])**2 for i, j in edges)
    rayleigh = energy/variance_sum
    q_lower, q_upper = lap_lower/denominator, rayleigh/denominator
    assert 0 < q_lower <= q_upper < 1
    gap_lower = minus_log_bounds(1-q_lower)[0]
    gap_upper = minus_log_bounds(1-q_upper)[1]
    functions = [[F((x % 3)-1) for x in d], [F(i % 2) for i in range(n)],
                 [F(i == 0)-F(1, n) for i in range(n)]]
    site, bond = reflection_forms(rows, functions)
    assert site >= 0 and bond >= floor*site
    # A conditionally contracted cylinder is only a diagnostic witness. The
    # proof covers every complex cylinder via the complete PSD transfer.
    serial = [[[j, str(w)] for j, w in sorted(row.items())] for row in rows]
    return {'states': n, 'undirected_edges': len(edges),
            'transfer': 'P=I-L/denominator', 'denominator': denominator,
            'actual_kernel_sha256': sha256(json.dumps(serial, separators=(',', ':')).encode()).hexdigest(),
            'complete_population_columns_checked': n,
            'strict_transfer_spectrum_floor_exact': str(floor),
            'connected_diameter_upper_bound': diameter_bound,
            'laplacian_gap_lower_exact': str(lap_lower),
            'laplacian_gap_rayleigh_upper_exact': str(rayleigh),
            'event_hamiltonian_gap_lower_exact': str(gap_lower),
            'event_hamiltonian_gap_upper_exact': str(gap_upper),
            'event_hamiltonian_definition': 'H_delta=-log(P)/delta',
            'event_hamiltonian_bound_units': '1/delta; delta is a declared positive clock conversion',
            'cylinder_site_reflection_form_exact': str(site),
            'cylinder_bond_reflection_form_exact': str(bond),
            'cylinder_bond_floor_margin_exact': str(bond-floor*site),
            'stationary_input': 'uniform distribution on the stated fixed-load sector',
            'vacuum': 'unique constant vector of this finite population transfer'}


def run():
    graphs = []
    for level in (0, 1, 2):
        graph = native.graph(level)
        rows = [population_column(graph, i) for i in range(graph.ports)]
        result = finite_certificate(graph.neighbours, rows, 2*len(graph.seams))
        result['level'] = level
        graphs.append(result)
    rows, actual = population_and_record_check([])
    neighbours = [tuple(j for j in row if j != i) for i, row in enumerate(rows)]
    two_load = finite_certificate(neighbours, rows, 60)
    two_load['complete_record_kernel_control'] = actual
    # A stationary occupation at one of twelve sites has P(n=1)=2/12.
    p = F(1, 6)
    kappa4 = p*(1-p)*(1-6*p+6*p*p)
    assert kappa4 == F(5, 216)
    # Positivity of probabilities alone does not imply bond reflection positivity.
    flip = [{1: F(1)}, {0: F(1)}]
    negative_site, negative_bond = reflection_forms(flip, [[F(1), F(-1)]])
    assert (negative_site, negative_bond) == (F(1), F(-1))
    return {'schema': 'oph.native_euclidean_transfer.v1',
            'verdict': 'PASS_NATIVE_FINITE_POPULATION_REFLECTION_POSITIVITY_AND_GAP_BOUNDS',
            'one_load_spatial_members': graphs, 'two_load_twelve_port_member': two_load,
            'non_gaussian_stationary_occupation_fourth_cumulant_exact': str(kappa4),
            'negative_control': {'source': 'two-state nonlazy flip',
                                 'site_form_exact': str(negative_site),
                                 'bond_form_exact': str(negative_bond)},
            'same_construction_relation': {
                'source': 'literal fair native C/E instrument with all raw records',
                'transfer': 'exact population marginal; committed histories retain their full instrument',
                'coherent_walk_relation': 'the same one-load P may be a compressed coherent continuation; real-time coherent evolution is a separate declared protocol',
                'clock': 'event ordinal; positive delta must be supplied for physical units'},
            'scope': {'finite_reflection_positive_stationary_path_measure': True,
                      'finite_logarithmic_transfer_hamiltonian_exists': True,
                      'native_record_history_replaced_by_classical_kernel': False,
                      'finite_occupation_non_gaussian': True,
                      'four_dimensional_field_limit': False,
                      'uniform_physical_continuum_gap': False,
                      'pure_yang_mills_identification': False,
                      'new_GeV_prediction': False, 'new_nature_match': False,
                      'OpenWave_criterion_awarded': False}}
