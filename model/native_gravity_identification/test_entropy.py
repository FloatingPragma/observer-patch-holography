"""Independent source, full-matrix and falsification controls for the finite bridge."""
from fractions import Fraction as F
from itertools import combinations
from math import comb

import numpy as np
import pytest
from scipy.linalg import expm, logm

from model.native_capacity_restriction import restriction
from model.native_program import interpreter as native
from model.native_fair_local_preparation import preparation
from tools import strict_json
from . import entropy as e, certificate as c


def density(rows, ports):
    return np.diag([float(rows[mask.bit_count()]['eigenvalue']) for mask in range(1 << ports)])


def entropy(matrix):
    eigenvalues = np.linalg.eigvalsh(matrix)
    values = eigenvalues[eigenvalues > 1e-15]
    return float(-values @ np.log(values))


@pytest.mark.parametrize('population,occupied,retained', [(6, 0, 3), (6, 1, 3), (6, 3, 3), (6, 5, 3), (6, 6, 3), (8, 3, 4)])
def test_actual_occupation_partial_trace_and_zero_sectors(population, occupied, retained):
    states = [sum(1 << bit for bit in chosen) for chosen in combinations(range(population), occupied)]
    rho = {(mask, 0, mask, 0): F(1, len(states)) for mask in states}
    actual = native.occupation_partial_trace(rho, population, tuple(range(retained)))
    rows = e.blocks(population, occupied, retained)
    expected = {(mask, 0, mask, 0): rows[mask.bit_count()]['eigenvalue']
                for mask in range(1 << retained) if rows[mask.bit_count()]['eigenvalue']}
    assert actual == expected
    assert len(rows) == retained + 1
    assert sum(row['probability'] for row in rows) == 1


@pytest.mark.parametrize('parameters', [(240, 80, 12), (240, 2, 12), (3840, 1280, 12)])
def test_exact_multiplicity_entropy_split(parameters):
    rows, bulk, edge, full = e.finite_entropy_data(*parameters)
    assert e.combine((1, bulk), (1, edge)) == full
    total = comb(parameters[0], parameters[1])
    assert sum(row['dimension'] * row['exterior_dimension'] for row in rows) == total
    for row in e.supported(rows):
        assert row['eigenvalue'] == F(row['exterior_dimension'], total)


def test_full_entropy_identity_and_sector_reweighted_state():
    rows, _, _, _ = e.finite_entropy_data(8, 3, 3)
    original = np.array([float(row['probability']) for row in rows])
    modified = original.copy()
    modified[0] -= .01
    modified[1] += .01
    weights = np.concatenate([np.full(row['dimension'] * row['exterior_dimension'], p / (row['dimension'] * row['exterior_dimension']))
                              for row, p in zip(rows, modified, strict=True)])
    assert len(weights) == comb(8, 3)
    direct_entropy = -weights @ np.log(weights)
    relative = modified @ np.log(modified / original)
    assert abs(direct_entropy - (np.log(comb(8, 3)) - relative)) < 2e-14
    assert relative > 0


def test_regional_modular_action_on_full_number_algebra():
    ports = 4
    rho = density(e.blocks(10, 4, ports), ports)
    generator = -logm(rho)
    for a in range(1 << ports):
        for b in range(1 << ports):
            if a.bit_count() == b.bit_count():
                unit = np.zeros_like(rho)
                unit[a, b] = 1
                assert np.array_equal(generator @ unit, unit @ generator)
    observable = np.diag([float(mask & 1) for mask in range(1 << ports)])
    unitary = expm(1.7j * generator)
    assert np.linalg.norm(unitary @ observable @ unitary.conj().T - observable, 2) < 1e-14
    # The full algebra has cross-number phases; centrality is not asserted there.
    unit = np.zeros_like(rho)
    unit[0, 1] = 1
    assert np.linalg.norm(generator @ unit - unit @ generator, 2) > .1


@pytest.mark.parametrize('endpoint_bits', range(4))
def test_complete_original_E_instrument_is_not_trivial_modular_time(endpoint_bits):
    a, b = preparation.SEAMS[0]
    local = ((endpoint_bits & 1) << a) | ((endpoint_bits >> 1) << b)
    spectators = ((1 << (80 - local.bit_count())) - 1) << 12
    mask = local | spectators
    command = {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': [a, b]}, 'law': 'E', 'repeat': 1}
    program = native.program(federation={'level': 0, 'total': 80}, biased=False, commands=[command])
    rho = {(mask, 0, mask, 0): F(1)}
    native._validate_density(rho, program['carrier'])
    outcomes = native._outcomes(program, command, rho)
    assert len(outcomes) == 25
    assert sum(native.trace(outcome.density) for outcome in outcomes) == 1
    contrast = sum(value * (((x >> a) & 1) - ((x >> b) & 1))
                   for outcome in outcomes for (x, ref, y, ref2), value in outcome.density.items()
                   if x == y and ref == ref2)
    assert contrast == 0
    if endpoint_bits in (1, 2):
        assert ((mask >> a) & 1) - ((mask >> b) & 1) in (-1, 1)
    assert sum(not outcome.density for outcome in outcomes) == 23


def test_native_reference_channel_keeps_actual_weights_and_joint_correlations():
    # All terms have total occupation two when the reference label 1 carries
    # one exterior particle. The local source has no cross-number coherence.
    rho = {(1, 1, 1, 1): F(1, 2), (3, 0, 3, 0): F(1, 2),
           (1, 1, 3, 0): F(1, 2), (3, 0, 1, 1): F(1, 2)}
    target = preparation.reference_channel(rho)
    one = sum(value for (x, a, y, b), value in target.items() if x == y and a == b and x.bit_count() == 1)
    two = sum(value for (x, a, y, b), value in target.items() if x == y and a == b and x.bit_count() == 2)
    assert one == two == F(1, 2)
    assert all(x.bit_count() == 1 if a == 1 else x.bit_count() == 2 for x, a, y, b in target)
    assert target[(1, 1, 1, 1)] == F(1, 24)
    assert target[(3, 0, 3, 0)] == F(1, 132)


def test_first_law_with_noncentral_tangent_and_full_matrix_entropy():
    rows = e.blocks(8, 3, 3)
    sigma = density(rows, 3)
    m = min(row['eigenvalue'] for row in rows)
    tangent = np.zeros_like(sigma)
    tangent[0, 3] = tangent[3, 0] = 1
    # This tests the matrix-analysis extension; fixed-number source marginals
    # do not prepare this cross-number coherence.
    t = m / 4
    rho = sigma + float(t) * tangent
    slope = float(np.trace(tangent @ (-logm(sigma))))
    relative = float(np.trace(rho @ (logm(rho) - logm(sigma))))
    bound = e.first_law_remainder(m, F(2), t, F(1))
    assert np.linalg.norm(sigma @ tangent - tangent @ sigma, 2) > 0
    assert slope == 0
    assert 0 < relative < float(bound)
    assert abs(entropy(rho) - entropy(sigma) - float(t) * slope + relative) < 2e-14


def test_modular_stability_for_same_sector_preparation_error():
    rows = e.blocks(8, 3, 3)
    sigma = density(rows, 3)
    m = min(row['eigenvalue'] for row in rows)
    epsilon = m / 5
    rho = sigma.copy()
    rho[1, 1] += float(epsilon)
    rho[2, 2] -= float(epsilon)
    observable = np.zeros_like(rho)
    observable[1, 2] = observable[2, 1] = 1
    time = F(7, 5)
    unitary = expm(1j * float(time) * logm(rho))
    error = np.linalg.norm(unitary @ observable @ unitary.conj().T - observable, 2)
    assert 0 < error <= float(e.modular_flow_allowance(m, epsilon, time))
    assert abs(np.sum(np.abs(np.linalg.eigvalsh(rho - sigma))) / 2 - float(epsilon)) < 1e-16


def test_fixed_cut_counterexample_has_a_positive_two_sided_state_family():
    value = e.sector_redistribution(240, 80, 12)
    rows = e.blocks(240, 80, 12)
    assert value['edge_ratio'] == F(80, 149)
    assert value['edge_derivative'] != ()
    assert value['static_cut_area_derivative'] == 0
    assert e.combine((1, value['edge_derivative']), (1, value['bulk_derivative'])) == ()
    for sign in (-1, 1):
        probabilities = [row['probability'] + sign * value['positive_parameter_radius'] * change
                         for row, change in zip(rows, value['tangent'], strict=True)]
        assert min(probabilities) > 0
        assert sum(probabilities) == 1


def test_no_capacity_uniform_logarithm_bound_from_fixed_trace_precision():
    population = restriction.federation_population(200)
    dilute = e.supported(e.blocks(population, 2, 12))
    dense = e.supported(e.blocks(population, population // 3, 12))
    epsilon = F(1, 2**64)
    assert min(row['eigenvalue'] for row in dilute) < epsilon
    floor = min(row['eigenvalue'] for row in dense)
    assert floor > epsilon
    assert e.modular_flow_allowance(floor, epsilon, F(1)) < F(1, 10**12)
    with pytest.raises(ValueError):
        e.modular_flow_allowance(min(row['eigenvalue'] for row in dilute), epsilon, F(1))


@pytest.mark.parametrize('parameters', [(True, 1, 1), (4, True, 2), (4, 2, True), (4, -1, 2), (4, 5, 2), (4, 2, 5)])
def test_source_integer_domain_fails_closed(parameters):
    with pytest.raises(ValueError):
        e.blocks(*parameters)


@pytest.mark.parametrize('values', [(F(0), F(0), F(1)), (F(1, 4), F(1, 4), F(1)), (F(1, 4), F(-1), F(1)), (F(1, 4), .01, F(1))])
def test_modular_bound_rejects_missing_spectral_premise(values):
    with pytest.raises(ValueError):
        e.modular_flow_allowance(*values)


def test_first_law_rejects_support_changing_direction():
    rows = e.blocks(240, 2, 12)
    with pytest.raises(ValueError):
        e.first_law_expression(rows, tuple(F(-1 if k == 0 else 1 if k == 12 else 0) for k in range(13)))
    with pytest.raises(ValueError):
        e.first_law_remainder(F(1, 4), F(2), F(1, 4), F(1))


def test_explicit_macroscopic_entropy_dimensions_are_not_materialized():
    with pytest.raises(ValueError):
        e.finite_entropy_data(10000, 4000, 12)


def test_retained_receipt_and_scope():
    value = c.report()
    assert c.canonical(strict_json.load(c.OUTPUT)) == c.canonical(value)
    assert value['scope']['finite_entropy_split_and_first_law']
    assert not value['scope']['same_stress_G2_and_G3_derived']
    assert not value['scope']['geometric_G4_area_and_stationarity_derived']
    assert not value['scope']['Einstein_equation_or_OPH_program_no_go']
