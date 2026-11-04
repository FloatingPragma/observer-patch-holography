"""Independent complete-Fock, signed-algebra and modular-flow controls."""
import copy
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest
from scipy import sparse

from model.field_thermal.thermal import car_restriction
from model.full_fock_field.field import FieldGeometry, FullFock
from model.finite_modular_inclusion import inclusion as inc
from model.finite_modular_inclusion import packet
from tools import strict_json


REGIONS = tuple(s for n in (1, 2, 3) for s in combinations(range(4), n))


def test_exact_source_matches_independent_all_number_assembly():
    for interaction in (0, 1):
        h, scale, _ = inc.source(interaction)
        independent = FullFock(8).hamiltonian(FieldGeometry((2, 2, 1), interaction=interaction))
        np.testing.assert_array_equal(h.matrix(scale).toarray(), independent.toarray())
    h, scale, _ = inc.source()
    assert np.max(np.linalg.eigvalsh(h.matrix(scale).toarray()-np.eye(256))) <= 3
    assert np.min(np.linalg.eigvalsh(h.matrix(scale).toarray()-np.eye(256))) >= -3


def test_all_proper_source_regions_and_complete_algebra_separation():
    rows = inc.certificate()['all_proper_site_regions']
    assert {tuple(r['sites']) for r in rows} == set(REGIONS)
    assert len(rows) == 14
    space = FullFock(8)
    h = space.hamiltonian(FieldGeometry((2, 2, 1))).toarray()
    for row in rows:
        p, q = row['probe_mode'], row['exterior_mode']
        probe = space.number(p).toarray()
        exterior = space.number(q).toarray()
        derivative = h@probe-probe@h
        separated = exterior@derivative-derivative@exterior
        assert np.linalg.norm(separated, 'fro')**2 == pytest.approx(
            float(F(row['double_commutator_Hilbert_Schmidt_squared'])), abs=1e-14)
        assert abs(derivative[1 << p, 1 << q]) == 1/16
    # These nonconsecutive regions exercise the signed CAR factorization.
    # An arbitrary matrix includes odd and number-changing local operators.
    rng = np.random.default_rng(7301)
    for sites in ((0, 2), (1, 3)):
        kept = inc.checked_region(sites)
        dimension = 1 << len(kept)
        local = rng.normal(size=(dimension, dimension))+1j*rng.normal(size=(dimension, dimension))
        operators = car_restriction(8, kept)
        embedded = sum((k.conj().T@local@k for k in operators), np.zeros((256, 256), complex))
        for q in set(range(8))-set(kept):
            outside = space.number(q).toarray()
            np.testing.assert_array_equal(outside@embedded, embedded@outside)
            for p in kept:
                assert embedded[1 << p, 1 << q] == 0


def test_full_spectrum_positive_time_bound_for_every_cut():
    space = FullFock(8)
    h = space.hamiltonian(FieldGeometry((2, 2, 1))).toarray()
    energies, basis = np.linalg.eigh(h)
    t, beta = float(inc.TIME), float(inc.BETA)
    unitary = (basis*np.exp(-1j*t*beta*energies))@basis.conj().T
    evolved = {}
    for sites in REGIONS:
        p, q, _, _ = inc.crossing_modes(sites)[0]
        if p not in evolved:
            evolved[p] = unitary@space.number(p).toarray()@unitary.conj().T
        assert abs(evolved[p][1 << p, 1 << q]) >= float(inc.leakage_bound())
    assert inc.leakage_bound() == F(23, 524288)
    assert inc.leakage_bound(beta=0) == 0
    assert inc.leakage_bound(time=0) == 0


def test_independent_boundary_deletion_and_noninteracting_controls():
    space = FullFock(8)
    geometry = FieldGeometry((2, 2, 1))
    h = space.hamiltonian(geometry)
    annihilators = [space.annihilator(p) for p in range(8)]
    for sites in REGIONS:
        kept = inc.checked_region(sites)
        deleted = h.copy()
        for group in geometry.groups:
            for p, q, coefficient in group:
                if (p in kept) != (q in kept):
                    term = coefficient*(annihilators[p].conj().T@annihilators[q])
                    deleted -= term+term.conj().T
        np.testing.assert_array_equal(inc.deleted_boundary(sites).matrix(16).toarray(), deleted.toarray())
        assert inc.tensor_sum_certificate(inc.deleted_boundary(sites), sites)
        assert not inc.tensor_sum_certificate(inc.source()[0], sites)
        assert inc.crossing_modes(sites, 0) == inc.crossing_modes(sites, 1)
    # A single restored crossing term destroys the exact product control.
    sites = (0, 2)
    p, q, _, _ = inc.crossing_modes(sites)[0]
    restored = inc.deleted_boundary(sites)+inc.GI.from_exact(
        annihilators[p].conj().T@annihilators[q]+annihilators[q].conj().T@annihilators[p])
    assert not inc.tensor_sum_certificate(restored, sites)


def test_literal_preparation_error_transfer_and_adversarial_certificate():
    prior = strict_json.loads((packet.ROOT/packet.PRIOR).read_text())
    result = packet.prepared_reference_bound(prior)
    error = F(result['full_modular_generator_error'])
    lower = F(result['distance_from_every_proper_regional_algebra_lower'])
    assert error < F(11, 10**10)
    assert lower == F(23, 524288)-error/16
    assert lower > F(1, 25000)
    assert F(result['full_success_reference_minimum_weight']) > 0
    corrupted = copy.deepcopy(prior)
    corrupted['modular_certificate']['full_modular_generator_error'] = '0'
    with pytest.raises(ValueError, match='recomputation'):
        packet.prepared_reference_bound(corrupted)
    corrupted = copy.deepcopy(prior)
    corrupted['preparation_certificate']['complete_preparation_output_error'] = '1/100'
    with pytest.raises(ValueError, match='bound changed'):
        packet.prepared_reference_bound(corrupted)
    corrupted = copy.deepcopy(prior)
    name = next(iter(corrupted['source_sha256']))
    corrupted['source_sha256'][name] = '0'*64
    with pytest.raises(ValueError, match='custody'):
        packet.prepared_reference_bound(corrupted)


@pytest.mark.parametrize('sites', [(), (0, 1, 2, 3), (1, 0), (0, 0), (True,), (4,), [0]])
def test_full_region_domain_is_enforced(sites):
    with pytest.raises(ValueError):
        inc.checked_region(sites)


@pytest.mark.parametrize('kwargs', [{'beta': .1}, {'time': True}, {'time': -1}, {'generator_error': -1}])
def test_bound_requires_exact_admitted_parameters(kwargs):
    with pytest.raises(ValueError):
        inc.leakage_bound(**kwargs)


def test_fresh_process_receipt_and_promotion_controls():
    actual = packet.report()
    retained = strict_json.loads(packet.OUTPUT.read_text())
    assert packet.canonical(actual) == packet.canonical(retained)
    for key in ('geometric_null_translation_constructed',
                'continuum_half_sided_inclusion_excluded',
                'Einstein_or_physical_particle_gate_accepted'):
        changed = copy.deepcopy(actual)
        changed['finite_source'][key] = True
        with pytest.raises(ValueError, match='promotion'):
            packet.validate(changed)
    changed = copy.deepcopy(actual)
    changed['finite_source']['all_proper_site_regions'].pop()
    with pytest.raises(ValueError, match='fourteen'):
        packet.validate(changed)
