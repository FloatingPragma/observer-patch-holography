from fractions import Fraction as F
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from . import ward as w
from ..full_fock_field.field import FieldGeometry, FullFock
from ..full_fock_energy.energy import local_energies


@pytest.mark.parametrize('sides', [(2,1,1),(1,2,1),(1,1,2),(3,1,1),(2,2,1)])
def test_actual_full_fock_matrix_binding(sides):
    geometry = FieldGeometry(sides)
    symbolic = w.energies(geometry)
    independent = local_energies(geometry)
    for a, b in zip(symbolic, independent.local, strict=True):
        np.testing.assert_array_equal(a.matrix(geometry.modes), b.matrix(independent.scale).toarray())
    np.testing.assert_array_equal(w.total(symbolic).matrix(geometry.modes),
                                  FullFock(geometry.modes).hamiltonian(geometry).toarray())


def test_clifford_multiplication_independently_on_all_four_mode_monomials():
    # Exhausts every product of 256 Majorana monomials in the complete dimer.
    matrices = [w.Clifford({i: 1}).matrix(4) for i in range(256)]
    for i in range(256):
        for j in range(256):
            mask, sign = w.multiply_masks(i,j)
            np.testing.assert_array_equal(matrices[i]@matrices[j], sign*matrices[mask])


def test_range_two_routing_and_nonzero_interaction():
    geometry = FieldGeometry((3,1,1))
    data = w.certify(geometry)
    assert (0,2) in data['currents']
    free = w.certify(FieldGeometry((3,1,1), interaction=0.))
    assert data['currents'][0,1] != free['currents'][0,1]
    assert w.divergence(data['currents'],3) == w.divergence(data['routed'],3)
    assert data['routed'][0,1] == data['currents'][0,1]+data['currents'][0,2]


def test_all_shortest_routes_and_incidence():
    sites = FieldGeometry((3,3,3)).sites
    for i, x in enumerate(sites):
        for j, y in enumerate(sites):
            distance = sum(abs(a-b) for a,b in zip(x,y))
            if not 1 <= distance <= 2:
                continue
            weights = w.route_pair(sites,i,j)
            incidence = [F() for _ in sites]
            for (a,b),weight in weights.items():
                incidence[a] += weight
                incidence[b] -= weight
                assert sum(abs(c-d) for c,d in zip(sites[a],sites[b])) == 1
            assert incidence == [F(k == i)-F(k == j) for k in range(len(sites))]
    assert len(w.shortest_paths((0,0,0),(1,1,0))) == 2
    assert len(w.shortest_paths((0,0,0),(2,0,0))) == 1


def test_interior_routing_weight_and_support_bounds():
    geometry = FieldGeometry((7,7,7))
    sites = geometry.sites
    edge = (sites.index((3,3,3)), sites.index((4,3,3)))
    mass = F()
    for i,x in enumerate(sites):
        for j in range(i+1,len(sites)):
            y=sites[j]
            if 1 <= sum(abs(a-b) for a,b in zip(x,y)) <= 2:
                mass += abs(w.route_pair(sites,i,j).get(edge,0))
    assert mass == 7


def test_arbitrary_block_maps_compose_and_cancel():
    data = w.certify(FieldGeometry((3,2,1)))
    labels = (2,0,1,0,2,1)  # Non-geometric, non-monotone: stronger algebra check.
    outer = (1,0,1)
    direct = tuple(outer[a] for a in labels)
    assert w.aggregate(w.aggregate(data['local'],labels),outer) == w.aggregate(data['local'],direct)
    assert w.aggregate_flows(w.aggregate_flows(data['routed'],labels),outer) == w.aggregate_flows(data['routed'],direct)
    assert w.divergence(w.aggregate_flows(data['routed'],labels),3) == w.aggregate(w.divergence(data['routed'],6),labels)
    assert not w.aggregate_flows(data['routed'],(0,)*6)


def test_full_three_dimensional_certificate():
    data = w.certify(FieldGeometry((3,3,3)))
    assert len(data['currents']) == 81
    assert len(data['routed']) == 54
    assert data['rectangular_partitions'] == 27
    assert any(sum(abs(a-b) for a,b in zip(FieldGeometry((3,3,3)).sites[x],
                                         FieldGeometry((3,3,3)).sites[y])) == 2
               for x,y in data['currents'])


def test_time_cell_uses_integrated_flux_and_composes():
    from .packet import time_cell_diagnostic
    result = time_cell_diagnostic()
    assert result['cell_identity_operator_residual'] < 1e-12
    assert result['time_composition_operator_residual'] < 1e-12
    assert result['endpoint_sample_counterexample_norm'] > 1e-4


@pytest.mark.parametrize('kwargs', [dict(radius=3),dict(duration=-1),dict(duration=.1),
                                    dict(steps=0),dict(speed=.125),dict(interaction=-1)])
def test_reject_invalid_error_domain(kwargs):
    args = dict(radius=12,duration=F(1,8),steps=10000)
    args.update(kwargs)
    with pytest.raises(ValueError):
        w.readout_error(**args)


def test_error_budget_improves_with_source_steps_and_buffer():
    a=w.readout_error(12,F(1,8),100000)
    b=w.readout_error(18,F(1,8),10000000)
    assert b['cell_Ward_residual_bound'] < a['cell_Ward_residual_bound']
    assert b['cell_Ward_residual_bound'] < F(1,10**8)
    assert w.readout_error(12,0,10)['cell_Ward_residual_bound'] == 0


def test_clifford_rejects_float_coefficients_and_outside_carrier():
    with pytest.raises(ValueError):
        w.Clifford({1:.125})
    with pytest.raises(ValueError):
        w.Clifford({256:1}).matrix(4)


def test_malformed_receipt_path_is_checked(tmp_path):
    receipt=tmp_path/'bad.json'
    receipt.write_text('{"all_pass":true}')
    proc=subprocess.run([sys.executable,'-m','model.coarse_energy_ward','--check',str(receipt)],
                        capture_output=True,text=True)
    assert proc.returncode != 0


def test_midpoint_and_readout_errors_are_explicit():
    exact=w.readout_error(12,F(1,8),10000)
    sampled=w.readout_error(12,F(1,8),10000,samples=1000)
    assert sampled['midpoint_quadrature_error'] > 0
    assert sampled['cell_Ward_residual_bound']-exact['cell_Ward_residual_bound'] == 6*sampled['midpoint_quadrature_error']
    noisy=w.readout_error(12,F(1,8),10000,samples=1000,
                           endpoint_readout=F(1,1000),face_readout=F(1,100))
    assert noisy['cell_Ward_residual_bound']-sampled['cell_Ward_residual_bound'] == F(2,1000)+6*F(1,8)*F(1,100)


@pytest.mark.parametrize('kwargs',[dict(samples=0),dict(samples=1.5),
                                  dict(endpoint_readout=-1),dict(face_readout=.1),
                                  dict(duration=100)])
def test_reject_invalid_sampling_or_inherited_tail_domain(kwargs):
    args=dict(radius=12,duration=F(1,8),steps=10000)
    args.update(kwargs)
    with pytest.raises(ValueError):
        w.readout_error(**args)
