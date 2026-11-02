from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
from . import algebra,criteria,experiment,packet


def test_exact_current_force_and_pauli_orthogonality():
    value=algebra.local_certificate()
    assert value['exact_graph_axis_checks']==12
    assert value['force_trace_orthogonal_to_all_linear_spins'] is True


def test_neighbor_pairs_and_cubic_trace_coefficient():
    assert algebra.regular_graph_pair_identity(6)['coefficient']==9
    value=algebra.cubic_trace_certificate()
    assert value['normalized_plane_wave_force_norm_squared_over_J4']=='E**2/8 + 9*E/8'
    assert value['tracial_state_is_not_the_ordered_ground'] is True


@pytest.mark.parametrize('bad',[True,1,13,3.0])
def test_invalid_degree_rejected(bad):
    with pytest.raises(ValueError):algebra.regular_graph_pair_identity(bad)


def test_exact_ratio_identity_retains_nonzero_two_ray_width():
    value=criteria.ratio_moments([1,2],[F(1,2),F(1,2)])
    assert value=={'mean':F(3,2),'inverse_mean':F(3,4),'variance':F(1,4),'defect':F(1,8)}


def test_radial_energy_spread_does_not_destroy_a_true_ratio_shell():
    radii=[F(1,2),F(1)];energies=[2*r for r in radii]
    raw=criteria.ratio_moments(energies,[1,1])
    ratio=criteria.ratio_moments([e/r for e,r in zip(energies,radii)],[1,1])
    assert raw['variance']>0
    assert ratio['mean']==2 and ratio['variance']==0 and ratio['defect']==0


def test_first_and_second_residual_norm_comparison_on_broad_band():
    alpha,beta,c=F(1,4),F(8),F(2)
    pairs=[(F(1,2),F(1,3)),(F(4),F(3,4)),(F(6),F(1))]
    first=sum((e-c*s)**2 for e,s in pairs)
    second=sum((e*e-c*c*s*s)**2 for e,s in pairs)
    assert alpha*alpha*first<=second<=(beta+c)**2*first
    constants=criteria.broad_band_constants(alpha,beta,c)
    assert constants['lattice_to_euclidean_error_coefficient']==F(1,6)


@pytest.mark.parametrize('values,weights',[([True],[1]),([0],[1]),([1],[0]),([1,2],[1]),([] ,[])])
def test_invalid_ratio_measures_rejected(values,weights):
    with pytest.raises(ValueError):criteria.ratio_moments(values,weights)


def test_continuous_ward_compatible_measure_is_not_a_shell():
    value=criteria.uniform_continuum_control()
    assert value['ratio_defect_lower']=='1/27'
    assert value['scaled_energy_variance']=='1/12'
    assert value['minimum_squared_energy_shell_residual']=='34/45'
    assert value['is_an_actual_Heisenberg_ground_state'] is False
    assert value['spectral_atoms'] is False


def test_complete_finite_source_and_current_control():
    value=experiment.finite_control()
    assert len(value['rows'])==3
    assert value['momenta_in_p_le1_regime'] is False
    assert all(row['current_moment_crosscheck'] for row in value['rows'])


def test_corrupted_eigensystem_is_not_admitted(monkeypatch):
    eig=np.linalg.eigh
    def broken(h):
        values,vectors=eig(h);vectors[:,0]=vectors[:,1];return values,vectors
    monkeypatch.setattr(np.linalg,'eigh',broken)
    with pytest.raises(ArithmeticError,match='eigensystem'):experiment.finite_control()


@pytest.fixture(scope='module')
def report():return packet.report()


def test_no_missing_mass_shell_is_promoted(report):
    assert report['mass_shell_established'] is False
    assert report['sharp_lattice_pole_established'] is False
    assert report['whole_response_concentration_necessary_for_any_sharp_branch'] is False
    assert report['physical_species_identified'] is False
    assert report['new_finite_attachment_claimed'] is False
    packet.validate(deepcopy(report),report)


@pytest.mark.parametrize('key,value',[('mass_shell_established',True),('mass_shell_established',0),('new_checks_pass',1),('whole_response_concentration_necessary_for_any_sharp_branch',True)])
def test_bool_numeric_scope_laundering_rejected(report,key,value):
    bad=deepcopy(report);bad[key]=value
    with pytest.raises(ValueError):packet.validate(bad,report)


def test_source_hash_change_rejected(report):
    bad=deepcopy(report);key=next(iter(bad['source_bindings']));bad['source_bindings'][key]='0'*64
    with pytest.raises(ValueError):packet.validate(bad,report)


def test_source_change_during_computation_is_rejected(monkeypatch):
    original=packet.bindings()
    changed=dict(original)
    changed[next(iter(changed))]='0'*64
    values=iter([original,changed])
    monkeypatch.setattr(packet,'bindings',lambda:next(values))
    with pytest.raises(ArithmeticError,match='source changed during'):
        packet.report()
