from fractions import Fraction as F
import copy
import numpy as np
import pytest
import sympy as s
from . import defect_field_join as d


def test_actual_source_parity_and_shared_core_Gram():
    value=d.spectral.source();h=value['h'];points=value['points'];core=value['core_indices']
    parity=np.repeat([(-1)**sum(p)for p in points],2)
    assert len(core)==4 and all(parity[j]==-1 for j in core)
    outside=[j for j in range(h.shape[0])if j not in core]
    b=h[outside][:,core]
    support=np.flatnonzero(np.asarray(abs(b).sum(axis=1)).ravel())
    assert all(parity[outside[j]]==1 for j in support)
    assert np.max(abs(np.linalg.eigvalsh((b.conj().T@b).toarray())-[1.25,1.25,1.75,1.75]))<1e-12
    # Delete all onsite blocks before checking the one-hop Schur term structurally.
    t=h.copy().tolil()
    for j in range(len(points)):t[2*j:2*j+2,2*j:2*j+2]=0
    t=t.tocsr()[outside][:,outside]
    assert (b.conj().T@t@b).nnz==0
    assert np.max(abs((h@value['C']+value['C']@h.conjugate()).data),initial=0)<1e-12


@pytest.mark.parametrize('m',[8,16,64,128])
def test_exact_parity_bound_controls_omitted_even_tail(m):
    r=d.parity_schur(m);b2=F(7,4);radius=b2/(m-4);q=2/(m-radius)
    partial=sum(q**(2*j)for j in range(1,20))
    assert b2*radius/(m*(m-radius))+b2/(m-radius)*partial<F(r['parity_resolvent_error'])
    assert F(r['parity_resolvent_error'])<1/F(4*m)
    assert r['cluster_count_each_sign']==2


def test_leading_Schur_uses_actual_texture_not_spectral_fit():
    r=d.topology.leading_schur(1)
    matrix=s.Matrix([[s.sympify(x)for x in row]for row in r['coefficient_of_inverse_M']])
    assert matrix.eigenvals()=={-s.Rational(1,4):2,s.Rational(1,4):2}
    assert r['this_term_is_the_full_spectrum']is False
    refined=d.parity_schur()
    assert F(refined['parity_resolvent_error'])==F(2365097,60283957504)<F(1,25000)


def test_complete_static_branches_and_actual_positive_gap():
    r=d.phase_zero_branches()
    assert len(r['rows'])==4
    assert [x['positive_Nambu_rows']for x in r['rows']]==[8,8,4,4]
    assert sum(x['positive_Nambu_rows']for x in r['rows'])==24
    assert all(F(x['positive_frequency_interval'][0])>F(1,40000)for x in r['rows'])
    assert r['three_observed_particle_species_assigned']is False
    assert r['within_group_distinct_splittings_proved']is False


def test_charge_pairing_and_uniform_free_action_are_distinct():
    r=d.local_pairing_identities()
    assert r['exact_complex_symmetric_Nambu_control_dimension']==16
    assert r['uniform_pairing_defect_Frobenius_squared']=='32768'
    assert r['same_uniform_free_action_claimed']is False


def test_localization_keeps_original_complement():
    r=d.localization_readout();p=F(r['normalized_core_probe_original_band_weight_lower'])
    assert F(99,100)<p<1
    assert r['other_modes_complement_retained']is True
    assert r['probability_is_conditioned_on_band']is False
    assert r['separate_addition_or_removal_weight_bound']is False


def test_full_Fock_ground_cannot_be_silently_reused():
    r=d.ground_countercontrol()
    assert r['complete_energies']==['-1/100','-1/256','1/256','1/100']
    assert r['old_filled_sea_stationary']is True
    assert r['old_filled_sea_is_paired_ground']is False
    assert r['electric_charge_conserved']is True


@pytest.mark.parametrize('bad',[True,7,0.5,'-1'])
def test_undeclared_source_regime_rejected(bad):
    with pytest.raises(ValueError):d.parity_schur(bad)


def test_retained_owner_is_authenticated_without_cold_eigensolve():
    r=d.source_receipt()
    assert r['retained_complete_spectrum_dimension']==686
    assert r['cold_spectrum_recomputed']is False
    assert r['floating_splittings_used_as_proof']is False


def test_same_native_Higgs_currents_and_detached_output():
    r=d.certificate();before=copy.deepcopy(r)
    assert r['same_current_Higgs']['charged_internal_family_supports']=={'-1':3,'1':3}
    assert r['same_current_Higgs']['localized_Nambu_rows']==360
    r['same_current_Higgs']['charged_internal_family_supports']['-1']=99
    assert d.certificate()==before
    assert r['physical_particle_assignments']==[]


def test_original_energy_instrument_keeps_both_charge_complements():
    r=d.ground_countercontrol()
    assert r['all_input_effect_sum_identity']is True
    assert r['complete_energy_records']==['1/2','0','0','1/2']
    assert sum(F(x)for x in r['complete_energy_records'])==1
