import sympy as s
from . import symmetry as m

def test_all_actual_link_coordinates_and_currents():
    row=m.link_checks()
    assert row['PW119_conjugation_coordinate_checks']==59
    assert row['PW119_antiunitary_endpoint_checks']==24

def test_full_local_polynomial_and_paths_with_seam_failure():
    assert m.local_checks()['full_F30_scalar5_Yukawa_polynomial_identities']==3
    from model.coherent_defect_wave.chiral_joint.oriented_links.rotations import rotate_words,hamiltonian_words
    from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration
    for n in (2,4):assert m.path_checks(n)['spatial_inversion_anticommutes_with_kinetic']
    words=hamiltonian_words(Configuration(side=4))
    assert rotate_words(words,s.ImmutableMatrix(-s.eye(3)),4,retain_seam=False)!=words

def test_occupied_full_H_gauss_and_antiunitary_square():
    row=m.fixture_checks()
    assert row['full_H_occupied_basis_controls']==13
    assert row['full_Gauss_occupied_basis_controls']==72

def test_native_ground_character_does_not_determine_continuous_spin():
    row=m.ground_character()
    assert row['native_cube_ground_quarter_turn_character']==1
    assert row['selected_origin_odd_channel_type']=='G1'
    assert not row['continuous_physical_spin_of_full_interacting_channel_proved']
    counts=(1,6,9,8);fund=(2,s.sqrt(2),0,1);higher=(8,0,0,1)
    assert sum(n*a*b for n,a,b in zip(counts,fund,higher))/24==1

def test_reference_cost_and_false_scalar_or_link_conjugation():
    assert 4*(30-0)==120 and 4*(30-1)!=4
    from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
    from model.coherent_defect_wave.chiral_joint.neutral_local_spectrum.symmetries import conjugate_polynomial
    parts,_=car.quantum_yukawa()
    conjugated=conjugate_polynomial(parts['charged_lepton'])
    # Undoing only the scalar exchange fails on the actual native Yukawa.
    wrong={(word,m.SCALAR_SWAP[row],m.SCALAR_SWAP[col]):v for (word,row,col),v in conjugated.items()}
    assert wrong!=parts['charged_lepton']
    from model.coherent_defect_wave.chiral_joint.peter_weyl_links.clebsch import coordinate,native_charges
    U=coordinate(native_charges()[0],0,0)
    assert U.conjugate()!=U.H


def test_public_fresh_comparison_and_antiunitary_clock():
    import copy
    from .certificate import raw_certificate,compare_research,certificate
    value=raw_certificate()
    assert compare_research(value)['all_inherited_scientific_fields_equal']
    bad=copy.deepcopy(value)
    bad['same_source_spectral_consequences']['every_origin_mode_mean_occupation']='0'
    import pytest
    with pytest.raises(ValueError):compare_research(bad)
    public=certificate()
    row=public['antiunitary_time_and_interpretation']
    assert row['source_time_evolution']=='Theta U(t) Theta^-1=U(-t)'
    assert not row['unitary_physical_CP_symmetry_claimed']
    assert not row['absence_of_measured_weak_CP_violation_inferred']
