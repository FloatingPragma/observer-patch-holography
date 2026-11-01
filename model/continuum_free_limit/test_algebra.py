from . import algebra


def test_cubic_field_creates_three_excitations_without_a_linear_component():
    row=algebra.cubic_control()
    assert row['cubic_addition_and_removal_energy']=='3'
    assert row['linear_addition_and_removal_energy']=='1'
    assert row['each_charge_weight']=='1/2'
    assert row['linear_cubic_overlap']=='0'


def test_vacuum_only_agreement_does_not_replace_outside_cyclicity():
    row=algebra.vacuum_only_countercontrol()
    assert row['vacuum_motion_alone_proves_full_free_dynamics'] is False
    assert row['two_particle_energy']=='3 rather than 2'
