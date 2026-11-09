from fractions import Fraction as F

import pytest
import sympy as s

from model.particle_observer_outputs import readout_contract as r


def small_menu():
    return {'first': s.ImmutableSparseMatrix(s.diag(1, 0, 0)),
            'degenerate': s.ImmutableSparseMatrix(s.diag(0, 1, 1))}


def test_complete_pointer_is_unitary_and_has_the_claimed_all_input_maps():
    menu = r.pad_menu(small_menu())
    blocks = r.pointer_blocks(menu)
    unitary = s.BlockMatrix(blocks).as_explicit()
    assert unitary.H*unitary == s.eye(16)
    for i, projection in enumerate(menu.values()):
        assert blocks[i][0] == projection
    assert blocks[3][0] == s.zeros(4)
    assert r.verify_projectors(menu)['retained_outcome_ranks'] == {
        'first': 1, 'degenerate': 2, 'outside_field_encoding': 1}


def test_degenerate_coherence_and_untouched_reference_are_preserved():
    menu = small_menu()
    coherent = s.Matrix([0, 1, s.I])/s.sqrt(2)
    assert menu['degenerate']*coherent == coherent
    assert menu['first']*coherent == s.zeros(3, 1)
    # Columns are reference labels, including a nonzero imaginary reference
    # off-diagonal. No inverse branch probability is used.
    columns = s.Matrix([[1/s.sqrt(2), s.I/2], [0, s.Rational(1, 2)], [0, 0]])
    reference = columns.T*columns.conjugate()
    after = s.zeros(2)
    total = 0
    for projection in menu.values():
        branch = projection*columns
        after += branch.T*branch.conjugate()
        total += s.trace(branch.H*branch)
    assert after == reference
    assert total == 1
    assert reference[0, 1] != 0


def test_incomplete_or_nonprojective_apparatus_is_rejected():
    with pytest.raises(ValueError):
        r.verify_projectors({'missing': s.diag(1, 0)})
    with pytest.raises(ValueError):
        r.verify_projectors({'half': s.eye(2)/2, 'other': s.eye(2)/2})
    with pytest.raises(ValueError):
        r.verify_projectors({'float': s.eye(2)*1.0})


def test_complete_error_contract_and_reused_compiler():
    error = F(1, 10**6)
    assert r.composition_error(ticks=2, detector_uses=2, preparation_error=error,
                               step_error=error, detector_error=error) == F(1, 200000)
    assert r.composition_error(ticks=100, detector_uses=0, step_error=F(1, 2)) == 1
    with pytest.raises(ValueError):
        r.composition_error(ticks=1, detector_uses=1, detector_error=0.0)
    plan = r.compiled_menu_plan(small_menu(), error)
    assert plan['coherent_record_bits'] == 2
    assert plan['unused_pointer_labels_retained_as_zero'] == 1
    assert plan['arbitrary_untouched_reference']
    assert not plan['compiled_tape_executed_by_this_packet']
    clock = r.macro_clock_plan()
    assert clock['unitary_compile_plan']['native_pulse_upper_bound'] == 5245018
    assert clock['identical_compiled_word_repeated']
    assert not clock['physical_seconds_or_GeV_assigned']


def test_actual_higgs_spectrum_and_number_instrument_use_the_same_source():
    from model.particle_observer_outputs import higgs_records
    data = higgs_records.sector_data(s.diag(1, 2, 3), 1)
    menu, energies = r.energy_menu(data)
    assert list(energies.values()) == [-7, 8, 1]
    assert r.verify_projectors(menu)['retained_outcome_ranks'] == {
        'lower_energy': 1, 'upper_energy': 1, 'dark_energy': 305}
    assert r.verify_projectors(r.scalar_number_menu(data))['retained_outcome_ranks'] == {
        'empty_scalar': 1, 'occupied_scalar': 306}
    changed = dict(data, H=data['H']+s.eye(307))
    with pytest.raises(ValueError):
        r.energy_menu(changed)


def test_actual_higgs_encoding_retains_dark_and_outside_columns():
    from model.particle_observer_outputs import higgs_records
    data = higgs_records.sector_data(s.diag(1, 2, 3), 1)
    menu, _ = r.energy_menu(data)
    plan = r.compiled_menu_plan(menu)
    assert plan['retained_outcome_ranks'] == {
        'lower_energy': 1, 'upper_energy': 1, 'dark_energy': 305,
        'outside_field_encoding': 205}
    assert plan['native_pulse_upper_bound'] == 83894362
    assert plan['unused_pointer_labels_retained_as_zero'] == 0
    number_plan = r.compiled_menu_plan(r.scalar_number_menu(data))
    assert number_plan['unused_pointer_labels_retained_as_zero'] == 1
