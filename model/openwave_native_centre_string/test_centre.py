import numpy as np
from .centre import native_quotient, spatial_labels, execute
from model.openwave_strong_step_scaling.source import chain


def test_actual_native_quotient_on_every_coordinate():
    c = native_quotient()
    rows = c['complete_native_one_generation_coordinates']
    assert len(rows) == 15
    assert all((x['electric_integer'] + x['colour_triality']) % 3 == 0 for x in rows)
    assert c['native_up_down_electric_integers'] == [2, -1]


def test_every_native_edge_term_gets_its_required_link_charge():
    for cells in (2, 3):
        d = chain(cells)
        labels = spatial_labels(d)
        shifts = labels['all_hamiltonian_link_shifts']
        assert any(3 in row or -3 in row for row in shifts)
        assert any(2 in row or -2 in row for row in shifts)
        assert any(1 in row or -1 in row for row in shifts)
        assert np.all((labels['flux'] + labels['triality']) % 3 == 0)


def test_complete_dynamical_phase_records_and_probe_energy():
    r = execute(2)
    assert r['complete_singlet_dimension'] == 120
    assert r['complete_time_records'][0]['nontrivial_string_probability_diagnostic'] == 0
    assert r['complete_time_records'][-1]['nontrivial_string_probability_diagnostic'] > .01
    from fractions import Fraction as F
    success = float(F(r['preparation']['success_exact']))
    for row in r['complete_time_records']:
        records = row['complete_joint_records']
        assert abs(sum(x['conditional_probability_diagnostic'] for x in records) - 1) < 2e-12
        assert abs(sum(x['original_input_probability_diagnostic'] for x in records) - success) < 2e-12
        assert abs(row['energy_after_complete_nonselective_probe_diagnostic'] - row['energy_before_probe_diagnostic'] - row['probe_work_diagnostic']) < 1e-12
    assert r['continuous_order_parameter_texture_winding_proved'] is False
