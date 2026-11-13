"""Input custody, exact declaration and complete calibration-menu regression."""
import copy
import json
import pytest
from . import source as s


def test_declared_actions_do_not_run_a_new_calibration(monkeypatch):
    def no_new_fit(*args,**kwargs):
        raise AssertionError('A declared source must not refit on execution')
    monkeypatch.setattr(s,'brentq',no_new_fit)
    s.parameters.cache_clear()
    for gradient in s.INPUTS['architecture_selection']['diagnostic_gradient_menu']:
        p=s.parameters(gradient)
        assert str(p['spin_gradient'])==gradient
        assert p['mbar']>p['delta']/2>0
    s.parameters.cache_clear()


def test_changed_inputs_recipe_or_incomplete_control_menu_are_rejected():
    record=json.loads((s.HERE/'calibration_realizations.json').read_text())
    changes=[lambda r:r.update(calibration_inputs_sha256='0'*64),
             lambda r:r['recipe_functions_sha256'].update(energy='0'*64),
             lambda r:r['actions'].pop('4'),
             lambda r:r['actions']['5/2'].update(spin_gradient='3'),
             lambda r:r['actions']['5/2'].update(U=True),
             lambda r:r['actions']['5/2'].update(spin_gradient='10/4')]
    for mutate in changes:
        bad=copy.deepcopy(record);mutate(bad)
        with pytest.raises(AssertionError):s.validate_declared_calibrations(bad)


def test_producer_pins_both_declared_data_and_exact_witnesses():
    from .__main__ import pins
    paths=pins()
    for name in('calibration_realizations.json','spectral_witnesses.json'):
        assert 'model/openwave_spin_gradient_nuclei/'+name in paths
