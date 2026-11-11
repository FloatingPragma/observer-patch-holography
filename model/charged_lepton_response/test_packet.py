from copy import deepcopy
from functools import lru_cache
import json
import math
import sys
from pathlib import Path

import pytest

from . import packet


@lru_cache(maxsize=1)
def built(): return packet.build()


def test_export_retains_calibration_ancestry_and_physical_boundary():
    value = built()
    assert value['source_selected_phase'] is None
    assert value['physical_particle_map'] == value['physical_mass_predictions'] == []
    source = value['certificates']['calibrated_sanity_test']
    assert source['phase']['tau_used_to_choose_phase'] is False
    assert source['comparison']['held_out_tau_ratio_inside_reported_one_sigma_window'] is True
    assert source['phase_calibration_is_source_derivation'] is False


@pytest.mark.parametrize('kind', ['phase', 'particle', 'mass', 'fit', 'source'])
def test_unproved_promotions_and_stale_bindings_rejected(kind):
    value = deepcopy(built())
    if kind == 'phase': value['source_selected_phase'] = '2/9'
    elif kind == 'particle': value['physical_particle_map'] = ['electron']
    elif kind == 'mass': value['physical_mass_predictions'] = [1776.969]
    elif kind == 'fit': value['calibrated_action_comparison_is_independent_prediction'] = True
    else: value['source_bindings'].pop(next(iter(value['source_bindings'])))
    with pytest.raises(ValueError): packet.validate(value)


def _same(stored, fresh):
    """Identical structure, text and flags; binary64 values agree to round-off.

    Eigen and singular-value diagnostics near 1e-15 differ between BLAS builds
    and runner CPUs, so they are compared with the round-off tolerance only.
    """
    if isinstance(stored, dict) and isinstance(fresh, dict):
        return stored.keys() == fresh.keys() and all(_same(stored[k], fresh[k]) for k in stored)
    if isinstance(stored, list) and isinstance(fresh, list):
        return len(stored) == len(fresh) and all(_same(a, b) for a, b in zip(stored, fresh))
    numbers = (int, float)
    if (isinstance(stored, numbers) and isinstance(fresh, numbers)
            and not isinstance(stored, bool) and not isinstance(fresh, bool)):
        return math.isclose(stored, fresh, rel_tol=1e-9, abs_tol=1e-12)
    return type(stored) is type(fresh) and stored == fresh


def test_generated_results_and_runtime_custody():
    value = built()
    stored = packet.OUTPUT.read_text()
    assert stored == packet.canonical(json.loads(stored))
    assert _same(json.loads(stored), json.loads(packet.canonical(value)))
    assert packet.MARKDOWN.read_text() == packet.render(value)
    assert packet.CANDIDATES.read_text() == packet.canonical(value['certificates']['candidate_particle_map'])
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.is_relative_to(packet.ROOT) and path.suffix == '.py':
                assert path.relative_to(packet.ROOT).as_posix() in value['source_bindings']


def test_cache_returns_independent_mutable_exports():
    value = packet.build(); value['certificates'].clear()
    assert len(packet.build()['certificates']) == 8


def _rendering_control(tau_inside, passing):
    # Synthetic display-only input; no physical certificate is asserted.
    ratios=[['1','1'],['2','2'],['3','3']]
    rows=[{'branch':branch,'source':{'charged_action_ratios':ratios},
           'comparison':{name:{'source_interval_inside_reference_one_sigma':branch in passing}
                         for name in ('muon_electron_ratio','tau_electron_ratio')}}
          for branch in ('forward','measured_endpoint')]
    return {'certificates':{
        'calibrated_sanity_test':{'comparison':{'source_action_ratios':ratios,
            'held_out_tau_ratio_inside_reported_one_sigma_window':tau_inside},
            'reference':{'edition':'display control','source':'https://example.invalid/reference',
                'muon_electron_ratio':{'value':'20','standard_uncertainty':'1'},
                'tau_electron_ratio':{'value':'30','standard_uncertainty':'2'}}},
        'conditional_P_dependence':{'branches':rows}}}


def test_rendering_uses_actual_comparison_and_reference_fields():
    text=packet.render(_rendering_control(False,[]))
    assert 'does not lie wholly inside' in text
    assert 'Branches wholly inside both quoted individual ratio intervals: none.' in text
    assert '| 20 ± 1 | Calibration |' in text and '| 30 ± 2 | Conditional comparison |' in text
    assert '[display control](https://example.invalid/reference)' in text
    assert 'Only the measured-alpha endpoint lies inside' not in text


def test_rendering_positive_conclusions_follow_the_reported_branches():
    text=packet.render(_rendering_control(True,['measured_endpoint']))
    assert 'The held-out tau ratio lies inside' in text
    assert 'Only the measured-alpha endpoint lies inside both' in text
    text=packet.render(_rendering_control(True,['forward','measured_endpoint']))
    assert 'intervals: forward, measured_endpoint.' in text
    assert 'Only the measured-alpha endpoint lies inside' not in text
