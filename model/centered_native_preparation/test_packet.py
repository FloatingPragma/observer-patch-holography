from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
from . import packet,bounds,tree,target


@pytest.fixture(scope='module')
def retained():return packet.read_pair(packet.OUTPUT)


def test_complete_native_data_matches_certified_centered_spectrum(retained):
    value,arrays=retained
    assert value['execution']['complete_output_error_diagnostic']<float(F(bounds.certificate()['complete_state_and_total_variation_error']))
    assert value['execution']['total_variation_diagnostic']<.03
    assert np.allclose(arrays['pointer_probabilities'].sum(axis=0),[.5,.5],atol=1e-8,rtol=0)
    assert value['full_source_programs']['tree_failure_weights']['native_pulses']==262137
    assert target.summary()['number']==8
    assert F(bounds.preparation_certificate()['complete_preparation_error'])<F(1,10**9)


@pytest.mark.parametrize('key',['physical_species_or_masses_identified','all228_lines_resolved_by_native_bins',
                                'continuing_bare_field_repair_included','interacting_continuum_or_Einstein_requirements_completed'])
def test_scope_promotions_rejected(retained,key):
    value,arrays=retained;changed=deepcopy(value);changed['scope'][key]=True
    with pytest.raises(ValueError):packet.validate(changed,arrays)


def test_bool_integer_and_resource_type_aliases_rejected(retained):
    value,arrays=retained
    for section,key,new in [('scope','actual_native_preparation_executed',1),('scope','physical_species_or_masses_identified',0)]:
        changed=deepcopy(value);changed[section][key]=new
        with pytest.raises(ValueError):packet.validate(changed,arrays)
    changed=deepcopy(value);changed['execution']['clock']['classical_charge_columns']=2.
    with pytest.raises(ValueError):packet.validate(changed,arrays)


def test_missing_rows_failed_probability_and_angle_changes_rejected(retained):
    value,arrays=retained
    changed=dict(arrays);changed['clock_success']=arrays['clock_success'][:-1]
    with pytest.raises(ValueError):packet.validate(value,changed)
    changed=dict(arrays);changed['tree_failure_weights']=arrays['tree_failure_weights'].copy();changed['tree_failure_weights'][0,0]=1e-3
    with pytest.raises(ValueError):packet.validate(value,changed)
    changed=dict(arrays);changed['literal_Gray_coefficients']=arrays['literal_Gray_coefficients'].copy();changed['literal_Gray_coefficients'][0]+=.001
    with pytest.raises(ValueError):packet.validate(value,changed)


def test_current_custody_cannot_be_fabricated(retained):
    value,arrays=retained;changed=deepcopy(value)
    changed['prior_custody'][packet.PRIORS[0]]['historical_expensive_experiment_reexecuted']=True
    with pytest.raises(ValueError):packet.validate(changed,arrays)
