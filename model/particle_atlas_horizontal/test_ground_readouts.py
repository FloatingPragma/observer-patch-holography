"""Prevent state mixing and interpretation changes in the ground panel."""
from copy import deepcopy
import pytest

from . import build, ground_readouts


@pytest.fixture(scope='module')
def report():
    return build.load()['observer_outputs']


def test_shared_preparation_and_complete_signed_records(report):
    ground_readouts.validate_inputs(report)


@pytest.mark.parametrize('mutation',[
    'temperature', 'ground', 'missing_branch', 'reference', 'matched_speed',
    'physical_mass', 'negative_phase_bin', 'continuum',
])
def test_display_rejects_incompatible_source_or_promoted_claim(report, mutation):
    value=deepcopy(report)
    native=value['spatial_low_energy_preparation']['native208']
    records=value['ground_energy_spectral_records']
    clock=records['cooled_source_records']['complete_clock_apparatus']
    gauge=value['gauge_matter_characteristics']
    if mutation=='temperature': native['beta']='1/64'
    elif mutation=='ground': native['ground_energy_interval']=['0','1']
    elif mutation=='missing_branch': clock['records'].pop()
    elif mutation=='reference': clock['all_quantum_references_retained']=False
    elif mutation=='matched_speed': gauge['declared_common_clock_example']['electric_scale']='1'
    elif mutation=='physical_mass': records['observed_particle_identification']=True
    elif mutation=='negative_phase_bin': clock['negative_bin_is_negative_energy_measurement']=True
    elif mutation=='continuum': value['spatial_low_energy_preparation']['thermodynamic_scope']['a_to_zero_limit_proved']=True
    with pytest.raises(ValueError): ground_readouts.validate_inputs(value)
