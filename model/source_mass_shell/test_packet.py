"""Strict receipt scope, interval and source-custody mutation controls."""
from copy import deepcopy

import pytest

from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.report()


def test_complete_receipt_replays(fresh):
    packet.validate(fresh)
    assert fresh['dispersion']['same_radial_source_exterior_witness']['full_packet_rows'] == 16*10**18
    assert fresh['physical_units']['accepted_same_generator_physical_clock_inputs'] == []


def test_numerical_admission_rejects_nonfinite_negative_and_wrong_type(monkeypatch):
    data = packet.exterior.diagnostic(mass=64)
    for key in ('norm_error', 'fiber_eigenvector_error', 'kinetic_envelope_residual',
                'texture_residual', 'actual_radial_residual'):
        for value in (float('nan'), float('inf'), -1., True):
            damaged = deepcopy(data); damaged[key] = value
            monkeypatch.setattr(packet.exterior, 'diagnostic', lambda *args, **kwargs: damaged)
            with pytest.raises(ArithmeticError, match='diagnostic'):
                packet.numerical_admission()


@pytest.mark.parametrize('mutation', ['mass', 'physical_scope', 'type', 'source', 'clock', 'halo'])
def test_receipt_rejects_changed_energy_scope_custody_and_clock(fresh, monkeypatch, mutation):
    # Freeze the independently computed reference inside these mutation tests;
    # the preceding replay exercises the actual calculation and source hashes.
    monkeypatch.setattr(packet, 'report', lambda: fresh)
    damaged = deepcopy(fresh)
    if mutation == 'mass':
        damaged['dispersion']['bulk_samples'][0]['rows'][0]['rest_energy_interval'][0] = '2'
    elif mutation == 'physical_scope':
        damaged['scope']['interacting_relativistic_continuum'] = True
    elif mutation == 'type':
        damaged['scope']['all_eight_valleys_and_sixteen_fiber_rows'] = 1
    elif mutation == 'source':
        damaged['source_bindings'].pop(next(iter(damaged['source_bindings'])))
    elif mutation == 'clock':
        damaged['physical_units']['accepted_same_generator_physical_clock_inputs'] = ['host timer']
    else:
        damaged['source_packet_control']['all_source_packet_and_halo_rows_retained'] = False
    with pytest.raises(ValueError, match='fresh calculation'):
        packet.validate(damaged)
