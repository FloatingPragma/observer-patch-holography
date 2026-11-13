"""Full native spectra and preparations survive without fresh proof frames."""
from copy import deepcopy
from fractions import Fraction as F
import json
from unittest.mock import patch
import pytest
from . import certificate as c
from . import source as s
from model.openwave_larmor_dynamics.certificate import equivalent


def test_complete_native_certificate_recomputed_from_fixed_proof_coordinates():
    retained=json.loads((s.HERE/'results.json').read_text())
    # All source construction and diagnostic energy calculations remain live.
    # The certificate itself cannot ask an eigensolver for new proof columns.
    with patch.object(c,'eigh',side_effect=AssertionError('fresh proof frame')):
        result=c.certificate()
    retained.pop('source_sha256')
    assert equivalent(retained,result)
    for key in ('parameters_rational','three_hadron_mass_calibrations',
                'native_axial_ground','complete_baryon_multiplets',
                'complete_meson_charge_spectra','complete_nuclear_spectrum',
                'selected_nuclear_record','all_fragment_comparisons',
                'unassigned_bound_native_channels'):
        assert result[key]==retained[key]
    assert result['complete_nuclear_spectrum']['complete_singlet_dimension']==2520
    assert F(result['selected_nuclear_record']['trial_residual_squared'])>0


def test_zero_native_trial_is_rejected_before_normalization():
    data=deepcopy(c.WITNESSES)
    data['axial']['trial_basis_integer_coefficients']=['0']*len(data['axial']['trial_basis_integer_coefficients'])
    c.axial_ground.cache_clear()
    try:
        with patch.object(c,'WITNESSES',data),pytest.raises(AssertionError):
            c.axial_ground()
    finally:
        c.axial_ground.cache_clear()


def test_nuclear_proof_coordinates_do_not_use_floating_pivot_selection():
    from . import fixed_nuclear_basis as fixed
    retained = json.loads((s.HERE/'results.json').read_text())
    fixed.channel.cache_clear()
    c.nuclear_exact.cache_clear()
    c.nuclear_spectra.cache_clear()
    c.deuteron_trial.cache_clear()
    # The complete native projector range and its ordered basis are rebuilt
    # even if the inherited floating QR chooser is unavailable or changed.
    with patch.object(c.nuclear, 'channel', side_effect=AssertionError('floating channel basis')):
        spectra = c.nuclear_spectra()
        trial = c.deuteron_trial()
    assert spectra == retained['complete_nuclear_spectrum']
    assert trial == retained['selected_nuclear_record']
    assert len(spectra['channels']) == 31
    assert spectra['complete_fixed_Ms0_Mi0_rank'] == 240
    assert spectra['complete_singlet_dimension'] == 2520
    assert fixed.channel(3, 3, 1) == {'rank': 0}


@pytest.mark.parametrize('bad', ('duplicate', 'out_of_range', 'bool', 'missing'))
def test_corrupt_native_projector_column_indices_are_rejected(bad):
    from . import fixed_nuclear_basis as fixed
    record = deepcopy(fixed.WITNESSES['nuclear_1_0_1'])
    indices = record['native_projector_column_indices']
    if bad == 'duplicate':
        indices[0] = indices[1]
    elif bad == 'out_of_range':
        indices[0] = fixed.projector(1, 0, 1)[0].shape[1]
    elif bad == 'bool':
        indices[0] = True
    else:
        indices.pop()
    with pytest.raises(AssertionError):
        fixed.channel_from_record(1, 0, 1, record)
