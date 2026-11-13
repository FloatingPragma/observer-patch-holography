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
