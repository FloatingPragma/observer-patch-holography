"""Complete fixed proof coordinates retain their exact native attachment."""
import json
from fractions import Fraction as F
from unittest.mock import patch
from . import certificate as c
from . import source as s


def test_all_complete_spectral_frames_and_native_trials_reproduce_certificates():
    # The source construction is independent of the proof-frame evaluation.
    c.reduced_line()
    for S, I in ((1, 1), (1, 3), (3, 1), (3, 3)):
        c.baryon_exact(S, I)
    for S in range(4):
        for I in range(4):
            for parity in (-1, 1):
                c.nuclear_exact(S, I, parity)
    retained = json.loads((s.HERE/'results.json').read_text())
    with patch.object(c, 'eigh', side_effect=AssertionError('new trial frame')):
        baryons = c.baryon_spectra()
        nuclei = c.nuclear_spectra()
        axial = c.axial_certificate()
        nuclear = c.nuclear_trial()
    assert baryons == retained['complete_baryon_multiplets']
    assert nuclei == retained['complete_nuclear_channels']
    assert axial == retained['native_axial_ground']
    assert nuclear == retained['nuclear_ground_record']
    assert sum(x['full_magnetic_dimension'] for x in baryons) == 120
    assert sum(x['full_magnetic_dimension'] for x in nuclei) == 2520
    assert F(axial['trial_variance']) > 0
