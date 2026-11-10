import pytest
import sympy as s

from .model import source_certificate,port_embed
from ..model import JointObserver


def test_actual_source_bands_and_complete_spectral_weights():
    value=source_certificate()
    rows=value['blank_triplet_source_bands']
    assert [r['carrier_multiplicity_per_triplet_component'] for r in rows]==[1,3,5,3]
    assert sum(s.Rational(r['localized_port_spectral_probability']) for r in rows)==1
    assert value['three_nonzero_source_energy_ratios_exact']==['1','3*sqrt(5)/10 + 3/2','sqrt(5)/2 + 3/2']
    assert not value['same_speed_for_three_positive_branches']
    assert not value['scope']['energies_are_physical_rest_masses']


def test_preparation_domain_rejects_wrong_vector_and_normalization():
    observer=JointObserver()
    with pytest.raises(ValueError):port_embed(observer,[1],1)
    with pytest.raises(ValueError):port_embed(observer,[1]*12,0)
