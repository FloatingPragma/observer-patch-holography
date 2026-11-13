from fractions import Fraction as F
import pytest
from .blocking import run,geometry,coefficient,pairwise_haar,graph


def test_actual_spatial_census_and_all_six_exact_haar_coefficients():
    r=geometry();assert(r['fine_vertices'],r['fine_links'],r['fine_faces'])==(27,54,36)
    assert len(r['all_six_tiling_integrals'])==6
    assert all(x['coarse_character_four_fine_Wilson_Haar_integral_exact']=='1/128'
               for x in r['all_six_tiling_integrals'])
    assert r['central_link_flip_fixes_coarse_links_and_reverses_all_fine_faces']


def test_positive_induced_coupling_and_retained_inconclusive_control():
    r=run();rows=r['declared_magnetic_steps']
    assert [x['positive_coarse_Wilson_coefficient_certified']for x in rows]==[False,True,True]
    for row in rows[1:]:assert F(row['coarse_Wilson_coefficient_interval_exact'][0])>0
    assert not r['scope']['autonomous_full_Hamiltonian_RG_derived']
    assert not r['physical_clock']['matching_electric_term_alone_matches_full_source']


def test_haar_and_domain_controls():
    _,edges,faces=graph(3)
    assert pairwise_haar((faces[0],),edges)==0
    assert pairwise_haar((faces[0],faces[0]),edges)==1
    with pytest.raises(ValueError):coefficient(0)
    with pytest.raises(ValueError):coefficient(F(1,100))
