from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace
import json
import pytest
from tools import native_physical_results as summary

ROOT=Path(__file__).resolve().parents[1]


def test_projection_preserves_calibrations_failed_comparisons_and_domains():
    result=summary.collect(ROOT)
    assert len(result['receipts_sha256'])==10
    assert not result['all_sources_identified_as_one_physical_theory']
    assert {r['name'] for r in result['strange_masses'] if r['comparison_role']=='calibration'}=={'Lambda','K+'}
    assert not result['nuclear_masses']['helion']['matches']
    assert not result['nuclear_masses']['alpha']['matches']
    assert not result['independent_unit_attachment']['physical_G_strength_match']
    assert all(r['rows']==15 and r['fitted_chronometer_parameters']==0
               for r in result['conditional_expansion_comparisons'])
    text='\n'.join(summary.render(result))
    assert 'Both unit comparisons fail' in text and 'virtual response' in text
    assert 'not independent likelihoods' in text
    # Measured input ancestry is projected from its source, not a second literal.
    result['conditional_expansion_comparisons'][0]['anchors']['H0_km_s_Mpc']='68.125'
    result['conditional_expansion_comparisons'][0]['anchors']['Omega_m']='0.3125'
    changed='\n'.join(summary.render(result))
    assert 'H0=68.125 and Omega_m=0.3125' in changed
    assert 'H0=67.4' not in changed
    for name in ('neutron','proton'):
        line=next(x for x in text.splitlines() if x.startswith(f'| {name} band |'))
        display=list(map(F,line.split('[',1)[1].split(']',1)[0].split(', ')))
        source=list(map(F,result['bulk_baryons'][name+'_GeV_interval']))
        assert display[0]<=source[0]<=source[1]<=display[1]
    gravity=result['transverse_photon_and_clocks']
    for label,centre_key,error_key in (
        ('mean arrival parameter','conditional_arrival_delay','arrival_delay_error_upper'),
        ('transverse position / cell','transverse_position_difference','position_difference_error_upper'),
        ('Fourier direction / rad','separate_Fourier_direction_deflection','direction_error_upper'),
    ):
        line=next(x for x in text.splitlines() if x.startswith(f'| {label} |'))
        centre,error=map(F,line.split('|')[2:4])
        source_centre,source_error=F(gravity[centre_key]),F(gravity[error_key])
        assert centre-error<=source_centre-source_error
        assert source_centre+source_error<=centre+error


def test_projected_errors_round_outwards():
    for x in (0.0005076474164068462,F(1,3),F(-1,3)):
        assert F(summary.decimal(x))<=F(x)<=F(summary.decimal(x,upper=True))


@pytest.mark.parametrize('value,error',[
    (F(1,3),F(1,10**14)),(-F(1,3),F(1,10**14)),
    (F(1,10**200),F(1,10**250)),(-F(1,10**200),F(1,10**250)),
    (F(1,10**400),F(0)),(-F(1,10**400),F(0)),(F(0),F(0)),
])
def test_printed_centre_error_contains_exact_signed_and_tiny_intervals(value,error):
    centre,radius=map(F,summary.difference_with_error(value,error))
    assert centre-radius<=value-error<=value+error<=centre+radius
    with pytest.raises(ValueError,match='nonnegative'):
        summary.difference_with_error(value,-F(1,10**200))


def test_incomplete_custody_inventory_fails_before_projection(tmp_path,monkeypatch):
    monkeypatch.setattr(summary,'PACKAGES',('probe',))
    monkeypatch.setattr(summary.importlib,'import_module',lambda _:SimpleNamespace(pins=lambda:{'a':'one','b':'two'}))
    p=tmp_path/'model/openwave_probe/results.json';p.parent.mkdir(parents=True)
    p.write_text(json.dumps({'source_sha256':{'a':'one'},'verdict':'PASS'}))
    with pytest.raises(ValueError,match='inventory changed'):
        summary.collect(tmp_path)
