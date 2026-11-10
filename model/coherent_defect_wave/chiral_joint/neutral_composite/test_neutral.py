"""Exact complete-domain instrument and same-ground neutral spectral checks."""
import sympy as s
from .model import instrument_certificate,local_instrument
from .readout import (exact_transfer_checks,star_translation_certificate,
                      complete_readout,translation_transfer)


def test_complete_quantum_scalar_car_instrument():
    result=instrument_certificate()
    assert all(row['effect_multiplicities']==[14,4,2] for row in result['branches'])
    assert all(row['all_twelve_gauge_currents_preserved'] for row in result['branches'])
    assert all(row['explicit_unitary_dilation_dimension']==40 for row in result['branches'])
    assert all(row['two_sided_unitarity_checked'] for row in result['branches'])


def test_complete_gauss_transfers_and_retained_failure_image():
    rows=exact_transfer_checks()
    assert all(row['unresolved_target_dimension']==298 for row in rows)
    assert all(row['complete_resolved_target_dimension']==298 for row in rows)
    assert all(row['failure_target_dimension']>row['source_dimension'] for row in rows)
    assert all(row['failure_has_nonzero_opposite_translation_image'] for row in rows)


def test_star_selection_and_phase_pairing():
    result=star_translation_certificate()
    assert all(row['star_anticommutes_translation_involution'] for row in result['odd_sector_checks'])
    assert not result['finite_translation_phase_is_particle_dispersion']


def test_complete_original_record_law_and_first_neutral_overlap():
    result=complete_readout()
    lo,hi=map(s.Rational,result['complete_original_probability_enclosure'])
    assert lo<=1<=hi
    count=0
    for branch in result['branches']:
        assert branch['first_energy_overlap_certifiably_positive']
        assert s.Rational(branch['retained_failure_probability_interval'][0])>0
        for phase in branch['complete_translation_branches']:
            assert sum(row['rank'] for row in phase['complete_energy_records'])==149
            count+=sum(row['rank'] for row in phase['complete_energy_records'])
    assert count==596


def test_invalid_translation_sign_rejected():
    import pytest
    with pytest.raises(ValueError):translation_transfer(sign=0)


def test_nearest_interval_bound_equals_full_complement_scan():
    from .spectral_bounds import projector_error
    from ..periodic_spectrum.readout import clusters,projector_error as inherited
    block=(0,1,-1,1);groups,intervals=clusters(block)
    selected=[groups[0],groups[-1],next(group for group in groups if len(group[0])>1)]
    for indices,_,_ in selected:assert projector_error(block,indices)==inherited(block,indices,intervals)


def test_nearest_interval_bound_rejects_hostile_spectral_menus(monkeypatch):
    import pytest
    from . import spectral_bounds as bounds
    cases=[
        ([(2,3),(0,1)],[s.Rational(5,2),s.Rational(1,2)],[([0],2,3),([1],0,1)]),
        ([(0,2),(1,3)],[s.Rational(1,2),s.Rational(5,2)],[([0],0,2),([1],1,3)]),
        ([(0,1),(2,3)],[s.Rational(1,2),s.Rational(7,2)],[([0],0,1),([1],2,3)]),
        ([(0,1),(2,3),(4,5)],[s.Rational(1,2),s.Rational(5,2),s.Rational(9,2)],[([0,2],0,5),([1],2,3)]),
        ([(0,1),(2,3)],[s.Rational(1,2),s.Rational(5,2)],[([0],0,1)])]
    for intervals,energies,groups in cases:
        bounds.projector_errors.cache_clear()
        monkeypatch.setattr(bounds,'approximation',lambda block:(None,energies,[s.Rational(1,100)]*len(energies),s.Integer(0)))
        monkeypatch.setattr(bounds,'clusters',lambda block:(groups,intervals))
        with pytest.raises(AssertionError):bounds.projector_errors(('hostile',))
    bounds.projector_errors.cache_clear()
