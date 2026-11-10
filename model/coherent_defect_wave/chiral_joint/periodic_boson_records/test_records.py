"""Full-register operator and complete same-source spectral checks."""
import sympy as s
from .model import SETTINGS,operator_certificate,operators
from .readout import transfer_checks,complete_readout


def test_exact_loop_norm_gauge_and_complete_instruments():
    result=operator_certificate()
    assert len(result['settings'])==3
    assert [row['exact_native_current_commutators'] for row in result['settings']]==[12,24,24]
    assert all(row['Hermitian_operator_norm']=='1' for row in result['settings'])


def test_full_gauss_transfers_and_translation_selection():
    rows=transfer_checks()
    assert all(row['complete_source_Kraus_identity'] for row in rows)
    assert [row['negative_translation_success_exactly_zero'] for row in rows]==[False,True,True]


def test_complete_original_probability_elastic_and_failure_records():
    result=complete_readout();lo,hi=map(s.Rational,result['complete_original_probability_enclosure'])
    assert lo<=1<=hi
    assert [row['setting'] for row in result['settings']]==list(SETTINGS)
    for setting in result['settings']:
        assert s.Rational(setting['retained_failure_probability_interval'][0])>0
        positive=setting['lowest_certified_positive_inelastic_record'];assert positive is not None
        for phase in setting['complete_translation_branches']:
            assert sum(row['rank'] for row in phase['complete_energy_records'])==phase['target_dimension']
            if phase['translation_sign']==1:assert phase['complete_energy_records'][0]['elastic_ground_projector']
            if phase['success_map_exactly_zero']:
                assert all(row['setting_success_energy_phase_probability_interval']==['0','0'] for row in phase['complete_energy_records'])


def test_wrong_loop_normalization_is_detected():
    a,_,_=operators('doublet_loop')
    assert (2*a)**3!=2*a
