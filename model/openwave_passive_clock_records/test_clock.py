import numpy as np
from .clock import pair_operators,records,native_attachment

def test_complete_CAR_instrument_and_source_mass_not_double_mass():
    N,Nr,A,X,Y=pair_operators()
    assert np.array_equal(N@A-A@N,A)
    assert np.linalg.matrix_rank(A)==1
    for phase in (.17,.9,2.4):
        r=records(phase);assert abs(sum(x['probability']for x in r['complete_eight_records'])-1)<1e-14
        xx=sum(x['first_sign']*x['final_sign']*x['probability']for x in r['complete_eight_records']if x['quadrature']=='X')*2
        assert abs(xx-np.cos(phase))<1e-14
        assert abs(xx-np.cos(2*phase))>1e-3

def test_dephasing_and_differential_drive_are_actual_failure_controls():
    plain=records(1.3);driven=records(1.3,differential_drive=.2);ref=records(1.3,reference_gap=.2)
    assert max(abs(a['probability']-b['probability'])for a,b in zip(plain['complete_eight_records'],driven['complete_eight_records']))>.01
    assert driven['complete_eight_records']==ref['complete_eight_records']
    assert all(abs(x['probability']-.125)<1e-14 for x in records(1.3,dephase=True)['complete_eight_records'])

def test_native_attachment_retains_all_probe_outputs():
    r=native_attachment();assert len(r)==3
    for x in r:assert len(x['native_probe_weights_96'])==96 and abs(sum(x['native_probe_weights_96'])-1)<1e-14
    assert r[0]['clock_cycles_per_second_diagnostic']<r[1]['clock_cycles_per_second_diagnostic']<r[2]['clock_cycles_per_second_diagnostic']
