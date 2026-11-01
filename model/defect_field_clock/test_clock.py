from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import clock


def probabilities(phase):
    c = -2/3
    return {(x,y): (1+(-1)**x*c*math.cos(phase))*(1+(-1)**y*c*math.sin(phase))/4
            for x in (0,1) for y in (0,1)}


@pytest.mark.parametrize('currents', [(), (-1,), (0,), (1,), (-1,0,1)])
def test_complete_source_clock_reconstructs_common_action_time(currents):
    result = clock.history(currents, F(1))
    expected = sum((1+F(j,4) for j in currents), F())
    interval = list(map(F, clock.infer_time(result['probabilities'])['effective_time_GeV_inverse_interval']))
    assert interval[0] <= expected <= interval[1]
    assert float(interval[1]-interval[0]) < 2e-6
    assert abs(sum(result['probabilities'].values())-1) < 1e-12
    assert result['source_full_input_columns'] == result['retained_purifier_columns'] == 4
    assert result['source_cost']['fresh_pointer_blanks'] == 2
    assert result['source_cost']['maximum_work_qubits'] == 4
    assert result['source_cost']['total_clock_qubits_including_purifiers_and_fresh_pointers'] == 6
    assert result['clock_resets'] == 0 and result['preparations'] == 1


def test_clock_is_complete_on_unknown_entangled_reference():
    _, words, _ = clock.advance(np.eye(4), 1, F(2))
    U = clock.compiler.execute(words,2)
    outputs, _, _ = clock.read(U)
    effect = sum((K.conj().T@K for K in outputs.values()), np.zeros((4,4),complex))
    assert np.linalg.norm(effect-np.eye(4)) < 1e-12
    # A fixed full-rank complex reference does not share the calibration basis.
    reference = np.arange(1,13).reshape(4,3)+1j*np.arange(12,0,-1).reshape(4,3)
    raw = sum(np.linalg.norm(K@reference)**2 for K in outputs.values())
    assert abs(raw-np.linalg.norm(reference)**2) < 1e-9


@pytest.mark.parametrize('phase', [0., .001, .3, 1., 1.5])
def test_exact_inverse_enclosure_including_phase_window_endpoints(phase):
    result = clock.infer_time(probabilities(phase))
    lo,hi = map(F,result['phase_interval'])
    assert float(lo) <= phase <= float(hi)
    assert not result['inverse_uses_matter_energies_or_program_duration']


@pytest.mark.parametrize('value', [0, -1, True, float('nan'), float('inf'), '1'])
def test_invalid_clock_calibration_rejected(value):
    with pytest.raises(ValueError): clock.calibration(value)


def test_unphysical_or_inconsistent_record_rejected():
    with pytest.raises(ValueError): clock.infer_time({(0,0):1})
    with pytest.raises(ValueError): clock.infer_time({k:2*v for k,v in probabilities(.3).items()})
    p=probabilities(.3)
    # Flip the X outcome, preserving the same Y marginal and total mass.
    p={(x,y):p[1-x,y] for x,y in p}
    with pytest.raises(ValueError,match='incompatible'): clock.infer_time(p)
    with pytest.raises(ValueError,match='phase window'): clock.history((1,), F(8))


def test_source_wrong_sign_mutation_is_not_hidden_by_inverse(monkeypatch):
    original=clock.advance
    def wrong(state,current,duration,energy=clock.REFERENCE_ENERGY):
        return original(state,current,duration,2*energy)
    monkeypatch.setattr(clock,'advance',wrong)
    with pytest.raises(ArithmeticError): clock.history((1,),F(1))


@pytest.mark.parametrize('gap,time', [(F(1,10),F(1)),(F(1),F(1)),(F(2),F(1,2)),(F(0),F(1))])
def test_one_energy_scale_and_record_only_gap_inference(gap,time):
    result=clock.infer_energy_gap((1+math.cos(float(gap*time)))/2,
                                 probabilities(float(clock.REFERENCE_ENERGY*time)))
    lo,hi=map(F,result['energy_gap_GeV_interval'])
    assert lo <= gap <= hi
    assert not result['source_energies_or_program_duration_used_by_inverse']
    assert not result['species_specific_scale_used']


def test_clock_energy_calibration_cannot_be_silently_changed():
    p=probabilities(.25)
    a=list(map(F,clock.infer_time(p)['effective_time_GeV_inverse_interval']))
    b=list(map(F,clock.infer_time(p,energy=F(1,2))['effective_time_GeV_inverse_interval']))
    assert a == [2*x for x in b]
    assert clock.calibration(F(1,2))['reference_energy_GeV']=='1/2'


def test_shot_error_enclosure_and_si_units():
    a=clock.sampling_envelope(10**7)
    b=clock.sampling_envelope(2*10**7)
    assert F(b['simultaneous_failure_probability_upper']) < F(a['simultaneous_failure_probability_upper']) < F(1,10**6)
    assert F(a['clock_binary_marginal_error'])==F(1,500)
    assert a['observed_shots_in_this_receipt']==0
    for n in (0,-1,True):
        with pytest.raises(ValueError):clock.sampling_envelope(n)
    units=clock.si_calibration()
    lo,hi=map(F,units['hbar_GeV_s_interval'])
    # This broad independent reference distinguishes hbar from h by 2*pi.
    assert F('6.582119569e-25') < lo < hi < F('6.582119570e-25')
    assert hi-lo < F(1,10**55)
