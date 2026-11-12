import sympy as s
from . import observables as o


def test_full_native_neutral_pairing_and_current():
    data=o.source()
    assert len(data['indices'])==6
    assert o.zero(data['current']-data['u'].H)
    assert data['masses'][0]==0
    assert o.zero(data['h'].H-data['h'])


def test_complete_instrument_retains_failure_and_spectators():
    r=o.instrument()
    assert r['complete_record_count']==8
    assert r['spectator_dimension']==174
    probs=[s.sympify(x['probability']['exact'])for x in r['electron_input_original_probabilities']]
    assert s.simplify(sum(probs)-1)==0
    assert probs[6]>0
    assert probs[:3]==[0,0,0]and probs[-1]==0


def test_phase_endpoints_change_LNV_without_changing_mass_norm():
    rows=o.phase_controls()
    assert rows[0]['electron_original_success_probability']['numeric']>rows[1]['electron_original_success_probability']['numeric']
    assert all(r['beta_mass_and_oscillation_kernel_unchanged']for r in rows)
    r=o.mass_observables();low,high=[x['numeric']for x in r['m_ee_phase_envelope_GeV']]
    assert 0<low<r['m_ee_at_packet_phases_GeV']['numeric']<high
    assert r['m_beta_GeV']['numeric']>high


def test_norm_certificate_rejects_discarded_failure():
    masses=o.source()['masses']
    success_trace=2*sum((x/masses[2])**2 for x in masses)
    assert success_trace<6
    assert o.instrument()['literal_native_pulse_compilation_claimed']is False
