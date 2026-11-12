import sympy as s
from .decay import interaction,matrix_element,prepared_states,act,dot,certificate

def test_full_native_quark_and_lepton_current_polynomial():
 poly,nu,e,sign,quads=interaction();assert(nu,e,sign,quads)==(10,11,1,19)
 assert len(poly)==31464 and sum(abs(v)for v in poly.values())==41952
 assert all(poly[c,d,a,b]==v for(a,b,c,d),v in poly.items())
 assert len({i for key in poly for i in key})==144

def test_actual_charged_weights_and_three_body_amplitude():
 amplitude,words=matrix_element();assert amplitude==-1 and len(words)==1 and words[0]['coefficient']=='-2'
 reference,initial,target,_=prepared_states();assert dot(initial,target)==0
 assert dot(act(reference,8,False),act(reference,8,False))==1
 assert dot(act(reference,9,False),act(reference,9,False))==0

def test_complete_probability_and_energy_accounting():
 r=certificate();lo,hi=r['original_receiver_probability_decimal_interval'];assert 0<lo<hi<1
 assert r['free_G_zero_control_receiver_probability']=='0'
 assert r['final_projector_measurement_absolute_energy_work_upper']=='1311'
 assert r['observed_muon_lifetime_weak_scale_or_continuum_S_matrix_derived']is False
