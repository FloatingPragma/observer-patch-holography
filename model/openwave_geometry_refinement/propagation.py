"""Complete finite histories with a certified Chebyshev source propagator."""
from functools import lru_cache
from fractions import Fraction as F
import math
import numpy as np
from .source import gamma,UNIT,DEN,upper,compression,initial,indices,block,vacuum
DEGREE=64
STEP=F(1,2)

def j_interval(n):
 # J_n(16), evaluated as an exact rational series. Beyond k128 the
 # alternating terms decrease; the first omitted term encloses its tail.
 terms=[(-1)**k*F(8**(n+2*k),math.factorial(k)*math.factorial(n+k))for k in range(129)]
 middle=sum(terms,F(0));tail=F(8**(n+258),math.factorial(129)*math.factorial(n+129))
 assert F(64,130*(n+130))<1
 return middle-tail,middle+tail

@lru_cache(None)
def coefficients():
 values=[];intervals=[];error=F(0);norm=F(0)
 for n in range(DEGREE+1):
  lo,hi=j_interval(n);weight=1 if n==0 else 2;lo*=weight;hi*=weight;c=float((lo+hi)/2);values.append(c*(-1j)**n);intervals.append([str(lo),str(hi)]);error+=max(abs(F(c)-lo),abs(F(c)-hi));norm+=max(abs(lo),abs(hi))
 # |J_n(16)| <= 8^n/n! exp(64/(n+1)); exp(64/66)<3.
 tail=6*F(8**65,math.factorial(65))/(1-F(8,66));assert tail<F(1,10**30)
 return {'values':np.array(values),'coefficient_intervals':intervals,'coefficient_error':upper(error),'coefficient_l1_upper':upper(norm),'tail':upper(tail,160)}

def step_error(source):
 d=len(source['H']);integer=source['integer_H'];rows=max(sum(abs(int(x))for x in row)for row in integer);a=F(rows,32*DEN);assert a<1
 # For the real-matrix/complex-vector dot products, gamma_(2d) bounds
 # each real component. The symmetric absolute matrix has norm<=a.
 zeta=2*gamma(2*d)*a+5*UNIT;growth=F(DEGREE*(DEGREE+1),2)*zeta;assert growth<F(1,1000)
 recurrence=growth/(1-growth);c=coefficients();arithmetic=c['coefficient_l1_upper']*recurrence+gamma(4*(DEGREE+1))*c['coefficient_l1_upper']*(1+recurrence)
 return upper(STEP*source['error']+c['tail']+c['coefficient_error']*(1+recurrence)+arithmetic)

def advance(source,psi):
 B=source['H']/32;c=coefficients()['values'];previous=psi;out=c[0]*psi;present=B@psi;out=out+c[1]*present
 for n in range(2,DEGREE+1):
  following=2*(B@present)-previous;out=out+c[n]*following;previous,present=present,following
 return out

def accumulated(error,steps):return upper(steps*error*(1+error)**(steps-1)+F(1,2**40))
def probability_error(delta):return upper(2*delta+delta*delta+F(1,2**40))
def interval(value,error):return [float(np.nextafter(float(value-error),-np.inf)),float(np.nextafter(float(value+error),np.inf))]
def lift(v,cap,wave=True):
 full=np.zeros(1225 if wave else 25,complex);full[indices(cap,wave)]=v;return full

def observations(left,right=0,cap=4):
 source=compression(left,right,cap,True);vs=compression(left,right,cap,False);v=initial(cap);vac=initial(cap,False);arrival=v.copy();mask=np.tile(np.repeat(np.arange(7)>=5,7),(cap+1)**2);records=[];hits=[];history=[];states={};dc=step_error(source);dv=step_error(vs)
 for step in range(1,17):
  v=advance(source,v);vac=advance(vs,vac);arrival=advance(source,arrival);hit=np.where(mask,arrival,0);history.append(lift(hit,cap));hits.append({'event':step,'all_geometry_wave_site_probabilities':(abs(hit.reshape((cap+1)**2,49))**2/2).tolist(),'unconditional_hit_probability':float(np.vdot(hit,hit).real/2)});arrival[mask]=0
  if step in(8,16):
   states[step]=(lift(v,cap),lift(vac,cap,False));p=abs(v.reshape((cap+1)**2,49))**2;delta=accumulated(dc,step);law=probability_error(delta);mean=float(p.sum(axis=0)@np.tile(np.arange(7),7))
   records.append({'event':step,'source_parameter':step/2,'all_wave_site_geometry_number_probabilities':(p/2).tolist(),'all_vacuum_geometry_probabilities':(abs(vac)**2/2).tolist(),'omitted_master_geometry_occupations_each_probability':'0 in this compressed comparison member',
    'remaining63_native_sectors_each_probability':'0','conditional_transverse_mean':mean,'conditional_transverse_mean_interval':interval(mean,6*law),'state_vector_error_upper':str(delta),'complete_source_dimension':64*25*50})
 history.append(lift(arrival,cap));success=2*sum(x['unconditional_hit_probability']for x in hits);first=2*sum(x['event']/2*x['unconditional_hit_probability']for x in hits);timeout=float(np.vdot(arrival,arrival).real);law_exact=probability_error(accumulated(dc,16));law=float(np.nextafter(float(law_exact),np.inf));assert abs(success+timeout-1)<1e-10
 mean=first/success;meaniv=[(first-8*law)/(success+law)-2**-40,(first+8*law)/(success-law)+2**-40]
 result={'native_levels':[left,right],'nested_cap':cap,'complete_single_time_records':records,'propagator_step_error_upper':str(dc),'vacuum_step_error_upper':str(dv),
  'complete_arrival_records':hits,'all_wave_timeout_geometry_site_probabilities':(abs(arrival.reshape((cap+1)**2,49))**2/2).tolist(),'all_vacuum_timeout_geometry_probabilities':(abs(vac)**2/2).tolist(),
  'conditional_wave_hit_probability':success,'conditional_wave_timeout_probability':timeout,'conditional_mean_arrival_source_parameter':mean,'conditional_mean_arrival_source_parameter_interval':meaniv,
  'complete_history_state_error_upper':str(accumulated(dc,16)),'complete_record_TV_error_upper':str(law_exact/2),'complete_native_wave_vacuum_and_reference_complements_retained':True}
 return result,{'history':np.concatenate(history),'states':states,'vacuum_final':lift(vac,cap,False)}

def refined_records():
 rows=[];internal={}
 for left in(0,1):
  for cap in range(5):
   row,data=observations(left,cap=cap);rows.append(row);internal[left,cap]=data
 controls=[]
 for cap in range(4):
  histories=[];singles={8:[],16:[]}
  for left in(0,1):
   r=next(x for x in rows if x['native_levels']==[left,0]and x['nested_cap']==cap);master=next(x for x in rows if x['native_levels']==[left,0]and x['nested_cap']==4)
   numerical=float(np.linalg.norm(internal[left,cap]['history']-internal[left,4]['history']));err=F(r['complete_history_state_error_upper'])+F(master['complete_history_state_error_upper'])+F(1,2**35)
   vacnum=float(np.linalg.norm(internal[left,cap]['vacuum_final']-internal[left,4]['vacuum_final']));vacerr=accumulated(F(r['vacuum_step_error_upper']),16)+accumulated(F(master['vacuum_step_error_upper']),16)+F(1,2**35)
   histories.append(max(numerical+float(err),vacnum+float(vacerr)))
   for step in(8,16):
    for index,key in((0,'propagator_step_error_upper'),(1,'vacuum_step_error_upper')):
     d=float(np.linalg.norm(internal[left,cap]['states'][step][index]-internal[left,4]['states'][step][index]));bound=d+float(accumulated(F(r[key]),step)+accumulated(F(master[key]),step)+F(1,2**35));singles[step].append(bound)
  controls.append({'nested_cap':cap,'full_arrival_instrument_reference_stable_TV_upper':max(histories),'single_time_reference_stable_TV_upper':{str(k/2):max(v)for k,v in singles.items()},
   'preparation_domain':'arbitrary input/reference state on span{native00 or10} tensor geometry vacuum tensor span{field vacuum, exact packet}; all outputs retained',
   'bound_method':'complete Stinespring history vectors; native and vacuum/wave branches are orthogonal, so their maximum norm difference bounds the isometry on the entire declared input span'})
 master=[x for x in rows if x['nested_cap']==4];a,b=master[1],master[0];delay=[a['conditional_mean_arrival_source_parameter_interval'][0]-b['conditional_mean_arrival_source_parameter_interval'][1],a['conditional_mean_arrival_source_parameter_interval'][1]-b['conditional_mean_arrival_source_parameter_interval'][0]]
 return {'all_nested_complete_records':rows,'certified_compression_comparisons':controls,'master_same_source_arrival_delay_interval':delay,
  'master_time4_bending_interval':[a['complete_single_time_records'][0]['conditional_transverse_mean_interval'][0]-b['complete_single_time_records'][0]['conditional_transverse_mean_interval'][1],a['complete_single_time_records'][0]['conditional_transverse_mean_interval'][1]-b['complete_single_time_records'][0]['conditional_transverse_mean_interval'][0]]}
