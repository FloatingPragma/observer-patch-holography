"""Static two-cycle controller with a normwise Taylor enclosure.

Certificates assume IEEE binary64 elementary arithmetic with correct rounding;
no claim of a formally verified numerical backend is made.
"""
from fractions import Fraction as F
from math import factorial, ceil, sqrt, atan2
import numpy as np
from scipy.sparse import diags, csr_matrix, kron, eye
from model.openwave_collective_finite_strength.bounds import exp_bounds
SIDE=512;INTERNAL=64;SCALE=4096;J=16;DT=F(1,16);DEGREE=28;STEPS=216
REGIONS={'prepare_one':(48,64),'read_one':(160,176),'reset':(208,224),'prepare_two':(256,272),'read_two':(400,416)}
FLIGHTS=[{'corridor':(64,160),'steps':(40,48)},{'corridor':(272,400),'steps':(148,164)}]
# A declared dyadic packet, not an exact Gaussian. These integers are the inputs.
WEIGHTS=(100263,186387,332883,571410,942077,1491865,2269767,3318706,4664652,6303480,8196567,10249203,12305649,14193031,15759223,16777216,17116198,16777216,15759223,14193031,12305649,10249203,8196567,6303480,4664652,3318706,2269767,1491865,942077,571410,332883,186387,100263)

def up_dyadic(x,bits=40):return F(ceil(float(x)*2**bits)+1,2**bits)
def upper_norm(v):return float(np.linalg.norm(v))*(1+2**-28)+2**-48

def atan_bounds(x):
    x=F(x)
    if abs(x)>1:
        plo,phi=pi_bounds()
        lo,hi=atan_bounds(1/abs(x))
        out=(plo/2-hi,phi/2-lo)
        return out if x>0 else (-out[1],-out[0])
    if abs(x)>F(1,2):
        plo,phi=pi_bounds();lo,hi=atan_bounds((abs(x)-1)/(abs(x)+1))
        out=(plo/4+lo,phi/4+hi)
        return out if x>0 else (-out[1],-out[0])
    total=sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(60)),F(0))
    error=abs(x)**121/121
    return total-error,total+error

def pi_bounds():
    a,b=atan_bounds(F(1,5));c,d=atan_bounds(F(1,239))
    return 16*a-4*d,16*b-4*c

def phase_frequency_enclosure(z0,z1,error0,error1,duration):
    z=z1*z0.conjugate()
    error=up_dyadic(error0*abs(z1)+error1*abs(z0)+error0*error1+2**-44,50)
    x,y=F(float(z.real)),F(float(z.imag))
    if x<=error:return None
    corners=[(y+dy*error)/(x+dx*error) for dx in (-1,1) for dy in (-1,1)]
    lo=atan_bounds(min(corners))[0];hi=atan_bounds(max(corners))[1]
    dt=F(duration)/SCALE
    grid=2**32
    return [str(F((lo/dt*grid).__floor__(),grid)),str(F((hi/dt*grid).__ceil__(),grid))]

def propagation_constants():
    z=40*DT;eu=exp_bounds(z/4)[1]**4
    gamma=F(1,2**47);k=2*DEGREE+2
    rounding=k*gamma/(1-k*gamma)*eu*(1+F(1,2**20))
    tail=eu*z**(DEGREE+1)/factorial(DEGREE+1)
    beta=F(1,2**46)
    delta=up_dyadic(rounding+tail+DT*beta+F(1,2**900),60)
    e0=F(1,2**43)
    e=(e0+STEPS*delta)/(1-STEPS*delta)
    return delta,e0,e,{'centered_exact_H_norm_upper_exact':'40','entry_rounding_operator_norm_upper_exact':str(beta),
        'matvec_and_axpy_relative_error_budget_exact':str(gamma),'Taylor_degree':DEGREE,'step_exact':str(DT),
        'one_step_state_error_upper_exact':str(delta),'initial_state_error_upper_exact':str(e0),
        'final_state_error_upper_exact':str(e),'complete_POVM_l1_error_upper_exact':str((2+e)*e),
        'arithmetic_assumption':'IEEE binary64, correctly rounded elementary operations, no overflow; absolute underflow contribution bounded by 2^-900. No formal backend verification is asserted.'}

def exchange_matrix(bit):
    out=np.zeros((INTERNAL,INTERNAL))
    for j in range(16):
        if j&bit:out[j,48+(j^bit)]=out[48+(j^bit),j]=1
    return out

def source(mass,coupling=2,twirl=False):
    mass=F(mass)
    if not 0<mass<=4096 or coupling not in (2,3):raise ValueError('certified native mass/strength domain')
    x=np.arange(SIDE);n=(np.arange(INTERNAL)//16>=2).astype(int)
    bits={name:((np.arange(INTERNAL)&bit)>0).astype(int) for name,bit in [('preparation',8),('first_readout',4),('battery',2),('second_readout',1)]}
    active={'preparation':(x<64)|((x>=256)&(x<272)),'first_readout':x>=160,'battery':np.ones(SIDE,bool),'second_readout':x>=400}
    counts=np.broadcast_to(n,(SIDE,INTERNAL)).copy()
    for name,v in bits.items():counts+=active[name][:,None]*v[None,:]
    h=mass/SCALE
    diag=np.array([float(h*int(c)) for c in counts.ravel()]).reshape(SIDE,INTERNAL)
    pieces={'native_source':np.broadcast_to(np.array([float(h*int(v)) for v in n]),diag.shape).copy()}
    for name,v in bits.items():pieces[name]=active[name][:,None]*np.array([float(h*int(k)) for k in v])[None,:]
    kinetic=diags([-J*np.ones(SIDE-1),-J*np.ones(SIDE-1)],[-1,1],format='csr')
    H=kron(kinetic,eye(INTERNAL),format='csr')+diags(diag.ravel(),format='csr')
    exchanges=[]
    for name,bit,strength in [('prepare_one',8,F(coupling)),('read_one',4,F(2)),('reset',2,F(22,7)),('prepare_two',8,F(coupling)),('read_two',1,F(2))]:
        lo,hi=REGIONS[name];mask=(x>=lo)&(x<hi);B=csr_matrix(exchange_matrix(bit))
        op=kron(diags(mask*float(strength)),B,format='csr');H+=op
        exchanges.append((name,op,strength))
    assert max(np.asarray(abs(H).sum(axis=1)).ravel())<40 and (H-H.T).nnz==0
    clock=np.zeros(SIDE,complex)
    for j,w in enumerate(WEIGHTS,4):clock[j]=w*(1j)**j
    clock/=sqrt(sum(w*w for w in WEIGHTS))
    seed=np.zeros((INTERNAL,8 if twirl else 1),complex)
    for col,j in enumerate((0,1,4,5,8,9,12,13)):seed[j,col if twirl else 0]=1/sqrt(8)
    initial=np.kron(clock[:,None],seed)
    return H,pieces,exchanges,initial

def effect(v,lo=0,hi=SIDE,bit=None,excited=False):
    a=v.reshape(SIDE,INTERNAL,-1)[lo:hi]
    if bit is not None:a=a[:,(np.arange(INTERNAL)&bit)>0]
    if excited:a=a[:,32:]
    return float(np.sum(abs(a)**2))

def coherence(v,lo,hi):
    a=v.reshape(SIDE,INTERNAL,-1)[lo:hi]
    return np.sum(a[:,48:].conj()*a[:,:16])

def bound_probability(p,error):return [max(0.,p-error),min(1.,p+error)]

def run(mass,coupling=2,twirl=False):
    H,pieces,exchanges,v=source(mass,coupling,twirl);initial=v.copy()
    delta,e0,E,cert=propagation_constants();deltaf=float(delta)
    records={};flight_integrals=[F(0),F(0)];overwrite=F(0)
    checkpoints={0,84,108,STEPS}
    for f in FLIGHTS:checkpoints.update(f['steps'])
    allnorm=[]
    for step in range(STEPS+1):
        err=float((e0+step*delta)/(1-step*delta));prob_error=(2+err)*err+2**-28
        if step in checkpoints:
            records[step]={'state':v.copy(),'error':prob_error}
        allnorm.append(float(np.sum(abs(v)**2)))
        if step==STEPS:break
        term=v.copy();total=v.copy()
        selectors=[]
        for k,f in enumerate(FLIGHTS):
            if f['steps'][0]<=step<f['steps'][1]:
                lo,hi=f['corridor'];selectors.append((k,[lo-1,lo,hi-1,hi]))
        need_read=step>=84
        bounds={k:np.array([upper_norm(v.reshape(SIDE,INTERNAL,-1)[x]) for x in cells]) for k,cells in selectors}
        readbound=upper_norm(v.reshape(SIDE,INTERNAL,-1)[160:176]) if need_read else 0.
        for degree in range(1,DEGREE+1):
            term=(-1j*float(DT)/degree)*(H@term);total+=term
            for k,cells in selectors:bounds[k]+=np.array([upper_norm(term.reshape(SIDE,INTERNAL,-1)[x]) for x in cells])
            if need_read:readbound+=upper_norm(term.reshape(SIDE,INTERNAL,-1)[160:176])
        localerr=err+deltaf*2
        for k,b in bounds.items():
            b+=localerr
            flight_integrals[k]+=DT*2*J*up_dyadic(b[0]*b[1]+b[2]*b[3])
        if need_read:overwrite+=DT*2*up_dyadic((readbound+localerr)**2)
        v=total
    phase=[];h=float(F(mass)/SCALE)
    for k,f in enumerate(FLIGHTS):
        lo,hi=f['corridor'];s0,s1=f['steps'];r0,r1=records[s0],records[s1]
        z0,z1=coherence(r0['state'],lo,hi),coherence(r1['state'],lo,hi)
        a=max(0.,abs(z0)-r0['error']);D=float(flight_integrals[k]);duration=float((s1-s0)*DT)
        p0=bound_probability(effect(r0['state'],lo,hi),r0['error'])
        drift=abs(z1*np.exp(-1j*h*duration)-z0)
        assert drift<=D+r0['error']+r1['error']
        frequency_error=(D/(a-D)/duration*SCALE if a>D else None)
        branch_safe=(a>D and F(mass)*(s1-s0)*DT/SCALE+flight_integrals[k]/(F(a)-flight_integrals[k])<pi_bounds()[0]/2)
        if not twirl:assert branch_safe
        if twirl:assert abs(z0)<r0['error'] and abs(z1)<r1['error']
        phase.append({'clock_corridor':list(f['corridor']),'scaled_time_window_exact':[str(s0*DT),str(s1*DT)],
            'native_time_window_exact':[str(s0*DT/SCALE),str(s1*DT/SCALE)],
            'complete_corridor_success_probability_interval':p0,
            'all_clock_outside_corridor_outcomes_retained':True,
            'coherence_start_real_diagnostic':float(z0.real),'coherence_start_imag_diagnostic':float(z0.imag),
            'coherence_amplitude_lower':a,'integrated_boundary_flux_upper_exact':str(flight_integrals[k]),
            'native_frequency_error_upper':frequency_error,
            'continuous_phase_branch_without_winding_certified':bool(branch_safe),
            'native_frequency_interval':None if frequency_error is None else [float(mass)-frequency_error,float(mass)+frequency_error],
            'certified_actual_phase_frequency_interval_exact':phase_frequency_enclosure(z0,z1,r0['error'],r1['error'],(s1-s0)*DT),
            'phase_difference_frequency_diagnostic':atan2((z1/z0).imag,(z1/z0).real)/duration*SCALE if abs(z0)>r0['error'] else None,
            'demodulated_drift_diagnostic':float(drift),'visibility_lower':max(0.,2*(a-D)/p0[1]) if p0[1] else 0.})
    probability=np.sum(abs(v.reshape(SIDE,INTERNAL,-1))**2,axis=2)
    energy=[]
    kinetic=H.copy()
    for name,op,strength in exchanges:kinetic-=op
    for diag in pieces.values():kinetic-=diags(diag.ravel())
    kinetic+=2*J*eye(SIDE*INTERNAL,format='csr')
    operators=[('finite_clock_kinetic',kinetic)]+[(n,diags(d.ravel())) for n,d in pieces.items()]+[(n,o) for n,o,_ in exchanges]
    for name,op in operators:
        first=float(np.vdot(initial,op@initial).real)*SCALE;last=float(np.vdot(v,op@v).real)*SCALE
        measured=float(np.dot(op.diagonal().real,probability.ravel()))*SCALE
        energy.append({'component':name,'initial_native_energy_diagnostic':first,'final_native_energy_diagnostic':last,
            'complete_terminal_pinching_energy_diagnostic':measured})
    e_bound=((2+float(E))*float(E)+2**-28)*144*SCALE
    balance=sum(x['final_native_energy_diagnostic']-x['initial_native_energy_diagnostic'] for x in energy)
    assert abs(balance)<=e_bound
    reset=records[108];rp=bound_probability(effect(reset['state'],228,252),reset['error'])
    re=bound_probability(effect(reset['state'],228,252,excited=True),reset['error'])
    read0=records[84];p1start=effect(read0['state'],bit=4);p1end=effect(v,bit=4)
    assert abs(p1end-p1start)<=float(overwrite)+read0['error']+records[STEPS]['error']
    if twirl:assert max(abs(coherence(r['state'],0,SIDE)) for r in records.values())<1e-10
    return {'native_mass_exact':str(mass),'preparation_exchange_strength_over_scale_exact':str(coupling),
        'all_three_phase_references_twirled':twirl,'propagation_certificate':cert,
        'two_free_flight_certificates':phase,'complete_terminal_probabilities':probability.tolist(),
        'probability_axes':['512 clock positions','4 native occupations times preparation/readout1/battery/readout2 qubits'],
        'complete_probability_sum_diagnostic':float(probability.sum()),
        'complete_probability_l1_error_upper':(2+float(E))*float(E)+2**-28,
        'orthogonal_native_Fock_complement_probability_exact':'0',
        'reset_clock_corridor_success_probability_interval':rp,
        'reset_excitation_joint_probability_interval':re,
        'reset_source_excitation_conditioned_upper':min(1.,re[1]/rp[0]) if rp[0] else 1.,
        'first_readout_excited_probability_before_reset_diagnostic':p1start,
        'first_readout_excited_probability_terminal_diagnostic':p1end,
        'second_readout_excited_probability_terminal_diagnostic':effect(v,bit=1),
        'first_readout_overwrite_probability_bound':min(1.,float(overwrite)),
        'clock_not_beyond_second_readout_probability_diagnostic':effect(v,0,416),
        'energy_components':energy,'closed_total_energy_change_diagnostic':balance,
        'closed_total_energy_change_absolute_bound':e_bound,
        'complete_terminal_measurement_work_diagnostic':sum(x['complete_terminal_pinching_energy_diagnostic']-x['final_native_energy_diagnostic'] for x in energy),
        'maximum_norm_drift_diagnostic':max(abs(x-1) for x in allnorm)}
