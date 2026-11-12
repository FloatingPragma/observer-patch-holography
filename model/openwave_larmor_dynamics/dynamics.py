"""Complete positive-kernel evolution with a quasistatic transverse self-field."""
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh
from model.openwave_intrinsic_magnetic.solver import dirac,geometry,GAMMA,BETA,SIGMA,moment,native_reduction
from model.openwave_field_observables.lattice import divergence_adjoint,curl,curl_adjoint


@lru_cache(maxsize=16)
def symbols(n,a):
    angle=2*np.pi*np.fft.fftfreq(n)
    d=np.stack(np.broadcast_arrays(*[np.expm1(1j*angle.reshape(tuple(n if j==axis else 1 for j in range(3))))/a for axis in range(3)]),axis=-1)
    return d,np.sum(abs(d)**2,axis=-1)


def transverse(current,a):
    d,lap=symbols(len(current),a);j=np.fft.fftn(current,axes=(0,1,2))
    v=np.sum(d.conj()*j,axis=-1)
    j-=d*np.divide(v,lap,out=np.zeros_like(v),where=lap>0)[...,None]
    j[0,0,0]=0
    return np.fft.ifftn(j,axes=(0,1,2)).real


def field_minimizer(current,a,kappa):
    j=transverse(current,a);_,lap=symbols(len(j),a)
    ft=np.fft.fftn(j,axes=(0,1,2))
    value=np.divide(kappa*ft/a**3,lap[...,None],out=np.zeros_like(ft),where=lap[...,None]>0)
    return np.fft.ifftn(value,axes=(0,1,2)).real


def spin_operators(n):
    zero=np.zeros((2,2))
    return tuple(np.kron(np.block([[s,zero],[zero,s]])/2,np.eye(n**3)) for s in SIGMA)


def diagonalize(field,params,external_B):
    h,phases=dirac(params['n'],params['a'],params['field_n'],field,
        params['mass'],params['trap'],params['wilson'],external_B)
    w,v=eigh(h.toarray(),check_finite=False)
    return w,v,phases


def kernel_apply(state,w,v):return v@(abs(w)*(v.conj().T@state))


def evolve(state,w,v,time,mass=0):
    return v@(np.exp(-1j*time*(abs(w)-mass))*(v.conj().T@state))


def current(state,w,v,phases,params):
    """J=-Tr[F_D(|psi><psi|) delta D], including both signed branches."""
    amplitudes=v.conj().T@state
    density=amplitudes[:,None]*amplitudes.conj()[None,:]
    factor=(w[:,None]+w[None,:])/(abs(w[:,None])+abs(w[None,:]))
    response=v@(factor*density)@v.conj().T
    n=params['n'];nn=n**3;N=params['field_n'];r=params['wilson']
    j=np.zeros((N,)*3+(3,))
    for left,right,axis,grid,phase in phases:
        block=response[np.arange(4)*nn+right][:,np.arange(4)*nn+left]
        j[grid+(axis,)]=np.real(np.sum(block.T*(GAMMA[axis]+1j*r*BETA))*phase)
    return j


def quasistatic(state,field,params,external_B,tolerance=2e-12):
    field=field.copy()
    for count in range(30):
        w,v,ph=diagonalize(field,params,external_B)
        j=current(state,w,v,ph,params)
        target=field_minimizer(j,params['a'],params['kappa'])
        residual=np.linalg.norm(target-field)
        if residual<tolerance:return field,w,v,ph,j,float(residual),count+1
        field=target
    raise ArithmeticError('quasistatic self-field iteration did not converge')


def stationary(params,axis=0):
    field=np.zeros((params['field_n'],)*3+(3,));spin=spin_operators(params['n'])[axis]
    for count in range(30):
        w,v,ph=diagonalize(field,params,0.)
        ids=np.flatnonzero(w>0)[:2];basis=v[:,ids]
        sw,su=eigh(basis.conj().T@spin@basis)
        state=basis@su[:,-1]
        j=current(state,w,v,ph,params);target=field_minimizer(j,params['a'],params['kappa'])
        residual=np.linalg.norm(target-field)
        if residual<2e-12:break
        field=target
    else:raise ArithmeticError('initial stationary spin branch did not converge')
    energy=float(np.vdot(state,kernel_apply(state,w,v)).real)
    eigen_residual=float(np.linalg.norm(kernel_apply(state,w,v)-energy*state))
    assert eigen_residual<2e-10
    return state,field,{'energy':energy,'field_residual':float(residual),'eigen_residual':eigen_residual,
        'spin':[float(np.vdot(state,s@state).real) for s in spin_operators(params['n'])],
        'moment':moment(j,params['a']).tolist()}


def midpoint_step(state,field,params,external_B,step):
    guess=field.copy()
    for count in range(30):
        w,v,ph=diagonalize(guess,params,external_B)
        mid=evolve(state,w,v,step/2,params['mass'])
        j=current(mid,w,v,ph,params)
        target=field_minimizer(j,params['a'],params['kappa'])
        residual=np.linalg.norm(target-guess)
        if residual<2e-12:
            return evolve(state,w,v,step,params['mass']),guess,float(residual),count+1
        guess=target
    raise ArithmeticError('implicit exponential midpoint did not converge')


def observe(state,field,params,B,complete=False):
    field,w,v,ph,j,resid,_=quasistatic(state,field,params,B)
    a=params['a'];k=params['kappa'];N=params['field_n'];n=params['n']
    spin=np.array([np.vdot(state,s@state).real for s in spin_operators(n)])
    b=curl(field,a);field_energy=a**3*np.sum(b*b)/(2*k)
    matter_energy=np.vdot(state,kernel_apply(state,w,v)).real
    density_rate=2*np.imag(state.conj()*kernel_apply(state,w,v)).reshape(4,-1).sum(axis=0)
    offset=(N-n)//2;div=divergence_adjoint(j,a)[offset:offset+n,offset:offset+n,offset:offset+n].ravel()
    output={'spin':spin.tolist(),'moment':moment(j,a).tolist(),
        'matter_energy':float(matter_energy),'field_energy':float(field_energy),'total_energy':float(matter_energy+field_energy),
        'continuity_residual':float(np.linalg.norm(density_rate+div)),
        'transverse_Maxwell_residual':float(np.linalg.norm(curl_adjoint(b,a)-k*transverse(j,a)/a**3)),
        'self_field_residual':resid,'norm_error':float(abs(np.vdot(state,state)-1)),
        'negative_D_branch_population':float(np.sum(abs(v[:,w<0].conj().T@state)**2)),
        'spin_measurement_probabilities':[[float(.5+s),float(.5-s)] for s in spin]}
    if complete:
        probability=abs(v.conj().T@state)**2
        output['complete_matter_energy_records']=[{'energy':float(abs(energy)),'probability':float(p),'D_sign':int(np.sign(energy))} for energy,p in zip(w,probability)]
        output['energy_record_total_probability']=float(np.sum(probability))
        output['magnetic_moment_records']=magnetic_records(state,w,v,ph,params)
    return output,field


def magnetic_records(state,w,v,phases,params):
    """Separate projective instruments of the full instantaneous moment."""
    _,positions,_=geometry(params['n'],params['a'],params['field_n'])
    nn=params['n']**3;result=[]
    factor=(w[:,None]+w[None,:])/(abs(w[:,None])+abs(w[None,:]))
    for axis in range(3):
        delta=np.zeros((4*nn,4*nn),complex)
        for left,right,direction,grid,phase in phases:
            midpoint=(positions[left]+positions[right])/2
            weight=np.cross(np.eye(3)[axis],midpoint)[direction]/2
            block=weight*phase*(-GAMMA[direction]-1j*params['wilson']*BETA)/2
            for i in range(4):
                for j in range(4):
                    delta[i*nn+left,j*nn+right]+=block[i,j]
                    delta[j*nn+right,i*nn+left]+=block[i,j].conjugate()
        observable=-v@(factor*(v.conj().T@delta@v))@v.conj().T
        energies,modes=eigh(observable,check_finite=False)
        probabilities=abs(modes.conj().T@state)**2
        result.append({'axis':axis,'records':[{'moment':float(e),'probability':float(p)} for e,p in zip(energies,probabilities)],
            'expectation':float(np.dot(energies,probabilities)),'total_probability':float(sum(probabilities))})
    return result


def trajectory(params,B,step,duration,record_stride=10):
    state,field,initial=stationary(params)
    field,*_=quasistatic(state,field,params,B)
    steps=int(round(duration/step));step=duration/steps
    rows=[];max_iter_residual=0.;maximum_iterations=0
    for index in range(steps+1):
        if index%record_stride==0 or index==steps:
            obs,field=observe(state,field,params,B,index in (0,steps));rows.append({'time':index*step,**obs})
        if index==steps:break
        state,field,resid,it=midpoint_step(state,field,params,B,step)
        max_iter_residual=max(max_iter_residual,resid);maximum_iterations=max(maximum_iterations,it)
    phase=np.unwrap([np.angle(r['spin'][0]+1j*r['spin'][1]) for r in rows]);times=np.array([r['time'] for r in rows])
    fit=np.polyfit(times,phase,1);moment_phase=np.unwrap([np.angle(-(r['moment'][0]+1j*r['moment'][1])) for r in rows])
    mfit=np.polyfit(times,moment_phase,1)
    return {'params':params,'external_B':B,'step':step,'steps':steps,'duration':duration,'initial_stationary':initial,
        'rows':rows,'spin_angular_frequency':float(fit[0]),'moment_angular_frequency':float(mfit[0]),
        'phase_line_maximum_residual':float(np.max(abs(phase-np.polyval(fit,times)))),
        'total_energy_drift':float(max(r['total_energy'] for r in rows)-min(r['total_energy'] for r in rows)),
        'maximum_midpoint_self_field_residual':max_iter_residual,'maximum_midpoint_iterations':maximum_iterations,
        'final_state':state,'final_field':field}
