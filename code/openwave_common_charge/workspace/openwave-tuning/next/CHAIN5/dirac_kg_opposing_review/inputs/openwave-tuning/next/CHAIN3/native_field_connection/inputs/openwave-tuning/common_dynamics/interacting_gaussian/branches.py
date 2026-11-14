"""Routed Dirac matter coupled to one complete quadratic Gaussian field.

The finite branch representation is an analytic realization of the declared
continuous field. It does not assert a native finite-register implementation.
"""
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from pathlib import Path
import importlib.util
import math
import sys

import numpy as np
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PROVIDER=HERE.parent/'quadratic_field/source.py'
spec=importlib.util.spec_from_file_location('dyn16_quadratic_source',PROVIDER)
field=importlib.util.module_from_spec(spec);sys.modules[spec.name]=field;spec.loader.exec_module(field)
matter=field.matter
MODEL=field.MODEL_ID
A=float(field.ELECTRIC)
V_H=.5


@lru_cache(None)
def data():
    g=field.geometry();ref=field.reference(harmonic_electric_variance=V_H)
    frame=np.column_stack((g['U'],g['harmonic_frame']))
    omega=ref['frequencies']
    return dict(g=g,frame=frame,omega=omega,
                vq=np.r_[A/(2*omega),np.full(3,1/(4*V_H))],
                vp=np.r_[omega/(2*A),np.full(3,V_H)],
                E0=float(omega.sum()/2))


def covariance(time):
    d=data();vq=d['vq'].copy();cross=np.zeros(193)
    vq[190:]+=A*A*time*time*V_H;cross[190:]=A*time*V_H
    return vq,d['vp'],cross


def coulomb(cells):
    return _coulomb(tuple(int(x) for x in cells))


@lru_cache(None)
def _coulomb(cells):
    q=np.zeros(96);q[0]=len(cells)
    for cell in cells:q[12*cell]-=1
    g=data()['g'];p=g['B'].T@g['G']@q
    return q,p,A*float(q@g['G']@q)/2


@dataclass
class State:
    cells: np.ndarray
    amps: np.ndarray
    q: np.ndarray
    p: np.ndarray
    time: float=0.
    macros: int=0

    def copy(self):
        return State(self.cells.copy(),self.amps.copy(),self.q.copy(),self.p.copy(),self.time,self.macros)


def seed(cells=(0,),spinors=None,plane=False):
    if len(cells) not in (1,2) or len(set(cells))!=len(cells) or any(type(c)is not int or not 0<=c<8 for c in cells):
        raise ValueError('one or two distinct physical root cells required')
    if plane and len(cells)!=1:raise ValueError('plane reference is one particle')
    if spinors is None:spinors=[np.array([1,0,0,0],complex) for _ in cells]
    if len(spinors)!=len(cells):raise ValueError('one spinor per occupied cell required')
    amp=np.ones(1,complex)
    for spinor in spinors:
        spinor=np.asarray(spinor,complex)
        if spinor.shape!=(4,) or not np.all(np.isfinite(spinor)) or abs(np.vdot(spinor,spinor)-1)>1e-12:
            raise ValueError('normalized spinor required')
        amp=np.kron(amp,spinor)
    if len(cells)==2 and cells[0]>cells[1]:
        cells=tuple(reversed(cells));amp=-amp.reshape(4,4).T.reshape(16)
    positions=np.array([[cell] for cell in range(8)] if plane else [cells],int)
    amplitudes=np.tile(amp,(len(positions),1))/math.sqrt(len(positions))
    return State(positions,amplitudes,np.zeros((len(positions),193)),np.zeros((len(positions),193)))


def axis_step(state,axis,*,coupling=True,omit_bch=False):
    if type(axis)is not int or axis not in (0,1,2):raise ValueError('physical axis required')
    particles=state.cells.shape[1];d=data();rows=[]
    projectors={s:(np.eye(4)+s*matter.ALPHA[axis])/2 for s in (-1,1)}
    # Every branch is retained, including exact-zero spin amplitudes.
    for signs in product((-1,1),repeat=particles):
        O=projectors[signs[0]]
        for sign in signs[1:]:O=np.kron(O,projectors[sign])
        amps=state.amps@O.T
        cells=np.empty_like(state.cells);j=np.zeros_like(state.p)
        for k,old in enumerate(state.cells):
            J=np.zeros(288)
            for slot,(cell,sign) in enumerate(zip(old,signs)):
                other,path,current=matter.route(int(cell),axis,sign)
                cells[k,slot]=other;J+=current
            j[k]=d['frame'].T@J if coupling else 0.
            q0,_,_=coulomb(old);q1,_,_=coulomb(cells[k])
            assert np.linalg.norm(q0+d['g']['B']@J-q1)<1e-12
            if particles==2 and cells[k,0]>cells[k,1]:
                cells[k]=cells[k,::-1];amps[k]=-amps[k].reshape(4,4).T.reshape(16)
        phase=np.ones(len(amps),complex) if omit_bch else np.exp(.5j*np.sum(j*state.q,axis=1))
        rows.append((cells,amps*phase[:,None],state.q.copy(),state.p+j))
    return State(*(np.concatenate([row[i] for row in rows]) for i in range(4)),state.time,state.macros)


def field_step(state,time=1.,*,electric_charge=True,free=True):
    if not math.isfinite(time) or time<0:raise ValueError('finite nonnegative field time required')
    out=state.copy();d=data()
    # exp(-i H0 t) W(z) = W(S(t)z) exp(-i H0 t), with no
    # unrecorded displacement phase in this Weyl convention.
    if free:
        c=np.cos(d['omega']*time);s=np.sin(d['omega']*time)
        out.q[:,:190]=state.q[:,:190]*c+A*state.p[:,:190]*s/d['omega']
        out.p[:,:190]=state.p[:,:190]*c-state.q[:,:190]*s*d['omega']/A
        out.q[:,190:]=state.q[:,190:]+A*time*state.p[:,190:]
        out.amps*=np.exp(-1j*d['E0']*time)
        out.time+=time
    if electric_charge:
        energies=np.array([coulomb(cells)[2] for cells in state.cells])
        out.amps*=np.exp(-1j*energies*time)[:,None]
    return out


def step(state,h=.125,mass=1.,*,coupling=True,free=True,electric_charge=True,omit_bch=False):
    if not math.isfinite(h) or h<=0 or not math.isfinite(mass):raise ValueError('finite mass and positive spacing required')
    if len(state.amps)*8**state.cells.shape[1]>8192:raise ValueError('exact branch budget exceeded; no pruning performed')
    out=state.copy();M=expm(-1j*h*mass*matter.BETA);O=M
    for _ in range(state.cells.shape[1]-1):O=np.kron(O,M)
    out.amps=out.amps@O.T
    for axis in range(3):out=axis_step(out,axis,coupling=coupling,omit_bch=omit_bch)
    out=field_step(out,1.,electric_charge=electric_charge,free=free)
    out.macros=state.macros+1
    return out


def gram(q,p,time):
    """G_ij=<W(q_i,p_i)phi_t|W(q_j,p_j)phi_t>."""
    return gram_cross(q,p,q,p,time)


def gram_cross(q,p,r,s,time):
    vq,vp,cross=covariance(time)
    left=np.sum(p*p*vq+q*q*vp-2*p*q*cross,axis=1)
    right=np.sum(s*s*vq+r*r*vp-2*s*r*cross,axis=1)
    bil=(p*vq)@s.T+(q*vp)@r.T-(p*cross)@r.T-(q*cross)@s.T
    variance=left[:,None]+right[None,:]-2*bil
    phase=(q@s.T-p@r.T)/2
    return np.exp(-variance/2+1j*phase)


def records(state,*,electric_moments=True):
    d=data();dim=state.amps.shape[1];particles=state.cells.shape[1]
    rho=np.zeros((dim,dim),complex);position=[];electric_mean=np.zeros(288,complex);electric_energy=0j
    currents=[];norm=0.
    vq,vp,cross=covariance(state.time)
    for cells in sorted(set(map(tuple,state.cells.tolist()))):
        ids=np.flatnonzero(np.all(state.cells==cells,axis=1));Q=state.q[ids];P=state.p[ids];amp=state.amps[ids]
        local=np.zeros((dim,dim),complex)
        mode_mean=np.zeros(193,complex);second_sum=0j
        left=P/2-1j*cross*P+1j*vp*Q
        right=P/2+1j*cross*P-1j*vp*Q
        for start in range(0,len(amp),128):
            sl=slice(start,start+128)
            G=gram_cross(Q[sl],P[sl],Q,P,state.time)
            local+=(amp[sl].conj().T@G@amp).T
            if electric_moments:
                W=(amp[sl].conj()@amp.T)*G
                rows=W.sum(axis=1);cols=W.sum(axis=0)
                mode_mean+=rows@left[sl]+cols@right
                second_sum+=(rows@np.sum(left[sl]**2,axis=1)
                    +cols@np.sum(right**2,axis=1)
                    +2*np.sum(left[sl]*(W@right))+W.sum()*sum(vp))
        weight=float(np.trace(local).real);norm+=weight;rho+=local
        position.append({'occupied_root_cells':list(cells),'probability':weight})
        # The last axis' spin rail is unchanged by the following field step.
        for signs in product((-1,1),repeat=particles):
            O=(np.eye(4)+signs[0]*matter.ALPHA[2])/2
            for sign in signs[1:]:O=np.kron(O,(np.eye(4)+sign*matter.ALPHA[2])/2)
            prob=float(np.trace(O@local).real)
            currents.append({'occupied_root_cells':list(cells),'last_axis_signs':list(signs),'probability':prob})
        if electric_moments:
            _,longitudinal,_=coulomb(cells)
            electric_mean+=weight*longitudinal+d['frame']@mode_mean
            electric_energy+=A*(weight*float(longitudinal@longitudinal)+second_sum)/2
    assert abs(norm-1)<5e-10
    assert np.max(abs(rho-rho.conj().T))<1e-10
    assert min(np.linalg.eigvalsh(rho))>-1e-10
    spin=[]
    for slot in range(particles):
        values=[]
        for S in matter.SIGMA:
            O=S/2
            if particles==2:O=np.kron(O,np.eye(4)) if slot==0 else np.kron(np.eye(4),O)
            values.append(float(np.trace(O@rho).real))
        spin.append(values)
    loads=np.zeros(96)
    for row in position:
        for cell in row['occupied_root_cells']:loads[12*cell]+=row['probability']
    Gsource=np.zeros(96);Gsource[0]=particles
    gauss=float(np.linalg.norm(loads+d['g']['B']@electric_mean-Gsource)) if electric_moments else None
    if electric_moments:assert gauss<2e-9 and abs(electric_energy.imag)<2e-9
    return {'norm_squared':norm,'branches':len(state.amps),'field_cell_time':state.time,'macros':state.macros,
            'positions':position,'spin_expectation':spin,'all_spin_basis_probabilities':np.diag(rho).real.tolist(),
            'last_axis_current_records':currents,'electric_mean_real':electric_mean.real.tolist() if electric_moments else None,
            'electric_energy':float(electric_energy.real) if electric_moments else None,'Gauss_mean_residual':gauss,
            'spin_density_real':rho.real.tolist(),'spin_density_imag':rho.imag.tolist(),
            'discarded_probability':0.,'success_normalized':False}


def state_moments(state):
    """Uniform free-interval bounds for all faces and electric coordinates.

    These are numerical evaluations of the stated exact Gaussian formulas.
    The theorem supplies the bound; interval certification is separate.
    """
    d=data();omega=d['omega'];g=d['g'];frame=d['frame']
    lam=np.sum(state.q[:,:190]**2*omega/A+state.p[:,:190]**2*A/omega,axis=1)/2
    W=float(np.sum(np.linalg.norm(state.amps,axis=1)));lm=float(max(lam))
    f=frame[:,:190].T@g['A'];vf=np.sum(f*f*(A/(2*omega))[:,None],axis=0)
    polynomial=256*lm**4+1792*lm**3+3360*lm**2+1680*lm+105
    face_norm_upper=W*vf*vf*math.sqrt(polynomial)
    edge_vp=np.sum(frame*frame*d['vp'],axis=1)
    co_vp=np.sum(g['U']**2*(omega/(2*A)),axis=1)
    mean_bound=np.zeros(288)
    for cells,q,p,l in zip(state.cells,state.q,state.p,lam):
        _,pl,_=coulomb(cells)
        base=pl+g['harmonic_frame']@p[190:]
        mean_bound=np.maximum(mean_bound,abs(base)+2*np.sqrt(co_vp*l))
    return {'branch_l1_weight':W,'all_coexact_occupation_means':lam.tolist(),
            'max_coexact_occupation_mean':lm,'all208_face_variances':vf.tolist(),
            'all208_face_fourth_power_norm_upper':face_norm_upper.tolist(),
            'sum_face_fourth_power_norm_upper':float(sum(face_norm_upper)),
            'all288_electric_variances':edge_vp.tolist(),'all288_uniform_electric_mean_absolute_upper':mean_bound.tolist(),
            'max_electric_variance':float(max(edge_vp)),'max_electric_mean_absolute_upper':float(max(mean_bound)),
            'harmonic_electric_variance':V_H,'harmonic_theta_variance':1/(4*V_H)+A*A*state.time**2*V_H,
            'harmonic_cross_covariance':A*state.time*V_H,
            'scope':'All modes and all branches; free-interval face/electric bounds. Values use floating spectral diagnostics, not outward-certified eigenvectors.'}


def photon_tail_bound(occupation_upper,branch_l1,cap):
    if type(cap)is not int or cap<0 or occupation_upper<0 or branch_l1<0:raise ValueError('nonnegative state domain required')
    if cap+2<=occupation_upper:return 1.
    if occupation_upper==0:return 0.
    logtail=(cap+1)*math.log(occupation_upper)-math.lgamma(cap+2)-math.log1p(-occupation_upper/(cap+2))
    return min(1.,branch_l1*math.exp(logtail/2))


def electric_characteristic(state,edge_test,spin_operator=None):
    """Complete joint spin/electric characteristic function, no sampling.

    With a spin projector O, this is Tr[O exp(i l dot E) rho]. Every
    physical electric coordinate commutes with every spin readout here.
    """
    edge_test=np.asarray(edge_test,float)
    if edge_test.shape!=(288,) or not np.all(np.isfinite(edge_test)):
        raise ValueError('finite288edge Fourier test vector required')
    d=data();l=d['frame'].T@edge_test;vq,vp,cross=covariance(state.time)
    O=np.eye(state.amps.shape[1]) if spin_operator is None else np.asarray(spin_operator,complex)
    if O.shape!=(state.amps.shape[1],)*2:raise ValueError('spin operator dimension mismatch')
    result=0j
    for cells in sorted(set(map(tuple,state.cells.tolist()))):
        ids=np.flatnonzero(np.all(state.cells==cells,axis=1));q=state.q[ids];p=state.p[ids];amp=state.amps[ids]
        left=(p/2-1j*cross*p+1j*vp*q)@l
        right=(p/2+1j*cross*p-1j*vp*q)@l
        _,pl,_=coulomb(cells)
        scalar=np.exp(1j*float(edge_test@pl)-float(vp@(l*l))/2)
        for start in range(0,len(amp),128):
            sl=slice(start,start+128);G=gram_cross(q[sl],p[sl],q,p,state.time)
            W=(amp[sl].conj()@O@amp.T)*G
            result+=scalar*np.sum(W*np.exp(1j*(left[sl,None]+right[None,:])))
    return result


def electric_density(state,cells,spin_index,electric_coordinates):
    """Density of the complete193coordinate electric and spin-basis law.

    Each sector's actual288edge field is P_(G-n)+frame electric_coordinates.
    The Gaussian momentum wave function includes free harmonic phases.
    """
    e=np.asarray(electric_coordinates,float)
    if e.shape!=(193,) or not np.all(np.isfinite(e)) or type(spin_index)is not int or not 0<=spin_index<state.amps.shape[1]:
        raise ValueError('all193real coordinates and physical spin index required')
    ids=np.flatnonzero(np.all(state.cells==tuple(cells),axis=1))
    vp=data()['vp'];answer=0j
    lognorm=-float(np.log(2*np.pi*vp).sum())/4
    for idx in ids:
        q=state.q[idx];p=state.p[idx];centered=e-p
        phase=-float(q@(e-p/2))-A*state.time*float(centered[190:]@centered[190:])/2
        wave=np.exp(lognorm-float(np.sum(centered*centered/vp))/4+1j*phase)
        answer+=state.amps[idx,spin_index]*wave
    return float(abs(answer)**2)


def electric_tail_bound(mean_upper,variance_upper,branch_l1,radius):
    if variance_upper<=0 or branch_l1<0 or mean_upper<0 or not all(math.isfinite(x) for x in (mean_upper,variance_upper,branch_l1,radius)):
        raise ValueError('finite Gaussian bounds required')
    if radius<=mean_upper:return 1.
    return min(1.,branch_l1*math.sqrt(576)*math.exp(-(radius-mean_upper)**2/(4*variance_upper)))


def free_reference_overlap(state,initial,h=.125,mass=1.):
    """Overlap with freely evolved field and the free rest spin phase.

    This comparison cancels a specified common field evolution; it is not
    stationarity of the absolute interacting state or its zero modes.
    """
    if initial.macros or initial.time:raise ValueError('time-zero reference required')
    particles=initial.cells.shape[1]
    O=matter.BETA if particles==1 else np.kron(matter.BETA,np.eye(4))+np.kron(np.eye(4),matter.BETA)
    if np.linalg.norm(initial.amps@O.T-particles*initial.amps)>1e-12:
        raise ValueError('positive rest-spin reference required')
    result=0j;vq,vp,cross=covariance(state.time)
    for cells,amp,q,p in zip(state.cells,state.amps,state.q,state.p):
        reference=sum((value for old,value in zip(initial.cells,initial.amps) if np.array_equal(cells,old)),start=np.zeros(initial.amps.shape[1],complex))
        variance=float(np.sum(p*p*vq+q*q*vp-2*p*q*cross))
        result+=np.vdot(reference,amp)*math.exp(-variance/2)
    # The field vacuum phase is explicitly matched to the reference.
    return result*np.exp(1j*(data()['E0']*state.time+particles*mass*h*state.macros))
