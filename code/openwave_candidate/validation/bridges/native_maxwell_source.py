"""One-way native coherent-current adapter into the existing Maxwell executor."""
from pathlib import Path
import sys,importlib.util,json,hashlib
import numpy as np
from scipy.integrate import quad
BASE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(BASE/'candidate'))
from oph_openwave.construction import Construction
from oph_openwave.uniform_native import UniformFiberReduction
PROVIDER=Path(__file__).resolve().parent/'providers/seam_maxwell_continuum.py'
spec=importlib.util.spec_from_file_location('existing_seam_maxwell',PROVIDER)
maxwell=importlib.util.module_from_spec(spec);spec.loader.exec_module(maxwell)


def source(model,psi,positions,k):
    """Effective port charge, exact half-hop current and straight-line Fourier map."""
    c=model.construction;psi=np.asarray(psi,complex);positions=np.asarray(positions,float);k=np.asarray(k,float)
    if positions.shape!=(c.ports,3) or k.shape!=(3,) or not np.isfinite(positions).all() or not np.isfinite(k).all():raise ValueError('finite declared embedding and 3-momentum required')
    if psi.shape!=(c.ports+len(c.seams),) or not np.isfinite(psi).all() or not np.isclose(np.vdot(psi,psi),1):raise ValueError('normalized finite native state required')
    rho=np.outer(psi,psi.conj());drho=-1j*(model.hamiltonian@rho-rho@model.hamiltonian)
    q=abs(psi[:c.ports])**2;dq=np.diag(drho).real[:c.ports].copy();currents=[];jhat=np.zeros(3,complex)
    for e,(a,b) in enumerate(c.seams):
        t=c.ports+e;q[a]+=.5*abs(psi[t])**2;q[b]+=.5*abs(psi[t])**2
        dq[a]+=.5*drho[t,t].real;dq[b]+=.5*drho[t,t].real
        ja=-2*model.hamiltonian[a,t]*rho[t,a].imag
        jb=-2*model.hamiltonian[b,t]*rho[t,b].imag
        current=(ja-jb)/2;currents.append(float(current))
        dx=positions[b]-positions[a];mid=(positions[a]+positions[b])/2
        jhat+=current*dx*np.exp(-1j*k@mid)*np.sinc((k@dx)/(2*np.pi))
    phase=np.exp(-1j*positions@k)
    return {'effective_port_charge':q,'effective_charge_derivative':dq,'whole_seam_current':np.array(currents),'charge_fourier':complex(q@phase),'charge_fourier_derivative':complex(dq@phase),'current_fourier':jhat}


def execute(steps=40):
    if type(steps) is not int or steps<1:raise ValueError('positive integer steps required')
    c=Construction(3,((0,1),(1,2),(2,0)),.23);model=UniformFiberReduction(c,6)
    positions=np.array([[0.,0,0],[1.,0,0],[.5,np.sqrt(3)/2,0]])
    k=np.array([.7,.4,.3]);omega,p,j=maxwell.blocks(k);pt=p[:3,:3]
    psi=np.array([1,2j,0,1-1j,0,0],complex);psi/=np.linalg.norm(psi)
    eig,v=np.linalg.eigh(model.hamiltonian);initial=psi.copy();field=np.zeros(6,complex)
    time=0.;dt=.4/steps;initial_source=source(model,psi,positions,k)
    held_charge=initial_source['charge_fourier'];energy_work=0.;continuity=0.;gauss=0.;records=[]
    initial_long=-1j*k*held_charge/(k@k);long=initial_long.copy();long_work=0.
    for n in range(steps):
        time=n*dt;psi=v@(np.exp(-1j*time*eig)*(v.T@initial));s=source(model,psi,positions,k)
        current=s['current_fourier'];jt=pt@current;jl=current-jt;forcing=np.r_[-jt,np.zeros(3)]
        continuity=max(continuity,float(abs(s['charge_fourier_derivative']+1j*k@current)))
        before=field.copy();u=maxwell.propagator(omega,dt,p,j);integ=maxwell.constant_forcing(omega,dt,p,j)
        field=u@before+integ@forcing
        work=quad(lambda t:float(np.vdot(maxwell.propagator(omega,t,p,j)@before+maxwell.constant_forcing(omega,t,p,j)@forcing,forcing).real),0,dt,epsabs=1e-13)[0]
        energy_work+=work
        previous_long=long.copy();long-=dt*jl
        long_work+=float(np.vdot(previous_long,-dt*jl).real+.5*dt*dt*np.vdot(jl,jl).real)
        held_charge-=dt*1j*k@current
        gauss=max(gauss,float(abs(1j*k@long-held_charge)))
        records.append({'time_end':float((n+1)*dt),'native_whole_seam_currents':s['whole_seam_current'].tolist(),'transverse_field_energy':float(np.vdot(field,field).real/2),'interval_source_work':work,'held_charge_real':float(held_charge.real),'held_charge_imag':float(held_charge.imag)})
    final=v@(np.exp(-1j*.4*eig)*(v.T@initial));actual=source(model,final,positions,k)
    energy=float(np.vdot(field,field).real/2);long_delta=float((np.vdot(long,long)-np.vdot(initial_long,initial_long)).real/2)
    return {'steps':steps,'source_hold_error':float(abs(held_charge-actual['charge_fourier'])),'fourier_continuity_error':continuity,'longitudinal_gauss_error':gauss,'transverse_energy':energy,'transverse_source_work':energy_work,'transverse_work_error':abs(energy-energy_work),'longitudinal_energy_change':long_delta,'longitudinal_source_work':long_work,'longitudinal_work_error':abs(long_delta-long_work),'records':records,'field_real':field.real.tolist(),'field_imag':field.imag.tolist()}


def run():
    coarse=execute(20);fine=execute(40)
    omega,p,j=maxwell.blocks(np.array([.7,.4,.3]));zero=maxwell.propagator(omega,.4,p,j)@np.zeros(6)
    paths=[Path(__file__).resolve(),BASE/'candidate/oph_openwave/construction.py',BASE/'candidate/oph_openwave/uniform_native.py',BASE/'candidate/oph_openwave/modular_native.py',PROVIDER]
    return {'schema':'oph.native.maxwell_source_bridge.v1','checks':{'exact_instantaneous_fourier_continuity':bool(fine['fourier_continuity_error']<1e-12),'held_source_longitudinal_gauss':bool(fine['longitudinal_gauss_error']<1e-12),'field_energy_source_work':bool(max(fine['transverse_work_error'],fine['longitudinal_work_error'])<1e-12),'nonzero_native_driven_field':bool(fine['transverse_energy']>1e-8),'zero_current_zero_field':bool(np.max(abs(zero))==0),'source_hold_refinement':bool(fine['source_hold_error']<coarse['source_hold_error'])},'parameters':{'native_graph':'triangle','mass':.23,'modulus':6,'preparation':'uniform modular Gauss fibers; declared coherent reduced state [1,2i,0,1-i,0,0]/sqrt7','commit_rate':0,'port_positions':[[0,0,0],[1,0,0],[.5,float(np.sqrt(3)/2),0]],'fourier_momentum':[.7,.4,.3],'field_interval':[0,.4],'source_scheme':'left-held exact native coherent current','field_dispersion':'existing continuum blocks omega=|k|; speed/embedding supplied'},'fine':fine,'coarse_source_hold_error':coarse['source_hold_error'],'assumptions':['one-way classical expectation current drives declared linear field; no quantum joint source/field state or backreaction','effective port charge gives each transit half to each endpoint, different from whole-transit regional apparatus','straight-edge Fourier line current matches exact effective point-charge continuity','longitudinal field initialized from same charge; source held each interval; native endpoint charge mismatch retained','single Fourier mode energy is an algebraic mode norm, not calibrated spatial distant flux or outgoing radiation'],'full_physical_criterion_met':False,'source_sha256':{str(path.relative_to(BASE)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}

if __name__=='__main__':
    r=run();Path(__file__).with_name('NATIVE_MAXWELL_SOURCE.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
