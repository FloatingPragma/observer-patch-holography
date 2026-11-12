"""Complete torque, field energy and angular records on a finite source."""
import numpy as np
from scipy.linalg import expm,eigh
from scipy.integrate import solve_ivp
from .source import oscillator,operators,evaluate,angular_expectation,field_energy,rhs,cross,native_reduction


def algebra(data,ops):
    J=data['J'];n=len(J[0]);zero=np.zeros((n,n),complex)
    h0,derivatives,probe,trap=ops
    lie=max(np.linalg.norm(J[i]@J[j]-J[j]@J[i]-1j*sum((cross(i,j,k)*J[k]for k in range(3)),zero))for i in range(3)for j in range(3))
    free=max(np.linalg.norm(j@h0-h0@j)+np.linalg.norm(j@trap-trap@j)for j in J)
    cov=max(np.linalg.norm(-1j*(J[i]@derivatives[j]-derivatives[j]@J[i])-sum((cross(k,i,j)*derivatives[k]for k in range(3)),zero))for i in range(3)for j in range(3))
    assert max(lie,free,cov)<2e-12
    theta=.37;U=expm(-1j*theta*J[2]);r=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
    b=np.array([.13,-.21,.08]);h=h0+sum((v*d for v,d in zip(b,derivatives)),zero)
    rotated=h0+sum((v*d for v,d in zip(r@b,derivatives)),zero)
    rotation_error=np.linalg.norm(U@h@U.conj().T-rotated)
    assert rotation_error<2e-12
    # Source spin alone does not rotate spatial modes or generate a symmetry.
    spatial=data['spatial_dimension'];from model.openwave_intrinsic_magnetic.solver import SIGMA
    spin=np.kron(np.kron(np.eye(2),SIGMA[2]/2),np.eye(spatial))
    spin_only=np.linalg.norm(spin@h0-h0@spin);assert spin_only>1
    # Replacing projected differential L by projected-coordinate cross-product
    # loses the shell-boundary term and fails its Lie relation.
    bad=tuple(sum((cross(i,j,k)*data['x'][j]@data['p'][k]for j in range(3)for k in range(3)),np.zeros((spatial,spatial),complex))for i in range(3))
    boundary=np.linalg.norm(bad[0]@bad[1]-bad[1]@bad[0]-1j*bad[2]);assert boundary>.1
    return {'su2_matrix_residual':float(lie),'free_rotation_residual':float(free),'localized_field_covariance_residual':float(cov),
        'noncubic_angle_rotation_residual':float(rotation_error),'spin_only_mutation_residual':float(spin_only),
        'projected_coordinate_orbital_mutation_residual':float(boundary)}


def run():
    native=native_reduction();data=oscillator();ops=operators(data,mass=4.,omega=.3);n=len(ops[0]);coupling=.2
    rng=np.random.default_rng(281013);psi=rng.normal(size=n)+1j*rng.normal(size=n);psi/=np.linalg.norm(psi)
    b=np.array([.04,-.03,.02]);pi=np.array([.03,.01,-.02]);h,g,_=evaluate(psi,b,ops)
    torque=np.array([(1j*np.vdot(psi,(h@j-j@h)@psi)).real for j in data['J']])
    identity=float(np.linalg.norm(torque-np.cross(b,g)));assert identity<2e-13
    derivative_errors=[]
    for axis in range(3):
        delta=np.eye(3)[axis]*1e-5
        ep=np.vdot(psi,evaluate(psi,b+delta,ops)[0]@psi).real
        em=np.vdot(psi,evaluate(psi,b-delta,ops)[0]@psi).real
        derivative_errors.append(abs((ep-em)/2e-5-g[axis]))
    assert max(derivative_errors)<3e-9
    initial=np.r_[psi,b,pi];times=np.linspace(0,8,33);runs=[]
    for tolerance in (2e-8,2e-10):
        sol=solve_ivp(rhs,(0,8),initial,args=(data,ops,coupling),t_eval=times,rtol=tolerance,atol=tolerance/30,method='DOP853')
        assert sol.success;energies=[];momenta=[];records=[];norms=[]
        for t,state in zip(sol.t,sol.y.T):
            v=state[:n];a=state[n:n+3].real;p=state[n+3:].real;H,_,_=evaluate(v,a,ops)
            energies.append(float(np.vdot(v,H@v).real+field_energy(a,p,data,coupling)))
            momenta.append(angular_expectation(v,data)+np.cross(a,p));norms.append(float(np.vdot(v,v).real))
            w,u=eigh(data['J'][2]);prob=np.abs(u.conj().T@v)**2
            records.append({'time':float(t),'all_Jz_eigenvalues':w.tolist(),'all_Jz_probabilities':prob.tolist(),
                'field_amplitude':a.tolist(),'field_canonical_momentum':p.tolist()})
        energy_error=max(abs(np.array(energies)-energies[0]));angular_error=max(np.linalg.norm(x-momenta[0])for x in momenta)
        assert energy_error<2e-6 and angular_error<2e-6
        runs.append({'tolerance':tolerance,'function_evaluations':sol.nfev,'maximum_norm_error':max(abs(np.array(norms)-1)),
            'maximum_energy_drift':float(energy_error),'maximum_total_angular_drift':float(angular_error),
            'complete_separate_copy_angular_records':records,'final_state':sol.y[:,-1]})
    final_difference=float(np.linalg.norm(runs[0].pop('final_state')-runs[1].pop('final_state')))
    static=np.array([np.linalg.norm(angular_expectation(expm(-1j*h*t)@psi,data)-angular_expectation(psi,data))for t in times])
    assert max(static)>1e-4
    return {'schema':'oph.localized_angular_torque.v1','verdict':'PASS','native_internal_reduction':native,
        'construction':{'spatial_shell':2,'one_block_rows':n,'complete_native_rows':4*n,'mass':4,'charge':-1,'trap':.3,'coupling':coupling,
            'positive_kernel':'absolute value of the complete signed Dirac matrix plus Omega times oscillator number',
            'field':'all three real Gaussian transverse dipole modes and all three conjugate momenta',
            'electric_metric':data['electric_metric'],'magnetic_metric':data['magnetic_metric']},
        'algebra':algebra(data,ops),'actual_signed_current_derivative_errors':derivative_errors,'instantaneous_total_torque_residual':identity,
        'time_refinement_runs':runs,'final_state_time_refinement_difference':final_difference,
        'frozen_field_missing_torque_mutation':float(max(static)),
        'scope':{'continuous_generator':'projected oscillator differential angular momentum plus native canonical spin; not a logarithm assigned to a cubic rotation',
            'field_angular_momentum':'b cross pi, derived from the Maxwell kinetic metric of the retained transverse mode space',
            'vacuum_and_fock':'empty vacuum, full dGamma frozen-field extension and all four native blocks; only one block is occupied by the declared input',
            'dynamics':'semiclassical Maxwell-mode mean field; not a linear quantum channel on arbitrary entangled inputs',
            'transfer':'same native Clifford algebra, declared rotational regulator and positive-kernel trap; no equality to the cubic or PW119 Hamiltonian',
            'omitted_physics':'longitudinal electrostatic backreaction, radiative continuum, physical electron identification and physical units'}}
