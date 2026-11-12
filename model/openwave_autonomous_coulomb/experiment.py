"""Autonomous arrivals, edge-current acceleration and retained loop records."""
import math
import numpy as np
import sympy as s
from scipy.linalg import eigh,expm
from .source import basis,matrices,column,native_certificate,ingredients,REFERENCE_MASK


def current_matrices():
    rows=basis();h,_=matrices();currents=[]
    for edge in range(4):
        j=s.zeros(len(rows))
        for a,(_,spin,y,f)in enumerate(rows):
            for b,(_,sp,z,g)in enumerate(rows):
                differing=[i for i in range(4)if f[i]!=g[i]]
                if a!=b and differing==[edge]:
                    j[a,b]=s.I*(int(y==edge)-int(z==edge))*h[a,b]
        assert j==j.H;currents.append(j)
    for cell in range(4):
        p=s.diag(*[int(y==cell)for _,_,y,_ in rows])
        assert s.I*(h*p-p*h)==currents[(cell-1)%4]-currents[cell]
    return tuple(currents)


def run():
    ancestry=native_certificate();rows=basis();hs,vs=matrices();js=current_matrices();dim=len(rows)
    candidates=[i for i,(_,spin,y,f)in enumerate(rows)if spin==0 and y==1]
    left=min(candidates,key=lambda i:sum(abs(v)for v in rows[i][3]))
    right=next(i for i,(_,spin,y,f)in enumerate(rows)if spin==0 and y==2 and hs[i,left]!=0)
    amplitude=hs[right,left];tau=s.sqrt(s.conjugate(amplitude)*amplitude)
    psi=s.zeros(dim,1);psi[left]=1/s.sqrt(2);psi[right]=-amplitude/(s.sqrt(2)*tau)
    V=vs;J=js[1];F=s.I*(V*J-J*V);force=s.simplify((psi.H*F*psi)[0])
    difference=V[right,right]-V[left,left]
    assert difference>0 and s.simplify(force+tau*difference)==0 and force<0
    # Every density operator is the literal native current of the full key.
    frozen=s.diag(*[int(y==1)for _,_,y,_ in rows]);assert hs*frozen!=frozen*hs
    theta=s.Symbol('theta',real=True);where={key:i for i,(key,_,_,_)in enumerate(rows)};derivative=s.zeros(dim)
    for b,(key,_,_,_)in enumerate(rows):
        for target,value in column(key,phase_edge=5,theta=theta).items():derivative[where[target],b]+=s.diff(value,theta).subs(theta,0)
    assert derivative==J
    identity=s.I*(hs*J-J*hs)-s.I*((hs-V)*J-J*(hs-V))
    assert identity==F
    p=np.asarray(psi,complex).ravel();position=np.array([y for _,_,y,_ in rows]);distance=np.minimum(position,4-position)
    runs=[];matrices_float={}
    for coupling in (0,1,4):
        # kappa=0 is a kinetic control; the full diagonal source is evaluated
        # by subtraction rather than violating its positive-coupling guard.
        h=np.asarray(hs+(coupling-1)*V,complex);h-=118*np.eye(dim);matrices_float[coupling]=h
        energy,u=eigh(h);weights=u.conj().T@p;initial_energy=np.vdot(p,h@p).real
        records=[];energy_error=0.;norm_error=0.
        for t in (0.,.1,.5,1.,3.,8.,20.,40.):
            state=u@(np.exp(-1j*energy*t)*weights);prob=abs(state)**2
            energy_error=max(energy_error,abs(np.vdot(state,h@state).real-initial_energy));norm_error=max(norm_error,abs(prob.sum()-1))
            records.append({'time':t,'all_joint_position_spin_flux_probabilities':prob.tolist(),
                'position_arrival_probabilities':[float(prob[position==y].sum())for y in range(4)],
                'mean_periodic_separation':float(prob@distance),
                'all_oriented_native_edge_currents':[float(np.vdot(state,np.asarray(j,complex)@state).real)for j in js],
                'retained_complement_probability':0.,'energy_including_reference_and_mass':float(118+np.vdot(state,h@state).real)})
        assert max(energy_error,norm_error)<1e-12
        runs.append({'coupling':coupling,'full_pair_sector_eigenvalues':(118+energy).tolist(),
            'all_joint_basis_labels':[{'spin':spin,'position':y,'native_flux':list(map(int,f))}for _,spin,y,f in rows],
            'initial_energy':118+initial_energy,'maximum_energy_drift':energy_error,'maximum_norm_drift':norm_error,'complete_records':records})
    h=matrices_float[1];j=np.asarray(J,complex);dt=1e-5
    expectation=lambda t:np.vdot(expm(-1j*h*t)@p,j@(expm(-1j*h*t)@p)).real
    numerical=(expectation(dt)-expectation(-dt))/(2*dt)
    exact=np.vdot(p,1j*(h@j-j@h)@p).real
    derivative_error=abs(numerical-exact);assert derivative_error<1e-10
    # A second complete evolution implementation with an explicit finite tail.
    t=1.;term=p.copy();series=term.copy()
    for k in range(1,49):term=(-1j*t/k)*(h@term);series+=term
    operator_bound=float(np.max(np.sum(abs(h),axis=1)))
    remainder=math.exp(t*operator_bound)*(t*operator_bound)**49/math.factorial(49)
    error=float(np.linalg.norm(series-expm(-1j*h*t)@p));assert error<remainder+2e-14
    tree=[i for i,(_,spin,y,f)in enumerate(rows)if spin==0 and f[3]==0]
    assert len(tree)==4 and left in tree and right in tree
    complement=[i for i in range(dim)if i not in tree]
    loop_leakage=float(np.linalg.norm(h[np.ix_(complement,tree)]))
    assert loop_leakage>.1
    full_state=expm(-20j*h)@p;cut_state=expm(-20j*h[np.ix_(tree,tree)])@p[tree]
    full_arrival=[float(sum(abs(full_state[i])**2 for i,(_,spin,y,f)in enumerate(rows)if y==cell))for cell in range(4)]
    cut_arrival=[float(sum(abs(cut_state[j])**2 for j,i in enumerate(tree)if rows[i][2]==cell))for cell in range(4)]
    loop_mutation=sum(abs(x-y)for x,y in zip(full_arrival,cut_arrival))/2
    assert loop_mutation>.01
    # Same native current pattern reproduces the pinned-minimum difference,
    # but the autonomous ground and arrivals are computed from H itself.
    potentials=[next(str(V[i,i])for i,(_,_,q,_)in enumerate(rows)if q==y)for y in range(4)]
    return {'schema':'oph.autonomous_native_ec_coulomb.v1','verdict':'PASS','source_ancestry':ancestry,
        'complete_clock':'exp(-it H) on the declared full carrier; the exact invariant physical sector is exponentiated without pinning or intermediate measurements',
        'initial_preparation':{'left_basis_row':left,'right_basis_row':right,'exact_full_physical_state':[str(v)for v in psi],
            'reference_mask':str(REFERENCE_MASK),'preparation':'declared normalized coherent source state; other reference and preparation outcomes belong to the retained complement'},
        'native_current_and_force':{'all_four_exact_continuity_identities':True,'actual_gauge_phase_derivative_equals_current':True,
            'complete_Ehrenfest_potential_force_identity':True,'selected_bond_potential_increase':str(difference),
            'exact_attractive_current_acceleration_contribution':str(force),'kinetic_amplitude':str(tau),
            'force_is_negative_coherence_times_actual_potential_difference':True,'independent_time_derivative_error':derivative_error},
        'same_Green_pinned_control_potentials':potentials,'autonomous_runs':runs,
        'independent_full_matrix_Taylor_control':{'terms':49,'time':t,'norm_bound':operator_bound,'analytic_tail':remainder,'measured_difference':error},
        'deleted_flux_loop_mutation':{'noninvariant_tree_subspace_leakage_norm':loop_leakage,'position_TV_at_time20':loop_mutation},
        'scope':{'source_is_the_unchanged_interacting_PW119_H_prime':False,'actual_native_PW_coordinate_and_all_twelve_Gauss_currents':True,
            'full_invariant_flux_loops_and_both_spin_branches_retained':True,'field_energy':'three-dimensional N=4 zero-mode-removed Green energy of the actual native current plus reference',
            'force_meaning':'autonomous edge-current acceleration from the Hamiltonian commutator; not a classical radial force substituted from pinned minima',
            'physical_alpha_GeV_or_hydrogen_identification':False,'free_loop_modes':'all selected circulation patterns retained; no transverse-loop energy is declared',
            'field_dynamics':'instantaneous longitudinal quadratic energy; no transverse radiation or full Maxwell kinetic field'}}
