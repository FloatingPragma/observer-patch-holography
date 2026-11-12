"""All-mode spectral and quasistatic error bounds for actual sparse eigenlines."""
from fractions import Fraction as F
from math import ceil,isqrt,factorial
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal
from model.openwave_intrinsic_magnetic.solver import GAMMA,BETA,SIGMA,dirac,geometry,source_current,moment
from model.openwave_larmor_dynamics.dynamics import field_minimizer
from model.openwave_field_observables.lattice import curl,curl_adjoint,divergence_adjoint
HERE=Path(__file__).resolve().parent


def cadd(a,b):return a[0]+b[0],a[1]+b[1]
def cmul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def cscale(a,t):return a[0]*t,a[1]*t
def conjugate(a):return a[0],-a[1]
def pair(z):return F(float(np.real(z))),F(float(np.imag(z)))


def phase_polynomial(theta,degree=18):
    term=(F(1),F(0));result=term
    for k in range(1,degree+1):
        term=cscale(cmul(term,(F(0),theta)),F(1,k));result=cadd(result,term)
    return result


def rational_residual(vector,energy,n,a,B,mass=F(8),trap=F(1,100),r=F(1,2)):
    nn=n**3;coords=np.indices((n,)*3).reshape(3,-1).T-(n-1)//2
    values=[[pair(vector[spin*nn+i]) for i in range(nn)] for spin in range(4)]
    residual=[]
    energy=F(float(energy));norm2=F(0)
    for spin in range(4):
        row=[];sign=1 if spin<2 else -1
        for i,x in enumerate(coords):
            coefficient=sign*(mass+trap*a*a*sum(int(t)**2 for t in x)+3*r/a)-energy
            row.append(cscale(values[spin][i],coefficient));norm2+=sum(q*q for q in values[spin][i])
        residual.append(row)
    phases={j:phase_polynomial(B*a*a*j/2) for j in range(-(n-1)//2,(n-1)//2+1)}
    blocks=[]
    for axis in range(3):
        block=[]
        for i in range(4):
            for j in range(4):
                gamma=pair(GAMMA[axis][i,j]);beta=F(float(BETA[i,j].real))
                coefficient=cscale(cadd(cmul((F(0),F(1)),gamma),(-r*beta,F(0))),1/(2*a))
                if coefficient!=(0,0):block.append((i,j,coefficient))
        blocks.append(block)
    for left,x in enumerate(coords):
        grid=x+(n-1)//2
        for axis in range(3):
            if grid[axis]>=n-1:continue
            y=grid.copy();y[axis]+=1;right=int((y[0]*n+y[1])*n+y[2])
            index=-int(x[1]) if axis==0 else int(x[0]) if axis==1 else 0
            phase=phases[index]
            for i,j,base in blocks[axis]:
                coefficient=cmul(base,phase)
                residual[i][left]=cadd(residual[i][left],cmul(coefficient,values[j][right]))
                residual[j][right]=cadd(residual[j][right],cmul(conjugate(coefficient),values[i][left]))
    infinity=max(abs(x)+abs(y) for row in residual for x,y in row)
    max_angle=abs(B)*a*a*((n-1)//2)/2
    assert max_angle<F(1,2) and F(999999999999,10**12)<norm2<F(1000000000001,10**12)
    tail=2*max_angle**19/F(factorial(19))
    operator_tail=3*(1+r)*tail/a
    raw_bound=F(isqrt(4*nn)+1)*infinity+2*operator_tail
    # ||v|| >= 1-10^-12 follows from the exact squared-norm check.
    # Cover both raw and normalized residuals for the preparation estimate
    # and the separate spectral root enclosures.
    bound=raw_bound/(1-F(1,10**12))
    scale=10**15;rounded=F(ceil(bound*scale),scale)
    return {'analytic_matrix_residual_upper':str(rounded),'decimal_upper':float(rounded),
        'matrix_exponential_Taylor_degree':18,'entrywise_tail_operator_upper':str(operator_tail),
        'exact_binary_vector_norm_squared':str(norm2),'all_internal_spatial_rows':4*nn,
        'bounds_both_raw_and_normalized_vector_residual':True,
        'method':'Exact rational matvec of degree18 link-phase polynomial, coordinate residual norm bound and analytic exponential tail; binary eigenvector and eigenvalue inputs interpreted exactly.'}


def scalar_mass_gap(n,a,B):
    r=F(1,2);eta=F(1,100)
    diagonal=[r/a+eta*(a*j)**2 for j in range(-(n-1)//2,(n-1)//2+1)]
    off=-r/(2*a)
    values=eigh_tridiagonal(np.array(diagonal,float),np.full(n-1,float(off)),eigvals_only=True)
    def negative_count(lam):
        count=0;previous=None
        for d in diagonal:
            pivot=d-lam if previous is None else d-lam-off*off/previous
            assert pivot!=0
            count+=pivot<0;previous=pivot
        return count
    bounds=[]
    for index in (0,1):
        scale=10**12;lo=F(int(np.floor(values[index]*scale))-2,scale);hi=F(int(np.ceil(values[index]*scale))+2,scale)
        assert negative_count(lo)==index and negative_count(hi)==index+1
        bounds.append((lo,hi))
    # Four internal copies of the scalar M enter beta D's Hermitian part.
    # The fifth singular value of D is at least scalar spatial eigenvalue2.
    phase_max=abs(B)*a*a*((n-1)//2)/2
    perturb=3*r*phase_max/a
    lower=8+2*bounds[0][0]+bounds[1][0]-perturb
    return {'one_dimensional_lowest_two_exact_Sturm_intervals':[[str(x) for x in row] for row in bounds],
        'uniform_probe_mass_operator_perturbation_upper':str(perturb),
        'fifth_absolute_D_eigenvalue_lower':str(lower),
        'proof':'Re<beta v,Dv>=<v,Mv> and singular-value min-max; the scalar M0 separates into three exact tridiagonal axes. Four internal copies give the fifth-eigenvalue bound.'}


def uniform_field_bound(n,a,N,kappa):
    r=F(1,2);m=F(8);trap=F(1,100)
    radius2=3*(a*((n-1)//2))**2
    M=m+trap*radius2+6*r/a+3/a
    first_squared=2*(1+r*r)
    C=((1+(M/m)**2)*first_squared+M*a*(1+r)/2)/m
    L=(M/m)*F(8,5) # sqrt(5/2)<8/5
    sin_lower=F(3,N)-F(22,7*N)**3/6
    gap=4*sin_lower**2/a**2
    q=kappa*C/(a**3*gap);Amax=kappa*L/(a**3*gap);delta=L*Amax
    assert q<1
    return {'kappa':str(kappa),'global_D_norm_upper':str(M),'current_dual_norm_upper':str(L),
        'transverse_Laplacian_gap_lower':str(gap),'contraction_upper':str(q),
        'uniform_self_field_norm_upper':str(Amax),'uniform_kernel_change_upper':str(delta),
        'finite_time_Duhamel_bound':'For every normalized trajectory, ||psi_qs(t)-exp(-itK0)psi(0)||<=|t|*uniform_kernel_change_upper. This uses a uniform perturbation of the fixed unitary, with no nonlinear Gronwall factor.'}


def projected_records(vectors,energies,n,a,N,B,times):
    field=np.zeros((N,)*3+(3,));h,phases=dirac(n,a,N,field,external_B=B)
    U=vectors.copy();nn=n**3
    spin=[np.block([[s,np.zeros((2,2))],[np.zeros((2,2)),s]])/2 for s in SIGMA]
    S=[np.einsum('ain,ab,bim->nm',U.reshape(4,nn,2).conj(),s,U.reshape(4,nn,2)) for s in spin]
    U[:,1]*=np.exp(-1j*np.angle(S[0][0,1]))
    blocks=U.reshape(4,nn,2)
    S=[np.einsum('ain,ab,bim->nm',blocks.conj(),s,blocks) for s in spin]
    E=np.asarray(energies,float)
    bonds=[];coords,position,_=geometry(n,a,N)
    moments=np.zeros((3,2,2),complex)
    weights=np.zeros(3)
    for left,right,axis,grid,phase in phases:
        half=blocks[:,left,:].conj().T@((GAMMA[axis]+.5j*BETA)@(phase*blocks[:,right,:]))/2
        op=half+half.conj().T
        midpoint=(position[left]+position[right])/2
        weight=np.cross(midpoint,np.eye(3)[axis])/2
        moments+=weight[:,None,None]*op
        weights+=weight*weight
        bonds.append((grid,axis,op))
    rows=[]
    # Full-space normalized input. The complete source acts on all modes;
    # the two columns specify the preparation and its certified residual.
    initial_norm=float(np.linalg.norm(U[:,0]+U[:,1]))
    for t in times:
        c=np.exp(-1j*(E-E[0])*t)/initial_norm
        spins=[float(np.vdot(c,s@c).real) for s in S]
        mu=[float(np.vdot(c,m@c).real) for m in moments]
        j=np.zeros_like(field)
        for grid,axis,op in bonds:j[grid+(axis,)]=np.vdot(c,op@c).real
        tiny_field=field_minimizer(j,a,1e-14)
        rows.append({'time':float(t),'spin':spins,'moment':mu,
            'complete_spin_probabilities':[[.5+x,.5-x] for x in spins],
            'instantaneous_Green_field_norm':float(np.linalg.norm(tiny_field))})
    return {'rows':rows,'energies':E.tolist(),'moment_weight_norms':np.sqrt(weights).tolist(),
        'initial_vector_normalization':initial_norm,'source_basis_Gram_residual':float(np.linalg.norm(U.conj().T@U-np.eye(2)))}


def spin_radius_certificate(vectors,n,normsum_upper):
    nn=n**3;v=[[pair(z) for z in vectors[:,j]] for j in range(2)]
    z={}
    for i in range(2):
        for j in range(2):
            value=(F(0),F(0))
            for offset in (0,2*nn):
                for x in range(nn):value=cadd(value,cmul(conjugate(v[i][offset+x]),v[j][offset+nn+x]))
            z[i,j]=value
    squares={key:sum(t*t for t in val) for key,val in z.items()}
    dominant=max(((0,1),(1,0)),key=lambda key:squares[key]);other=tuple(reversed(dominant))
    square=squares[dominant];scale=10**14
    lower=F(isqrt(square.numerator*scale*scale//square.denominator),scale)
    diagonal=cadd(z[0,0],z[1,1]);error=sum(abs(t) for t in diagonal)+sum(abs(t) for t in z[other])
    radius=(lower-error)/normsum_upper
    assert radius>F(49,100) and error/lower<F(1,2)
    return {'dominant_spin_transition_modulus_lower':str(lower),'counterrotation_and_offset_upper':str(error),
        'uniform_frozen_transverse_spin_radius_lower':str(radius),
        'phase_deviation_from_single_split_frequency_upper':str(4*error/lower),
        'method':'Exact binary-vector rational matrix elements of Sx+iSy; triangle inequality and integer-square-root bound, uniform over every time.'}


def same_member_stationary(zero,n,a,N,kappa):
    basis=np.linalg.qr(zero['vectors'][:,:2])[0];nn=n**3
    spin_diagonal=np.repeat((.5,-.5,.5,-.5),nn)
    sw,su=np.linalg.eigh(basis.conj().T@(spin_diagonal[:,None]*basis))
    state=basis@su[:,-1];state/=np.linalg.norm(state)
    blank=np.zeros((N,)*3+(3,));h0,ph0=dirac(n,a,N,blank)
    j0=source_current(state,ph0,N)
    field=field_minimizer(j0,a,kappa)
    h,ph=dirac(n,a,N,field);energy=float(np.vdot(state,h@state).real)
    residual=float(np.linalg.norm(h@state-energy*state))
    j=source_current(state,ph,N);b=curl(field,a)
    field_energy=a**3*np.sum(b*b)/(2*kappa)
    scaled_field_defect=float(np.linalg.norm(field_minimizer(j,a,kappa)-field)/kappa)
    scaled_Maxwell=float(np.linalg.norm(curl_adjoint(b,a)-kappa*j/a**3)/kappa)
    zero_field_mutation=float(np.linalg.norm(j0/a**3))
    assert residual<1e-11 and scaled_field_defect<1e-9 and scaled_Maxwell<1e-8
    assert zero_field_mutation>1e-3 and np.linalg.norm(field)>0
    spin=float(np.vdot(state,spin_diagonal*state).real);mu=moment(j,a)
    return {'source_kappa':kappa,'n':n,'spacing':a,'field_n':N,'external_B':0.,
        'preparation':'positive-spin line selected within the actual full zero-field ground doublet, with its computed transverse self-field',
        'eigenstate_residual':residual,'self_field_defect_divided_by_kappa':scaled_field_defect,
        'Maxwell_residual_divided_by_kappa':scaled_Maxwell,
        'zero_self_field_mutation_scaled_residual':zero_field_mutation,
        'current_divergence':float(np.linalg.norm(divergence_adjoint(j,a))),
        'nonzero_self_field_norm':float(np.linalg.norm(field)),
        'actual_intrinsic_current_moment':mu.tolist(),'canonical_spin':spin,
        'bound_excitation_plus_field_energy':float(energy+field_energy),
        'static_g_canonical_spin':float(-2*(energy+field_energy)*mu[2]/spin),
        'static_g_doublet_normalization':float(-4*(energy+field_energy)*mu[2]),
        'scope':'Numerical stationary-residual witness on the same refined weak-coupling member, including its nonzero self-field. These residuals are not an interval proof of a particular exact nonlinear stationary eigenvector.'}


def run():
    zero=np.load(HERE/'witnesses/probe_zero.npz',allow_pickle=False)
    probe=np.load(HERE/'witnesses/probe_B001.npz',allow_pickle=False)
    n=int(probe['n']);a=F(6,n+1);N=int(probe['field_n']);B=F(1,100);kappa=F(1,10**14)
    zero_cert=rational_residual(zero['vectors'][:,0],zero['energies'][0],n,a,F(0))
    certificates=[rational_residual(probe['vectors'][:,j],probe['energies'][j],n,a,B) for j in range(2)]
    zero_gap=scalar_mass_gap(n,a,F(0));probe_gap=scalar_mass_gap(n,a,B)
    zenergy=F(float(zero['energies'][0]));zres=F(zero_cert['analytic_matrix_residual_upper'])
    energies=[F(float(x)) for x in probe['energies'][:2]];residuals=[F(x['analytic_matrix_residual_upper']) for x in certificates]
    assert zenergy+zres<F(zero_gap['fifth_absolute_D_eigenvalue_lower'])
    assert energies[0]+residuals[0]<energies[1]-residuals[1]
    assert energies[1]+residuals[1]<F(probe_gap['fifth_absolute_D_eigenvalue_lower'])
    bounds=uniform_field_bound(n,a,N,kappa);delta=F(bounds['uniform_kernel_change_upper'])
    duration=F(6000);records=projected_records(probe['vectors'][:,:2],probe['energies'][:2],n,float(a),N,float(B),np.linspace(0,float(duration),25))
    # The sum vector has norm>1 by the exact near-unit norm certificates
    # and numerical Gram control guarded independently by the exact residual
    # energy separation: |<u,v>|<=(r_u+r_v)/|E_u-E_v|.
    overlap=(residuals[0]+residuals[1])/(energies[1]-energies[0])
    normsum_lower=sum(F(c['exact_binary_vector_norm_squared']) for c in certificates)-2*overlap
    normsum_upper=sum(F(c['exact_binary_vector_norm_squared']) for c in certificates)+2*overlap
    assert normsum_lower>1
    propagation_error=duration*(sum(residuals)+delta)
    assert propagation_error<F(1,1000)
    # Unit-normalized exact initial state has normsum>1, so the sum of the
    # two residual bounds safely bounds its linear preparation residual.
    difference=energies[1]-energies[0]
    mass_interval=[zenergy-zres-delta,zenergy+zres]
    frequency_interval=[difference-sum(residuals),difference+sum(residuals)]
    phase=np.unwrap([np.angle(x['spin'][0]+1j*x['spin'][1]) for x in records['rows']])
    omega=float((phase[-1]-phase[0])/float(duration))
    moment_phase=np.unwrap([np.angle(-(x['moment'][0]+1j*x['moment'][1])) for x in records['rows']])
    momega=float((moment_phase[-1]-moment_phase[0])/float(duration))
    # The approximation is initially normalized, not exactly normalized
    # at every time. A spin-component error is at most
    # (1+epsilon/2)epsilon; the phase bound covers both components.
    spin_certificate=spin_radius_certificate(probe['vectors'][:,:2],n,normsum_upper)
    spin_floor=F(spin_certificate['uniform_frozen_transverse_spin_radius_lower'])
    phase_error=4*propagation_error/spin_floor+F(spin_certificate['phase_deviation_from_single_split_frequency_upper'])
    omega_error_exact=2*phase_error/duration
    omega_error=float(omega_error_exact)
    assert abs(abs(omega)-float(difference))<1e-10
    assert abs(momega-omega)<1e-10
    g_center=2*float(zenergy)*abs(omega)/float(B)
    assert difference>omega_error_exact
    g_lower=2*mass_interval[0]*(difference-omega_error_exact)/B
    g_upper=2*mass_interval[1]*(difference+omega_error_exact)/B
    def decimal_endpoint(value,up=False):
        scaled=value*10**12;integer=ceil(scaled) if up else scaled.numerator//scaled.denominator
        whole,fraction=divmod(integer,10**12)
        return str(whole)+'.'+str(fraction).zfill(12)
    L=F(bounds['current_dual_norm_upper']);R=F(bounds['uniform_self_field_norm_upper'])
    M=F(bounds['global_D_norm_upper']);r=F(1,2)
    C=((1+(M/8)**2)*2*(1+r*r)+M*a*(1+r)/2)/8
    negative_component=sum(residuals)/(min(energies)+8)
    current_error=L*(2+propagation_error)*propagation_error+C*R+4*L*negative_component*(1+negative_component+propagation_error)
    weight_squared=F((n-1)*n,2)*a*a*sum(j*j for j in range(-(n-1)//2,(n-1)//2+1))
    weight_upper=F(isqrt(weight_squared.numerator//weight_squared.denominator)+1)
    moment_error=weight_upper*current_error
    moment_bound=weight_upper*L
    for row in records['rows']:
        row['complete_weak_moment_probabilities']=[[float((1+x/float(moment_bound))/2),float((1-x/float(moment_bound))/2)] for x in row['moment']]
    return {'schema':'oph.refined_all_mode_larmor_certificate.v1','n':n,'spacing':str(a),'field_n':N,'external_B':str(B),
        'source':'same finite |D_W| matter and unique quasistatic field at kappa=1e-14; this coupling is a declared weak-self-field accuracy choice, not a fitted physical electromagnetic strength',
        'zero_probe_exact_matrix_residual':zero_cert,'probe_exact_matrix_residuals':certificates,
        'zero_probe_all_mode_gap':zero_gap,'probe_all_mode_gap':probe_gap,
        'ground_identification':'The zero-field antiunitary with square-1 and exact chiral symmetry give a fourfold |D| ground quartet. At nonzeroB two isolated positive eigenlines plus chiral partners account for all four levels below the rigorous fifth-level bound. Four identical native blocks retain all16internal rows.',
        'uniform_quasistatic_bounds':bounds,'exact_initial_sum_norm_squared_lower':str(normsum_lower),
        'same_member_intrinsic_stationary_witness':same_member_stationary(zero,n,float(a),N,float(kappa)),
        'duration':str(duration),'complete_state_error_upper_at_duration':str(propagation_error),
        'state_error_decimal_upper':float(propagation_error),'rest_energy_interval':[str(x) for x in mass_interval],
        'uniform_spin_phase_certificate':spin_certificate,
        'moment_readout_error':{'uniform_component_error_upper':str(moment_error),'decimal_upper':float(moment_error),
            'full_moment_operator_norm_upper':str(moment_bound),
            'complete_binary_instrument':'effects(I+/-M_j/L_mu)/2 on the full signed kernel; probabilities are computed from the prepared pair with the stated state/field/negative-branch error',
            'probability_error_upper':str(moment_error/(2*moment_bound))},
        'isolated_energy_split_interval':[str(x) for x in frequency_interval],
        'records':records,'spin_cycle_average_frequency':omega,'dipole_cycle_average_frequency':momega,
        'cycle_average_frequency_error_upper':omega_error,
        'cycle_average_frequency_error_exact_upper':str(omega_error_exact),
        'g_from_actual_rest_energy_and_precession':g_center,
        'finite_g_interval_exact':[str(g_lower),str(g_upper)],
        'finite_g_interval_outward_decimal':[decimal_endpoint(g_lower),decimal_endpoint(g_upper,True)],
        'finite_g_definition':'2 E_star |Delta unwrapped arg<Sx+iSy>|/(B T), where E_star is the true zero-probe minimum of the same one-particle quasistatic construction and T=6000. This is a finite-duration precession normalization, distinct from dividing static moment by measured canonical spin.',
        'rigorous_inputs':'Exact rational source-matrix residuals, exact Sturm index intervals, uniform rational self-field/state bounds and exact spin-matrix phase bounds imply the rational finite-g enclosure. Individual plotted/printed phase and moment samples are numerical diagnostics.',
        'g_interval_is_a_continuum_or_physical_electron_prediction':False,
        'finite_grid_g_is_exactly_two':False,'unprepared_modes_removed_from_source':False,
        'nonlinear_evolution_is_a_linear_quantum_channel':False,
        'verdict':'PASS_REFINED_PRECESSION_WITH_ALL_MODE_STATE_ERROR_BOUND'}
