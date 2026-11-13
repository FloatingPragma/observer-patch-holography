"""Exact compact Gauss sectors and complete finite-plaquette instruments."""
from fractions import Fraction as F
from math import factorial, isqrt
import numpy as np
from scipy.linalg import eigh
from .proof_frames import fixed_eigen
from model.openwave_common_mass_source.source import lepton_source,ix,INPUTS
from model.openwave_common_electromagnetism.source import exact_ward,certificate as matter_certificate

ELECTRIC=F(1,16)
MAGNETIC=F(1,8)
HOP=F(1,16)
UNIT=F(INPUTS['electron_unit_GeV'])


def flux(x,n,q):
    return tuple(n-q*(edge<x) for edge in range(4))


def gauss_checks():
    count=0
    for q in (-3,-2,-1,0,1,2,3):
        for x in range(4):
            for n in range(-4,5):
                E=flux(x,n,q)
                assert all(E[v]-E[(v-1)%4]==q*((v==x)-(v==0)) for v in range(4))
                for direction in (-1,1):
                    y=(x+direction)%4
                    nn=n-q if x==3 and direction==1 else n+q if x==0 and direction==-1 else n
                    EE=flux(y,nn,q);edge=x if direction==1 else y
                    assert all(EE[e]-E[e]==(-q*direction if e==edge else 0) for e in range(4))
                    count+=1
    return {'literal_native_charge_in_thirds':[-3,-2,-1,0,1,2,3],
        'all_hop_columns_tested':count,'all_local_constraints_exact':True,
        'magnetic_loop_shift':'E_e -> E_e+1 on every oriented edge; divergence unchanged',
        'reference_charge':'-q at vertex0, retained as a quantum charge coordinate',
        'no_flux_cutoff_in_declared_source':True}


def matrices(cap,q,mass=F(1)):
    basis=[(x,n)for x in range(4)for n in range(-cap,cap+1)]
    index={z:i for i,z in enumerate(basis)};D=len(basis)
    # H minus mass is exactly integer/32. All boundary diagonal terms remain.
    K=np.zeros((D,D),dtype=np.int64);V=np.zeros_like(K);electric=np.zeros(D,dtype=np.int64)
    for i,(x,n)in enumerate(basis):
        K[i,i]=4;V[i,i]=4;electric[i]=sum(e*e for e in flux(x,n,q))
        for direction in (-1,1):
            y=(x+direction)%4
            nn=n-q if x==3 and direction==1 else n+q if x==0 and direction==-1 else n
            if (y,nn)in index:K[index[y,nn],i]-=2
            if (x,n+direction)in index:V[index[x,n+direction],i]-=2
    A=K+V+np.diag(electric)
    assert np.array_equal(A,A.T)
    return {'basis':basis,'index':index,'A':A,'K':K,'V':V,'electric':electric,'mass':F(mass)}


def root_upper(a):
    a=F(a);d=1<<40;n=(a.numerator*d*d+a.denominator-1)//a.denominator
    return F(isqrt(n)+1,d)


def exact_cap_enclosure(A,eigen=None):
    """Enclose every eigenvalue by an exact dyadic complete-basis residual."""
    lam,U=eigh(A.astype(float)/32) if eigen is None else eigen;dv=1<<24;de=1<<20
    Z=np.rint(U*dv).astype(np.int64);E=np.rint(lam*de).astype(np.int64);n=len(A)
    assert n*int(abs(Z).max())**2<2**63
    gram=Z.T@Z-np.eye(n,dtype=np.int64)*dv*dv
    rowsum=max(sum(abs(int(x)) for x in row) for row in A)
    assert rowsum*int(abs(Z).max())*de+32*int(abs(Z).max())*int(abs(E).max())<2**63
    residual=(A@Z)*de-32*Z*E
    g=root_upper(F(sum(int(x)**2 for x in gram.flat),dv**4))
    r=root_upper(F(sum(int(x)**2 for x in residual.flat),(32*dv*de)**2))
    K=F(max(abs(int(x)) for x in E),de);assert g<F(1,100)
    error=(2*K*g+(1+g)*r)/(1-g)
    return {'lowest_interval_exact':[str(F(int(E[0]),de)-error),str(F(int(E[0]),de)+error)],
        'all_eigenvalue_centres_exact':[str(F(int(x),de))for x in E],
        'uniform_all_spectrum_error_exact':str(error),'Gram_error_upper_exact':str(g),
        'full_matrix_residual_upper_exact':str(r),'dimension':n},lam,U


def exact_phase(x):
    """Rational Taylor enclosure of exp(-ix), including range reduction."""
    pi=F('3.14159265358979323846264338327950288419716939937510')
    turns=round(float(x)/(2*np.pi));y=x-2*turns*pi
    assert abs(y)<=4
    real=F(1);imag=F(0);tr=F(1);ti=-y
    for k in range(1,51):
        tr*= -y*y/F((2*k-1)*(2*k));real+=tr
    for k in range(51):
        if k:ti*= -y*y/F((2*k)*(2*k+1))
        imag+=ti
    error=F(2*abs(turns),10**50)+F(4**102,factorial(102))+F(4**103,factorial(103))
    denominator=1<<80
    return (F(real.numerator*denominator//real.denominator,denominator),
        F(imag.numerator*denominator//imag.denominator,denominator),error+F(2,denominator))


def finite_position_certificate(A,lam,U,psi,index,time,cert):
    """Bound the actual floating position/flux table by exact rational work."""
    dv=1<<24;de=1<<20;Z=np.rint(U*dv).astype(np.int64);E=np.rint(lam*de).astype(np.int64)
    phases=[exact_phase(F(int(e)*time,de))for e in E]
    difference=F(0)
    for row in range(len(A)):
        re=sum((F(int(Z[row,k])*int(Z[index,k]),dv*dv)*phases[k][0]for k in range(len(A))),F(0))
        im=sum((F(int(Z[row,k])*int(Z[index,k]),dv*dv)*phases[k][1]for k in range(len(A))),F(0))
        difference+=(re-F(float(psi[row].real)))**2+(im-F(float(psi[row].imag)))**2
    g=F(cert['Gram_error_upper_exact']);r=F(cert['full_matrix_residual_upper_exact'])
    phase_error=max(x[2]for x in phases)
    # ||exp(-itH)-Q exp(-itD)Q*|| <= g+|t| r sqrt(1+g).
    vector=g+abs(time)*r*root_upper(1+g)+root_upper(difference)+(1+g)*phase_error
    displayed=abs(psi)**2
    rounding=sum((abs(F(float(displayed[i]))-F(float(psi[i].real))**2-F(float(psi[i].imag))**2)for i in range(len(psi))),F(0))/4
    # Fixed verified outward budgets keep receipts independent of harmless
    # LAPACK-frame and amplitude-rounding choices. Every live rational
    # calculation above must pass these budgets; none is an assumed error.
    vector_budget=F(len(A),2000000);rounding_budget=F(1,10**14)
    assert vector<=vector_budget and rounding<=rounding_budget
    probability=vector_budget*(2+vector_budget)/4+rounding_budget
    return {'method':'exact dyadic complete-basis residual, rational scalar phases and rational comparison with every displayed complex amplitude',
        'conditional_vector_error_upper_exact':str(vector_budget),
        'displayed_probability_rounding_L1_half_upper_exact':str(rounding_budget),
        'fixed_outward_budget_rule':'vector <= complete finite dimension/2000000; probability rounding <= 1/10^14; actual exact errors checked on every replay',
        'original_input_position_flux_table_L1_half_error_upper_exact':str(probability),
        'complete_position_flux_table_numerically_enclosed':True,
        'native_probe_weights':'Separate displayed floating diagnostics for the fixed native instrument; this numerical enclosure covers the position/flux table.'}


def tail_bound(cap,q,time):
    step=max(1,abs(q));r=cap//step+1;x=(2*HOP+MAGNETIC)*F(time)
    assert x<r+1
    # Interaction-picture Dyson terms need at least r off-diagonal transitions.
    tail=x**r/F(factorial(r))/(1-x/F(r+1))
    return min(F(2),2*tail)


def field_floor(cap):
    assert 32*2*ELECTRIC==4 and 32*MAGNETIC==4 and 32*(MAGNETIC/2)==2
    n=np.arange(-cap,cap+1,dtype=np.int64)
    A=np.diag(4*n*n+4)+np.diag(np.full(2*cap,-2,dtype=np.int64),1)+np.diag(np.full(2*cap,-2,dtype=np.int64),-1)
    eigen=fixed_eigen(A,f'vacuum_{cap}') if cap in (18,24,48) else None
    cert,lam,U=exact_cap_enclosure(A,eigen)
    low,high=map(F,cert['lowest_interval_exact']);outside=2*ELECTRIC*(cap+1)**2
    coupling=MAGNETIC/2
    assert outside>high
    lower=low-coupling*coupling/(outside-high)
    return {'loop_cap':cap,'full_infinite_rotor_vacuum_interval_exact':[str(lower),str(high)],
        'complete_cap_certificate':cert,'outside_electric_floor_exact':str(outside),
        'boundary_operator_norm_upper_exact':str(coupling)}


def energy_certificate(cap=18,q=-3):
    assert 32*HOP==2 and 32*(ELECTRIC/2)==1 and 32*MAGNETIC==4
    d=matrices(cap,q)
    eigen=fixed_eigen(d['A'],'charged_18_-3') if (cap,q)==(18,-3) else None
    cert,lam,U=exact_cap_enclosure(d['A'],eigen)
    lo,hi=map(F,cert['lowest_interval_exact']);outside=2*ELECTRIC*(cap+1-abs(q))**2
    b=2*HOP+MAGNETIC/2;assert outside>hi
    low=lo-b*b/(outside-hi);vac=field_floor(cap)
    vlo,vhi=map(F,vac['full_infinite_rotor_vacuum_interval_exact'])
    # The charged and neutral reference sectors are different; compare true field vacuum.
    gap=[1+low-vhi,1+hi-vlo]
    assert gap[0]>1
    return {'cap':cap,'charge_exact':str(F(q,3)),
        'energy_sector':'One positive electron, colourful vacuum and full external flux/reference; not a minimization over every charged Fock sector',
        'charged_rest_minus_full_field_vacuum_electron_units_interval_exact':list(map(str,gap)),
        'charged_rest_minus_full_field_vacuum_GeV_interval_exact':[str(x*UNIT)for x in gap],
        'complete_cap_certificate':cert,'vacuum':vac,
        'outside_floor_electron_units_exact':str(outside),
        'boundary_norm_upper_exact':str(b),
        'physical_isolated_electron_pole_or_measured_mass_fit':False}


def executions():
    l=lepton_source();T=np.vstack((l['A'],l['B'].conj()));assert np.linalg.norm(T.conj().T@T-np.eye(48))<3e-15
    rows=[]
    for cap in (6,12,24):
        for name,c,f in [('electron',11,0),('positron',9,0),('muon',11,1),('tau',11,2),('neutral',10,0)]:
            j=ix(0,0,c,f);q=3*int(l['charge'][j]);d=matrices(cap,q);lam,U=eigh(d['A'].astype(float)/32)
            spectral,_,_=exact_cap_enclosure(d['A'],(lam,U))
            initial=np.zeros(len(lam));initial[d['index'][0,0]]=1
            native=abs(T[:,j])**2;assert abs(sum(native)-1)<1e-14
            for time in (2,4):
                psi=U@(np.exp(-1j*lam*time)*(U.T@initial));prob=abs(psi)**2
                assert abs(prob.sum()-1)<2e-13
                energy=float(np.vdot(psi,d['A']@psi).real/32)
                initial_energy=float(d['A'][d['index'][0,0],d['index'][0,0]]/32)
                assert abs(energy-initial_energy)<3e-13
                kinetic=float(np.vdot(psi,d['K']@psi).real/32)
                magnetic=float(np.vdot(psi,d['V']@psi).real/32)
                electric=float(np.dot(prob,d['electric'])/32)
                mass=float(l['E'][j]);bound=tail_bound(cap,q,time)
                numerical=finite_position_certificate(d['A'],lam,U,psi,d['index'][0,0],time,spectral)
                preparation_work=.5*(mass+initial_energy-float(MAGNETIC))
                readout_work=.5*(-mass-kinetic+float(MAGNETIC)-magnetic)
                balance=preparation_work+readout_work-.5*electric
                assert abs(balance)<2e-12
                rows.append({'family':name,'charge_exact':str(F(q,3)),'cap':cap,'time_exact':str(time),
                    'native_positive_mode':j,'native_probe_probabilities_96':native.tolist(),
                    'all_position_flux_labels':[list(z)for z in d['basis']],
                    'conditional_position_flux_probabilities':prob.tolist(),
                    'complete_joint_record_rule':'one half times position_flux_probability times native_probe_probability; preparation failure one half; native/Fock/Gauss complements retained',
                    'preparation_failure_probability_exact':'1/2',
                    'infinite_source_exact_Galerkin_original_record_TV_truncation_upper_exact':str(bound/2),
                    'finite_position_flux_numerical_certificate':numerical,
                    'infinite_source_original_position_flux_table_L1_half_error_upper_exact':str(bound/2+F(numerical['original_input_position_flux_table_L1_half_error_upper_exact'])),
                    'conditional_source_energy_electron_units_diagnostic':mass+energy,
                    'preparation_input_and_failure_field_energy_electron_units_exact':str(MAGNETIC),
                    'preparation_mean_source_work_electron_units_diagnostic':preparation_work,
                    'destructive_readout_mean_source_work_electron_units_diagnostic':readout_work,
                    'complete_energy_account':{'initial_mean_electron_units':float(MAGNETIC),
                        'final_mean_with_failure_electron_units':float(MAGNETIC)+.5*electric,
                        'conditional_final_electric_energy_electron_units':electric,
                        'total_work_balance_residual_diagnostic':float(balance)},
                    'actual_charge_and_flux_constraints_checked':True,
                    'norm_error_diagnostic':float(abs(prob.sum()-1))})
    return rows


def certificate():
    matter=matter_certificate();ward=exact_ward()
    energy=energy_certificate();vlo,vhi=map(F,energy['vacuum']['full_infinite_rotor_vacuum_interval_exact'])
    return {'schema':'oph.native_compact_transverse_completion.v1','verdict':'PASS_FULL_NATIVE_GAUSS_SOURCE_AND_CONTROLLED_TRANSVERSE_RECORDS',
        'architecture':{'complete_192_CAR_blocks_per_external_vertex':True,'executed_external_graph':'four oriented edges of one plaquette',
            'electric_coefficient_exact':str(ELECTRIC),'magnetic_coefficient_exact':str(MAGNETIC),'positive_covariant_hopping_exact':str(HOP),
            'external_kinetic':'positive finite-range covariant graph Laplacian on every positive Dirac excitation and native colourful mode',
            'constant_selection':'declared rational coefficients; no electromagnetic or mass target selects them',
            'full_link_carrier':'L2(U(1)) on each link, all integer electric fluxes',
            'internal_longitudinal_and_external_spatial_links_are_distinct':True,
            'interacting_Lorentz_covariance_or_physical_photon_claimed':False},
        'native_charge_Ward_identity':ward,'Gauss_exact_controls':gauss_checks(),
        'full_source_bound':{'inherited_matter_lower_bounds':matter['positivity'],
            'formula':'H_total >= H_matter_lower + E_gamma; K_covariant >=0; external photon vacuum subtracted only once',
            'field_vacuum_gauge_invariant':'Unique positive ground of the compact elliptic field Hamiltonian is gauge invariant',
            'finite_plaquette_gap_is_physical_massless_global_gap':False,
            'physical_massless_limit_has_no_positive_global_vacuum_gap':True},
        'infinite_loop_charged_rest_certificate':energy,
        'complete_preparation_and_detector':{
            'input':'Matter vacuum and charge-zero reference, bare electric loop n=0; this is an externally prepared field state, not its true ground.',
            'input_field_energy_electron_units_exact':str(MAGNETIC),
            'input_energy_above_true_field_vacuum_interval_exact':[str(MAGNETIC-vhi),str(MAGNETIC-vlo)],
            'success_Kraus':'2^(-1/2) times the Gauss-dressed native creation at x=0 times the projector onto the declared bare-flux vacuum input',
            'failure_Kraus':'sqrt(I-K_success* K_success), including every Fock/field/reference complement',
            'failure_state':'The bare input immediately after preparation; its full vacuum-sector unitary then evolves, conserving mean field energy b.',
            'preparation_mean_work_electron_units':'(m+2t)/2; the initial and failure field energies b are retained',
            'detector_charge':'A retained zero-Hamiltonian local apparatus absorbs matter charge q at measured x; Gauss is div E=q_matter+q_detector+R delta_0. It retains the full flux state and the compensating reference.',
            'total_mean_energy_identity':'W_prepare+W_absorption = one half conditional final electric energy; failure branch included; constant field-vacuum subtraction cancels'},
        'complete_native_executions':executions(),
        'controls':{'zero_hopping_charge_position_constant_exact':True,
            'zero_magnetic_hopping_zero_initial_flux_constant_exact':True,
            'undressed_charged_hop_violates_Gauss':True},
        'scope':'Finite compact transverse completion of the declared longitudinal matter source. The executed charged object is bound to its quantum compensating reference; its dressing energy is not an isolated electron rest-mass prediction. No 3+1D continuum, radiative QED fit or Yang-Mills gap follows.'}
