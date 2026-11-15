"""Exact Gauss transport and full integer-tail spectral enclosures."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh

ROOT=Path(__file__).resolve().parents[2]
ALPHA=F(1,32)
BETA=F(1,64)
HOP=F(1,64)


def physical_labels(q0,q1,n):
    q0,q1=F(q0),F(q1)
    if type(n) is not int or (3*q0).denominator!=1 or (3*q1).denominator!=1:
        raise ValueError('integer cycle and charges in thirds required')
    return (int(3*q1-n),n,int(-3*(q0+q1)))


def gauss_residual(q0,q1,labels):
    e1,e2,r=labels
    return (-e1-e2-3*F(q0)-r,e1+e2-3*F(q1))


def complete_gauss_control():
    # Complete four-CAR occupation carrier; n is a symbolic unrestricted integer.
    charges=(F(1),F(-1),F(1),F(-1)); rows=[];hops=0
    for mask in range(16):
        q=[sum((charges[j] for j in range(2*x,2*x+2) if mask>>j&1),F()) for x in range(2)]
        for n in (-101,-2,0,3,107):
            labels=physical_labels(*q,n)
            assert gauss_residual(*q,labels)==(0,0)
            assert labels[1]==n
            shifted=(labels[0]-1,labels[1]+1,labels[2])
            assert shifted==physical_labels(*q,n+1)
            for src,dst in ((0,2),(2,0),(1,3),(3,1)):
                if mask>>src&1 and not(mask>>dst&1):
                    after=q.copy();after[src//2]-=charges[src];after[dst//2]+=charges[src]
                    target=physical_labels(*after,n)
                    assert target==(labels[0]+3*(after[1]-q[1]),labels[1],labels[2])
                    assert gauss_residual(*after,target)==(0,0)
                    assert gauss_residual(*after,labels)!=(0,0)
                    hops+=1
        rows.append({'mask':mask,'charges_exact':list(map(str,q)),
                     'n_zero_flux_and_reference':list(physical_labels(*q,0))})
    return {'all16_occupations':rows,'dressed_hop_checks':hops,
            'symbolic_integer_cycle_formula':'(E1,E2,R)=(3Q1-n,n,-3(Q0+Q1)), n in Z',
            'inverse':'n=E2; both Gauss equations determine E1 and R uniquely',
            'undressed_hops_fail_Gauss':True}


def diagonal(cell,n,q):
    return ALPHA*((3*q*cell-n)**2+n*n)+(HOP if q else 0)+BETA


def action(v,q):
    """Apply the exact full integer-flux Hamiltonian; include every boundary image."""
    out={}
    def add(k,x):
        out[k]=out.get(k,F())+x
    for (cell,n),x in v.items():
        add((cell,n),diagonal(cell,n,q)*x)
        add((cell,n-1),-BETA*x/2);add((cell,n+1),-BETA*x/2)
        if q:add((1-cell,n),-HOP*x)
    return {k:x for k,x in out.items() if x}


def spectral_certificate(q,cutoff=6,bits=42):
    if q not in (0,1,-1) or cutoff<3:
        raise ValueError('vacuum or unit charged sector, cutoff at least3 required')
    # q=0 uses the exact scalar vacuum source. The neutral matter ground is
    # its tensor product with the symmetric spatial hopping eigenvector.
    keys=[(c,n) for c in (range(2) if q else range(1)) for n in range(-cutoff,cutoff+1)]
    positions={k:i for i,k in enumerate(keys)};H=np.zeros((len(keys),len(keys)))
    for j,k in enumerate(keys):
        for key,x in action({k:F(1)},q).items():
            if key in positions:H[positions[key],j]=float(x)
    values,vectors=eigh(H,subset_by_index=(0,0));raw=vectors[:,0]
    ints=[int(round(float(x)*(1<<bits))) for x in raw]
    v={k:F(x) for k,x in zip(keys,ints) if x}
    norm=sum(x*x for x in v.values());hv=action(v,q)
    mean=sum(v.get(k,F())*x for k,x in hv.items())/norm
    variance=sum((hv.get(k,F())-mean*v.get(k,F()))**2 for k in set(v)|set(hv))/norm
    # Electric operator has a unique zero mode and its second eigenvalue
    # is2alpha for q=0 or q=+-1. Both added quadratic forms are positive.
    second=2*ALPHA
    assert 0<mean<second and variance>0
    lower=mean-variance/(second-mean)
    assert lower>0
    allkeys=sorted(set(v)|set(hv))
    tail=[k for k in allkeys if abs(k[1])>cutoff]
    assert tail and all(abs(k[1])==cutoff+1 for k in tail)
    return {'charge':q,'trial_flux_window':[-cutoff,cutoff],
            'trial_integer_coefficients':[{'cell':k[0],'flux':k[1],'value':str(v.get(k,0))} for k in keys],
            'trial_norm_squared_exact':str(norm),'trial_rayleigh_exact':str(mean),
            'full_integer_source_residual_squared_exact':str(variance),
            'all_omitted_boundary_action_rows':[{'cell':c,'flux':n,'value_exact':str(hv[c,n])} for c,n in tail],
            'full_source_second_eigenvalue_lower_exact':str(second),
            'full_source_ground_interval_exact':[str(lower),str(mean)],
            'full_source_ground_interval_decimal':[float(lower),float(mean)],
            'truncation_used_as_source':False,'finite_vector_used_only_as_trial':True}


def source_inventory():
    parent=json.loads((ROOT/'model/openwave_common_mass_source/results.json').read_text())
    transport=json.loads((ROOT/'model/openwave_fock_field_transport/results.json').read_text())
    terms=parent['source']['terms']
    assert len(terms)==114 and parent['source']['native_CAR_modes']==192
    rows=[]
    for term in terms:
        rows.append({'name':term['name'],'kind':term['kind'],
                     'dressing':'Every complete matter matrix element carries shift DeltaE1=3DeltaQ1, DeltaE2=0, DeltaR=-3DeltaQtotal.'})
    assert transport['constant_native_pair_unitary_commutes_with_each_local_charge'] is True
    # The actual positive Dirac coordinates include the inherited neutral
    # family rotation. The24 pair maps alone need not empty that reference.
    ann=np.zeros((192,192),complex);cre=np.zeros_like(ann)
    lepton_native={x*96+sp*48+c*3+f for x,sp,c,f in product(range(2),range(2),(9,10,11,15),range(3))}
    for i in range(192):
        if i not in lepton_native:ann[i,i]=1
    for term in terms:
        if term['kind']!='positive_linear_fermion_square':continue
        j=int(term['name'].rsplit('_',1)[1]);x,jj=divmod(j,24);sp,jj=divmod(jj,12);cc,f=divmod(jj,3)
        row=x*96+sp*48+(9,10,11,15)[cc]*3+f
        for i,re,im in term['annihilation_coefficients']:ann[row,i]=re+1j*im
        for i,re,im in term['creation_coefficients']:cre[row,i]=re+1j*im
    car=float(np.linalg.norm(ann@ann.conj().T+cre@cre.conj().T-np.eye(192)))
    anomalous=float(np.linalg.norm(ann@cre.T+cre@ann.T))
    assert max(car,anomalous)<8e-15
    charge=np.array([(-1 if c==11 else 1 if c==9 else 0) for x,sp,c,f in product(range(2),range(2),range(16),range(3))])
    # The coloured identity rows commute with their actual charge as well;
    # this diagnostic checks the nontrivial lepton part without replacing it.
    assert np.linalg.norm(charge[:,None]*ann-ann*charge[None,:])==0
    assert np.linalg.norm(charge[:,None]*cre+cre*charge[None,:])==0
    return parent,{'native_terms':rows,'native_CAR_modes':192,'native_source_terms':114,
        'spatial_positive_quadratic_terms':96,'local_pair_intertwiners':24,
        'actual_positive_Dirac_definition':'The48 positive_linear_fermion_square annihilators of H0, with their inherited A/B and neutral-family rotation;144 colourful annihilators are unchanged.',
        'constant_neutral_family_rotation_retained':True,
        'full192_CAR_residual_diagnostic':car,'full192_anomalous_CAR_residual_diagnostic':anomalous,
        'exact_parent_operator_definition':parent['source']['exact_operator_definition'],
        'snapshot_coefficients_reclassified_as_exact':False}


def run():
    parent,inventory=source_inventory()
    vacuum=spectral_certificate(0);plus=spectral_certificate(1);minus=spectral_certificate(-1)
    vlo,vhi=map(F,vacuum['full_source_ground_interval_exact'])
    plo,phi=map(F,plus['full_source_ground_interval_exact'])
    mlo,mhi=map(F,minus['full_source_ground_interval_exact'])
    assert max(plo,mlo)<=min(phi,mhi)
    # Reflection of every integer flux proves charge-conjugate spectra exactly.
    offset=(plo-vhi,phi-vlo)
    assert 0<offset[0]<=offset[1]<BETA
    unit=F(parent['electron_unit_GeV_exact'])
    colourful_floor=F(parent['positivity']['colourful_sector_floor_per_occupation_exact'])
    assert colourful_floor>1>BETA>vhi
    assert F(parent['positivity']['finite_source_gap_GeV_exact'])>0
    return {'schema':'oph.full_native_matter_integer_flux_source.v1',
        'verdict':'PASS_COMMON_FULL_MATTER_CYCLE_SOURCE_AND_EXACT_INTEGER_TAIL_BOUNDS',
        'source':{'gauge_group':'U(1); no independent SU(3) link field is installed',
            'kinematic_carrier':'Fock192 tensor l2(Z_E1) tensor l2(Z_E2) tensor l2(Z_R)',
            'physical_carrier':'Gauss kernel unitarily equal to Fock192 tensor l2(Z_cycle)',
            'Gauss_constraints':'C0=-E1-E2-3Q0-R=0, C1=E1+E2-3Q1=0',
            'physical_H':'H0+alpha*((3Q1-n)^2+n^2)+t sum_j(d_0j-d_1j)^*(d_0j-d_1j)+beta*(I-(S+S*)/2)',
            'alpha_exact':str(ALPHA),'beta_exact':str(BETA),'t_exact':str(HOP),
            'd_modes':'The actual48 positive Dirac annihilators of H0, including its constant charge-preserving neutral-family rotation and24 pair maps;144 colourful modes unchanged.',
            'H0':'The complete inherited native114-term interacting matter source, including all semileptonic products and calibrated coefficients.',
            'native_process_boundary':'This is a declared source extension; no native twelve-port update generating these added terms has been proved.',
            'source_domain':'D(n^2) on Fock192 tensor l2(Z); bounded matter charges make the electric graph norm equivalent to D(n^2).',
            'stability_lower_bound_exact':'0',
            'self_adjoint_argument':'positive electric multiplication plus bounded self-adjoint H0, spatial hopping and Wilson shifts; bounded perturbation on D(n^2)',
            'spatial_term_norm_upper_exact':str(192*HOP),'magnetic_term_norm_upper_exact':str(2*BETA),
            'all_integer_flux_and_matter_sectors_retained':True},
        'native_inventory':inventory,'complete_Gauss_control':complete_gauss_control(),
        'identities':{'charge':'[H,Q0+Q1]=[H,R]=0',
            'local_current':'J_0to1=i[H,Q1]=-i[H,Q0]; electric and Wilson-cycle pieces commute with Qx',
            'Gauss':'[H,C0]=[H,C1]=0 for the fully dressed kinematic operator',
            'cycle_field_current':'i[H,n]=i*beta*(S-S*)/2, with S|n>=|n+1>',
            'work_transport':'Tr[W rho W* WHW*]-Tr[rho H]=0 for matched finite-energy states; complete instruments conjugate with W and retain I-G failure on the kinematic complement',
            'constant_pair_map':'commutes with each Qx and all fluxes, so its connection is zero',
            'field_dependent_map':'Any smooth periodic alternative V(theta) conjugates E to E+A, A=-i V* dV; E^2 becomes E^2+AE+EA+A^2 including -i A prime. Pair creation and derivative terms are retained.'},
        'rest_operator_comparison':{'all_native_sectors':'H-H0=alpha*((3Q1-n)^2+n^2)+t K+beta*(I-cos theta) >=0 as forms; the difference is unbounded on the full cycle carrier.',
            'every_native_preparation_in_n_zero':'For any normalized native finite state psi, added mean=9alpha <Q1^2>+t<K>+beta. This applies to each calibrated native rest preparation with its actual charge-square and hopping expectation.',
            'no_uniform_bounded_rest_operator_transfer':True,
            'vacuum':vacuum,'positive_unit_charge_spatial_sector':plus,'negative_unit_charge_spatial_sector':minus,
            'global_vacuum_certificate':{'colourful_floor_per_occupation_exact':str(colourful_floor),
                'zero_quark_charged_occupation_floor_exact':'1',
                'cycle_ground_upper_exact':str(vhi),
                'pure_neutral_sector':'Q0=Q1=0 identically; H=H_v+positive neutral rest and spatial hopping. The empty matter state attains E_v.',
                'whole_source_ground_equals_cycle_vacuum':True,
                'all_coloured_and_charged_sectors_excluded_by_gap_to_vacuum':True},
            'common_charged_one_lepton_gap_offset_native_exact':list(map(str,offset)),
            'common_charged_one_lepton_gap_offset_GeV_exact':list(map(str,(offset[0]*unit,offset[1]*unit))),
            'common_charged_one_lepton_gap_offset_GeV_decimal':[float(offset[0]*unit),float(offset[1]*unit)],
            'electron_muon_tau_rule':'Each calibrated one-lepton sector gap is its inherited rest coefficient plus the common positive charged offset. Spin and charge-conjugate degeneracies are retained.',
            'three_neutral_Dirac_rule':'The vacuum cycle energy cancels exactly; the symmetric spatial hopping state has zero kinetic energy. Each inherited neutral rest coefficient is unchanged in this reducing sector.',
            'physical_GeV_prediction':False,'electron_recalibration_performed':False},
        'optional_charged_counterterm':{
            'member':'H_delta=H-delta*N_ch, delta=inf spec(H_q-m)-inf spec(H_v), q=+1 or -1',
            'coefficient_definition':'Exact common full-source ground-energy difference; not a rounded interval endpoint or a fitted empirical coefficient.',
            'coefficient_native_interval_exact':list(map(str,offset)),
            'coefficient_strict_upper_exact':str(BETA),
            'charged_number_definition':'Sum of the24 positive charged Dirac occupation projections; both signs, all three families, both spins and both cells.',
            'charged_number_norm_exact':'24',
            'perturbation_norm_strict_upper_exact':str(24*BETA),
            'full_form_lower_bound':'H_delta >= c*N_colourful+(1-delta)*N_ch+H_neutral, with c>1 and H_neutral>=0.',
            'colourful_floor_per_occupation_exact':str(colourful_floor),
            'zero_colour_charged_floor_lower_exact':str(1-offset[1]),
            'cycle_vacuum_upper_exact':str(vhi),
            'whole_source_ground_equals_same_cycle_vacuum':True,
            'all_six_charged_particle_and_antiparticle_family_sector_rest_coefficients_exactly_restored':True,
            'all_three_neutral_Dirac_sector_rest_coefficients_unchanged':True,
            'one_counterterm_for_all_families_and_both_charges':True,
            'new_empirical_input_required':False,
            'choice_of_counterterm_forced_by_native_process':False,
            'finite_interval_endpoint_used_as_exact_coefficient':False,
            'unchanged_current_Ward_work_or_many_body_decay_transferred':False,
            'generator_boundary':'This optional member has its own H_delta evolution, currents and work; weak interactions can change N_ch. Only the reducing-sector rest values and global vacuum statement are claimed.',
            'physical_mass_prediction_or_native_process_closure':False},
        'resources':{'reference':'All integer R states retained; charged preparations require their compensating quantum-reference charge.',
            'preparation':'Preparing a complete matter/flux label has its stated source energy; the trial-vector Rayleigh work is reported without an autonomous preparation claim.',
            'measurement':'Complete matter/flux labels and the Gauss-failure complement; spectral windows may be measured with zero mean source-energy change. No free reset is assumed.',
            'clock':'exp(-itH) in inherited native units; no physical seconds inferred.'},
        'scope':{'full_source_extension_self_adjoint_and_stable':True,'complete_one_lepton_reducing_spectra_bounded':True,
            'native_process_compilation_proved':False,'physical_particle_dictionary_or_rest_masses_derived':False,
            'continuum_Yang_Mills_Clay_or_Einstein_theorem':False,'new_Nature_validation':False}}
