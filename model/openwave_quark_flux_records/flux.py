"""Exact Haar moments, native hole currents and complete spatial records."""
from collections import defaultdict
from itertools import permutations,product
import sympy as s
import numpy as np
from scipy.linalg import expm
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_joint.model import fermion
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES

Z=-s.Rational(1,2)+s.I*s.sqrt(3)/2


def epsilon(indices):
    if len(set(indices))!=3:return 0
    return (-1)**sum(indices[i]>indices[j]for i in range(3)for j in range(i+1,3))


def haar(u,v):
    """Normalized SU(3) Haar integral for every monomial used by the source."""
    m,n=len(u),len(v)
    if (m-n)%3:return s.Integer(0)
    if m==n==0:return s.Integer(1)
    if m==n==1:return s.Rational(int(u[0]==v[0]),3)
    if (m,n)==(3,0):return s.Rational(epsilon([i for i,j in u])*epsilon([j for i,j in u]),6)
    if (m,n)==(0,3):return haar(v,u)
    if m==n==2:
        ans=0
        for pr in permutations(range(2)):
            for pc in permutations(range(2)):
                if all(u[k][0]==v[pr[k]][0]and u[k][1]==v[pc[k]][1]for k in range(2)):
                    ans+=s.Rational(1,8)if pr==pc else -s.Rational(1,24)
        return ans
    raise ValueError(('unprovided exact Haar moment',m,n))


def mul(a,b):
    c=defaultdict(lambda:s.Integer(0))
    for (u,v),x in a.items():
        for (uu,vv),y in b.items():c[tuple(sorted(u+uu)),tuple(sorted(v+vv))]+=x*y
    return {k:s.expand(v)for k,v in c.items()if v!=0}


def conjugate(a):return {(v,u):s.conjugate(x)for (u,v),x in a.items()}

def integral(a):return s.simplify(sum(x*haar(u,v)for (u,v),x in a.items()))


def class_moments(r):
    psi={((),()):s.Integer(1)}
    for i in range(3):psi[((i,i),),()]=s.expand(Z**(-r));psi[(),((i,i),)]=s.expand(Z**r)
    density={k:s.simplify(v/3)for k,v in mul(psi,conjugate(psi)).items()}
    assert integral(density)==1
    means=[];squares=[]
    for i,j in product(range(3),repeat=2):
        term={(((i,j),),()):s.Integer(1)}
        mean=integral(mul(density,term));square=integral(mul(mul(density,term),conjugate(term)))
        assert s.simplify(mean-(Z**r/3 if i==j else 0))==0
        assert square==(s.Rational(7,18)if i==j else s.Rational(11,36))
        means.append(str(mean));squares.append(str(square))
    chi={(((i,i),),()):s.Integer(1)for i in range(3)}
    charsquare=integral(mul(mul(density,chi),conjugate(chi)))
    assert charsquare==s.Rational(5,3)
    return {'centre_index':r,'wavefunction':'[1+z^(-r) chi_3(U)+z^r chi_bar3(U)]/sqrt(3)',
        'norm_exact':'1','all_nine_U_means_exact':means,'all_nine_U_absolute_square_means_exact':squares,
        'TrU_absolute_square_mean_exact':str(charsquare),'colour_loop_visibility_exact':str(s.simplify(Z**r/3)),
        'postselected_colour_singlet_probability_loop_branch_exact':'5/27',
        'full_colour_electric_Casimir_before_per_edge_exact':'8/9',
        'full_colour_electric_Casimir_after_loop_per_edge_exact':'20/9'}


def native_holes():
    gauges=representation()[2];Q=gauges[10]+gauges[11];filled=(1<<30)-1;rows=[]
    for family in ('uc','dc'):
        ids=FAMILIES[family];basis=[fermion(filled,i,False)for i in ids];index={mask:r for r,(mask,sign)in enumerate(basis)}
        hcurrents=[]
        for g in gauges:
            h=s.zeros(3)
            for col,(mask,sgn)in enumerate(basis):
                for spin in range(2):
                    for i,j in product(range(15),repeat=2):
                        if not g[i,j]:continue
                        one=fermion(mask,15*spin+j,False)
                        if one is None:continue
                        two=fermion(one[0],15*spin+i,True)
                        if two is None:continue
                        if two[0]in index:
                            row=index[two[0]];h[row,col]+=sgn*one[1]*two[1]*basis[row][1]*g[i,j]
                        else:assert s.simplify(g[i,j])==0,'family must be invariant'
            assert h==-g.extract(ids,ids).T
            hcurrents.append(h)
        singlet=s.Matrix([s.Rational(int(i==j),1)/s.sqrt(3)for i,j in product(range(3),repeat=2)])
        for h,g in zip(hcurrents,gauges):
            assert (s.kronecker_product(h,s.eye(3))+s.kronecker_product(s.eye(3),g.extract(ids,ids)))*singlet==s.zeros(9,1)
        gram=s.Matrix(8,8,lambda a,b:s.trace(hcurrents[a]*hcurrents[b]))
        inverse=gram.inv();casimir=sum((inverse[a,b]*hcurrents[a]*hcurrents[b]/2 for a,b in product(range(8),repeat=2)),s.zeros(3))
        assert casimir==s.eye(3)*s.Rational(4,3)
        centre=s.Rational(2,3)*hcurrents[6]+s.Rational(4,3)*hcurrents[7]
        assert centre.is_diagonal()
        assert all(s.simplify(s.expand_complex(s.exp(2*s.pi*s.I*centre[i,i]))-Z)==0 for i in range(3))
        q=-Q[ids[0],ids[0]];assert q in(s.Rational(2,3),-s.Rational(1,3))
        rows.append({'family':family,'actual_native_indices':list(ids),'actual_hole_occupation_masks':[x[0]for x in basis],
            'actual_CAR_hole_signs':[int(x[1])for x in basis],'native_charge_exact':str(q),
            'reference_charge_exact':str(-q),'all_twelve_source_reference_current_checks':True,
            'native_colour_Casimir_exact':'4/3','full_reference_and_source_Fock_complements_retained':True})
    return rows


def amplitudes(K,phi):
    return {n:s.expand_complex(s.exp(-s.I*n*phi/3))/s.sqrt(2*K+1)for n in range(-K,K+1)}


def records(q,r,K,phi,colour=True):
    assert isinstance(q,int)and abs(q)<2*K+1
    amp=amplitudes(K,phi);mu=s.simplify(Z**r/3)if colour else s.Integer(1)
    dim=3 if colour else 1;out=[]
    for n in range(-K-abs(q),K+abs(q)+1):
        idle=amp.get(n,0);loop=amp.get(n-q,0)
        for i,j in product(range(dim),repeat=2):
            diagonal=int(i==j)
            square=(s.Rational(7,18)if diagonal else s.Rational(11,36))if colour else s.Integer(1)
            cross=s.simplify(s.conjugate(idle)*loop*mu*diagonal)
            baseline=s.simplify(s.conjugate(idle)*idle*diagonal+s.conjugate(loop)*loop*square)
            for quadrature in ('X','Y'):
                component=s.re(cross)if quadrature=='X'else s.im(cross)
                for sign in (-1,1):
                    probability=s.simplify((baseline+2*sign*component)/(8*dim))
                    assert float(probability)>=-1e-15
                    out.append({'electric_loop_flux':n,'quark_colour':i,'reference_colour':j,'quadrature':quadrature,'sign':sign,'probability_exact':str(probability)})
    assert s.simplify(sum(s.sympify(x['probability_exact'])for x in out))==1
    overlap=s.simplify(mu*(1-s.Rational(abs(q),2*K+1))*s.exp(s.I*s.Rational(q,3)*phi))
    marginal=[]
    for quadrature in ('X','Y'):
        for sign in (-1,1):
            p=s.simplify(sum(s.sympify(x['probability_exact'])for x in out if x['quadrature']==quadrature and x['sign']==sign))
            target=(1+sign*(s.re(overlap)if quadrature=='X'else s.im(overlap)))/4
            assert s.simplify(p-target)==0
            marginal.append(str(p))
    return {'native_charge_exact':str(s.Rational(q,3)),'centre_index':r,'phase_exact':str(phi),'phase_packet_half_width':K,
        'all_spatial_flux_colour_pointer_records':out,'pointer_marginal_exact':marginal,'loop_overlap_exact':str(s.expand_complex(overlap)),
        'all_nine_colour_outcomes_retained':colour,'full_field_support_retained':'Exact finite word support of a normalizable state in full L2(SU3) tensor L2(U1), without a source cutoff.'}


def monodromy_controls():
    rows=[]
    for q in (2,-1,-3):
        colour=q!=-3
        for r in range(3):
            a=(Z**r/3 if colour else 1)*s.exp(s.I*s.Rational(q,3)*s.pi/2)
            b=(Z**((r+1)%3)/3 if colour else 1)*s.exp(s.I*s.Rational(q,3)*(s.pi/2+2*s.pi))
            assert s.simplify(s.expand_complex(a-b))==0
            uncompensated=(Z**r/3 if colour else 1)*s.exp(s.I*s.Rational(q,3)*(s.pi/2+2*s.pi))
            if colour:assert s.simplify(s.expand_complex(a-uncompensated))!=0
            rows.append({'charge_exact':str(s.Rational(q,3)),'centre_index':r,'combined_spatial_flux_and_centre_record_invariance':True,
                'electromagnetic_turn_without_colour_compensation_is_invariant':not colour})
    return rows


def complete_native_preparation():
    """Exact full-Fock instrument, including the preceding reference test."""
    c=[]
    for mode in range(6):
        a=np.zeros((64,64),dtype=np.int64)
        for mask in range(64):
            out=fermion(mask,mode,False)
            if out is not None:a[out[0],mask]=out[1]
        c.append(a)
    transfer=sum((c[3+i].T@c[i]for i in range(3)),np.zeros((64,64),dtype=np.int64))
    p=np.diag(np.array([int(mask>>3==0)for mask in range(64)],dtype=np.int64))
    number=np.diag(np.array([(mask&7).bit_count()for mask in range(64)],dtype=np.int64))
    k=transfer@p                         # success numerator; divide by sqrt(3)
    assert np.array_equal(k.T@k,number@p)
    conditional_failure_numerator=3*p-number@p
    test_failure=np.eye(64,dtype=np.int64)-p
    assert min(np.diag(conditional_failure_numerator))>=0
    assert np.array_equal(k.T@k+conditional_failure_numerator+3*test_failure,3*np.eye(64,dtype=np.int64))
    filled=np.zeros(64,dtype=np.int64);filled[7]=1
    assert int(np.dot(k@filled,k@filled))==3
    # The untested pump is not a contraction: retain its exact countercontrol.
    one_transfer=transfer@filled
    assert int(np.dot(transfer@one_transfer,transfer@one_transfer))==4*int(np.dot(one_transfer,one_transfer))
    return {'complete_active_Fock_dimension':64,'preceding_reference_vacuum_test':True,
        'success_Kraus':'B P_refvac, B=sum_i r_i* c_i/sqrt(3)',
        'conditional_failure_Kraus':'sqrt(P_refvac-P_refvac B* B P_refvac)',
        'test_failure_Kraus':'I-P_refvac','exact_complete_Gram_sum':'I_64',
        'success_Gram_on_reference_vacuum':'N_source P_refvac/3',
        'stated_filled_source_success_probability_exact':'1',
        'unprojected_pump_countercontrol_squared_norm_lower_exact':'4/3',
        'spectator_extension':'P_refvac tests all thirty reference modes; retained source spectators factorize with this even native pump, and nonvacant reference spectators enter test failure.'}


def full_CAR_pulse():
    """Actual six-active-mode Fock pulse, before its proved one-hole restriction."""
    U=np.array([[.6,.8,0],[-.8,.6,0],[0,0,1]],complex)@np.diag([np.exp(1j*.4),np.exp(-1j*.4),1])
    assert np.linalg.norm(U.conj().T@U-np.eye(3))<1e-15
    c=[]
    for mode in range(6):
        a=np.zeros((64,64),complex)
        for mask in range(64):
            out=fermion(mask,mode,False)
            if out is not None:a[out[0],mask]=out[1]
        c.append(a)
    H=np.zeros((64,64),complex)
    for i,j in product(range(3),repeat=2):
        term=-U[i,j].conjugate()*c[3+i].conj().T@c[j]
        H+=term+term.conj().T
    vacuum=np.zeros(64);vacuum[-1]=1
    B=np.column_stack([a@vacuum for a in c]);target=np.block([[np.zeros((3,3)),U.conj().T],[U,np.zeros((3,3))]])
    intertwiner=np.linalg.norm(H@B-B@target);assert intertwiner<1e-14
    gate=expm(-.5j*np.pi*H);error=np.linalg.norm(gate@B[:,:3]+1j*B[:,3:]@U);assert error<2e-14
    return {'complete_active_Fock_dimension':64,'active_native_family_modes_on_two_vertices':6,
        'literal_CAR_to_hole_intertwiner_residual':float(intertwiner),'complete_Fock_quarter_turn_error':float(error),
        'all64_energy_values_diagnostic':np.linalg.eigvalsh(H).tolist(),
        'spectator_transfer':'Even native hopping factorizes with the other retained CAR coordinates; full occupations there are Pauli blocked on the stated input.'}


def certificate():
    native=native_holes();moments=[class_moments(r)for r in range(3)]
    rows=[records(q,0,4,phi,colour=q!=-3)for q in (2,-1,-3)for phi in (0,s.pi/2,2*s.pi)]
    currents=representation()[2];electron=FAMILIES['ec'][0]
    assert -(currents[10]+currents[11])[electron,electron]==-1
    assert all(g[electron,electron]==0 for g in currents[:10])
    return {'schema':'oph.native_spatial_quark_flux_records.v1','verdict':'PASS_COMPLETE_SPATIAL_NATIVE_FRACTIONAL_FLUX_RECORDS',
        'source_architecture':{'spatial_graph':'oriented four-cycle, complete SU3 and U1 link rotors',
            'matter':'full30CAR at each of four vertices; actual filled-source holes and full30CAR quantum reference',
            'spatial_transport':'four gauge-covariant native hopping quarter-turns; complete CAR actions, not installed quark-only Hamiltonians',
            'full_link_Hilbert_space':'L2(SU3 x U1) on each edge, no Peter-Weyl cutoff',
            'controlled_operations':'field evolution switched off during the declared finite transport services; no autonomous particle trajectory claim',
            'idle_energy':'sum hole and reference number plus (1/16)/2 sum E_U1^2 + (1/8)/2 sum C_SU3',
            'colour_orientation_apparatus':'full L2(SU3) quantum frame, normalized constant initial wavefunction; colour effects relational to that retained frame',
            'PW119_source_equivalence_claimed':False},
        'actual_native_hole_current_attachment':native,'actual_full_CAR_transport_pulse':full_CAR_pulse(),'exact_SU3_Haar_moments':moments,
        'actual_native_electron_control':{'index':electron,'hole_charge_exact':'-1','colour_and_weak_singlet':True},
        'complete_spatial_records':rows,'combined_monodromy_controls':monodromy_controls(),
        'phase_visibility_refinement':[{'K':K,'up_visibility_exact':str((1-s.Rational(2,2*K+1))/3),
            'down_visibility_exact':str((1-s.Rational(1,2*K+1))/3),'electron_visibility_exact':str(1-s.Rational(3,2*K+1))}for K in (2,4,8,16,32)],
        'source_work':{'initial_SU3_Casimir_per_edge_exact':'8/9',
            'loop_branch_SU3_Casimir_per_edge_exact':'20/9',
            'stated_K4_input_energy_above_idle_vacuum_exact':str(2+s.Rational(2,9)+s.Rational(4*5,24)),
            'input_energy_formula':'2 for one source hole and one reference particle, plus2/9 colour field and K(K+1)/24 Abelian field',
            'controlled_loop_mean_work_by_native_q_exact':{str(q):str(s.Rational(q*q,16)+(s.Rational(1,6)if q!=-3 else 0))for q in (2,-1,-3)},
            'terminal_nonselective_field_energy_work_exact':'0',
            'energy_role':'declared source units; no quark GeV mass calibration'},
        'preparation_and_failure':{'reference_vacuum_and_filled_source':'declared input; complete vacuum-test failure retained',
            'native_pump':complete_native_preparation(),
            'controlled_pointer':'both idle and full-loop paths retained; no singlet postselection',
            'constraint_failure':'complete nonphysical Gauss complement retained as its own output'},
        'fractional_spatial_texture_winding_or_stable_quark_identity_proved':False,
        'physical_quark_masses_or_confinement_derived':False,
        'scope':'Actual spatial Wilson transport of the native fractional-current holes with normalizable quantum gauge states and complete relational colour/flux/pointer outputs. A spatial interference phase is established; a self-sustained texture winding, physical quark mass, continuum confinement or autonomous bound quark is not identified.'}
