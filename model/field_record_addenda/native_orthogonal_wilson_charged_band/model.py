"""Exact bounds and complete charged records for the declared Gamma cube."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import argparse, hashlib, importlib.util, json, sys
import sympy as s

HERE=Path(__file__).resolve().parent; WORK=HERE.parent
RER=HERE.parents[2]
from model.native_wilson_maxwell import field
from model.native_wilson_maxwell_refinement.bounds import root,rounded
from model.native_block_dirac.algebra import annihilator
from model.affine_gauss_charged_clock.model import tail_bound
from model.particle_common_contract.provenance import import_closure
BRIDGE=WORK/'native_orthogonal_wilson_bridge'
BRIDGE_SHA='75590a8e62090717d5b91ada6ac472498952c3b4ace19b15de815419ecdadc6c'
SOURCE='model/native_wilson_maxwell/results.json'
SOURCE_SHA='eba4963148b62f6841a038119ed5b9baf4a6f142009915b31e49d72dacd12b72'
CLOCK='model/affine_gauss_charged_clock/results.json'
CLOCK_SHA='7c8f0214f1412e38bdc410c607d0a8f8765df323c3a33a778da923f3892e8613'

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def enc(x,upper=True):return str(rounded(x,upper))
def interval(lo,hi):return [enc(lo,False),enc(hi)]

@lru_cache(None)
def admit():
    assert digest(BRIDGE/'research_reference.json')==BRIDGE_SHA
    bridge=json.loads((BRIDGE/'research_reference.json').read_text()); public={};research={}
    research['native_orthogonal_wilson_bridge/results.json']=BRIDGE_SHA
    for f,h in bridge['owned_source_sha256'].items():
        assert digest(BRIDGE/'retained'/(f+'.raw'))==h;research['native_orthogonal_wilson_bridge/'+f]=h
    for f,h in bridge['research_parent_source_sha256'].items():
        p=WORK/'native_wilson_cone_audit'/('research_reference.json' if f=='results.json' else 'retained/'+f+'.raw');assert digest(p)==h
        research['native_wilson_cone_audit/'+f]=h
    for path,sha in ((SOURCE,SOURCE_SHA),(CLOCK,CLOCK_SHA)):
        assert digest(RER/path)==sha;v=json.loads((RER/path).read_text());public[path]=sha
        for f,h in v['input_sha256'].items():
            assert digest(RER/f)==h
            assert f not in public or public[f]==h;public[f]=h
    for f,h in bridge['input_sha256'].items():
        assert digest(RER/f)==h
        assert f not in public or public[f]==h;public[f]=h
    seeds=set(public)|{'model/native_wilson_maxwell/field.py','model/native_wilson_maxwell_refinement/bounds.py',
        'model/native_block_dirac/algebra.py','model/affine_gauss_charged_clock/model.py',
        'model/affine_gauss_charged_clock/PROOF.md','model/particle_common_contract/provenance.py'}
    for f in import_closure(seeds):
        h=digest(RER/f);assert f not in public or public[f]==h;public[f]=h
    return bridge,json.loads((RER/SOURCE).read_text()),public,research

@lru_cache(None)
def native_algebra():
    admit()
    from ..native_orthogonal_wilson_bridge import algebra
    return algebra

@lru_cache(None)
def coulomb():
    G=field.cube()['G'];rows=[]
    for fractional in range(8):
        for positive in combinations([j for j in range(8) if j!=fractional],4):
            rho=s.Matrix([-1 if j==fractional else 2 if j in positive else -2 for j in range(8)])
            assert sum(rho)==1;ref=rho-s.eye(8)[:,0];v=(ref.T*G*ref)[0]/2
            rows.append((tuple(rho),v))
    neutral=[]
    for positive in combinations(range(8),4):
        rho=s.Matrix([2 if j in positive else -2 for j in range(8)])
        neutral.append((rho.T*G*rho)[0]/2)
    assert len(rows)==280 and max(v for _,v in rows)==s.Rational(67,8)
    assert len(neutral)==70 and max(neutral)==8
    root0=s.eye(8)[:,0];adj=field.cube()['edges'][0][1];d=s.eye(8)[:,adj]-root0
    moved=(d.T*G*d)[0]/2;assert moved==s.Rational(7,24)
    return {'Q_plus_one_vertex_count':280,'Q_minus_one_by_sign_reversal':True,
        'Q_plus_minus_one_Coulomb_over_m_g_squared_max':'67/8',
        'Q_zero_vertex_count':70,'Q_zero_Coulomb_over_m_g_squared_max':'8',
        'charged_maximizing_patterns':[list(map(int,rho)) for rho,v in rows if v==s.Rational(67,8)],
        'fixed_root_rotation_counterexample_Coulomb_over_m_g_squared':str(moved)}

@lru_cache(None)
def exact_controls():
    a=native_algebra();H,Js,minimum,powers,edges=a.cube()
    aa,beta=a.current.matrices();Gamma=aa[0]*aa[1]*aa[2]*beta;C0=s.I*beta*aa[1]
    for j in range(3):
        B=-(Gamma+s.I*aa[j])/2
        assert C0*B.conjugate()*C0.H==-B
        assert B*B.H*B==B and s.trace(B.H*B)==2
    assert C0*beta.conjugate()*C0.H==-beta and C0*Gamma.conjugate()*C0.H==-Gamma
    cut=field.cutoff(3);parity=s.diag(*[(-1)**sum(n) for n in cut['states']])
    assert parity*cut['H_no_vacuum_at_a1']*parity==cut['H_no_vacuum_at_a1']
    assert all(parity*x*parity==-x for x in cut['annihilators'])
    x=s.Symbol('x');q=s.zeros(4)
    for j in range(3):q[j,j+1]=q[j+1,j]=s.sqrt(j+1)
    assert s.expand(q.charpoly(x).as_expr()-(x**4-6*x*x+3))==0
    return {'actual_32_mode_characteristic_and_minimal_polynomial_replayed':True,
        'all_12_literal_current_and_rank_two_partial_isometry_checks':True,
        'charge_conjugation_bond_coefficient_identity':True,
        'full_56_state_photon_parity':True,'full_finite_coordinate_norm_squared':'(13+5sqrt(6))/48',
        'paired_actual_Q_plus_minus_one_low_spectra':True,
        'native_high_complement_charge_conjugation_claimed':False,
        'central_two_pi_action_on_fixed_charged_cluster':'-I',
        'full_double_cubic_action_with_fixed_root_claimed':False,
        'eight_states_identified_as_species_or_spin_multiplicity':False}

@lru_cache(None)
def pump_controls():
    # Four source modes, followed by two actual fermionic pointer modes.
    c=[annihilator(6,j) for j in range(6)];n=[x.H*x for x in c]
    Q=sum(n[:4],s.zeros(64))-2*s.eye(64)
    parity=s.diag(*[(-1)**j.bit_count() for j in range(64)])
    dual=Q-n[4]+n[5];initial=[j for j in range(16) if j.bit_count()==2]
    for sign,D,f in ((1,c[2].H,c[4]),(-1,c[0],c[5])):
        X=D*f.H-f*D.H
        assert X.H==-X and X**3==-X
        U=s.eye(64)+X+X*X
        assert U.H*U==s.eye(64) and U*dual==dual*U and U*parity==parity*U
        assert Q*D-D*Q==sign*D
        mode=4 if sign==1 else 5
        for j in initial:
            ket=s.eye(64)[:,j];out=U*ket
            success=s.diag(*[int((k>>mode)&1) for k in range(64)])*out
            predicted=D*ket
            target=s.zeros(64,1)
            for k in range(16):target[k|(1<<mode),0]=predicted[k,0]
            assert success==target
            empty=s.diag(*[int(((k>>mode)&1)==0) for k in range(64)])*out
            assert empty==(s.eye(64)-D.H*D)*ket
    return {'source_modes':4,'pointer_modes':2,'literal_full_CAR_dimension':64,
        'graded_order':'source before pointer','full_even_unitarity_and_joint_parity':True,
        'conserved_generators':['Q+R via physical Gauss compensation','Q-n_f_plus+n_f_minus'],
        'all_four_pointer_patterns_retained':True,'neutral_even_input_Kraus':['D_q','I-D_q^dagger D_q']}


def atan_series(z,terms=64):
    z=F(z);assert abs(z)<1
    term=z;total=F(0)
    for k in range(terms):total+=term/(2*k+1);term*=-z*z
    error=abs(z)**(2*terms+1)/(2*terms+1)
    return rounded(total-error,False,65),rounded(total+error,True,65)

@lru_cache(None)
def pi_interval():
    a,b=atan_series(F(1,5)),atan_series(F(1,239))
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def phase_interval(xlo,xhi):
    # Monotonicity plus rational argument reduction, no floating transcendental.
    p=pi_interval();lo=rounded(xlo,False,20);hi=rounded(xhi,True,20)
    def atan(v):
        if abs(v)<=F(1,2):return atan_series(v)
        assert v>0;z=(v-1)/(v+1);a=atan_series(z)
        return p[0]/4+a[0],p[1]/4+a[1]
    al,ah=atan(lo),atan(hi)
    return al[0]/p[1],ah[1]/p[0]


def log_interval(x,terms=100):
    x=F(x);assert x>0;k=0
    while x>=2:x/=2;k+=1
    while x<1:x*=2;k-=1
    def local(y):
        z=(y-1)/(y+1);v=2*sum((z**(2*j+1)/(2*j+1) for j in range(terms)),F(0))
        return v,v+2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    a,b=local(x),local(F(2))
    return (a[0]+k*b[0],a[1]+k*b[1]) if k>=0 else (a[0]+k*b[1],a[1]+k*b[0])

@lru_cache(None)
def calculate():
    bridge,source,_,_=admit();g=F(1,64);c=F(23,2)
    r3=(root(3,False),root(3));r7=(root(7,False),root(7))
    levels=[(root(5,False)/2,root(5)/2),
        (root(17-4*r3[1],False)/2,root(17-4*r3[0])/2),
        (root(13,False)/2,root(13)/2),
        (root(33-4*r7[1],False)/2,root(33-4*r7[0])/2),
        (root(17+4*r3[0],False)/2,root(17+4*r3[1])/2),
        (root(33+4*r7[0],False)/2,root(33+4*r7[1])/2)]
    assert all(levels[j][1]<levels[j+1][0] for j in range(5))
    l1,l3,l4=levels[0],levels[2],levels[3]
    assert l1[0]+2>l4[1] and 3*l1[0]>l4[1]
    q=root((13+5*root(6))/48);nu0=24*g*q+8*g*g;nu1=24*g*q+F(67,8)*g*g
    deltaL,deltaU=map(F,bridge['energy_alignment']['signed_code_vacuum_shift_over_m'])
    gap=l4[0]-l3[1]-2*nu1;comp=c-deltaU-l3[1]-nu1
    assert gap>0 and comp>0 and levels[1][1]-l1[0]-2*nu1<0
    assert c-deltaU>l4[0]-nu1 and c-deltaU>2-nu0
    theta=nu0/(2-nu0);F0=1-theta*theta;eta0=root(2-2*root(F0,False))
    kappa=nu1/(l4[0]-nu1-l1[1]);A=root(1-kappa*kappa,False)-eta0
    odds=(1-F0)/F0;eps=F(3,128);r1=odds*odds+2*eps*eps*odds
    Fprep=1/(1+r1);etap=root(2-2*root(Fprep,False));a=A-etap
    pcode=F(source['actual_preparation']['global_code_probability_lower']);qprep=pcode*pcode*F0*F0
    assert a>0 and 4*(2-2*nu0)**2>F(16,9)*(1+nu0**2)*(1+(2-nu0)**2)
    assert 4>F(16,9)*(1+nu0**2)
    chargedlow=phase_interval(l1[0]-nu1,l1[1]-nu1)
    chargedhigh=phase_interval(l3[0]+nu1,l3[1]+nu1)
    highcomplement=phase_interval(l4[0]-nu1,l4[1]-nu1)
    neutralground=phase_interval(nu0,nu0);neutralexcited=phase_interval(2-nu0,2-nu0)
    assert chargedlow[0]>F(1,8)+F(1,4096)
    assert chargedhigh[1]<F(1464,4096) and highcomplement[0]>F(1466,4096)
    assert neutralground[1]<F(1,8)-F(1,64) and neutralexcited[0]>F(1,8)+F(1,64)
    beta=tail_bound(2**24,F(1,4096));assert beta==tail_bound(2**18,F(1,64))
    t=(1-beta)**2*a*a;success=qprep*t;fidelity=t/(t+beta*(1-a*a))
    massL,massU=map(F,source['clock_and_capacity']['source_selected_excitation_interval'])
    B=32*(600+360*F(1,2**34))+massU*(173+18*g+F(464,3)*g*g+8*c)+128
    coeffs={2:2*32*256+53,5:2*32*8,119:2*32*24,2**31+1:2*32,9:14,65:2,56:2}
    capacity=(sum(k*log_interval(x)[0] for x,k in coeffs.items()),sum(k*log_interval(x)[1] for x,k in coeffs.items()))
    lower={'charged_cluster_gap_over_m':gap,'charged_native_complement_margin_over_m':comp,
        'neutral_code_gap_over_m':2-2*nu0,'neutral_native_complement_margin_over_m':c-deltaU-nu0,
        'free_vacuum_overlap_squared':F0,'actual_ground_probe_cluster_amplitude':A,
        'actual_ground_probe_cluster_weight':A*A,'prepared_pair_ground_fidelity':Fprep,
        'original_preparation_success':qprep,'prepared_joint_cluster_ground_amplitude':a,
        'original_accepted_energy_pair_probability':success,'original_each_signed_accepted_probability':success/2,
        'accepted_joint_cluster_ground_fidelity':fidelity}
    upper={'finite_link_coordinate_norm':q,'neutral_perturbation_norm_over_m':nu0,
        'charged_perturbation_norm_over_m':nu1,'neutral_ground_free_vector_error':eta0,
        'free_charged_trial_outside_cluster_norm':kappa,'prepared_pair_purification_error':etap,
        'each_clock_false_positive_and_false_negative':beta,'full_raw_source_norm':B,
        'each_declared_control_norm':2*B+pi_interval()[1],'gross_switching_work':10*(2*B+pi_interval()[1])}
    exact={k+'_lower':enc(v,False) for k,v in lower.items()}
    exact.update({k+'_upper':enc(v) for k,v in upper.items()})
    phases={'charged_cluster_lower':interval(*chargedlow),'charged_cluster_upper':interval(*chargedhigh),
        'higher_charged_lower':interval(*highcomplement),'neutral_ground_max_circle_distance':interval(*neutralground),
        'neutral_excited_lower':interval(*neutralexcited),'pi':interval(*pi_interval())}
    return {'schema':'oph.research.orthogonal_wilson_charged_band.v1',
        'source':{'frozen_bridge_sha256':BRIDGE_SHA,'g':'1/64','photon_total_cap':3,
            'positive_complement_c':'23/2','only_changed_source_term':'m*(23/2-21/2)*sum_x(I-P_x)',
            'original_Psi_rho_Y_frame_retained':True,'old_beta_sea_raw_scalar_retained':True,
            'neutral_active_ground_unchanged_from_bridge':True,'full_native_and_outer_Gauss_complements_retained':True,
            'static_R_energy':'4|R|','fixed_charged_sector_R':'R=-Q',
            'active_cluster_rank_each_Q_plus_minus_one':8,'ranks_include_zero_energy_spectators':False},
        'exact_source_controls':exact_controls(),'fixed_charge_Coulomb':coulomb(),'even_pointer_pump':pump_controls(),
        'free_positive_levels_over_m':[interval(*v) for v in levels],
        'positive_level_multiplicities':[2,4,2,2,4,2],
        'rank_two_isolation_from_this_norm_bound':False,
        'reference_corrected_charged_transition_over_m':interval(l1[0]-nu1-nu0,l3[1]+nu1+nu0),
        'bounds_exact':exact,'diagnostic_binary64':{k:float(F(v)) for k,v in exact.items()},
        'phase_enclosures':phases,
        'complete_instrument':{'native_trial_copies':2,'global_code_flags':2,'all_code_patterns':4,
            'energy_agreement_size':64,'energy_agreement_success_label':0,
            'signed_setting_probability':'1/2','fermionic_pointer_dimension':4,
            'charged_pointer_size':2**24,'charged_accept_fraction':['1/8','1465/4096'],
            'neutral_pointer_size':2**18,'neutral_accept_circle_radius':'1/8',
            'all_record_labels':2**53,'all_failure_and_unused_outcomes_retained':True,
            'ideal_interacting_ground_projector_supplied':False,
            'accepted_good_subspace':'charged P_q on main source tensor actual neutral ground on retained partner'},
        'clock_and_capacity':{'source_time':'5','duration_one_driven_slots':5,
            'full_background_cancellation_in_every_slot_and_failure':True,
            'sources_and_partners_retained_for_all_slots':True,'finite_exact_matrix_services':True,
            'elementary_native_gate_synthesis_or_full_source_matrix_executed':False,
            'capacity_formula_nats':'2*[32*(256log2+8log5+24log119+log(2^31+1))+7log9+log65+log56]+53log2',
            'total_log_Hilbert_dimension_nats':interval(*capacity),
            'coefficient_and_synthesis_workspace_included':False},
        'boundaries':{'fixed_charge_sector_cluster_not_all_charge_isolation':True,
            'relativistic_rest_mass_or_SM_identity':False,'continuous_spin_or_Lorentz_invariance':False,
            'infinite_volume_or_photon_cutoff_limit':False,'exact_eight_state_eigenvectors_executed':False}}


def certificate():
    value=calculate();_,_,public,research=admit();value['input_sha256']=public;value['research_input_sha256']=research
    value['owned_source_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}
    return value

if __name__=='__main__':
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');args=parser.parse_args()
    value=certificate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    else:assert value==json.loads(path.read_text())
    print('PASS: full fixed-charge rank8 cluster, actual-ground residue, complete finite two-clock records')
