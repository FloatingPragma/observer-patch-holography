"""Native invariant tensors, retained neutral modes and charged-sector bounds."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import numpy as np
import sympy as s
from scipy.linalg import expm
from model.common_source import native_matter_attachment as matter
from model.common_source.chiral_candidate import kinetic
from model.lepton_mass_compiler import analyzer,compiler as native


@lru_cache(maxsize=1)
def native_tensors():
    rep=matter.representation();idx=list(rep['selected']);weak=rep['operators']['weak_C2'].extract(idx,idx)
    L=s.diag(*(int(weak[i,i]!=0)for i in range(15)));Q=(rep['operators']['T3']+rep['operators']['q']/6).extract(idx,idx)
    vertices={tuple(map(str,r['origin'])):r for r in matter.invariant_higgs_lines()};mapping=matter.action_intertwiner()[1]
    D=s.zeros(15);channels=[]
    for row in mapping:
        C=s.sympify(row['same_current_tensor_scale_exact'])*vertices[tuple(row['charge_origin'])]['components'][1]
        d=L*C*(s.eye(15)-L);assert d+d.T==C
        assert Q.T*C+C*Q==s.zeros(15)
        D+=d;channels.append({'name':row['existing_action_channel'],'D':d,'C':C,'origin':row['charge_origin'],'scale':row['same_current_tensor_scale_exact']})
    P=D.H*D+D*D.H;assert D**2==s.zeros(15) and P**2==P and P.rank()==14
    neutral=s.eye(15)-P;assert Q*neutral==s.zeros(15)
    assert P==s.diag(*(int(Q[i,i]!=0)for i in range(15)))
    assert all(abs(Q[i,i])<=1 for i in range(15))
    assert D.rank()==7
    # Every internal pair is disjoint. The actual native current attaches charge.
    pairs=[]
    for row in channels:
        for (i,j),v in row['D'].todok().items():
            assert v in(-1,1) and Q[i,i]==-Q[j,j]
            pairs.append({'channel':row['name'],'left':i,'right':j,'coefficient':str(v),
                'left_charge':str(Q[i,i]),'right_charge':str(Q[j,j])})
    assert len({i for r in pairs for i in(r['left'],r['right'])})==14
    return {'D':D,'Q':Q,'P':P,'channels':channels,'pairs':pairs,'neutral':neutral}


def exact_symbol():
    d=native_tensors();m=s.symbols('m0:3',positive=True);Y=s.diag(*m);C=d['D']+d['D'].T
    K=s.kronecker_product(C,Y);Q=s.kronecker_product(d['Q'],s.eye(3))
    assert K.H*K==s.kronecker_product(d['P'],Y**2)
    assert Q.T*K+K*Q==s.zeros(45)
    sig=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1));eps=s.Matrix([[0,1],[-1,0]])
    delta=s.kronecker_product(eps,K.conjugate());M=s.Matrix.vstack(s.Matrix.hstack(s.zeros(90),delta),s.Matrix.hstack(delta.H,s.zeros(90)))
    charge=s.diag(s.kronecker_product(s.eye(2),Q),-s.kronecker_product(s.eye(2),Q.T))
    assert M*charge==charge*M
    for sigma in sig:
        h=-s.kronecker_product(sigma,s.eye(45));A=s.diag(h,h.T)
        assert A*M+M*A==s.zeros(180)
    return {'native_internal_rank':14,'directed_native_pair_rank':7,'all_native_pairs':d['pairs'],
        'charged_internal_family_coefficients':42,'neutral_internal_family_coefficients':3,
        'exact_native_charge_conservation':True,'exact_three_axis_Clifford_mass_anticommutators':3,
        'exact_squared_mass':'P_charged tensor diag(e_k^2); neutral projector in kernel',
        'all_remaining_massless_coefficients_are_electrically_neutral':True,
        'complete_Nambu_dimension':180,'full_CAR_modes_per_cell':90}


def census():
    d=native_tensors();boxes=analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2'];values=np.array([float(sum(map(F,b),F(0))/2)for b in boxes])
    K=np.kron(np.asarray(d['D']+d['D'].T,float),np.diag(values));rows=[];ground=0.
    for address in product(range(4),repeat=3):
        p=np.array([(0,1,0,-1)[i]for i in address]);p2=int(p@p);H=kinetic.bdg_symbol(p,K,speed=1.)
        eig=np.linalg.eigvalsh(H);positive=np.r_[np.repeat(np.sqrt(p2+values**2),28),np.repeat(np.sqrt(p2),6)]
        # 42 internal charged entries, each contributes two positive Nambu energies.
        expected=np.sort(np.r_[-positive,positive]);assert len(expected)==180
        err=float(np.max(abs(eig-expected)));assert err<3e-13
        ground-=sum(positive)/2
        rows.append({'momentum_address':list(address),'kinetic_norm_squared':p2,'all_180_Nambu_energies':eig.tolist(),'spectral_error':err})
    return {'side':4,'sites':64,'physical_CAR_modes':5760,'rest_valleys':8,
        'charged_Dirac_channels_per_family_per_valley':7,'all_64_momentum_records':rows,
        'neutral_zero_energy_CAR_modes':48,'reference_degeneracy':'2^48, with every zero mode neutral',
        'chosen_reference':'negative bands filled, zero modes vacant; the same full quadratic action',
        'unshifted_source_ground_energy_diagnostic':float(ground),'charged_rest_coefficients_k0_k1_k2':values.tolist()}


def pair_records():
    d=native_tensors();boxes=analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2'];rows=[]
    c=[native._annihilator(2,i)for i in range(2)];N=[a.conj().T@a for a in c]
    for pair in d['pairs']:
        q=float(s.Rational(pair['left_charge']));sign=int(pair['coefficient']);Q=q*(N[0]-N[1]);V=sign*(c[0].conj().T@c[1].conj().T+c[1]@c[0])
        assert np.linalg.norm(V@Q-Q@V)<1e-14
        ev,vec=np.linalg.eigh(V);ground=vec[:,0];lo=(np.eye(4)-V@V);assert np.linalg.matrix_rank(lo)==2
        # Fixed one-leg creation is a contraction; its failure is retained.
        success=c[0].conj().T@ground;w=float(np.vdot(success,success).real);assert abs(w-.5)<1e-14
        for k,b in enumerate(boxes):
            e=float(sum(map(F,b),F(0))/2);t=.7/e;state=expm(-1j*e*V*t)@success
            energy=float(np.vdot(state,(e*V+e*np.eye(4))@state).real/w)
            charge=float(np.vdot(state,Q@state).real/w)
            assert abs(energy-e)<1e-14 and abs(charge-q)<1e-14
            probs=np.abs(state)**2
            assert abs(probs.sum()+1-w-1)<1e-14
            rows.append({'channel':pair['channel'],'native_left_index':pair['left'],'character':k,
                'source_gap':e,'charge':str(pair['left_charge']),'all_four_occupation_probabilities':probs.tolist(),
                'failure_probability':1-w,'compensating_reference_charge':str(-s.Rational(pair['left_charge'])),
                'idle_reference_energy':'0','source_time':t})
    return {'pair_orientation_records':rows,'spin_pair_count_per_native_pair':2,
        'other_spin_orientation':'same complete instrument under the native epsilon-spin CAR intertwiner',
        'remaining_carrier_and_failure_outputs_retained':True}


def countercontrols():
    d=native_tensors();off=next(r['D']for r in d['channels']if r['name']=='charged_lepton');old=off+off.T
    zeros=[i for i in range(15)if old[:,i]==s.zeros(15,1)and d['Q'][i,i]!=0]
    assert len(zeros)==12
    # Switching one actual channel off restores charged zero modes without
    # altering any current or charge convention.
    rows=[]
    for channel in ('up','down'):
        C=sum((r['C']for r in d['channels']if r['name']!=channel),s.zeros(15))
        charged_zero=[i for i in range(15)if C[:,i]==s.zeros(15,1)and d['Q'][i,i]!=0]
        assert len(charged_zero)==6;rows.append({'removed_channel':channel,'charged_zero_internal_indices':charged_zero})
    return {'charged_lepton_only_charged_zero_internal_count':len(zeros),'removed_channel_controls':rows,
        'shared_up_down_mass_spectrum_equal':True,'unequal_observed_quark_masses_derived':False}


def certificate():
    exact=exact_symbol();data=census();records=pair_records();masses=data['charged_rest_coefficients_k0_k1_k2'];light=masses[1];anchor=.00051099895069
    calibrated=[{'channel':name,'GeV_values_k0_k1_k2':[anchor*m/light for m in masses],
        'role':'conditional shared-family free-source energies; quark channels are not empirical mass fits'}for name in('charged_lepton','up','down')]
    return {'schema':'oph.full_native_Higgs_charged_floor.v1','exact_native_symbol':exact,'complete_spectrum':data,
        'complete_native_pair_instruments':records,'countercontrols':countercontrols(),
        'charged_sector_theorem':{'operator_bound':'H−E_ground >= e_min*abs(Q) on the complete Fock source',
            'charge_one_minimum':'e_min, attained at a charged-lepton rest mode',
            'colour_singlet_charge_one_quark_triple_energy':'3*e_min for three same-family rest quasiparticles',
            'quark_triple_is_bound_baryon':False,'all_neutral_zero_modes_retained':True},
        'calibrated_channel_values':calibrated,
        'finite_analyzer_implementation':{'physical_CAR_per_site':90,'enlarged_CAR_per_site':324000,
            'active_CAR_per_site':168,'disjoint_four_mode_factors_per_site':42,
            'complete_trace_error':'min(2,42*s*R^5*t^5/(10*n^4)), R=(7/3)(4/a+1)',
            'same_four_phase_and_inversion_weak_leg_parity_services':True},
        'architecture':{'native_Higgs_channels':['up','down','charged_lepton'],'same_existing_tensor_normalizations':True,
            'shared_family_matrix':'Y_u=Y_d=Y_e=diag(e_k), same declared analyzer',
            'neutral_background_norm':'1','coherent_phase':'2/9, retrospective input',
            'relative_channel_choice_source_forced':False,'physical_interacting_quark_or_lepton_identification':False},
        'scope':'An explicitly extended free source removes charged zero modes with native Higgs tensors, retaining neutral zeros and eight valleys. Shared quark/lepton spectra are conditional architecture outputs; they do not reproduce observed quark masses or select an interacting Higgs vacuum.',
        'verdict':'PASS_FULL_NATIVE_CHARGED_MASS_FLOOR'}
