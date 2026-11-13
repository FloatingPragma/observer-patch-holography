"""Exact source reductions, complete spectra and native preparation records."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import hashlib
import numpy as np
import scipy.sparse as sp
import sympy as sy
from scipy.linalg import eigh
from model.openwave_six_quark import source as native6
from model.openwave_six_quark.certificate import spectral_enclosure
from model.openwave_axial_state.certificate import exact_action,sqrt_upper,reduced_line
from model.openwave_calibrated_hadrons.source import angular4
from model.openwave_nuclear_distance import nuclear
from . import source as s

UNIT=F(s.INPUTS['electron_unit_GeV'])

def reduce_exact(B,Gram,parts,p):
    L,D=Gram.LDLdecomposition(hermitian=True);diag=list(D.diagonal());assert all(x>0 for x in diag)
    A=sy.zeros(B.shape[1])
    for k,w in s.weights(p).items():
        reduced=native6.exact_product(B.T,native6.exact_product(parts[k],B)).toarray()
        A+=sy.Rational(w.numerator,w.denominator)*sy.Matrix(reduced.tolist())
    Li=L.inv();T=Li*A*Li.T
    H=sy.Matrix(len(diag),len(diag),lambda i,j:T[i,j]/sy.sqrt(diag[i]*diag[j]))
    transform=np.asarray(Li.T.evalf(18),float)@np.diag([1/np.sqrt(float(x))for x in diag])
    return H,transform

@lru_cache(None)
def baryon_multiplet(S2,I2):
    d=s.components('baryon');eye=sp.eye(len(d['basis']),dtype=np.int64,format='csr')
    S4=angular4(d['modes'],d['basis'],(0,1),True);I4=angular4(d['modes'],d['basis'],(0,1),False)
    ids=[i for i,row in enumerate(d['labels'])if sum(f==0 for x,z,f in row)==1 and sum(z==0 for x,z,f in row)==2]
    a=(15*eye-S4)if S2==1 else(S4-3*eye)
    b=(15*eye-I4)if I2==1 else(I4-3*eye)
    B=native6.exact_product(a,native6.exact_product(b,d['Binteger'][:,ids]))
    G=sy.Matrix(native6.exact_product(B.T,B).toarray().tolist());_,cols=G.rref();B=B[:,list(cols)].tocsr()
    Gram=sy.Matrix(native6.exact_product(B.T,B).toarray().tolist())
    assert (S4@B-S2*(S2+2)*B).nnz==0 and(I4@B-I2*(I2+2)*B).nnz==0
    H,T=reduce_exact(B,Gram,d['parts'],s.parameters());cert=spectral_enclosure(H)
    delta=s.parameters()['delta'];center=[[str(F(a)-delta/2),str(F(b)-delta/2)]for a,b in cert['all_eigenvalue_intervals']]
    return {'twice_spin':S2,'twice_isospin':I2,'fixed_Ms_half_Mi_minus_half_dimension':B.shape[1],
            'full_magnetic_dimension':B.shape[1]*(S2+1)*(I2+1),'complete_reference_spectrum':cert,
            'all_isospin_center_energy_intervals':center,
            'isospin_magnetic_shift':'-delta*Mi; the reference spectrum uses Mi=-1/2'}

@lru_cache(None)
def baryon_spectra():
    rows=[baryon_multiplet(S,I)for S,I in product((1,3),repeat=2)]
    assert sum(x['full_magnetic_dimension']for x in rows)==120
    return rows

def two_baryon_threshold(S,I):
    choices=[]
    for a,b in product(baryon_spectra(),repeat=2):
        sa,sb=F(a['twice_spin'],2),F(b['twice_spin'],2);ia,ib=F(a['twice_isospin'],2),F(b['twice_isospin'],2)
        if not(abs(sa-sb)<=S<=sa+sb and abs(ia-ib)<=I<=ia+ib):continue
        ea=list(map(F,a['all_isospin_center_energy_intervals'][0]));eb=list(map(F,b['all_isospin_center_energy_intervals'][0]))
        choices.append({'fragment_twice_spin_isospin':[[int(2*sa),int(2*ia)],[int(2*sb),int(2*ib)]],
                        'interval':[str(ea[0]+eb[0]),str(ea[1]+eb[1])]})
    assert choices
    return {'complete_admitted_multiplet_pairs':choices,'threshold_interval':[str(min(F(x['interval'][0])for x in choices)),str(min(F(x['interval'][1])for x in choices))],
            'scope':'Two independent complete copies of this source, with compatible native total spin and isospin; no continuum scattering or orbital-parity identification.'}

@lru_cache(None)
def nuclear_exact(S,I,parity):
    d=nuclear.carrier();c=nuclear.channel(S,I,parity)
    if not c['rank']:return None
    G4=2*d['parts']['S4']-4*d['Sg'];parts={**d['parts'],'G4':G4}
    assert (native6.exact_product(G4,d['P'])-native6.exact_product(d['P'],G4)).nnz==0
    H,T=reduce_exact(c['Binteger'],c['Gram'],parts,s.parameters())
    return dict(H=H,transform=T,B=c['Binteger'],parts=parts,carrier=d,channel=c)

@lru_cache(None)
def nuclear_spectra():
    rows=[];rank=0;full=0
    for S,I,parity in product(range(4),range(4),(-1,1)):
        c=nuclear.channel(S,I,parity);rank+=c['rank'];full+=c['rank']*(2*S+1)*(2*I+1)
        if not c['rank']:continue
        d=nuclear_exact(S,I,parity);spec=spectral_enclosure(d['H']);threshold=two_baryon_threshold(S,I)
        lo,hi=map(F,spec['all_eigenvalue_intervals'][0]);tl,th=map(F,threshold['threshold_interval'])
        rows.append({'spin':S,'isospin':I,'cell_parity':parity,'fixed_Ms0_Mi0_dimension':c['rank'],
                     'full_magnetic_dimension':c['rank']*(2*S+1)*(2*I+1),'complete_spectrum':spec,
                     'ground_GeV_interval':[str(UNIT*lo),str(UNIT*hi)],'same_source_two_baryon_threshold':threshold,
                     'binding_GeV_interval':[str(UNIT*(tl-hi)),str(UNIT*(th-lo))],
                     'strictly_bound_in_this_native_channel':tl>hi,'above_this_channel_breakup_threshold':lo>th,
                     'all_Iz_energy_shifts':[{'twice_Iz':2*m,'shift':str(-s.parameters()['delta']*m)}for m in range(-I,I+1)]})
    assert rank==240 and full==2520
    return {'complete_fixed_Ms0_Mi0_rank':rank,'complete_singlet_dimension':full,'channels':rows}

@lru_cache(None)
def axial_ground():
    old=reduced_line();d=s.components('baryon');H,T=reduce_exact(old['B'],old['Gram'],d['parts'],s.parameters())
    e,U=eigh(np.asarray(H.evalf(18),float));coef=T@U[:,0];integers=[round(float(x)*2**45)for x in coef]
    v=exact_action(old['B'],integers);norm=sum(x*x for x in v);Hv=[F(0)for _ in v]
    for k,c in s.weights(s.parameters()).items():
        for i,x in enumerate(exact_action(d['parts'][k],v)):Hv[i]+=c*x
    mean=sum(x*y for x,y in zip(v,Hv))/norm;variance=sum((x-mean*y)**2 for x,y in zip(Hv,v))/norm
    spec=spectral_enclosure(H);gap=F(spec['all_eigenvalue_intervals'][1][0])-mean;assert gap>0
    sine=sqrt_upper(variance)/gap;current=F(sum(x*y for x,y in zip(v,exact_action(old['X'],v))),norm);error=12*sine
    assert abs(current-F(s.INPUTS['axial_calibration_target_magnitude']))+error<F(1,10**6)
    return {'exact_selected_dimension':H.rows,'complete_selected_spectrum':spec,'integer_native_trial':list(map(str,v)),
            'integer_trial_norm_squared':str(norm),'trial_mean':str(mean),'trial_variance':str(variance),'ground_vector_sine_bound':str(sine),
            'gV_exact':'1','gA_interval':[str(current-error),str(current+error)],'gA_interval_decimal':[float(current-error),float(current+error)],
            'current_unchanged':True,'scope':'Actual native current on the certified even spin/isospin-half ground preparation; axial magnitude is a disclosed input.'}

@lru_cache(None)
def deuteron_trial():
    d=nuclear_exact(1,0,1);H=d['H'];e,U=eigh(np.asarray(H.evalf(18),float));a=d['transform']@U[:,0]
    coeff=[round(float(x)*2**45)for x in a];v=exact_action(d['B'],coeff);norm=sum(x*x for x in v);hv=[F(0)for _ in v]
    for k,c in s.weights(s.parameters()).items():
        for i,x in enumerate(exact_action(d['parts'][k],v)):hv[i]+=c*x
    mean=sum(x*y for x,y in zip(v,hv))/norm;var=sum((x-mean*y)**2 for x,y in zip(hv,v))/norm
    spec=spectral_enclosure(H);gap=F(spec['all_eigenvalue_intervals'][1][0])-mean;assert gap>0
    sine=sqrt_upper(var)/gap;c=d['carrier'];N=len(c['src']);n0=c['n0'];quad=F(sum(x*x*int(n*(6-n))for x,n in zip(v,n0)),18*norm)
    nonzero=[{'full_input_index':int(c['keep'][i]),'native_mask':str(c['src'][int(c['keep'][i])]),'integer_amplitude':str(x)}for i,x in enumerate(v)if x]
    assert sum(int(x['integer_amplitude'])**2 for x in nonzero)==norm
    return {'exact_preparation_input_dimension':N,'quantum_numbers':{'spin':1,'isospin':0,'cell_parity':1,'Sz':0,'Iz':0},
            'nonzero_native_trial_amplitudes':nonzero,'integer_trial_norm_squared':str(norm),
            'trial_energy':str(mean),'trial_energy_GeV':str(mean*UNIT),'trial_residual_squared':str(var),'ground_vector_sine_bound':str(sine),
            'quadrupole_graph_squared_interval':[str(max(F(0),quad-sine)),str(quad+sine)],
            'complete_native_instrument':{'uniform_input_basis_rule':'all six occupied native modes among24 with exactly three spin-up and three up-flavour occupations; basis order is the inherited exact CAR carrier',
                'rank_one_filter_success':str(F(1,N)),'filter_failure':str(1-F(1,N)),
                'each_success_occupation_probability':'v_i^2/(norm*N), with v_i=0 at every unlisted full input index',
                'each_failure_occupation_probability':'1/N-v_i^2/(norm*N), including every unlisted occupation',
                'complete_probability_sum_exact':'1','other_Fock_sector_outcome':'identity on the full orthogonal complement, with zero probability for the declared input',
                'reference':'idle charge-minus-one reference, retained in every outcome','total_nonselective_filter_and_occupation_work_exact':'0 on the full uniform input',
                'ground_complete_channel_half_trace_error_upper':str(2*sine),'arbitrary_external_reference_transfer':True},
            'physical_length_calibration':None,'scope':'Spin-one native bound-channel trial; quadrupole uses graph separation squared, without a femtometre assignment.'}

def partitions(number,least=1):
    if number==0:yield ()
    for n in range(least,number+1):
        for rest in partitions(number-n,n):yield(n,)+rest

def fragment_bounds():
    p=s.parameters();trial=F(deuteron_trial()['trial_energy']);Cas=lambda n:F(4,3)if n%3 else F(0)
    controls=[]
    for kappa in(1,2,4):
        floor={n:min(p['A']*F(x*x+(n-x)**2,9)+p['U']*(Cas(x)+Cas(n-x))+kappa*Cas(n)for x in range(n+1))for n in(1,2,4,5)}
        rows=[]
        for part in partitions(6):
            if part in((6,),(3,3)):continue
            # The elementary number/triality floor remains valid for an N3 fragment too.
            f3=min(p['A']*F(x*x+(3-x)**2,9)+p['U']*(Cas(x)+Cas(3-x))for x in range(4))
            lower=6*p['mbar']+sum(f3 if n==3 else floor[n]for n in part)
            rows.append({'partition':list(part),'same_source_threshold_lower':str(lower),'margin_above_trial':str(lower-trial),'margin_GeV':str((lower-trial)*UNIT)})
        if kappa==4:assert all(F(x['margin_above_trial'])>0 for x in rows)
        controls.append({'global_colour_coefficient':kappa,'all_non_3plus3_partitions':rows,'minimum_margin':str(min(F(x['margin_above_trial'])for x in rows))})
    colour3=3*p['mbar']+min(p['A']+3*p['U'],F(5,9)*p['A']+F(8,3)*p['U'])+3*p['kappa_colour']
    singlet3=min(F(x['all_isospin_center_energy_intervals'][0][0])for x in baryon_spectra())
    assert min(2*colour3,colour3+singlet3)>trial
    return {'native_local_colour_minimum':'0 when local quark number is divisible by3; otherwise at least4/3. N<=5 is a finite native SU3 representation bound.',
            'positive_gradient_spin_weak_and_hopping_terms_dropped_for_lower_bounds':True,'all_nontrivial_partitions_retained':True,
            'global_colour_controls':controls,'coloured_three_quark_center_energy_lower':str(colour3),
            'every_3plus3_coloured_or_singlet_combination_margin_lower':str(min(2*colour3,colour3+singlet3,2*singlet3)-trial),
            'scope':'Every independent-copy quark-number partition of six on the identical two-cell action; no asymptotic scattering or spatial-confinement assertion.'}

def certificate():
    p=s.parameters();b=baryon_spectra();n=nuclear_spectra();a=axial_ground();trial=deuteron_trial();fragments=fragment_bounds()
    selected=baryon_multiplet(1,1)['complete_reference_spectrum']['all_eigenvalue_intervals'][0]
    neutron=list(map(F,selected));proton=[x-p['delta']for x in neutron];meson=[]
    H=s.exact_matrix('meson');d=s.components('meson');q=np.rint(np.diag(d['q']).real).astype(int)
    for charge in sorted(set(q)):
        ids=np.flatnonzero(q==charge).tolist();other=np.flatnonzero(q!=charge).tolist();assert all(H[i,j]==0 for i in ids for j in other)
        meson.append({'charge':int(charge),'complete_spectrum':spectral_enclosure(H.extract(ids,ids))})
    pion=list(map(F,next(x for x in meson if x['charge']==1)['complete_spectrum']['all_eigenvalue_intervals'][0]))
    fits=[]
    for name,interval in(('proton',proton),('neutron',neutron),('charged_pion',pion)):
        target=F(s.INPUTS['mass_calibration_targets_GeV'][name]);gev=[x*UNIT for x in interval];assert gev[0]<=target<=gev[1]
        fits.append({'state':name,'target_GeV':str(target),'energy_GeV_interval':list(map(str,gev)),'role':'empirical calibration input'})
    even=next(x for x in n['channels']if(x['spin'],x['isospin'],x['cell_parity'])==(1,0,1))
    odd=next(x for x in n['channels']if(x['spin'],x['isospin'],x['cell_parity'])==(0,0,-1))
    assert even['strictly_bound_in_this_native_channel'] and odd['above_this_channel_breakup_threshold']
    target=F(s.INPUTS['deuteron_binding_calibration_target_GeV']);assert F(even['binding_GeV_interval'][0])<=target<=F(even['binding_GeV_interval'][1])
    controls=[]
    for gradient in(0,1,2,4):
        pp=s.parameters(F(gradient));threshold=sum(float(s.INPUTS['mass_calibration_targets_GeV'][k])for k in('proton','neutron'))
        controls.append({'gradient':gradient,'parameters_rational':{k:str(v)for k,v in pp.items()},
                         'even_energy_GeV':s.nuclear_energy(pp)*float(UNIT),'odd_energy_GeV':s.nuclear_energy(pp,0,0,-1)*float(UNIT),
                         'odd_binding_to_pn_GeV':threshold-s.nuclear_energy(pp,0,0,-1)*float(UNIT),
                         'scope':'Five targets are recalibrated for each distinct diagnostic source; the architecture menu is retrospective.'})
    return {'schema':'oph.positive_native_spin_gradient_nuclei.v1','verdict':'PASS_COMPLETE_SPIN_GRADIENT_SOURCE_WITH_ODD_CHANNEL_ABOVE_BREAKUP',
            'input_ancestry':s.INPUTS,'parameters_rational':{k:str(v)for k,v in p.items()},'positive_gradient_CAR_identity_exact':True,
            'three_hadron_mass_calibrations':fits,'native_axial_ground':a,'complete_baryon_multiplets':b,'complete_meson_charge_spectra':meson,
            'complete_nuclear_spectrum':n,'selected_nuclear_record':trial,'all_fragment_comparisons':fragments,'retrospective_gradient_controls':controls,
            'unassigned_bound_native_channels':[{'spin':x['spin'],'isospin':x['isospin'],'cell_parity':x['cell_parity'],'binding_GeV_interval':x['binding_GeV_interval']}for x in n['channels']if x['strictly_bound_in_this_native_channel']and(x['spin'],x['isospin'],x['cell_parity'])!=(1,0,1)],
            'positive_full_Fock_floor':'H>=min(m_up,m_down)*N_colourful; every displayed current square and graph hopping is positive, and spectator/reference sectors are retained.',
            'scope':'Five empirical calibrations on one declared native action. The odd spin-zero channel is excluded below its spin/isospin-compatible two-baryon threshold. Other bound channels are retained; native cell parity is not continuum orbital parity, and no complete physical nuclear spectrum or continuum pole is asserted.'}
