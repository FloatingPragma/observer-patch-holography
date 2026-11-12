from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh
from model.openwave_six_quark.certificate import spectral_enclosure,sqrt_upper_fraction
from model.openwave_six_quark.source import exact_product
from model.openwave_calibrated_hadrons.calculation import spectra
from model.openwave_calibrated_hadrons.source import INPUTS,parameters
from .source import reduction,exact_matvec,rational_action,native


def partitions(n,least=1):
    if n==0:yield ();return
    for x in range(least,n+1):
        for rest in partitions(n-x,x):yield(x,)+rest


def trial(d,spectrum):
    e,V=eigh(d['h']);c=d['coefficient_transform']@V[:,0]
    integers=[int(round(x*2**52))for x in c];v=exact_matvec(d['Binteger'],integers);D=sum(x*x for x in v);Hv=rational_action(d,v)
    E=sum((x*y for x,y in zip(v,Hv)),F(0))/D
    rho2=sum(((y-E*x)**2 for x,y in zip(v,Hv)),F(0))/D;rho=sqrt_upper_fraction(rho2)
    distance=F(spectrum['all_eigenvalue_intervals'][1][0])-E;assert distance>0;trace=2*rho/distance
    probs=np.zeros(d['full_dimension']);quad=F(0)
    for i,j in enumerate(d['keep']):
        prob=F(v[i]**2,D);probs[j]=float(prob);mask=d['basis'][j]
        n0=sum(1 for k,m in enumerate(native.MODES)if m[0]==0 and(mask>>k)&1)
        quad+=prob*F(n0*(6-n0),18)
    assert 0<quad<=F(1,2)
    return {'integer_trial_coefficients':list(map(str,integers)),'exact_norm':str(D),'exact_energy':str(E),
       'energy':float(E),'exact_residual_squared':str(rho2),'ground_trace_distance_upper':str(trace),
       'quadrupole_trial':str(quad),'stationary_ground_quadrupole_interval':[str(quad-trace/2),str(min(F(1,2),quad+trace/2))],
       'complete_success_probabilities':(probs/d['full_dimension']).tolist(),
       'complete_failure_probabilities':((1-probs)/d['full_dimension']).tolist(),
       'success_probability':str(F(1,d['full_dimension'])),'full_Fock_complement_probability':0,
       'nonselective_filter_and_occupation_mean_work':0,'idle_reference_charge':-1,
       'record_boundary':'Exact rational trial filter on the complete11700-word uniform input. Trace comparison to the stationary ground is uniform under autonomous evolution; complete failures and every occupation are retained.'}


def baryon_threshold(total_charge):
    rows={x['charge']:x['certificate']['all_eigenvalue_intervals'][0]for x in spectra()['baryon']}
    candidates=[]
    for a in rows:
        b=total_charge-a
        if b in rows and a<=b:candidates.append((a,b,F(rows[a][0])+F(rows[b][0]),F(rows[a][1])+F(rows[b][1])))
    return candidates


def exchange_control():
    # Spin/weak exchange is an actual single-particle permutation. Its CAR lift
    # includes the canonical occupation-word sign, not a relabelled spectrum.
    permutation=[native.INDEX[x,f,s,c]for x,s,f,c in native.MODES]
    assert sorted(permutation)==list(range(24))
    for x in range(2):
        for axis in (1,2):
            diagonal=np.array([(1 if v[axis]==0 else -1)if v[0]==x else 0 for v in native.MODES])
            expected=np.array([(1 if v[3-axis]==0 else -1)if v[0]==x else 0 for v in native.MODES])
            assert np.array_equal(diagonal[permutation],expected)
    original=native.basis();target=native.basis(4,3);lookup={x:i for i,x in enumerate(target)};mapped=[];signs=[]
    for mask in original:
        images=[permutation[i]for i in range(24)if(mask>>i)&1]
        out=sum(1<<i for i in images);mapped.append(lookup[out]);signs.append((-1)**sum(a>b for i,a in enumerate(images)for b in images[i+1:]))
    assert len(set(mapped))==len(target)==11700
    return {'native_one_body_permutation':permutation,'full_sector_permutation':mapped,'canonical_CAR_signs':signs,
            'exact_source_identity':'On S1,I0,Nu3,Nspinup4, spin/weak exchange gives S0,I1,Nu4,Nspinup3. Equal local Js=Jw leaves the core unchanged; eta I_total² adds2eta and the number term adds-delta. I-lowering gives the other two charge components with the declared delta shifts.'}


@lru_cache(None)
def certificate():
    d=reduction();p=parameters();unit=F(INPUTS['electron_unit_GeV']);spectrum=spectral_enclosure(d['H_exact']);v=trial(d,spectrum)
    E=F(v['exact_energy']);thresholds=baryon_threshold(1);low=min(x[2]for x in thresholds);gap=low-E;assert gap>1
    other=[]
    for parts in partitions(6):
        if len(parts)<2 or parts==(3,3):continue
        colour=sum(n%3!=0 for n in parts)*F(4,3)*p['U']
        lower=6*p['mbar']+p['A']*sum((n*n+1)//2 for n in parts)/9+colour
        assert lower>E;other.append({'partition':list(parts),'threshold_lower':str(lower)})
    # N3 nontrivial colour sectors have a stronger floor than either singlet threshold.
    nonsinglet=6*p['mbar']+F(16,3)*p['U']+F(10,9)*p['A'];assert nonsinglet>E
    partners=[]
    for Iz in (-1,0,1):
        shift=2*p['eta']-Iz*p['delta'];intervals=[[str(F(a)+shift),str(F(b)+shift)]for a,b in spectrum['all_eigenvalue_intervals']]
        candidates=baryon_threshold(Iz+1);threshold_upper=min(x[3]for x in candidates);above=F(intervals[0][0])-threshold_upper
        assert above>F(2,5)
        partners.append({'charge':Iz+1,'internal_spin':0,'isospin':1,'complete_spectrum_intervals':intervals,
                         'ground_above_same_charge_two_baryon_threshold_lower':str(above),'bound_below_two_baryons':False})
    comparison=F(INPUTS['held_out_comparisons_GeV']['deuteron_binding']);binding=unit*gap
    return {'schema':'oph.native_nuclear_exchange_calibrated_source.v1','verdict':'PASS_NATIVE_BOUND_TRIPLET_AND_UNBOUND_EXCHANGED_PARTNER_WITH_HELD_OUT_BINDING_TEST',
      'same_source_parameters':{k:str(x)for k,x in p.items()},'complete_source_sector_dimensions':[11700,1800,27],
      'fixed_triplet_quantum_numbers':{'quark_number':6,'up_flavour_number':3,'spin':1,'spin_projection':1,'isospin':0,'isospin_projection':0},
      'all_spin_projections_triplet_dimension':81,'partner_spectrum_scope':'27 levels per fixed isospin component, spin zero; all three components retained.',
      'source':'The exact same positive baryon-density/native-edge-current action used for calibrated nucleons and charged pions; no nuclear value enters its parameters.',
      'complete_triplet_spectrum':spectrum,'complete_original_records':v,
      'same_charge_baryon_thresholds':[{'charges':[a,b],'lower':str(lo),'upper':str(hi)}for a,b,lo,hi in thresholds],
      'binding_lower':str(gap),'binding_lower_GeV':str(binding),'stationary_energy_GeV_interval':[str(unit*F(x))for x in spectrum['all_eigenvalue_intervals'][0]],
      'all_other_fragment_thresholds':other,'two_coloured_three_quark_fragment_floor':str(nonsinglet),
      'exact_spin_weak_exchange':exchange_control(),'complete_spin_zero_isotriplet_partners':partners,
      'held_out_deuteron_binding_GeV':str(comparison),'binding_fraction_of_observed_diagnostic':float(binding/comparison),
      'physical_deuteron_binding_reproduced':False,
      'zero_edge_current_control':'Removing eta I_total² leaves the S1,I0 spectrum unchanged and lowers each baryon threshold by3eta/4, while removing the partner splitting2eta. This is a different declared action, not a silent change of thresholds.',
      'stability':'Sector-preserving generator perturbations bounded by epsilon retain binding when3epsilon is below the certified margin; partner exclusion persists under the analogous3epsilon bound.',
      'scope':'Complete finite S1,I0 and exchanged S0,I1 spectra, same-source fragment comparisons and original records. Other spin/isospin channels are retained; the graph-refinement package tests additional low bands. No fm length, continuum scattering or complete observed nuclear spectrum follows.'}
