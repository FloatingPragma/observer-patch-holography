"""Nested full-native charged windows from complete local spectral roots."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
import json
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.local import positive_response
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import field_data

def frac(x):
    x=s.sympify(x);assert x.is_Rational
    return F(int(x.p),int(x.q))

def root(x,upper=False):
    x=F(x);assert x>=0
    scale=2**128;k=isqrt(x.numerator*scale**2//x.denominator)
    if upper and k*k*x.denominator!=x.numerator*scale**2:k+=1
    return F(k,scale)

def round_up(x):
    x=F(x);scale=2**128
    return F(-((-x.numerator*scale)//x.denominator),scale)

def load(path):return json.loads((ROOT/path).read_text())

@lru_cache(None)
def data():
    charged=load('model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json')
    atomic=load('model/coherent_defect_wave/chiral_joint/atomic_reference_control/results.json')
    row=charged['same_H_interacting_channels']['ec']['full_interacting_energy_readout']
    exact=positive_response('ec');e,f,_=field_data();q=exact['q']
    assert str(q.as_expr())==charged['exact_local_channels']['ec']['complete_target_annihilator']
    eroots=f.intervals(eps=s.Rational(1,10**100));elo,ehi=map(frac,eroots[0][0])
    assert eroots[0][1]==1 and elo<ehi<0
    oldelo,oldehi=map(frac,exact['e_interval']);assert oldelo<=elo<ehi<=oldehi
    qroots=q.intervals(eps=s.Rational(1,10**100))
    assert len(qroots)==9 and all(mult==1 for _,mult in qroots)
    energies=[];weights=[]
    for (bounds,_),witness in zip(qroots,exact['rows'],strict=True):
        lo,hi=map(frac,bounds)
        assert frac(witness['energy_interval'][0])<=lo<=hi<=frac(witness['energy_interval'][1])
        energies.append((lo+7*elo+4,hi+7*ehi+4))
        weights.append(tuple(map(frac,witness['weight_interval'])))
    v=F(row['full_hopping_perturbation_norm']);eta=F(row['actual_neutral_ground_trial_vector_error'])
    gap=F(atomic['atomic_lattice_reference']['full_physical_atomic_gap_lower_bound'])
    assert atomic['full_PW_hopping_controls']['exact_atomic_P0_T_P0']=='0'
    assert v==F(360,2**34) and 0<v<gap
    neutral_shift=v*v/(gap-v)
    return {'f':str(f.as_expr()),'q':str(q.as_expr()),'e':(elo,ehi),'energies':energies,
      'weights':weights,'v':v,'eta':eta,'gap':gap,'neutral_shift':neutral_shift,
      'old_absolute_bin':tuple(map(F,row['selected_total_absolute_energy_bin'])),
      'old_source_bin':tuple(map(F,row['source_excitation_energy_bin_after_subtracting_declared_reference_cost'])),
      'old_weight_lower':F(row['conditional_probabilities']['line'][0])}

def window(multiplier):
    d=data();v=d['v'];radius=multiplier*v
    lo,hi=d['energies'][0];center=(lo+hi)/2
    I=(center-radius,center+radius)
    assert d['old_absolute_bin'][0]<I[0]<I[1]<d['old_absolute_bin'][1]
    distances=[radius-(hi-lo)/2]
    for lo,hi in d['energies'][1:]:
        distance=max(I[0]-hi,lo-I[1]);assert distance>0;distances.append(distance)
    assert distances[0]>v
    selected=root(d['weights'][0][0]*(1-(v/distances[0])**2))
    other=round_up(sum((v*root(w[1],True)/dist for w,dist in zip(d['weights'][1:],distances[1:],strict=True)),F(0)))
    lower_amplitude=selected-other-d['eta'];assert lower_amplitude>0
    lower=lower_amplitude**2
    upper=min(F(1,2),(root(d['weights'][0][1],True)+other+d['eta'])**2)
    elo,ehi=d['e']
    source=(I[0]-4-8*ehi,I[1]-4-8*elo+d['neutral_shift'])
    assert 0<lower<=upper<F(1,2)
    return {'radius_in_units_of_full_hopping_norm':multiplier,'absolute_source_plus_reference_bin':list(map(str,I)),
      'source_excitation_interval_relative_actual_native_ground':list(map(str,source)),
      'source_excitation_width':str(source[1]-source[0]),
      'old_to_refined_width_ratio_lower':str((d['old_source_bin'][1]-d['old_source_bin'][0])/(source[1]-source[0])),
      'selected_atom_distance_to_bin_complement_lower':str(distances[0]),
      'other_visible_atom_distances_to_bin':list(map(str,distances[1:])),
      'complete_other_atom_leakage_norm_upper':str(other),
      'unnormalized_bin_weight_interval':[str(lower),str(upper)],
      'complete_outside_weight_interval':[str(F(1,2)-upper),str(F(1,2)-lower)],
      'literal_native_CAR_failure_probability':'1/2',
      'retains_original_certified_weight_lower':lower>=d['old_weight_lower'],
      'diagnostic_lower_weight':float(lower),'diagnostic_source_interval':list(map(float,source)),
      'diagnostic_width_improvement':float((d['old_source_bin'][1]-d['old_source_bin'][0])/(source[1]-source[0]))}

def calculate():
    d=data();rows=[window(c) for c in (2,4,8,16,32,64,128,256,512)]
    assert rows[6]['retains_original_certified_weight_lower']
    assert F(rows[6]['unnormalized_bin_weight_interval'][0])>F('.419706')
    assert F(rows[0]['unnormalized_bin_weight_interval'][0])>F('.31477')
    return {'schema':'oph.refined_native_charge_resolution.v1',
      'source':'unchanged full native H0+2^-34 T with complete physical Gauss carrier',
      'local_neutral_minimal_polynomial':d['f'],'complete_ec_local_target_annihilator':d['q'],
      'refined_neutral_local_root_interval':list(map(str,d['e'])),
      'full_hopping_norm':str(d['v']),'full_atomic_gap_lower':str(d['gap']),
      'full_neutral_ground_shift_interval_relative_atomic':['-'+str(d['neutral_shift']),'0'],
      'actual_ground_trial_vector_error':str(d['eta']),
      'old_bin_lower_weight':str(d['old_weight_lower']),
      'reference_cost_retained':'4','complete_record_law':['1/2','w_sharp','1/2-w_sharp'],
      'nested_window_frontier':rows,'recommended_same_residue_multiplier':128,
      'narrower_selection_is_explicit':True,'original_maximum_residue_eigenvalue_relocalized':False,
      'refined_Estar_selection':'maximum per-spin residue within refined bin, lower-energy tie break',
      'refined_Estar_weight_lower':'refined bin lower weight / full native source dimension D',
      'finite_group_representation':'G1, paired by exact native antiunitary with equal full opposite-reference costs',
      'affine_or_Dirac_couplers_automatically_unchanged_by_refined_selection':False,
      'physical_spin_mass_or_SM_identity_claimed':False}
