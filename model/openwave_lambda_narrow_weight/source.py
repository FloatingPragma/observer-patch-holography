"""Rational approximate inverses with complete source residuals."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement
from pathlib import Path
import json
import sympy as sy
from model.openwave_lambda_spectral_measure.source import native_current,variational_columns,full_column,packet_inner
from model.openwave_full_strong_current.native import action_vector,native_vector,add,clean,dot
from model.openwave_scattering_current_support.native import epsilon,normalized_rational_column,pion_field
from model.openwave_pole_weak_current.native import contraction
from model.openwave_baryon_poles.source import phase,vector,trial_moments,residue_bounds
from model.openwave_baryon_dispersion import exact_native as base
from model.openwave_cluster_weyl_sequence.source import wedge,shifted
from model.openwave_cluster_bandwidth.meson import matrices
from model.openwave_hadron_refinement_flow.certified_rest import eigen_intervals,sqrt_upper
from model.openwave_capacity_stress.source import enclosure

WITNESSES=Path(__file__).with_name('witnesses.json')

@lru_cache(None)
def clusters():
    recorded=json.loads(WITNESSES.read_text())
    labels=tuple(combinations_with_replacement(range(16),3))
    shapes=sorted({base.canonical(lab)[1]for lab in labels if max(base.ORBITALS[i][0]for i in lab)-min(base.ORBITALS[i][0]for i in lab)<=1})
    selected=[shape for shape in shapes if sum(w==0 for x,spin,w in shape)==2 and sum(spin==0 for x,spin,w in shape)==2]
    proton=phase(True,F(1),F(0));v=vector(proton['trial_source_coordinates'])
    weights=contraction()['native_proton_embedding']*sy.diag(sy.sqrt(6),1,1,1,1)*v
    raw=defaultdict(lambda:sy.S.Zero)
    for shape,c in zip(selected,weights,strict=True):
        for mask,a in epsilon(shape).items():raw[mask]+=c*a
    p=native_vector(normalized_rational_column(raw))
    pion=defaultdict(F);den=recorded['pion_denominator'];nums=recorded['pion_coordinates_integer']
    assert den==2**40 and nums[0]==den and nums[1]==nums[2]>0
    for (x,y),n in zip(((0,0),(0,1),(1,0)),nums,strict=True):
        add(pion,native_vector(pion_field(x,y)),F(n,den))
    A,C,params,unit=matrices(F(1),F(0),repair=False)
    v=sy.Matrix([sy.Rational(n,den)for n in nums]);n,mean,variance=trial_moments(A,sy.eye(3),v)
    bands=eigen_intervals(A,sy.eye(3));error=4*params['t']**2/(C-bands[0][1])
    certificate=residue_bounds(bands,C,4*params['t']**2,error,mean,variance)
    assert dot(p,p)>0 and dot(pion,pion)==6*n
    return p,clean(pion),{'proton_prepared_ground_weight_lower_exact':proton['pole']['prepared_trial_pole_weight_lower'],
        'pion_prepared_ground_weight_lower_exact':certificate['prepared_trial_pole_weight_lower'],
        'proton_native_occupations':len(p),'pion_native_occupations':len(pion),
        'pion_rational_coordinates':list(map(str,v)),'pion_compact_residual_variance_exact':str(variance),
        'role':'Compact rational cluster-ground trials; complete isolated-cluster complement errors are retained. Their products are variational vectors, not exact full-H5 scattering states.'}

@lru_cache(None)
def frame():
    p,m,certificate=clusters();recorded=json.loads(WITNESSES.read_text())
    offsets=recorded['relative_offsets'];assert offsets==[-3,-2,-1,0,1,2,3,4]
    columns=variational_columns()+[wedge(p,shifted(m,x))for x in offsets]
    assert len(columns)==32 and all(dot(x,x)>0 for x in columns)
    images=[action_vector(full_column,x)for x in columns]
    return columns,images,certificate

def window_upper(residual_square_upper,trial_square_upper,width,norm):
    if min(residual_square_upper,trial_square_upper,width)<0 or norm<=0:raise ValueError('invalid norm or width')
    return min(F(1),(sqrt_upper(residual_square_upper)+width*sqrt_upper(trial_square_upper))**2/norm)

def witness(row):
    d=native_current();columns,images,_=frame();x=defaultdict(F);hx=defaultdict(F)
    assert len(row['coordinates_integer'])==len(columns)
    for col,image,num in zip(columns,images,row['coordinates_integer'],strict=True):
        assert type(num)is int
        add(x,col,F(num,row['denominator']));add(hx,image,F(num,row['denominator']))
    x=clean(x);hx=clean(hx);r=defaultdict(F);add(r,d['u']);add(r,hx,-1);add(r,x,d['center']);r=clean(r)
    rr,returns=packet_inner(r,r);xx,_=packet_inner(x,x)
    width=F(row['half_width_GeV'])/d['unit']
    inverse=window_upper(rr[1],xx[1],width,d['norm'])
    distance=d['mean'][0]-d['center']-width;assert distance>0
    cantelli=d['variance'][1]/(d['variance'][1]+distance*distance)
    selected=min(inverse,cantelli)
    return {'half_width_GeV_exact':row['half_width_GeV'],
        'full_source_residual_squared_interval':list(map(str,rr)),
        'normalizable_trial_norm_squared_interval':list(map(str,xx)),
        'residual_occupation_count':len(r),'trial_occupation_count':len(x),'full_H_trial_successor_count':len(hx),
        'translation_return_path_orders':sorted(returns),'inverse_witness_upper_exact':str(inverse),
        'inherited_Cantelli_upper_exact':str(cantelli),'selected_conditional_window_interval':['0',str(selected)],
        'selected_conditional_window_display':enclosure((F(0),selected)),
        'original_input_window_upper_exact':str(d['success']*selected/7),
        'inverse_to_Cantelli_upper_ratio_display':float(inverse/cantelli),
        'positive_lower_bound_supplied':False,'full_complement_retained':True}

def run():
    d=native_current();_,_,certificate=frame();rows=[witness(row)for row in json.loads(WITNESSES.read_text())['witnesses']]
    assert F(rows[0]['selected_conditional_window_interval'][1])<F(27,1000)
    assert F(rows[0]['inverse_witness_upper_exact'])<F(3,5)*F(rows[0]['inherited_Cantelli_upper_exact'])
    return {'schema':'oph.native_lambda_narrow_current_weight.v1','verdict':'PASS_COMPLETE_H5_NARROW_WINDOW_UPPER_CERTIFICATES',
        'native_source':'Unchanged full infinite-chain H5, all72 native current words, normalized kappa=1/512 packet and exact inherited Lambda trial.',
        'source_carrier_or_coefficients_changed':False,'variational_frame':{'columns':32,'original_channel_and_current_columns':24,
            'cluster_relative_centres':[-3,-2,-1,0,1,2,3,4],'cluster_trial_certificate':certificate,
            'floating_search_defines_scientific_values':False,'fixed_rational_witnesses_checked_against_full_source':True},
        'center_GeV_exact':str(d['center']*d['unit']),'windows':rows,
        'projection_inequality':'For r=u-(H-E)x, ||1_I(H)u|| <= ||r||+w||x|| when I=[E-w,E+w]. All norms use the complete normalizable packet.',
        'complete_records':{'preparation_success_exact':'1/7','preparation_failure_exact':'6/7',
            'native_current_success_exact':str(d['success']),'current_failure_exact':str(1-d['success']),
            'spectral_outcomes':'1_I(H5) and I-1_I(H5) on every retained preparation/current branch; arbitrary idle reference.',
            'nonselective_energy_readout_work_exact':'0','finite_local_compilation_of_spectral_projector':False},
        'lower_witness_scope':'These upper certificates supply no positive narrow-window lower bound and prove no vanishing spectral measure.',
        'on_shell_density_or_outgoing_wave_operator_proved':False,'physical_decay_rate_or_lifetime_proved':False,
        'source_current_total_mean_is_rest_mass':False,'new_empirical_calibration':False,'new_match_with_nature':False}
