"""Exact full-H moments and variational spectral-window certificates."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import sympy as sy
from model.openwave_full_strong_current.native import action_vector,full_action,native_vector,dot,translate,clean,add,coefficients,colour_casimir
from model.openwave_strange_shell_current.native import contraction,frame
from model.openwave_strange_shell_current.source import trial,lambda_pole
from model.openwave_scattering_current_support.native import channel
from model.openwave_normalizable_particle_records.source import cosine_mean
from model.openwave_lambda_wavepacket.source import KAPPA,packet
from model.openwave_capacity_stress.source import enclosure
WITNESSES=Path(__file__).with_name('witnesses.json')


@lru_cache(None)
def full_column(state):
    """Identical occupation columns are shared by the four exact witnesses."""
    return full_action(state)


def interval_sum(coefficients):
    lo=hi=F(0)
    for shift,value in coefficients.items():
        a,b=cosine_mean(abs(shift)*KAPPA)if shift else(F(1),F(1))
        lo+=value*(a if value>=0 else b);hi+=value*(b if value>=0 else a)
    return lo,hi


def fibres(vector):
    rows=defaultdict(lambda:defaultdict(F))
    for state,value in vector.items():
        x=min(m[1]for m in state);rows[translate(state,-x)][x]+=value
    return {state:dict(row)for state,row in rows.items()}


def packet_inner(left,right):
    """Even packet integrals of every translation return path, before rounding."""
    result=defaultdict(F);left_rows=fibres(left);right_rows=fibres(right)
    for state,a in left_rows.items():
        for x,u in a.items():
            for y,v in right_rows.get(state,{}).items():result[x-y]+=u*v
    assert all(abs(k)*KAPPA<1 for k in result)
    return interval_sum(result),dict(result)


@lru_cache(None)
def native_current():
    d=contraction();_,P=frame();rest=lambda_pole();v=trial(rest['actual_trial_coordinates'])
    U=d['matrix']*P*v
    u=native_vector({mask:U[i,0]for i,mask in enumerate(d['raw_output_masks'])if U[i,0]})
    assert all(min(m[1]for m in state)==0 for state in u)
    assert not action_vector(colour_casimir,u,(0,1))
    # Actual integer native charge, including the antiquark sign.
    assert all(sum((1 if sp==0 else-1)*(2 if w==0 else-1)for sp,x,spin,w,c in state)==0 for state in u)
    n=F((v.T*(P.T*P)*v)[0]);h=action_vector(full_column,u)
    mean,mean_words=packet_inner(u,h);second,second_words=packet_inner(h,h)
    norm=dot(u,u);mean=tuple(x/norm for x in mean);second=tuple(x/norm for x in second)
    variance=(max(F(0),second[0]-mean[1]**2),second[1]-mean[0]**2)
    assert 0<variance[0]<=variance[1]
    scale=d['Cabibbo_squared']/(1296*6*n)
    success=scale*norm
    assert success==packet()['trial_total_success']
    center=sum(map(F,rest['full_ground_GeV_interval']))/2/coefficients()['unit']
    return {'u':u,'h':h,'norm':norm,'mean':mean,'second':second,'variance':variance,
        'mean_words':mean_words,'second_words':second_words,'normalization':n,'success':success,
        'current_scale_squared':scale,'center':center,'unit':coefficients()['unit'],'rest':rest}


@lru_cache(None)
def variational_columns():
    ch=channel();Y=ch['basis']
    return [native_vector({mask:Y[i,j]for i,mask in enumerate(ch['raw_masks'])if Y[i,j]})for j in range(Y.cols)]+[native_current()['u']]


def witness(row):
    d=native_current();cols=variational_columns();assert len(row['coordinates_integer'])==len(cols)==24
    x=defaultdict(F)
    for col,numerator in zip(cols,row['coordinates_integer'],strict=True):add(x,col,F(numerator,row['denominator']))
    x=clean(x);assert all(min(m[1]for m in state)==0 for state in x)
    hx=action_vector(full_column,x);r=defaultdict(F);add(r,hx);add(r,x,-d['center']);r=clean(r)
    overlap=dot(d['u'],x);nx=dot(x,x);residual,_=packet_inner(r,r)
    width=F(row['half_width_GeV'])/d['unit'];eta=width*F(row['eta_over_half_width'])
    denominator=residual[1]+eta*eta*nx
    effect=eta*eta*overlap*overlap/denominator
    outside=eta*eta/(width*width+eta*eta)
    window=max(F(0),(effect-outside*d['norm'])/(1-outside))/d['norm']
    centered_second_upper=d['second'][1]-2*d['center']*d['mean'][0]+d['center']**2
    moment_window=max(F(0),1-centered_second_upper/(width*width))
    combined=max(window,moment_window)
    assert 0<effect<d['norm'] and 0<=window<=1
    return {'half_width_GeV_exact':row['half_width_GeV'],'regularization_source_exact':str(eta),
       'trial_Gram_norm_exact':str(nx),'native_current_trial_overlap_exact':str(overlap),
       'full_source_residual_squared_interval':list(map(str,residual)),
       'full_H_trial_successors':len(hx),'all_translation_return_paths_retained':True,
       'normalized_Lorentzian_effect_lower_exact':str(effect/d['norm']),
       'normalized_Lorentzian_effect_lower_display':float(effect/d['norm']),
       'normalized_energy_window_probability_lower_exact':str(window),
       'normalized_energy_window_probability_lower_display':float(window),
       'normalized_moment_window_probability_lower_exact':str(moment_window),
       'normalized_best_window_probability_lower_exact':str(combined),
       'normalized_best_window_probability_lower_display':float(combined),
       'original_input_window_probability_lower_exact':str(d['success']*combined/7),
       'variational_resolvent_method_is_on_shell_density_proof':False}


def run():
    d=native_current();unit=d['unit'];mean=tuple(x*unit for x in d['mean']);variance=tuple(x*unit*unit for x in d['variance'])
    rows=[witness(row)for row in json.loads(WITNESSES.read_text())['witnesses']]
    assert rows[-1]['normalized_energy_window_probability_lower_display']>0
    assert F(rows[-1]['normalized_best_window_probability_lower_exact'])>F(1,2)
    upper=[]
    for width in(F('0.001'),F('0.01'),F('0.1')):
        endpoint=d['center']*unit+width;distance=mean[0]-endpoint;assert distance>0
        bound=variance[1]/(variance[1]+distance*distance)
        upper.append({'half_width_GeV_exact':str(width),'conditional_current_window_probability_interval':['0',str(bound)],
                      'conditional_probability_upper_display':float(bound),
                      'original_input_probability_upper_exact':str(d['success']*bound/7),
                      'positive_lower_bound_proved':False})
    initial=packet()['energy'];work=(mean[0]-initial[1],mean[1]-initial[0])
    assert work[0]>1
    assert (mean[0]-d['center']*unit-F('0.1'))**2>variance[1]
    return {'schema':'oph.native_lambda_current_spectral_measure.v1',
      'verdict':'PASS_COMPLETE_NATIVE_CURRENT_PACKET_MOMENTS_AND_VARIATIONAL_WINDOW_CERTIFICATES',
      'source_and_preparation':{'Hamiltonian':'Complete autonomous infinite-chain Nq4/Na1 H5, including both common inverse-rest correlated terms, colour strings and every CAR successor.',
        'input':'Exact inherited rational native Lambda trial times the normalized kappa=1/512 polynomial packet.',
        'native_vertex_words':72,'raw_current_output_count':len(d['u']),'full_H_current_successor_count':len(d['h']),
        'native_current_Gram_norm_exact':str(d['norm']),
        'actual_current_source_charge_exact':'0','full_colour_Gauss_annihilation_checked':True,
        'trial_preparation_success_exact':'1/7','trial_preparation_failure_exact':'6/7',
        'native_W_over_36_success_exact':str(d['success']),'probe_failure_exact':str(1-d['success']),
        'normalizable_packet_translation_return_paths_retained':True,'new_source_cutoff':False},
      'actual_successful_current_energy':{'mean_GeV_interval':list(map(str,mean)),'mean_GeV_display':enclosure(mean),
        'variance_GeV_squared_interval':list(map(str,variance)),'variance_GeV_squared_display':enclosure(variance),
        'source_mean_Fourier_coefficients':{str(k):str(v)for k,v in d['mean_words'].items()},
        'source_second_moment_Fourier_coefficients':{str(k):str(v)for k,v in d['second_words'].items()},
        'successful_branch_mean_energy_change_GeV_interval':list(map(str,work)),
        'successful_branch_mean_energy_change_GeV_display':enclosure(work),
        'conditioned_energy_change_is_nonselective_probe_work':False,
        'nonselective_probe_and_failure_work':'Tr[H(Phi(rho)-rho)] with both sector-changing success and full failure; no zero-work assertion.'},
      'energy_windows':{'center_GeV_exact':str(d['center']*unit),
        'center_role':'Midpoint of the inherited full Lambda rest enclosure; a finite-width source energy service, not the exact unknown pole.',
        'narrow_window_Cantelli_upper_bounds':upper,'exact_variational_resolvent_witnesses':rows,
        'full_spectral_instrument':'1_I(H5) and I-1_I(H5), both retained; zero nonselective energy work, arbitrary idle quantum reference.',
        'Lorentzian_effect':'eta^2[(H5-E)^2+eta^2]^-1 with its complete effect complement.',
        'finite_local_compilation_of_spectral_service_proved':False},
      'controls':{'zero_current_spectral_measure_exact':'0',
        'wrong_total_charge_projected_current_exact':'0, by the checked literal CAR charge of every raw output.',
        'finite_moment_nonimplication':'A two-atom measure at mean +/- sqrt(variance) has the same first two normalized moments and no weight near the Lambda energy; exact moments alone cannot imply an on-shell overlap.',
        'disconnected_channel_countercontrol':'On C direct-sum L2([0,3]), let H=4 direct-sum multiplication by E and u=(1,0). Every E in (0,3) belongs to essential spectrum, but the current spectral measure there is exactly zero. This is a countercontrol to the implication, not a replacement for the native H5.',
        'compact_resolvent_method':'All 23 actual native proton/pion channel columns and the native current are retained. Their finite variational witnesses use full H5 residuals; zero certified lower bounds do not prove a zero spectral measure.'},
      'essential_channel_relation':'The separate complete-source Weyl sequence places the Lambda rest energy in essential spectrum. It supplies no lower overlap for this particular native current.',
      'shrinking_window_density_obligation':'A density claim requires a controlled limit of measure([E-w,E+w])/(2w), or an outgoing wave-operator overlap. The positive broad-window and finite-eta effects do not establish that limit.',
      'on_shell_native_weak_density_proved':False,'physical_Lambda_decay_rate_or_lifetime_proved':False,
      'new_calibration':False,'new_match_with_nature':False,
      'scope':'Exact energy moments and complete finite-resolution spectral bounds for the actual normalizable native-current experiment. The current success branch is energetically far from an energy-conserving weak decay; the full on-shell scattering overlap is a separate calculation.'}
