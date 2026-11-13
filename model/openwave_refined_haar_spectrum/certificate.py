"""Complete-complement spectrum and local common-clock spectral moments."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import sympy as sy
from model.openwave_refined_haar_source import gauss
from model.openwave_cubic_haar_source import cube
from model.openwave_coupled_haar_spectrum.certificate import algebraic_bounds,source_digest,interval_product
from model.openwave_colour_string.fixed_witness import spectral_enclosure
from model.openwave_su2_spectral.spectrum import sqrt_upper
from model.openwave_euclidean_refinement.response import exponential_decay_interval,dyadic_outer

HERE=Path(__file__).resolve().parent
D=2**48

def strings(interval):return [str(v)for v in dyadic_outer(interval)]
def display(interval):return [float(v)for v in interval]


@lru_cache(None)
def data(which):
    if which=='refined':
        states=gauss.basis();actions=gauss.total_action();energies=tuple(gauss.electric(s)for s in states)
        alpha,beta,faces,floor=F(16),F(1,4),36,F(80)
        operators=(('interior_fine_plaquette',gauss.graph()[2][16]),('embedded_coarse_boundary',gauss.coarse_boundary_loop()))
        apply=lambda state,loop:gauss.successors(state,loop)
    elif which=='coarse':
        states=tuple(s for _,s,_ in cube.basis());energies=tuple(F(sum(k*(k+2)for k in s),4)for s in states)
        all_columns=[]
        for state in states:
            vector={}
            for face in range(6):
                for target,value in cube.successors(state,face):vector[target]=vector.get(target,0)+value
            all_columns.append(tuple((s,sy.simplify(v))for s,v in vector.items()if v))
        actions=tuple(all_columns);alpha,beta,faces,floor=F(32),F(1,32),6,F(208)
        operators=(('coarse_boundary_plaquette',0),)
        apply=lambda state,face:tuple(cube.successors(state,face))
    else:raise ValueError('unknown declared spatial source')
    index={s:i for i,s in enumerate(states)};entries={}
    for j,column in enumerate(actions):
        entries[j,j]=sy.Rational(alpha*energies[j]+faces*beta)
        for target,value in column:
            if target in index:entries[index[target],j]=entries.get((index[target],j),0)-sy.Rational(beta)*value
    H=sy.SparseMatrix(len(states),len(states),entries);assert H==H.T
    return {'states':states,'actions':actions,'energies':energies,'alpha':alpha,'beta':beta,'faces':faces,'floor':floor,'operators':operators,'apply':apply,'H':H}


def rounded_action(states,coordinates,columns):
    vector={};row_error={};col_error=[]
    for amplitude,column in zip(coordinates,columns):
        error_sum=F(0)
        for target,value in column:
            lo,hi=algebraic_bounds(value);integer=((lo+hi)*D/2+F(1,2)).__floor__()
            error=max(abs(F(integer,D)-lo),abs(F(integer,D)-hi))
            row_error[target]=row_error.get(target,F(0))+error;error_sum+=error
            vector[target]=vector.get(target,0)+integer*amplitude
        col_error.append(error_sum)
    norm_error=sqrt_upper(max(row_error.values(),default=F(0))*max(col_error,default=F(0)))
    return vector,norm_error


def offblock_bound(states,columns):
    retained=set(states);rows={};cols=[]
    for column in columns:
        total=F(0)
        for target,value in column:
            if target in retained:continue
            lo,hi=algebraic_bounds(value);bound=max(abs(lo),abs(hi))
            total+=bound;rows[target]=rows.get(target,F(0))+bound
        cols.append(total)
    return max(rows.values(),default=F(0))*max(cols,default=F(0)),len(rows)


@lru_cache(None)
def calculate(which):
    source=data(which);states=source['states'];n=len(states)
    packet=json.loads((HERE/'numerical_witnesses.json').read_text())[which]
    H=source['H'];assert packet['exact_matrix_sha256']==source_digest(H)
    assert packet['dimension']==n and packet['ground_denominator']==D
    frame=spectral_enclosure(H,packet['complete_finite_frame'])
    intervals=[tuple(map(F,p))for p in frame['all_eigenvalue_intervals']]
    beta,alpha,floor=source['beta'],source['alpha'],source['floor']
    bound,external=offblock_bound(states,source['actions']);bound*=beta**2
    assert floor>intervals[1][1]
    shift=bound/(floor-intervals[1][1]);excited=(intervals[1][0]-shift,intervals[1][1])
    coordinates=tuple(int(v)for v in packet['ground_integer_coordinates'])
    assert len(coordinates)==n and all(str(v)==raw and abs(v)<=D for v,raw in zip(coordinates,packet['ground_integer_coordinates']))
    norm=sum(v*v for v in coordinates);assert norm>0
    vector,error=rounded_action(states,coordinates,source['actions']);u=dict(zip(states,coordinates))
    all_states=set(u)|set(vector);electric=dict(zip(states,source['energies']))
    hu={s:(alpha*electric.get(s,F(0))+source['faces']*beta)*u.get(s,0)-beta*F(vector.get(s,0),D)for s in all_states}
    operator_error=beta*error
    mu=sum(F(u.get(s,0))*value for s,value in hu.items())/norm
    residual=sqrt_upper(sum((value-mu*u.get(s,0))**2 for s,value in hu.items())/norm)+operator_error
    denominator=excited[0]-mu-operator_error;assert denominator>0
    ground=(max(F(0),mu-operator_error-residual**2/denominator),mu+operator_error)
    distance=residual/denominator;gap=(excited[0]-ground[1],excited[1]-ground[0]);assert gap[0]>0
    records=[]
    for name,loop in source['operators']:
        length=len(loop)if which=='refined' else 4
        columns=tuple(source['apply'](state,loop)for state in states)
        w,err=rounded_action(states,coordinates,columns)
        mean=sum(F(u.get(s,0)*v,D)for s,v in w.items())/norm
        mean_interval=(mean-err-2*distance,mean+err+2*distance)
        second=sum(F(v*v,D*D)for v in w.values())/norm
        second_error=2*err+err**2+2*distance
        second_interval=(max(F(0),second-second_error),min(F(1),second+second_error))
        sq=interval_product(mean_interval,mean_interval)
        variance=(max(F(0),second_interval[0]-max(F(0),sq[1])),min(F(1),second_interval[1]-max(F(0),sq[0])))
        assert variance[0]>0
        moment=(alpha*F(length,4)*(1-second_interval[1]),alpha*F(length,4)*(1-second_interval[0]))
        normalized=(moment[0]/variance[1],moment[1]/variance[0])
        correlations=[]
        for time in(F(0),F(1,128),F(1,64),F(1,32)):
            lower=variance[0]*exponential_decay_interval(time*moment[1]/variance[0])[0]
            upper=variance[1]*exponential_decay_interval(time*gap[0])[1]
            assert 0<lower<=upper
            correlations.append({'time_exact':str(time),'centered_Euclidean_correlation_interval_exact':strings((lower,upper)),'interval':display((lower,upper)),'proof':'positive spectral measure Jensen lower bound and full-source spectral-gap upper bound'})
        records.append({'observable':name,'loop_length':length,'full_Wilson_trial_output_states':len(w),
            'mean_interval_exact':strings(mean_interval),'second_moment_interval_exact':strings(second_interval),'variance_interval_exact':strings(variance),
            'source_clock_first_spectral_moment_interval_exact':strings(moment),'normalized_first_spectral_moment_interval_exact':strings(normalized),'normalized_first_spectral_moment_interval':display(normalized),
            'exact_double_commutator_identity':'m1=alpha*L/4*(1-<W^2>)','correlations':correlations})
    report={'source':which,'spacing_exact':'1/2'if which=='refined'else'1','g_squared_exact':'16'if which=='refined'else'64',
        'alpha_exact':str(alpha),'beta_exact':str(beta),'complete_frame_dimension':n,'complete_spectral_certificate':frame,
        'full_complement_floor_exact':str(floor),'all_omitted_action_states':external,'full_offblock_norm_squared_upper_exact':str(bound),'Schur_shift_upper_exact':str(shift),
        'full_ground_energy_interval_exact':strings(ground),'full_gap_interval_exact':strings(gap),'full_gap_interval':display(gap),
        'full_trial_residual_upper_exact':str(residual),'ground_projector_distance_upper_exact':str(distance),'observables':records,
        'all_Haar_representations_and_intertwiner_channels_retained':True,'physical_GeV_mass':None}
    return report,{'source':source,'bound':bound,'ground':ground,'gap':gap,'distance':distance,'coordinates':coordinates,'norm':norm,'residual':residual}


def run():
    coarse,_=calculate('coarse');fine,_=calculate('refined')
    left=tuple(map(F,coarse['observables'][0]['normalized_first_spectral_moment_interval_exact']))
    right=tuple(map(F,fine['observables'][1]['normalized_first_spectral_moment_interval_exact']))
    return {'schema':'oph.refined_full_Haar_common_clock.v1','verdict':'PASS_INTERIOR_RESOLVED_FULL_SOURCE_SPECTRUM_AND_COMMON_CLOCK_MOMENTS',
        'records':[coarse,fine],'common_clock_electric_matching':'coarse alpha=32 equals twice refined alpha=16 under the Haar edge subdivision',
        'embedded_boundary_minus_coarse_normalized_moment_interval_exact':strings((right[0]-left[1],right[1]-left[0])),
        'boundary_moment_intervals_overlap':max(left[0],right[0])<=min(left[1],right[1]),
        'fine_interior_modes_are_retained':True,'derived_running_coupling_or_continuum_limit':False,
        'physical_GeV_calibration_or_new_Nature_validation':False,'new_empirical_fit_inputs':[]}
