from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from itertools import combinations
import json
import sympy as sy
from .cube import recoupled_matrices,spin_census,successors,trial_state_set
from model.openwave_coupled_haar_spectrum.certificate import source_digest,algebraic_bounds,strings,display,interval_product
from model.openwave_colour_string.fixed_witness import spectral_enclosure
from model.openwave_su2_spectral.spectrum import sqrt_upper

HERE=Path(__file__).resolve().parent
D=2**48


@lru_cache(None)
def action(states,face):
    """Every nonzero output is retained, including all outside spin labels."""
    columns=[]
    for spins in states:columns.append(tuple(successors(spins,face)))
    return tuple(columns)


def rounded_action(states,coordinates,face):
    vector={};row_error={};column_errors=[]
    for amplitude,column in zip(coordinates,action(states,face)):
        column_error=F(0)
        for target,value in column:
            lo,hi=algebraic_bounds(value);integer=((lo+hi)*D/2+F(1,2)).__floor__()
            error=max(abs(F(integer,D)-lo),abs(F(integer,D)-hi))
            row_error[target]=row_error.get(target,F(0))+error;column_error+=error
            vector[target]=vector.get(target,0)+integer*amplitude
        column_errors.append(column_error)
    error=sqrt_upper(max(row_error.values(),default=F(0))*max(column_errors,default=F(0)))
    return vector,error


@lru_cache(None)
def calculate(g2):
    packet=json.loads((HERE/'numerical_witnesses.json').read_text())[str(g2)]
    g2=F(g2);alpha,beta=g2/2,2/g2
    states,E,W=recoupled_matrices(packet['complete_base_quarter_energy'])
    assert packet['ordered_complete_base_spins']==[list(s)for s in states]
    H=sy.Rational(alpha)*E+sy.Rational(beta)*(6*sy.eye(E.rows)-sum(W,sy.zeros(E.rows)))
    assert source_digest(H)==packet['exact_complete_base_matrix_sha256']
    frame=spectral_enclosure(H,packet['complete_finite_frame'])
    finite=[tuple(map(F,row))for row in frame['all_eigenvalue_intervals']]
    floor=alpha*F(packet['complete_base_quarter_energy']+1,4)
    assert floor>finite[1][1]
    shift=36*beta**2/(floor-finite[1][1])
    excited=(finite[1][0]-shift,finite[1][1])
    assert packet['trial_quarter_energy']==56
    trial_states=trial_state_set()
    assert packet['ordered_trial_spins']==[list(s)for s in trial_states]
    assert packet['ground_denominator']==D
    raw=packet['ground_integer_coordinates']
    assert len(raw)==len(trial_states)and all(type(x)is str and str(int(x))==x and abs(int(x))<=D for x in raw)
    coordinates=tuple(int(x)for x in raw);norm=sum(x*x for x in coordinates);assert norm>0
    vectors=[];errors=[]
    for face in range(6):
        vector,error=rounded_action(trial_states,coordinates,face);vectors.append(vector);errors.append(error)
    u=dict(zip(trial_states,coordinates));all_labels=set(u).union(*(set(v)for v in vectors))
    hu={s:(alpha*F(sum(k*(k+2)for k in s),4)+6*beta)*u.get(s,0)
             -beta*sum(F(v.get(s,0),D)for v in vectors)for s in all_labels}
    operator_error=beta*sum(errors)
    mu=sum(F(u.get(s,0))*value for s,value in hu.items())/norm
    residual=sqrt_upper(sum((hu[s]-mu*u.get(s,0))**2 for s in all_labels)/norm)+operator_error
    assert mu+operator_error<excited[0]
    distance=residual/(excited[0]-mu)
    ground=(max(F(0),mu-operator_error-residual**2/(excited[0]-mu-operator_error)),mu+operator_error)
    gap=(excited[0]-ground[1],excited[1]-ground[0]);assert gap[0]>0
    means=[];variances=[];nominal=[]
    for vector,error in zip(vectors,errors):
        mean=sum(F(u.get(s,0)*v,D)for s,v in vector.items())/norm
        error_mean=error+2*distance
        means.append((mean-error_mean,mean+error_mean));nominal.append(mean)
        square=sum(F(v*v,D*D)for v in vector.values())/norm
        square_error=2*error+error**2+2*distance
        mean_square=interval_product(means[-1],means[-1])
        variances.append((max(F(0),square-square_error-mean_square[1]),min(F(1),square+square_error-mean_square[0])))
    covariances=[]
    for a,b in combinations(range(6),2):
        mean=sum(F(value*vectors[b].get(s,0),D*D)for s,value in vectors[a].items())/norm
        error=errors[a]+errors[b]+errors[a]*errors[b]+2*distance
        product=interval_product(means[a],means[b]);covariance=(mean-error-product[1],mean+error-product[0])
        covariances.append({'faces':[a,b],'opposite':a//2==b//2,
            'connected_covariance_interval_exact':strings(covariance),'connected_covariance_interval':display(covariance)})
    return {'g_squared_exact':str(g2),'spacing_exact':'1','complete_base_dimension':len(states),
        'trial_dimension':len(trial_states),'complete_trial_output_dimension':len(all_labels),
        'all_trial_outputs_retained':True,'complete_finite_spectrum':frame,
        'full_Q_floor_exact':str(floor),'full_offblock_norm_squared_upper_exact':str(36*beta**2),
        'Schur_error_upper_exact':str(shift),'full_source_ground_energy_interval_exact':strings(ground),
        'full_source_gap_interval_exact':strings(gap),'full_source_gap_interval':display(gap),
        'full_native_trial_residual_upper_exact':str(residual),
        'ground_projector_distance_upper_exact':str(distance),
        'all_six_Wilson_means_intervals_exact':[strings(x)for x in means],
        'all_six_Wilson_variances_intervals_exact':[strings(x)for x in variances],
        'all_fifteen_connected_covariances':covariances,
        'all_source_representations_retained':True,
        'physical_4D_Yang_Mills_or_GeV_mass_identified':False}
