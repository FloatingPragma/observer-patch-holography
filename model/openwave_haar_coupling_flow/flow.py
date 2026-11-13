from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from model.openwave_su2_spectral.spectrum import sqrt_upper
from .certificate import calculate, witnesses, strings, display, ROUND_DENOMINATOR


def subtract(a,b):
    return a[0]-b[1],a[1]-b[0]


def quotient(a,b):
    assert b[0]>0
    values=[x/y for x in a for y in b]
    return min(values),max(values)


def excited_mean_interval(record,data):
    raw=data['packet']['low_trial_integer_coordinates'][1]
    u=np.array([int(x)for x in raw],dtype=object);norm=sum(int(x)**2 for x in u)
    energies=[tuple(map(F,row))for row in record['full_Haar_first_four_energy_intervals_exact']]
    center=F(int(data['packet']['complete_finite_frame']['dyadic_eigenvalues'][1]),2**40)
    window=tuple(map(F,record['full_native_residual_spectral_windows_exact'][1]))
    residual=max(center-window[0],window[1]-center)
    separation=min(center-energies[0][1],energies[2][0]-center)
    assert separation>0
    distance=residual/separation;assert distance<F(1,100)
    out=[]
    for operator,error in ((data['w1'],data['error1']),(data['w2'],data['error2'])):
        mean=F(sum(int(x)*int(y)for x,y in zip(u,operator@u)),ROUND_DENOMINATOR*norm)
        delta=error+2*distance;out.append((mean-delta,mean+delta))
    return out,distance


def derivative_record(key):
    record,data=calculate(key);x=F(record['g_squared_exact'])
    energy0=tuple(map(F,record['ground_energy_interval_exact']))
    energy1=tuple(map(F,record['full_Haar_first_four_energy_intervals_exact'][1]))
    means0=[tuple(map(F,row))for row in record['Wilson_means_intervals_exact']]
    means1,distance=excited_mean_interval(record,data)
    derivatives=[]
    for energy,means in ((energy0,means0),(energy1,means1)):
        magnetic=(2-sum(row[1]for row in means),2-sum(row[0]for row in means))
        assert 0<=magnetic[0]<=magnetic[1]<=4
        # H'(x)=H(x)/x-4M/x². Only M=2-W1-W2 is needed as
        # a bounded observable; unbounded electric expectations are not
        # estimated from a trace-distance inequality.
        derivatives.append((energy[0]/x-4*magnetic[1]/x**2,
                            energy[1]/x-4*magnetic[0]/x**2))
    gap=tuple(map(F,record['full_Haar_gap_interval_exact']))
    gap_derivative=subtract(derivatives[1],derivatives[0])
    assert gap_derivative[0]>0
    # Holding Delta(x(a))/a fixed imposes dx/dlog(a)=Delta/Delta'.
    # This necessary source-level tangent is not an accepted RG beta function.
    tangent=quotient(gap,gap_derivative)
    return {'proof_coordinate_key':key,'g_squared_exact':str(x),
        'ground_energy_derivative_interval_exact':strings(derivatives[0]),
        'first_excitation_derivative_interval_exact':strings(derivatives[1]),
        'gap_derivative_interval_exact':strings(gap_derivative),'gap_derivative_interval':display(gap_derivative),
        'excited_projector_distance_upper_exact':str(distance),
        'constant_physical_gap_required_dx_dlog_a_interval_exact':strings(tangent),
        'constant_physical_gap_required_dx_dlog_a_interval':display(tangent),
        'source_gap_interval_exact':record['full_Haar_gap_interval_exact'],
        'source_covariance_interval_exact':record['connected_Wilson_covariance_interval_exact'],
        'slope_is_a_renormalization_group_beta_function':False}


def run():
    keys=sorted(witnesses(),key=lambda k:F(witnesses()[k]['g_squared_exact']))
    rows=[derivative_record(k)for k in keys]
    steps=[]
    for left,right in zip(rows,rows[1:]):
        x,y=F(left['g_squared_exact']),F(right['g_squared_exact'])
        a,b=tuple(map(F,left['source_gap_interval_exact'])),tuple(map(F,right['source_gap_interval_exact']))
        secant=tuple(v/(y-x)for v in subtract(b,a));assert secant[0]>0
        ratio=quotient(b,a)
        steps.append({'from_g_squared_exact':str(x),'to_g_squared_exact':str(y),
            'finite_step_gap_ratio_interval_exact':strings(ratio),
            'finite_step_gap_slope_interval_exact':strings(secant),
            'physical_spacing_ratio_required_for_equal_gap_interval_exact':strings(ratio)})
    return {'schema':'oph.full_Haar_coupling_family.v1',
        'verdict':'PASS_FULL_SOURCE_COUPLING_RESPONSE_AND_CLOCK_TANGENTS',
        'declared_couplings_exact':[r['g_squared_exact']for r in rows],
        'records':rows,'finite_coupling_steps':steps,
        'all_full_spectrum_receipts':[calculate(k)[0]for k in keys],
        'exact_derivative_identity':'dE_k/dx=E_k/x-4<W_magnetic>_k/x², W_magnetic=2-W1-W2, x=g², a=1',
        'all_tested_gap_slopes_positive':True,
        'input_nature_data_used':False,'a_GeV_conversion_fitted':False,
        'physical_scope':'A common full-Haar finite graph with declared source clock. The constant-mass tangent is a necessary level-set condition, not a source-derived running coupling or a four-dimensional continuum trajectory.',
        'spatial_refinement_performed':False,'Yang_Mills_continuum_constructed':False}
