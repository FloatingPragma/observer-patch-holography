"""Independent covariance, exact cochain and complete-instrument checks."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import json
import numpy as np
import mpmath as mp
from scipy.linalg import qr

ROOT=Path(__file__).resolve().parents[2]


def covariance_review(receipt):
    directory=ROOT/'model/openwave_capacity_stress/inputs'
    systematic=[[F(x)for x in line.split()]for line in(directory/'data_MM20.dat').read_text().splitlines()if line and not line.startswith('#')]
    def interpolate(z,column):
        if z<=systematic[0][0]:return systematic[0][column]
        if z>=systematic[-1][0]:return systematic[-1][column]
        a,b=next((a,b)for a,b in zip(systematic,systematic[1:])if a[0]<=z<=b[0])
        return ((b[0]-z)*a[column]+(z-a[0])*b[column])/(b[0]-a[0])
    out=[]
    with mp.workdps(70):
        cv=lambda x:mp.mpf(x.numerator)/x.denominator
        for record in receipt['comparisons']:
            rows=[tuple(F(x)for x in line.split(',')[:3])for line in(directory/record['table']).read_text().splitlines()if line and not line.startswith('#')]
            assert len(rows)==record['rows']==15
            assert record['fitted_chronometer_parameters']==0
            u=[[h*interpolate(z,c)/100 for z,h,s in rows]for c in(1,record['sps_column'])]
            C=[[((rows[i][2]**2)if i==j else F(0))+sum(col[i]*col[j]for col in u)for j in range(15)]for i in range(15)]
            assert C==[[F(x)for x in row]for row in record['full_covariance_rational']]
            matrix=mp.matrix([[cv(x)for x in row]for row in C])
            residual=mp.matrix([mp.mpf('67.4')*mp.sqrt(mp.mpf('.315')*(1+cv(z))**3+mp.mpf('.685'))-cv(h)for z,h,s in rows])
            chi=(residual.T*mp.lu_solve(matrix,residual))[0]
            lo,hi=map(mp.mpf,record['chi_square']);assert lo<=chi<=hi
            diagonal=sum(residual[i]**2/matrix[i,i]for i in range(15))
            assert abs(diagonal-chi)>mp.mpf('.01')
            out.append({'table':record['table'],'systematic_column':record['sps_column'],
                'complete_correlated_chi_square_70digit_diagnostic':mp.nstr(chi,60),
                'inside_exact_owner_interval':True,'offdiagonal_omission_changes_result':True})
    assert receipt['openwave_row28_promotion'] is False
    for cycle in receipt['complete_cycles']:
        assert sum(map(F,cycle['work_strokes']))==0
        assert sum(F(r['probability'])for r in cycle['complete_appended_reference_records'])==1
    return out


def cube(cells):
    vertices=list(product(range(cells+1),repeat=3));vi={x:i for i,x in enumerate(vertices)}
    edges=[]
    for a in range(3):
        for x in vertices:
            if x[a]<cells:edges.append((a,x))
    ei={x:i for i,x in enumerate(edges)};G=np.zeros((len(edges),len(vertices)),dtype=np.int64)
    for i,(a,x) in enumerate(edges):
        y=list(x);y[a]+=1;G[i,vi[x]]=-1;G[i,vi[tuple(y)]]=1
    rows=[]
    for a,b in combinations(range(3),2):
        for x in vertices:
            if x[a]>=cells or x[b]>=cells:continue
            row=np.zeros(len(edges),dtype=np.int64);xa=list(x);xa[a]+=1;xb=list(x);xb[b]+=1
            for edge,sgn in(((a,x),1),((b,tuple(xa)),1),((a,tuple(xb)),-1),((b,x),-1)):row[ei[edge]]=sgn
            rows.append(row)
    C=np.array(rows);assert not np.any(C@G)
    return G,C


def packet_guard(G,physical,allowed=1e-9):
    assert abs(np.vdot(physical,physical).real-1)<1e-11
    assert np.linalg.norm(G.T@physical)<allowed


def packet_review(receipt):
    cells=receipt['source']['complete_spatial_cube_cells'];G,C=cube(cells)
    assert G.shape==(882,343) and C.shape==(756,882)
    Q,_=qr(G[:,1:].astype(float),mode='full');frame=Q[:,342:]
    assert frame.shape==(882,540)
    proposed=np.array([complex(a,b)for a,b in receipt['source']['packet_integer_coordinates']])
    coefficients=frame.T@proposed;ratio=np.linalg.norm(proposed)/np.linalg.norm(coefficients)
    physical=frame@coefficients/np.linalg.norm(coefficients);packet_guard(G,physical)
    for trial in receipt['full_source_executions']:
        c=trial['packet_certificate'];assert abs(c['accepted_to_projected_norm_ratio_diagnostic']-ratio)<1e-11
        assert 2*1e-9*ratio< c['normalized_packet_error_upper']<=1e-8
        allocation=trial['extra_vector_error_allocation'];assert sum(allocation.values())<=5.000001e-8
    bad=physical+G[:,0]*.01;bad/=np.linalg.norm(bad)
    try:packet_guard(G,bad)
    except AssertionError:pass
    else:raise AssertionError('longitudinal packet accepted')
    return {'vertices':343,'edges':882,'faces':756,'transverse_dimension':540,
        'integer_curl_gradient_zero':True,'all_boundary_Gauss_constraints':True,
        'independent_packet_norm_ratio_diagnostic':float(ratio),'longitudinal_packet_rejected':True}


def record_probability(record):
    if 'all4_setting4_native2_geometry_vacuum_and_every_edge_probabilities'in record:
        p=np.array(record['all4_setting4_native2_geometry_vacuum_and_every_edge_probabilities'])
        assert p.shape==(4,4,2,883)
    else:p=np.array(record['complete_native_geometry_marginal']);assert p.shape==(4,4,2)
    assert p.min()>=0
    assert abs(p.sum()-record['total_probability'])<2e-10
    assert abs(record['postmeasurement_source_energy']-record['premeasurement_source_energy']-record['complete_absorption_and_clock_readout_work'])<2e-10
    return float(p.sum())


def trial_review(trial):
    for terminal in trial['terminal_joint_records']:assert abs(record_probability(terminal)-1)<2e-10
    arrival=sum(record_probability(x)for x in trial['complete_first_arrival_histories'])
    failure=record_probability(trial['retained_no_arrival_and_photon_vacuum'])
    assert failure>=.5-2e-10
    assert abs(arrival+failure-1)<2e-10
    assert abs(arrival-trial['arrival_success_probability'])<2e-10
    work=trial['complete_history_work']
    coarse=sum(x['coarse_detector_work']for x in trial['complete_first_arrival_histories'])
    fine=sum(x['complete_absorption_and_clock_readout_work']for x in trial['complete_first_arrival_histories'])+trial['retained_no_arrival_and_photon_vacuum']['complete_absorption_and_clock_readout_work']
    assert abs(coarse-work['coarse_detector'])<2e-10
    assert abs(fine-work['terminal_and_native_readout'])<2e-10
    assert abs(work['final']-work['initial']-coarse-fine)<2e-9
    prep=trial['complete_preparation_and_work'];assert F(prep['input_energy_exact'])==0
    assert abs(prep['mean_preparation_work']-work['initial'])<2e-10
    assert abs(prep['mean_total_work']-work['initial']-coarse-fine)<2e-9
    assert abs(prep['mean_total_work']-work['final'])<2e-9
    return {'source_level':trial['source_native_level'],'gamma':trial['gamma_declared'],
        'all_terminal_and_twelve_arrival_tables_checked':True,'failure_probability':failure,
        'complete_energy_budget_checked':True}


def gravity_review(receipt):
    packets=packet_review(receipt)
    trials=list(receipt['full_source_executions'])
    for control in receipt['spatial_metric_controls']:trials+=control['complete_executions']
    trials+=receipt['second_physical_polarization_control']['complete_executions']
    out=[trial_review(t)for t in trials]
    c=receipt['same_source_photon_comparison']
    assert c['conditional_arrival_delay']>c['arrival_delay_error_upper']
    assert c['separate_Fourier_direction_deflection']+c['direction_error_upper']<0
    assert c['transverse_position_difference']+c['position_difference_error_upper']<0
    assert all(0<lo<hi<1 for lo,hi in c['same_joint_native_clock_phase_ratio_intervals'])
    u=receipt['independent_energy_length_attachment']
    assert u['physical_G_strength_match']is False and u['physical_common_light_speed_match']is False
    assert u['conditional_to_observed_strength_ratio']>1e39
    assert receipt['source']['canonical_metric_QED_Bogoliubov_transfer_claimed']is False
    return {'exact_carrier_and_packet':packets,'complete_executions':out,
        'positive_delay_inward_bending_and_both_joint_clock_redshifts_certified':True,
        'failed_physical_energy_length_attachment_retained':True}


def run():
    a=json.loads((ROOT/'model/openwave_capacity_stress/results.json').read_text())
    b=json.loads((ROOT/'model/openwave_transverse_gravity/results.json').read_text())
    return {'schema':'oph.independent_quantum_observable_review.v1',
        'verdict':'PASS_INDEPENDENT_COVARIANCE_GAUSS_PACKET_AND_COMPLETE_RECORD_REVIEW',
        'capacity_comparisons':covariance_review(a),'transverse_photon':gravity_review(b),
        'scope':'independent arithmetic and instrument review; no physical G, Lambda amplitude, continuum Einstein or official scoreboard award'}
