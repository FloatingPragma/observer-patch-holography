"""Controlled, computed Feshbach correction from every omitted Haar channel."""
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter
import sympy as sy
from model.openwave_refined_haar_source.gauss import basis,electric,total_action,compressed_sum
from model.openwave_refined_haar_source.source import exact_checks
from model.openwave_refined_haar_spectrum.certificate import offblock_bound
from model.openwave_coupled_haar_spectrum.certificate import source_digest
from model.openwave_euclidean_refinement.response import exponential_decay_interval,dyadic_outer


def strings(x):return [str(v)for v in dyadic_outer(x)]


@lru_cache(None)
def correction():
    states=basis();retained=set(states);rows={};alpha,beta=F(16),F(1,4)
    for j,column in enumerate(total_action()):
        for target,value in column:
            if target not in retained:rows.setdefault(target,[]).append((j,value))
    entries={}
    for target,row in rows.items():
        energy=electric(target);assert energy>=5
        denominator=alpha*energy+36*beta
        for j,x in row:
            for k,y in row:
                entries[j,k]=entries.get((j,k),0)+sy.Rational(beta**2/denominator)*x*y
    entries={k:sy.simplify(v)for k,v in entries.items()if v}
    matrix=sy.SparseMatrix(len(states),len(states),entries);assert matrix==matrix.T
    assert all(matrix[i,i]>=0 for i in range(len(states)))
    assert all(matrix[0,j]==0 for j in range(len(states)))
    return matrix,rows


def probe_defects(time):
    time=F(time)
    if time<0:raise ValueError('positive Euclidean time required')
    # The exact norm is calculated from the full-source Haar moments, not
    # a spectral cutoff or a Taylor approximation asserted to converge.
    alpha,beta=F(16),F(1,4)
    norm2=81*alpha**2*beta**2+3888*alpha*beta**3+F(6999363,4)*beta**4
    assert norm2==F(9321795,1024)<96**2
    boundary_lower=max(F(0),3*beta*time-48*time*time)**2
    boundary_upper=(3*beta*time+48*time*time)**2
    retained_upper=norm2*time**4/4
    return {'time_exact':str(time),'boundary_only_vacuum_semigroup_defect_interval_exact':strings((boundary_lower,boundary_upper)),
        'retained_201_mode_vacuum_semigroup_defect_interval_exact':strings((F(0),retained_upper)),
        'retained_upper_divided_by_boundary_lower_exact':str(retained_upper/boundary_lower)if boundary_lower else None,
        'comparison_is_same_Haar_vacuum_probe':True,'interacting_vacuum_or_operator_norm_improvement_claimed':False}


def run():
    checks=exact_checks();matrix,rows=correction();alpha,beta=F(16),F(1,4)
    S=compressed_sum();v=S[:,0];assert (v.T*S*v)[0]==0
    s_bound,_=offblock_bound(basis(),total_action());b2=beta**2*s_bound
    d=F(80);d0=F(89);magnetic=F(9)
    error=b2*magnetic/(d0*(d0-magnetic))
    nonzero=matrix.todok();trace=sy.trace(matrix)
    assert trace>0 and sum(i!=j for i,j in nonzero)>0
    assert checks['discarded_S_squared_vacuum_norm_squared_exact']=='303/2'
    controls=[probe_defects(F(1,128)),probe_defects(F(1,1024))]
    assert F(controls[1]['retained_upper_divided_by_boundary_lower_exact'])<F(1,225)
    memory=[]
    for t in(F(1,128),F(1,64),F(1,32),F(1,8)):
        decay=exponential_decay_interval(d*t)
        # Integral of the norm-bound kernel is bounded without a time grid.
        tail=b2*decay[1]/d
        norm_difference=min(b2*t*t/2,b2*(t/d-(1-decay[1])/d**2))
        memory.append({'source_time_exact':str(t),'memory_kernel_norm_upper_exact':str(dyadic_outer((F(0),b2*decay[1]))[1]),
            'integrated_memory_tail_after_time_upper_exact':str(dyadic_outer((F(0),tail))[1]),
            'compressed_transfer_minus_autonomous_PHP_norm_upper_exact':str(dyadic_outer((F(0),norm_difference))[1])})
    return {'schema':'oph.full_Haar_retained_memory_closure.v1','verdict':'PASS_COMPUTED_RETAINED_SELF_ENERGY_WITH_FULL_COMPLEMENT_ERROR',
        'retained_dimension':201,'all_full_offblock_states':len(rows),'full_QHQ_lower_exact':str(d),
        'full_offblock_norm_squared_upper_exact':str(b2),
        'exact_Feshbach_operator':'F(z)=PHP-z-PHQ(QHQ-z)^(-1)QHP for z<80',
        'actual_computed_self_energy':{'definition':'Sigma0_leading=PHQ (alpha E_Q+36 beta)^(-1) QHP',
            'exact_sparse_matrix_sha256':source_digest(matrix),'nonzero_entries':len(nonzero),
            'offdiagonal_entries':sum(i!=j for i,j in nonzero),'trace_exact':str(trace),
            'positive_semidefinite_by_complete_weighted_Gram':True,'vacuum_row_exactly_zero':True,
            'first_four_face_diagonal_values_exact':[str(matrix[i,i])for i in range(1,5)],
            'positive_nonzero_correction_executed':True,
            'infinite_complement_Neumann_ratio_upper_exact':str(magnetic/d0),
            'full_Sigma0_minus_computed_matrix_norm_upper_exact':str(error),
            'full_Sigma_z_minus_computed_matrix_norm_bound':'error + b2*z/(80*(80-z)) for 0<=z<80',
            'all_Q_representations_and_multiplicities_in_remainder':True},
        'exact_memory_equation':'Uprime(t)=-A U(t)+integral_0^t Bstar exp(-(t-s)D) B U(s) ds; A=PHP,D=QHQ,B=QHP',
        'kernel_bound':'||Bstar exp(-tD) B|| <= b2 exp(-80t)',
        'memory_tail_records':memory,'same_clock_semigroup_defect_controls':controls,
        'exact_retained_vacuum_H_first_step':True,
        'retained_vacuum_second_step_leakage_norm_squared_exact':str(beta**4*F(303,2)),
        'boundary_only_autonomous_semigroup_rejected':True,
        'retained_201_modes_exactly_autonomous':False,
        'closed_finite_premise':'An explicit computed leading energy-dependent coarse correction with an all-representation norm remainder, rather than discarded interior dynamics',
        'uniform_refinement_or_physical_4D_limit_proved':False,'new_empirical_fit_inputs':[]}
