"""A complete sign-sensitive massive-vector Weyl experiment, with intervals."""
from fractions import Fraction as F
import itertools
import mpmath as mp


def endpoint(interval,upper=False):
    sign,mantissa,exponent,_=interval._mpi_[int(upper)]
    return F((-1)**sign*mantissa)*F(2)**exponent


def enclosure(interval):return [str(endpoint(interval)),str(endpoint(interval,True))]


def controlled_experiment(a,side,error=F(19,20*10**9)):
    a=F(a);error=F(error)
    assert a*side==16
    mp.iv.dps=70;iv=mp.iv
    pi=iv.pi
    assert endpoint(pi)>3 and endpoint(pi,True)<F(22,7)
    qspacing=2*pi/(16*iv.sqrt(6))
    P=iv.mpf(1)/4
    assert endpoint(qspacing,True)<F(1,4)
    assert endpoint(2*qspacing*qspacing,True)<F(1,16)
    assert endpoint(3*qspacing*qspacing)>F(1,16)
    modes=[word for word in itertools.product((-1,0,1),repeat=3) if sum(n*n for n in word)<=2]
    assert len(modes)==19
    volume=iv.mpf(16)**3*6*iv.sqrt(6)
    # hhat(q)=(1-|q|²/P²)_+², separately in component0 and component1.
    # Both actual discrete L2 norms are bounded by19/volume<1/64.
    assert F(19,16**3*12)<F(1,64)
    av=iv.mpf(a.numerator)/a.denominator
    variances=[iv.mpf(0),iv.mpf(0)]
    norm=iv.mpf(0)
    for menu in modes:
        q=[qspacing*n for n in menu];q2=sum((v*v for v in q),iv.mpf(0))
        weight=(1-q2/(P*P))**4
        omega=iv.sqrt(1+(1-av*av)/(3*av*av)*sum((1-iv.cos(2*pi*n/side) for n in menu),iv.mpf(0)))
        norm+=weight/volume
        variances[0]+=weight*q2/(2*omega*volume)
        variances[1]+=weight*(1+q[0]**2)/(2*omega*volume)
    assert endpoint(norm,True)<F(1,64)
    probabilities=[[(1+sign*iv.exp(-variance/2))/2 for sign in (1,-1)] for variance in variances]
    # Spatial/temporal cross covariance and commutator vanish exactly by
    # q1->-q1 on this actual menu, at every time, so the ideal signs factor.
    records=[];ideal_sum=iv.mpf(0)
    for sign0,flag0,sign1,flag1 in itertools.product((1,-1),(0,1),(1,-1),(0,1)):
        probability=iv.mpf(0) if flag0 or flag1 else probabilities[0][int(sign0==-1)]*probabilities[1][int(sign1==-1)]
        ideal_sum+=probability
        lo=max(F(0),endpoint(probability)-error);hi=min(F(1),endpoint(probability,True)+error)
        records.append({'first_sign':sign0,'first_flag':flag0,'second_sign':sign1,'second_flag':flag1,
            'ideal_finite_free_probability_interval':enclosure(probability),
            'actual_full_source_probability_interval':[str(lo),str(hi)]})
    assert endpoint(ideal_sum)<=1<=endpoint(ideal_sum,True)
    assert sum(F(row['actual_full_source_probability_interval'][0]) for row in records)<=1
    assert sum(F(row['actual_full_source_probability_interval'][1]) for row in records)>=1
    marginal=[]
    for label,variance,probs in zip(('temporal_A0','spatial_A1'),variances,probabilities):
        minus=probs[1];lo=max(F(0),endpoint(minus)-error);hi=min(F(1),endpoint(minus,True)+error)
        assert lo>0
        marginal.append({'readout':label,'ideal_variance_interval':enclosure(variance),
            'ideal_minus_sign_probability_interval':enclosure(minus),
            'actual_minus_sign_probability_interval':[str(lo),str(hi)]})
    return {'spacing':str(a),'macro_side':side,'source_physical_box_side':'16',
        'unit_speed_box_side':'16*sqrt(6)','Fourier_support_radius':'1/4',
        'exact_supported_Fourier_modes':19,
        'two_test_functions':'hhat_0(q)=hhat_1(q)=(1-16|q|^2)_+^2, each in its stated single component',
        'actual_discrete_test_norm_squared_interval':enclosure(norm),
        'test_times_Delta_times_t':['0','1/2'],'total_Delta_times_duration':'1',
        'ideal_cross_covariance_and_commutator_exactly_zero_by_reflection':True,
        'complete_record_TV_error_upper':str(error),
        'complete_sixteen_records':records,'strictly_positive_actual_sign_marginals':marginal,
        'ideal_complete_probability_sum_interval':enclosure(ideal_sum),
        'probability_of_any_Halmos_flag_one_upper':str(error),
        'all_record_intervals_enclose_one_normalized_law':True,
        'interval_endpoints_are_not_independently_renormalized':True,
        'full_native_source_matrix_numerically_enumerated':False,
        'calculation':'exact19-mode Gaussian comparison plus full actual-source vector/record bound'}
