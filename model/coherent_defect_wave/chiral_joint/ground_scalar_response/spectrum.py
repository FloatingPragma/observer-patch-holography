"""Exact two-level scalar weights with every zero-energy record retained."""
from functools import lru_cache
import sympy as s

from .model import matrices
from ..ground_composite_response.spectrum import projector_coefficients
from ..family_cell.quotient import add,mul,divide,evaluate,outward
from ..local_reference.spectrum import certificate as ground_spectrum


@lru_cache(maxsize=1)
def certificate():
    data=matrices();h=data['H'];p=data['minimal_polynomial'];x=p.gen
    factor=s.Poly((x*(x-1))**3-313*(x*(x-1))**2+31708*x*(x-1)-1041468,x)
    assert p.rem(factor).is_zero and factor.is_irreducible
    coefficients=projector_coefficients(p,factor);denominator=p.diff().rem(factor)
    numerators={name:s.Integer(0) for name in ('I','N','Y','S','current')}
    power=s.eye(32)
    for coefficient in coefficients:
        for name in numerators:
            op=s.eye(32) if name=='I' else data[name]
            numerators[name]+=coefficient.as_expr()*s.trace(power*op)
        power=h*power
    numerators={name:s.Poly(value,x).rem(factor) for name,value in numerators.items()}
    assert numerators['I']==2*denominator
    assert (numerators['N']*s.Poly(2*x-1,x)-2*s.Poly(x,x)*denominator).rem(factor).is_zero
    assert numerators['Y']==(2*s.Poly(x,x)*denominator-numerators['N']).rem(factor)
    assert numerators['S']==(2*s.Poly(x*(x-1),x)*denominator).rem(factor)
    assert numerators['current'].is_zero
    roots=s.polys.polytools.intervals(factor,eps=s.Rational(1,10**110))
    assert len(roots)==6 and all(m==1 for interval,m in roots)
    ground=tuple(map(s.Rational,roots[0][0]));declared=ground_spectrum()['global_physical_ground_energy_interval']
    assert s.Rational(declared[0])<=ground[0]<=ground[1]<=s.Rational(declared[1])
    xi=mul(ground,add(ground,(-1,-1)))
    omega=add((1,1),mul((-2,-2),ground))
    mean_n=divide(ground,add(mul((2,2),ground),(-1,-1)))
    norm=add((1,1),mul((4,4),xi))
    variance=divide(xi,norm)
    y_elastic=divide(mul((4,4),mul(xi,xi)),norm)
    assert 0<mean_n[0]<mean_n[1]<1 and 0<variance[0]<variance[1]<s.Rational(1,4)
    full=s.Poly(h.charpoly(x).as_expr(),x)
    factors=s.factor_list(full)[1]
    records=[]
    for target_factor,rank in factors:
        target_factor=target_factor.monic()
        for root_index,(interval,multiplicity) in enumerate(s.polys.polytools.intervals(target_factor,eps=s.Rational(1,10**110))):
            assert multiplicity==1
            interval=tuple(map(s.Rational,interval))
            same_factor=target_factor.as_expr()==factor.as_expr()
            role=('elastic' if same_factor and root_index==0 else
                  'inelastic' if same_factor and root_index==5 else 'zero')
            if role=='elastic':excitation=(s.Rational(0),s.Rational(0))
            elif role=='inelastic':
                excitation=omega
                assert max(interval[0],1-ground[1])<=min(interval[1],1-ground[0])
            else:
                excitation=(interval[0]-ground[1],interval[1]-ground[0])
                assert excitation[0]>0
            records.append({'target_energy_factor':str(target_factor.as_expr()),'target_root_index':root_index,
                'target_energy_interval':outward(interval,48),'excitation_energy_interval':outward(excitation,48),
                'complete_spectral_rank':rank,'role':role})
    records.sort(key=lambda row:s.Rational(row['target_energy_interval'][0]))
    assert len(records)==9 and sum(row['complete_spectral_rank'] for row in records)==32
    assert [row['complete_spectral_rank'] for row in records if row['role']!='zero']==[2,2]
    settings=[]
    for name in ('scalar_number','total_Yukawa'):
        scale=s.Rational(1,8) if name=='scalar_number' else s.Rational(1,512)
        elastic=mul(mul(mean_n,mean_n) if name=='scalar_number' else y_elastic,(scale,scale))
        inelastic=mul(variance,(scale,scale))
        total=mul(mean_n if name=='scalar_number' else xi,(scale,scale))
        assert max(elastic[0]+inelastic[0],total[0])<=min(elastic[1]+inelastic[1],total[1])
        failure=(s.Rational(1,2)-total[1],s.Rational(1,2)-total[0])
        assert 0<elastic[0] and 0<inelastic[0] and 0<failure[0]
        output=[]
        for record in records:
            probability={'elastic':elastic,'inelastic':inelastic,'zero':(s.Rational(0),s.Rational(0))}[record['role']]
            output.append({**record,'joint_probability_interval':outward(probability,48)})
        settings.append({'setting':name,'setting_probability':'1/2',
            'success_operator':'N/2' if name=='scalar_number' else 'Y/16',
            'complete_energy_records':output,'joint_success_probability_interval':outward(total,48),
            'joint_failure_probability_interval':outward(failure,48),
            'outside_spin_half_maximal_module_record':{'probability_interval':['0','0'],'retained':True},
            'all_zero_energy_outcomes_retained':True,'conditioned_success_renormalization_applied':False})
    probabilities=[record['joint_probability_interval'] for setting in settings for record in setting['complete_energy_records']]
    probabilities.extend(setting['joint_failure_probability_interval'] for setting in settings)
    assert sum(s.Rational(row[0]) for row in probabilities)<=1<=sum(s.Rational(row[1]) for row in probabilities)
    return {'reference':'P_ground/2','complete_target_spin_half_dimension':32,'settings':settings,
        'complete_joint_records':22,'positive_probability_records':6,'exact_zero_probability_records':16,
        'ground_mean_scalar_excitation_interval':outward(mean_n,48),
        'scalar_inelastic_transition_interval':outward(omega,48),
        'unscaled_scalar_inelastic_residue_interval':outward(variance,48),
        'unscaled_scalar_current_inelastic_residue_interval':outward(xi,48),
        'exact_formula_variables':{'e':'global physical ground energy','xi':'e(e-1)','omega':'1-2e'},
        'exact_ground_scalar_mean':'e/(2e-1)',
        'exact_unscaled_scalar_residue':'xi/(1+4xi)',
        'exact_N_elastic_weight':'(e/(2e-1))^2',
        'exact_Y_elastic_weight':'4xi^2/(1+4xi)',
        'exact_raw_retarded_NN_and_YY_response':'2 omega xi / ((1+4xi)(z^2-omega^2)), Im(z)>0',
        'exact_raw_retarded_NY_response':'minus the NN response',
        'retarded_zero_frequency_pole':False,'positive_and_negative_transition_poles':['omega','-omega'],
        'transition_frequency_identified_with_relativistic_rest_mass':False,
        'source_projector_trace_and_means_checked_exactly':True,
        'scalar_current_stationary_mean':'0','numerical_eigensolver_used':False}
