"""Complete positive elementary-removal measures and compensated local bounds."""
from functools import lru_cache
import sympy as s
from ..charged_ground_exchange.model import transfer,weight_basis
from ..charged_ground_exchange.resolvent import field_data,field_value,polynomial,encoded
from ..ground_composite_response.model import source,reduced
from ..local_reference import physical
from ..local_reference.reduction import quark_module,lepton_operators,coordinate_kernel
from ..local_reference.characters import certificate as character_certificate
from ..atomic_reference_control.ground_exchange import actual_spin_coefficients
from ..family_cell.quotient import add,mul,divide,evaluate,outward
from ...chiral_higgs.model import higgs_memory


@lru_cache(maxsize=None)
def positive_response(family):
    if family not in ('ec','L'):raise ValueError('actual ec or L elementary-removal channel required')
    e,f,field=field_data();sb=source(1);h=sb['H'];x=s.Symbol('x')
    p=s.Poly(h.charpoly(x).as_expr(),x).sqf_part()
    blocks=[transfer(family,spin) for spin in range(2)]
    q=s.Poly(1,x)
    for block in blocks:q=s.lcm(q,s.Poly(block['H'].charpoly(x).as_expr(),x).sqf_part())
    coefficients=[];powers=[s.eye(h.rows)]
    for i in range(p.degree()):
        c=s.Poly(sum(p.nth(j)*x**(j-i-1) for j in range(i+1,p.degree()+1)),x)
        coefficients.append(field_value(c)/field_value(p.diff()))
        if i+1<p.degree():powers.append(h*powers[-1])
    images=[block['map'] for block in blocks];moments=[]
    for _ in range(q.degree()):
        pullback=sum((block['adjoint']*image for block,image in zip(blocks,images)),s.zeros(h.rows))
        moments.append(sum((coefficient*field.convert(s.trace(power*pullback))
            for coefficient,power in zip(coefficients,powers)),field.zero))
        images=[block['H']*image for block,image in zip(blocks,images)]
    assert moments[0]==field.one
    entries={}
    for i,moment in enumerate(moments):
        for j in range(i+1,q.degree()+1):
            entries[j-i-1]=entries.get(j-i-1,field.zero)+field.convert(q.nth(j))*moment
    z=s.Symbol('z')
    numerator=s.Poly.from_dict({(k,):v for k,v in entries.items()},z,domain=field)
    denominator=s.Poly.from_dict({(k,):field.convert(q.nth(k)) for k in range(q.degree()+1)},z,domain=field)
    roots=q.intervals(eps=s.Rational(1,10**70))
    e_interval=tuple(map(s.Rational,actual_spin_coefficients()['ground_energy_isolating_interval']))
    rows=[]
    for interval,multiplicity in roots:
        assert multiplicity==1
        interval=tuple(map(s.Rational,interval));value=(s.Integer(0),s.Integer(0))
        for k in range(numerator.degree(),-1,-1):
            coefficient=evaluate(polynomial(numerator.rep.nth(k)),e_interval)
            value=add(mul(value,interval),coefficient)
        weight=divide(value,evaluate(q.diff(),interval))
        assert 0<weight[0]<=weight[1]<=1
        instrument=(weight[0]/2,weight[1]/2)
        rows.append({'energy_interval':interval,'weight_interval':instrument,
            'source_energy_difference_interval':(interval[0]-e_interval[1],interval[1]-e_interval[0])})
    assert len(rows)==q.degree() and all(row['source_energy_difference_interval'][0]>0 for row in rows)
    assert sum(row['weight_interval'][0] for row in rows)<=s.Rational(1,2)<=sum(row['weight_interval'][1] for row in rows)
    return {'q':q,'numerator':numerator,'denominator':denominator,'moments':moments,
        'rows':rows,'e_interval':e_interval,'blocks':blocks}


@lru_cache(maxsize=None)
def compensated_blocks(family):
    if family not in ('ec','L'):raise ValueError('actual color-singlet reference irrep required')
    weak,hyper=(0,-2) if family=='ec' else (1,1)
    q=quark_module();l=lepton_operators();scalar=higgs_memory()[1];records=[]
    x=s.Symbol('x')
    for number in range(7):
        for spin in range(4):
            keys=weight_basis(number,spin,weak,hyper)
            if not keys:continue
            wm=physical.tensor_map(keys,weight_basis(number,spin,weak+2,hyper),
                quark=q['currents'][0]+s.I*q['currents'][1],
                lepton=l['currents'][0]+s.I*l['currents'][1],scalar=scalar[8]+s.I*scalar[9])
            sm=physical.tensor_map(keys,weight_basis(number,spin+2,weak,hyper),
                lepton=l['spin'][0]+s.I*l['spin'][1])
            basis,_=coordinate_kernel(wm.col_join(sm))
            expected=(len(keys)-len(weight_basis(number,spin,weak+2,hyper))
                -len(weight_basis(number,spin+2,weak,hyper))+len(weight_basis(number,spin+2,weak+2,hyper)))
            assert basis.cols==expected
            if not basis.cols:continue
            h=reduced(keys,basis)['H'];p=s.Poly(h.charpoly(x).as_expr(),x)
            interval=tuple(map(s.Rational,p.sqf_part().intervals(eps=s.Rational(1,10**60))[0][0]))
            records.append({'transformed_lepton_number':number,'raw_L':number-2,
                'twice_source_spin':spin,'complete_highest_spin_dimension':h.rows,
                'complete_characteristic_polynomial':str(p.as_expr()),'ground_interval':interval})
    minimum=min(records,key=lambda row:row['ground_interval'][0])['ground_interval']
    excluded=character_certificate()
    assert excluded['all_other_color_admissible_quark_Casimir_upper_bound']=='42'
    lower=s.Rational(excluded['excluded_physical_sectors_energy_lower_bound'])
    ground=positive_response(family)['e_interval']
    gap=min(lower,minimum[0])-ground[1]
    assert gap>0
    return {'blocks':records,'maximal_module_ground_interval':minimum,
        'every_other_color_singlet_module_lower':lower,'complete_compensated_local_gap_lower':gap,
        'global_compensated_ground_identified':bool(minimum[1]<lower)}


def certificate(family):
    data=positive_response(family);sector=compensated_blocks(family)
    return {'family':family,'native_internal_dimension':1 if family=='ec' else 2,
        'probe':'B_r=sum_(spin,a in r) r^dagger_(spin,a) c_(origin,spin,a)/sqrt(2 d_r), on reference vacuum',
        'reference_orientation':'same native representation for removal; opposite orientation would violate hypercharge',
        'source_reference':'actual local physical ground doublet, with reduced density I/2',
        'actual_target_dimensions_by_removed_spin':[block['H'].rows for block in data['blocks']],
        'complete_target_annihilator':str(data['q'].as_expr()),
        'positive_total_raw_moments_coefficients_in_E0':[encoded(value) for value in data['moments']],
        'complete_local_energy_records':[{'energy_interval':outward(row['energy_interval'],60),
            'source_energy_difference_interval':outward(row['source_energy_difference_interval'],60),
            'pump_probability_interval':outward(row['weight_interval'],80)} for row in data['rows']],
        'all_visible_weights_strictly_positive':True,'exact_pump_success_probability':'1/2',
        'exact_pump_failure_probability':'1/2',
        'outside_complete_local_response_support':{'probability':'0','retained':True},
        'complete_compensated_local_maximal_module_blocks':[{**row,'ground_interval':outward(row['ground_interval'],60)} for row in sector['blocks']],
        'every_other_color_admissible_quark_module_lower':str(sector['every_other_color_singlet_module_lower']),
        'full_compensated_local_gap_above_neutral_ground_lower':str(sector['complete_compensated_local_gap_lower']),
        'lowest_compensated_local_energy_exactly_identified':sector['global_compensated_ground_identified'],
        'local_static_reference_with_zero_energy_has_lower_state_than_neutral_ground':False,
        'all_raw_lepton_charges_and_source_spins_compared':True,
        'higher_reference_occupations_covered_by_this_local_single_reference_theorem':False}
