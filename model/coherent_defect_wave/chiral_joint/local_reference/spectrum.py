"""Exact complete maximal-module levels and the physical global ground gap."""
from collections import defaultdict
from functools import lru_cache
import sympy as s

from .physical import certificate as physical_certificate
from .characters import certificate as character_certificate
from ..family_cell.quotient import add,mul,sqrt_interval,outward


def ordered_records(records,dimension):
    """Reject incomplete counts and intervals which do not prove strict ordering."""
    records=sorted(records,key=lambda row:s.Rational(row['energy_interval'][0]))
    if sum(row['physical_rank'] for row in records)!=dimension:
        raise ValueError('incomplete physical spectral ranks')
    if any(row['physical_rank']<=0 for row in records):
        raise ValueError('nonpositive physical spectral rank')
    for left,right in zip(records,records[1:]):
        if s.Rational(left['energy_interval'][1])>=s.Rational(right['energy_interval'][0]):
            raise ValueError('energy intervals do not prove strict spectral separation')
    return records


@lru_cache(maxsize=1)
def certificate():
    physical=physical_certificate();excluded=character_certificate();x=s.Symbol('x')
    roots={};ranks=defaultdict(int);sectors=defaultdict(list);dark=0;dark_sectors=[]
    for block in physical['complete_spin_blocks']:
        polynomial=s.Poly(block['squared_transfer_polynomial'],x)
        sector={name:block[name] for name in ('raw_L','twice_spin','spin_multiplicity')}
        for factor,multiplicity in s.factor_list(polynomial)[1]:
            factor=factor.monic();label=str(factor.as_expr())
            intervals=s.polys.polytools.intervals(factor,eps=s.Rational(1,10**80))
            assert sum(m for _,m in intervals)==factor.degree()
            for index,(interval,m) in enumerate(intervals):
                assert m==1 and interval[0]>0
                roots[label,index]=tuple(map(s.Rational,interval))
                for sign in (-1,1):
                    key=label,index,sign
                    rank=multiplicity*block['spin_multiplicity']
                    ranks[key]+=rank;sectors[key].append({**sector,'physical_rank':rank})
        rank=block['dark_energy_one_dimension']*block['spin_multiplicity']
        dark+=rank
        if rank:dark_sectors.append({**sector,'physical_rank':rank})
    records=[]
    for (label,index,sign),rank in ranks.items():
        root=roots[label,index]
        discriminant=sqrt_interval(add((1,1),mul((4,4),root)))
        energy=mul((s.Rational(1,2),)*2,add((1,1),mul((sign,sign),discriminant)))
        records.append({'squared_factor':label,'squared_root_index':index,'energy_sign':sign,
            'squared_root_interval':outward(root,64),'energy_interval':outward(energy,48),
            'physical_rank':rank,'sectors':sectors[label,index,sign]})
    records.append({'exact_energy':'1','energy_interval':['1','1'],
                    'physical_rank':dark,'sectors':dark_sectors})
    records=ordered_records(records,physical['full_maximal_module_physical_dimension'])
    ground,first=records[:2]
    assert ground['squared_factor']=='x**3 - 313*x**2 + 31708*x - 1041468'
    assert ground['squared_root_index']==2 and ground['energy_sign']==-1
    assert ground['physical_rank']==2 and ground['sectors']==[
        {'raw_L':1,'twice_spin':1,'spin_multiplicity':2,'physical_rank':2}]
    assert first['squared_factor']=='x**3 - 313*x**2 + 31738*x - 1041088'
    assert first['squared_root_index']==2 and first['energy_sign']==-1
    assert first['physical_rank']==2 and {r['raw_L'] for r in first['sectors']}=={0,2}
    lower=s.Rational(excluded['excluded_physical_sectors_energy_lower_bound'])
    assert s.Rational(first['energy_interval'][1])<lower
    gap=(s.Rational(first['energy_interval'][0])-s.Rational(ground['energy_interval'][1]),
         s.Rational(first['energy_interval'][1])-s.Rational(ground['energy_interval'][0]))
    assert gap[0]>0
    return {'full_maximal_module_physical_dimension':104,'complete_maximal_module_energy_records':records,
        'dark_energy_one_rank':dark,'excluded_physical_sectors_strict_energy_lower_bound':str(lower),
        'global_physical_ground_energy_interval':ground['energy_interval'],
        'global_physical_ground_rank':2,'global_physical_ground_raw_L':1,
        'global_physical_ground_raw_B3':0,'global_physical_ground_twice_spin':1,
        'first_excited_energy_interval':first['energy_interval'],'first_excited_rank':2,
        'global_physical_gap_interval':outward(gap,48),
        'stationary_rotation_invariant_ground_density':'P_ground/2',
        'reference_selection':'normalized projector on the complete local physical ground doublet',
        'all_zero_boundary_charge_local_Gauss_sectors_compared':True,
        'full_physical_complement_excluded_by_strict_bound':True,
        'numerical_eigensolver_used':False,
        'complete_spectrum_outside_maximal_quark_module_computed':False,
        'unconstrained_or_boundary_compensated_local_ground_proved':False,
        'interacting_many_cell_vacuum_proved':False}
