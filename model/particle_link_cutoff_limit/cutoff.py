"""Exact finite-volume cutoff bounds with an unchanged native rotor control."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s

from model.particle_observer_outputs import field_route
from model.particle_ground_spectral_records import records

SOURCE_INPUTS=(
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_ground_spectral_records/records.py',
    'model/particle_spatial_ground_state/PROOF.md',
    'Lean/QFT/ElectricCutoffFloor.lean',
)


def rational(value):
    if type(value) not in (int,str,F):
        raise ValueError('exact rational input required')
    return F(value)


def block_floor(retained_lower, omitted_lower, coupling_norm):
    """H_PP>=L, H_QQ>=b>L and ||H_PQ||<=c imply H>=L-c²/(b-L)."""
    lo,b,c=map(rational,(retained_lower,omitted_lower,coupling_norm))
    if b<=lo or c<0:
        raise ValueError('omitted floor must exceed retained floor; norm is nonnegative')
    return lo-c*c/(b-lo)


@lru_cache(maxsize=8)
def rotor_blocks(cutoff):
    """Full native CAR and every flux -K..K, in its real phase chart."""
    if type(cutoff) is not int or cutoff<6:
        raise ValueError('native flux cutoff must be an integer at least six')
    c=field_route.rotor_coefficients(); n=2*cutoff+1
    phase=s.diag(*[s.I**(mask>>2).bit_count() for mask in range(16)])
    onsite=phase.H*c['onsite']*phase
    hop=s.simplify(phase.H*c['hop']*phase)
    entries={}
    def add(i,j,value): entries[i,j]=entries.get((i,j),F(0))+F(value)
    for mask in range(16):
        for j in range(n): add(mask*n+j,mask*n+j,F((j-cutoff)**2,2))
    for (i,j),v in onsite.todok().items():
        for k in range(n): add(i*n+k,j*n+k,v)
    charge=c['charge']
    for (i,j),v in hop.todok().items():
        for k in range(n-charge):
            add(i*n+k+charge,j*n+k,v)
            add(j*n+k,i*n+k+charge,v)
    entries={key:v for key,v in entries.items() if v}
    neighbors=[set() for _ in range(16*n)]
    for i,j in entries:
        if entries[i,j]!=entries.get((j,i),0):raise ArithmeticError('lost native Hermiticity')
        if i!=j:neighbors[i].add(j)
    unseen=set(range(16*n)); blocks=[]
    while unseen:
        todo=[min(unseen)]; component=set()
        while todo:
            i=todo.pop()
            if i not in component:component.add(i);todo.extend(neighbors[i]-component)
        indices=tuple(sorted(component));unseen-=component
        occupations={(i//n).bit_count() for i in indices}
        if len(occupations)!=1:raise ArithmeticError('occupation is not conserved')
        matrix=tuple(tuple(entries.get((i,j),F(0)) for j in indices) for i in indices)
        blocks.append((indices,occupations.pop(),matrix))
    return tuple(blocks)


@lru_cache(maxsize=16)
def sector_intervals(cutoff=6,bits=24):
    if type(bits) is not int or not 4<=bits<=64:raise ValueError('bits outside4..64')
    result={}
    for number in range(5):
        matrices=tuple(set(m for _,q,m in rotor_blocks(cutoff) if q==number))
        lo,hi=F(-4),F(0)
        if not all(records._shift_is_positive(m,lo) for m in matrices):
            raise ArithmeticError('native lower bracket fails')
        if all(records._shift_is_positive(m,hi) for m in matrices):
            raise ArithmeticError('native upper bracket fails')
        while hi-lo>F(1,1<<bits):
            mid=(lo+hi)/2
            if all(records._shift_is_positive(m,mid) for m in matrices):lo=mid
            else:hi=mid
        result[number]=(lo,hi)
    return result


def rotor_limit_control(cutoff=6,bits=24):
    retained=sector_intervals(cutoff,bits)
    omitted=F((cutoff+1)**2,2)-4
    parent={q:(block_floor(lo,omitted,4),hi) for q,(lo,hi) in retained.items()}
    ground_lo=min(v[0] for v in parent.values());ground_hi=min(v[1] for v in parent.values())
    gaps={q:(parent[q][0]-ground_hi,parent[q][1]-ground_lo) for q in (1,3)}
    return {'cutoff':cutoff,'retained_dimension':16*(2*cutoff+1),
            'all_occupations_retained':True,'retained_sector_intervals':retained,
            'omitted_electric_energy_lower':F((cutoff+1)**2,2),
            'omitted_full_action_lower':omitted,'off_diagonal_norm_upper':F(4),
            'uncut_sector_intervals':parent,'uncut_ground_interval':(ground_lo,ground_hi),
            'uncut_adjacent_charge_threshold_intervals':gaps,
            'source_time_rescaled':False,'spatial_particle_mass_inferred':False}


def serial(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serial(v) for v in value]
    return value


def certificate():
    return serial({'schema':'oph.finite_volume_link_cutoff_limit.v1',
        'source_inputs':list(SOURCE_INPUTS),
        'native_rotor_controls':[rotor_limit_control(k) for k in (6,12,18)],
        'fixed_volume_Wilson_word_form_core':True,
        'fixed_volume_ground_energy_convergence':True,
        'fixed_volume_ground_state_trace_norm_subsequence':True,
        'same_clock_complete_bounded_history_subsequence':True,
        'joint_volume_spacing_cutoff_limit_proved':False,
        'physical_particle_or_Einstein_identification':False,
        'empirical_parameters_used':[]})
