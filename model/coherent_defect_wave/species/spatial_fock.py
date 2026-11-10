"""Two-cell execution of a declared spatial Fock circuit on common local links.

Exterior one-body hopping differs from partial SWAP of two occupied tensor
registers. No tensor-wave reference or dispersion is imported by this module.
"""
from dataclasses import dataclass
from fractions import Fraction as F

import sympy as s

from .fock import FockObserver, FockState, MODES, coordinates, wedge_insert

A=(1-s.I)/2
B=(1+s.I)/2
STAR=6553634


@dataclass(frozen=True)
class PairFockState:
    amplitudes: tuple
    norm_denominator: int=1

    def __post_init__(self):
        if type(self.norm_denominator) is not int or self.norm_denominator<=0:
            raise ValueError('positive integer normalization required')
        entries=self.amplitudes.items() if isinstance(self.amplitudes,dict) else self.amplitudes
        allowed=set(FockObserver().base.link_words)
        out={}
        for (links,occupied),value in entries:
            links,occupied=tuple(links),tuple(occupied)
            if len(links)!=2 or any(type(g) is not int or g not in allowed for g in links):
                raise ValueError('two owned local link words required')
            if tuple(sorted(set(occupied)))!=occupied or any(type(i) is not int or i not in range(2*MODES) for i in occupied):
                raise ValueError('ordered two-cell occupation set required')
            value=s.sympify(value)
            if not value.is_number or not value.is_finite or value.has(s.Float):
                raise ValueError('exact finite amplitudes required')
            key=links,occupied
            out[key]=s.expand(out.get(key,0)+value)
        object.__setattr__(self,'amplitudes',tuple(sorted((k,v) for k,v in out.items() if v!=0)))

    def as_dict(self):
        return dict(self.amplitudes)

    def norm_squared(self):
        return s.expand(sum((s.conjugate(v)*v for _,v in self.amplitudes),s.Integer(0))/self.norm_denominator)


def frame_sign(word,index):
    p,j,_,_=coordinates(index)
    return -1 if word>(word^STAR) and p==5 and j>=3 else 1


def exchange(state,inverse=False):
    """Exterior power of a dressed one-body root-SWAP, with both links fixed."""
    aa,bb=(s.conjugate(A),s.conjugate(B)) if inverse else (A,B)
    output={}
    for (links,occupied),value in state.amplitudes:
        terms={():value}
        for index in occupied:
            cell,local=divmod(index,MODES); other=1-cell
            sign=frame_sign(links[cell],local)*frame_sign(links[other],local)
            column=((index,aa),(other*MODES+local,bb*sign))
            following={}
            for wedge,amplitude in terms.items():
                for target,weight in column:
                    inserted,parity=wedge_insert(wedge,target)
                    if parity:
                        following[inserted]=s.expand(following.get(inserted,0)+parity*amplitude*weight)
            terms={k:v for k,v in following.items() if v!=0}
        for wedge,amplitude in terms.items():
            key=links,wedge
            output[key]=s.expand(output.get(key,0)+amplitude)
    return PairFockState(output,state.norm_denominator)


def gamma(state,cell):
    if cell not in (0,1):
        raise ValueError('two-cell index required')
    out={}
    for (links,occupied),value in state.amplitudes:
        changed=list(links);changed[cell]^=STAR
        power=sum(index//MODES==cell and coordinates(index%MODES)[0]==5
                  and coordinates(index%MODES)[1]>=3 for index in occupied)
        out[tuple(changed),occupied]=(-1)**power*value
    return PairFockState(out,state.norm_denominator)


def positive_basis(links,occupied):
    """Normalized simultaneous positive-Gauss preparation of a definite word."""
    state=PairFockState({(tuple(links),tuple(occupied)):1})
    for cell in (0,1):
        values=state.as_dict()
        for key,value in gamma(state,cell).amplitudes:
            values[key]=values.get(key,0)+value
        state=PairFockState(values,2*state.norm_denominator)
    return state


def onsite(state,cell,inverse=False):
    """Same local U_F with one shared field register per cell."""
    if cell not in (0,1):
        raise ValueError('two-cell index required')
    observer=FockObserver(); out={}
    for (links,occupied),value in state.amplitudes:
        active=tuple(i%MODES for i in occupied if i//MODES==cell)
        spectators=tuple(i for i in occupied if i//MODES!=cell)
        local=FockState({(links[cell],active):value})
        evolved=observer.inverse_step(local) if inverse else observer.step(local)
        for (word,wedge),weight in evolved.amplitudes:
            words=list(links);words[cell]=word
            joined=tuple(sorted((*spectators,*(cell*MODES+i for i in wedge))))
            key=tuple(words),joined
            out[key]=s.expand(out.get(key,0)+weight)
    return PairFockState(out,state.norm_denominator)


def number_measurement(state,cell):
    if cell not in (0,1):
        raise ValueError('two-cell index required')
    maximum=max((len(occupied) for (_,occupied),_ in state.amplitudes),default=0)
    return tuple((n,PairFockState({key:value for key,value in state.amplitudes
                                if sum(i//MODES==cell for i in key[1])==n},state.norm_denominator))
                 for n in range(maximum+1))


def certificate():
    from .fock import mode
    i=mode(5,3,0,0);j=mode(6,0,0,0)
    state=positive_basis((1,3),(i,MODES+j))
    result=exchange(state)
    assert state.norm_squared()==result.norm_squared()==1
    assert exchange(result,inverse=True)==state
    for cell in (0,1):
        assert gamma(state,cell)==state
        assert gamma(result,cell)==result
        assert exchange(gamma(state,cell))==gamma(result,cell)
    counts=number_measurement(result,0)
    assert [branch.norm_squared() for _,branch in counts]==[s.Rational(1,4),s.Rational(1,2),s.Rational(1,4)]
    assert all(gamma(branch,cell)==branch for _,branch in counts for cell in (0,1))
    assert {words for (words,_),_ in result.amplitudes}=={words for (words,_),_ in state.amplitudes}
    # Pauli blocking and the determinant phase distinguish this gate from tensor partial SWAP.
    identical=PairFockState({((1,1),(j,MODES+j)):1})
    expected=PairFockState({((1,1),(j,MODES+j)):-s.I})
    assert exchange(identical)==expected
    # Vacuum is fixed by the actual spatial and onsite gates, in both constraints.
    vacuum=positive_basis((1,3),())
    assert exchange(vacuum)==vacuum
    assert onsite(onsite(vacuum,0),1)==vacuum
    # One traveling wave then executes the literal owned receiver U.
    one=positive_basis((1,3),(i,))
    hopped=exchange(one)
    evolved=onsite(hopped,1)
    assert evolved.norm_squared()==1
    assert onsite(evolved,1,inverse=True)==hopped
    assert all(gamma(evolved,cell)==evolved for cell in (0,1))
    return {'one_body_gate':'a I + b R_g SWAP R_g on h_left direct_sum h_right',
            'many_body_gate':'Gamma(one_body_gate)',
            'coefficients':{'a':str(A),'b':str(B)},
            'two_wave_input_mode_coordinates':[list(coordinates(i)),list(coordinates(j))],
            'two_wave_output_support':len(result.amplitudes),
            'complete_receiver_number_probabilities':['1/4','1/2','1/4'],
            'receiver_number_outcomes':[0,1,2],
            'each_branch_preserves_both_Gauss_constraints':True,
            'links_retained_during_hopping':True,'total_occupation_conserved':True,
            'exact_inverse':True,'same_mode_pair_exchange_phase':'-I',
            'equals_tensor_wave_partial_SWAP':False,
            'two_cell_positive_Gauss_matter_vacuum_fixed':True,
            'single_wave_hop_then_literal_receiver_U_checked':True,
            'one_wave_receiver_U_output_support':len(evolved.amplitudes),
            'spatially_extended_Fock_circuit':'onsite U_F and six disjoint nearest-neighbor hopping layers',
            'graded_locality':'disjoint even gate algebras commute; disjoint odd CAR fields anticommute',
            'physical_vacuum_selected':False,'relativistic_mass_shell_derived':False,
            'tensor_wave_dispersion_transferred':False}
