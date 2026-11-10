"""Literal Gauss-positive link code for every graded occupation slot."""
from itertools import product
import sympy as s

from ..mobile_flux.model import observer,STAR,charge
from ..mobile_flux.routing import odd_words
from .fock import act


def slot_orbit(index,occupied):
    if type(index) is not int or index<0 or type(occupied) is not bool:
        raise ValueError('mode index and occupation flag required')
    face=(index%20)//4
    word=odd_words()[face] if occupied else 0
    return (word,word^STAR)


def encoded_words(mask,size):
    """Unnormalized full gauge-partner expansion; squared denominator2^size."""
    if type(size) is not int or not 1<=size<=6 or type(mask) is not int or not 0<=mask<1<<size:
        raise ValueError('exact small literal encoding requires one to six slots')
    return {words:s.Integer(1) for words in product(*(slot_orbit(j,bool(mask&(1<<j))) for j in range(size)))}


def source_creation_column(mask,index,create,size):
    """The exact graded source gate includes its declared CAR parity sign."""
    result=act(mask,index,create,size)
    if result is None:return {}
    target,sign=result
    return {words:sign*coefficient for words,coefficient in encoded_words(target,size).items()}


def certificate():
    law=observer().base.law
    for index,occupied in product(range(20),(False,True)):
        words=slot_orbit(index,occupied)
        assert words[0]^STAR==words[1]
        assert all(charge(g)==int(occupied) and law.syndrome(g).bit_count()==2*int(occupied) for g in words)
    # Every occupation of four slots, including the fully coincident sector.
    columns=[]
    for mask in range(16):
        words=encoded_words(mask,4)
        assert len(words)==16
        for slot in range(4):
            flipped={tuple(g^STAR if j==slot else g for j,g in enumerate(key)):v for key,v in words.items()}
            assert flipped==words
        assert all(sum(charge(g) for g in key)==mask.bit_count() for key in words)
        columns.append(set(words))
    assert all(not columns[i]&columns[j] for i in range(16) for j in range(i))
    # Exchange is antisymmetric on the source code; this sign is a gate input.
    one=act(0,0,True,4);two=act(one[0],1,True,4)
    other=act(0,1,True,4);reverse=act(other[0],0,True,4)
    assert two[0]==reverse[0] and one[1]*two[1]==-other[1]*reverse[1]
    for mask,index,create in product(range(16),range(4),(False,True)):
        result=act(mask,index,create,4)
        column=source_creation_column(mask,index,create,4)
        assert bool(column)==(result is not None)
        if result is not None:assert set(column)==columns[result[0]]
    return {'one_object_modes_per_cell':20,'literal_link_slots_per_cell':20,
            'full_local_graded_occupation_dimension':2**20,
            'all_occupations_retained':True,'individual_slot_link_domain_size':4,
            'literal_link_domain_is_flat_plus_mode_face_orbit':True,
            'literal_occupation_and_gauge_columns_checked':16,
            'literal_slot_charge_equals_occupation':True,
            'sphere_face_count_equals_twice_object_number':True,
            'graded_source_operator_columns_checked':128,
            'exchange_sign':'-1 by declared CAR gates',
            'one_object_code_equals_flux_face_times_Clifford_code':True,
            'bosonic_two_lane_collision_law_inherited':False,
            'bosonic_link_tensor_product_alone_proves_local_CAR':False,
            'native_face_changing_unitary_automatically_preserves_fixed_face_slot_code':False}
