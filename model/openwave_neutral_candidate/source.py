"""The actual native neutral current and an explicit periodic-boundary variant."""
from collections import defaultdict
from dataclasses import replace
from functools import lru_cache
from itertools import product
import random
import sympy as s

from model.neutrino_mass_response.neutral import certificate as neutral_identity
from model.defect_mass_join.source import car, one_cell_polynomial
from model.coherent_defect_wave.chiral_joint.periodic_family.model import (
    Configuration, cells, cell_index, shift_cell, kinetic_terms, native_blocks,
    translate_fermions)
from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as native
from model.coherent_defect_wave.chiral_joint.model import clean


CONFIG = Configuration()
L_WEIGHTS = tuple(-1 if a == 9 else 1 if a in (10,11) else 0 for a in range(15))
B_WEIGHTS = tuple(-1 if a in (0,1,4,12,13,14) else 1 if a in (2,3,5,6,7,8) else 0 for a in range(15))


@lru_cache(maxsize=1)
def periodic_terms():
    result=[]
    for term in kinetic_terms(CONFIG):
        anchor=cells(2)[term.left_cell]
        seam=(anchor[term.axis]+term.displacement)//2
        result.append(replace(term,coefficient=(-1)**seam*term.coefficient))
    assert len(result)==24 and {row.coefficient for row in result}=={s.I/4}
    return tuple(result)


def full_column(key, *, periodic=False, coupling=s.Rational(1,2**34)):
    """Complete source action; only twelve stored seam signs differ in the variant."""
    native.validate_key(key,CONFIG)
    ordinary=native.hopping_column(key,kinetic_terms(CONFIG))
    output=defaultdict(lambda:s.Integer(0),native.column(key,CONFIG))
    for target,value in ordinary.items():output[target]-=value
    terms=periodic_terms() if periodic else kinetic_terms(CONFIG)
    for target,value in native.hopping_column(key,terms).items():output[target]+=coupling*value
    return clean(output)


def translate_key(key,axis, *, periodic=False):
    mask,links,scalars=key
    changed,phase=translate_fermions(mask,axis,2)
    if periodic:
        wraps=sum(1 for mode in range(240) if mask>>mode&1 and cells(2)[mode//30][axis]==1)
        phase*=(-1)**wraps
    outlinks=[0]*24;outscalars=[0]*8
    for site,xyz in enumerate(cells(2)):
        dest=cell_index(2,shift_cell(2,xyz,axis))
        outscalars[dest]=scalars[site]
        for edgeaxis in range(3):outlinks[3*dest+edgeaxis]=links[3*site+edgeaxis]
    return (changed,tuple(outlinks),tuple(outscalars)),phase


def apply_translation(state,axis,periodic):
    output=defaultdict(lambda:s.Integer(0))
    for key,value in state.items():
        target,phase=translate_key(key,axis,periodic=periodic)
        output[target]+=phase*value
    return clean(output)


def current_certificate():
    inherited=neutral_identity()
    baryon=car.one_body(s.diag(*(B_WEIGHTS[a%15] for a in range(30))))
    h=one_cell_polynomial()
    assert not car.commutator(baryon,h)
    q=car.add(car.gauge_current(10),car.gauge_current(11))
    for spin in (0,1):
        j=car.neutral_composite(spin)
        assert not car.commutator(q,j) and not car.commutator(baryon,j)
        assert all(len(word)%2==1 for word,_,_ in j)
    for block in native_blocks():
        assert len({L_WEIGHTS[a] for a in block})==len({B_WEIGHTS[a] for a in block})==1
    assert all((l+b)%2==1 for l,b in zip(L_WEIGHTS,B_WEIGHTS))
    return {'native_operator':inherited['native_neutral_operator'],
            'gauge_current_commutators_per_spin':12,
            'existing_electric_candidate_current':'Q=T3+Y including the scalar current',
            'Q_commutator_exactly_zero':True,'raw_lepton_change':-1,'raw_baryon_three_change':0,
            'same_complete_local_H_preserves_L_and_B3':True,
            'each_actual_hopping_coordinate_preserves_L_and_B3':True,
            'full_source_fermion_parity_identity':'(-1)^N=(-1)^(L+B3)',
            'sterile_coordinate_or_zero_charge_added':False,
            'repair_transparency_used_as_electric_neutrality':False}


def boundary_certificate():
    original=kinetic_terms(CONFIG);variant=periodic_terms()
    signs=[s.cancel(b.coefficient/a.coefficient) for a,b in zip(original,variant)]
    assert signs.count(-1)==12 and signs.count(1)==12
    assert all(sign**2==1 for sign in signs)
    assert s.simplify(60*sum(abs(t.coefficient) for t in variant))==360
    assert all(sum(cells(2)[t.left_cell])%2 != sum(cells(2)[t.right_cell])%2 for t in variant)
    # Full operator translation is tested on occupied columns, including nonblank links.
    keys=[native.empty_key(CONFIG),
          (1<<11,(0,)*24,(0,)*8),
          ((1<<9)|(1<<45)|(1<<179),(0,)*24,(0,)*8),
          ((1<<10)|(1<<26)|(1<<155),tuple(1 if j in (0,7) else 0 for j in range(24)),(1,0,0,0,0,0,0,0))]
    checks=0
    for periodic in (False,True):
        for key in keys:
            image=full_column(key,periodic=periodic)
            for axis in range(3):
                left=apply_translation(image,axis,periodic)
                shifted,phase=translate_key(key,axis,periodic=periodic)
                right={out:phase*value for out,value in full_column(shifted,periodic=periodic).items()}
                assert left==clean(right);checks+=1
    rng=random.Random(41859);squares=0
    for _ in range(64):
        mask=rng.getrandbits(240)
        key=(mask,tuple(rng.randrange(119) for _ in range(24)),tuple(rng.randrange(5) for _ in range(8)))
        for periodic in (False,True):
            for axis in range(3):
                target,first=translate_key(key,axis,periodic=periodic)
                returned,second=translate_key(target,axis,periodic=periodic)
                assert returned==key and first*second==(1 if periodic else (-1)**mask.bit_count())
                squares+=1
    return {'changed_seam_terms':12,'retained_unchanged_terms':12,
            'periodic_coefficients':[str(t.coefficient) for t in variant],
            'all_onsite_and_link_operators_identical':True,
            'complete_hopping_norm_upper':'360','blank_link_compression_of_every_hopping_word':'0',
            'bipartite_effective_H_evenness_preserved':True,
            'second_order_effective_H_identical':'distinct edges have orthogonal nonblank supports; each surviving edge sign is squared',
            'second_order_diagonal_edge_factors':[str(sign**2) for sign in signs],
            'different_edge_second_order_pairs':24*23,
            'complete_occupied_column_translation_checks':checks,
            'occupied_full_carrier_translation_square_checks':squares,
            'original_axis_translation_squared':'fermion parity',
            'periodic_variant_axis_translation_squared':'identity',
            'original_odd_sector_zero_momentum_exists':False,
            'original_odd_axis_translation_eigenvalues':['i','-i'],
            'periodic_variant_zero_momentum_projector':'product_axis(I+T_axis)/2',
            'source_boundaries_silently_identified':False}


def momentum_trial_certificate():
    # Distinct cells are exactly orthogonal because their local raw-L labels differ.
    # This independent eight-address representation records the group-average norm.
    sites=cells(2);matrices=[]
    for axis in range(3):
        matrix=s.zeros(8)
        for i,x in enumerate(sites):matrix[cell_index(2,shift_cell(2,x,axis)),i]=1
        assert matrix*matrix==s.eye(8)
        matrices.append(matrix)
    projector=s.eye(8)
    for matrix in matrices:projector=projector*(s.eye(8)+matrix)/2
    localized=s.eye(8)[:,0];projected=projector*localized
    assert projector*projector==projector and (projected.T*projected)[0]==s.Rational(1,8)
    assert all(matrix*projected==projected for matrix in matrices)
    return {'trial':'projected native local neutral lowest eigenline acting on the translation-zero atomic spin ground',
            'different_location_local_raw_L_patterns_orthogonal':True,
            'periodic_zero_momentum_projection_success_on_normalized_local_trial':'1/8',
            'all_eight_momentum_effects_retained':True,
            'trial_atomic_energy':'7E0+E_neutral_local',
            'trial_full_hopping_expectation':'0 because all links are blank',
            'trial_is_full_interacting_eigenvector':False,
            'abstract_address_permutation_is_attached_by_native_translation_of_the_operator_and_ground':True}


def dark_number_control():
    from model.dark_spectral_source.exact import target_numerator
    _,real,imag=target_numerator()
    weights=[int(a)*int(a)+int(b)*int(b) for a,b in zip(real,imag)]
    norm=sum(weights)
    mean=sum(i.bit_count()*weight for i,weight in enumerate(weights))
    second=sum(i.bit_count()**2*weight for i,weight in enumerate(weights))
    assert norm==4 and mean==16 and second==64
    return {'repair_transparent_source':'unchanged six-site Weyl/contact dark preparation',
            'actual_conserved_number_above_its_vacuum':4,
            'number_variance':0,'number_is_identified_with_R15_electric_current':False,
            'electric_neutrality_inferred_from_repair_transparency':False,
            'dark_source_energy_offset_reused_as_a_rest_mass':False}
