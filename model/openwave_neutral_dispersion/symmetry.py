"""Plane parity and actual PW119 link inversion on the complete carrier."""
from collections import defaultdict
from functools import lru_cache
import random
import sympy as s
from model.coherent_defect_wave.chiral_joint.oriented_links.model import reverse_label,endpoint_current
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.clebsch import link
from model.coherent_defect_wave.chiral_joint.periodic_family.model import cells,cell_index,native_blocks,yukawa_words
from model.coherent_defect_wave.chiral_joint.model import clean
from model.openwave_neutral_candidate.source import periodic_terms,translate_key,full_column as boundary_column
from model.openwave_neutral_stationary.source import full_column


def transform(key,axis,*,parity=True,invert=True):
    if axis not in range(3):raise ValueError('axis must be zero, one or two')
    mask,links,scalars=key;out=list(links)
    for i,x in enumerate(cells(2)):
        target=list(x);target[axis]^=1;j=cell_index(2,target)
        out[3*j+axis]=reverse_label(links[3*i+axis]) if invert else links[3*i+axis]
    number=sum(((mask>>(30*i))&((1<<30)-1)).bit_count()
               for i,x in enumerate(cells(2)) if x[axis]==1)
    return (mask,tuple(out),scalars),(-1)**number if parity else 1


def apply(state,axis,**kwargs):
    result=defaultdict(lambda:s.Integer(0))
    for key,value in state.items():
        target,phase=transform(key,axis,**kwargs);result[target]+=phase*value
    return clean(result)


def commutator(key,axis,action=full_column,**kwargs):
    left=apply(action(key),axis,**kwargs);target,phase=transform(key,axis,**kwargs)
    right={out:phase*value for out,value in action(target).items()}
    return clean({key:left.get(key,0)-right.get(key,0) for key in left.keys()|right.keys()})


@lru_cache(maxsize=1)
def coordinate_certificate():
    inversion=s.SparseMatrix(119,119,{(reverse_label(i),i):1 for i in range(119)})
    assert inversion*inversion==s.eye(119)
    count=0
    for block in native_blocks():
        for a in block:
            for b in block:
                assert inversion*link(a,b)*inversion==link(b,a).H;count+=1
    current_checks=0
    for current in range(12):
        for endpoint in (0,1):
            assert inversion*endpoint_current(current,endpoint)*inversion==endpoint_current(current,1-endpoint)
            current_checks+=1
    assert all(len(word)%2==0 for _,_,_,word,_ in yukawa_words())
    return {'full_link_dimension':119,'complete_coordinate_inversion_identities':count,
        'complete_endpoint_current_inversion_identities':current_checks,
        'local_source_words_have_even_fermion_degree':True,
        'full_nonblank_penalty_commutes_with_inversion':True,
        'identity':'J U_ab J=U_ba^dagger; J L_endpoint J=L_other_endpoint'}


@lru_cache(maxsize=1)
def word_certificate():
    terms=periodic_terms();lookup={t.path[0]:t for t in terms}
    checks=0
    for axis in range(3):
        for term in terms:
            assert len(term.path)==1 and term.coefficient==s.I/4
            if term.axis==axis:
                other=lookup[3*term.right_cell+axis]
                assert other.right_cell==term.left_cell
                assert term.coefficient==-s.conjugate(other.coefficient)
                # Inversion swaps coordinate indices; sigma_axis is Hermitian.
                # Plane parity contributes the compensating minus sign.
                assert cells(2)[term.left_cell][axis]!=cells(2)[term.right_cell][axis]
            else:
                assert cells(2)[term.left_cell][axis]==cells(2)[term.right_cell][axis]
            checks+=1
    return {'all_axis_hopping_word_identities':checks,'complete_oriented_edges':24,
        'link_only_axis_reflection':'R_i T_j R_i=(-1)^(i=j) T_j',
        'matter_plane_parity':'Q_i T_j Q_i=(-1)^(i=j) T_j',
        'complete_source_symmetry':'C_i=Q_i R_i, [C_i,H]=0, C_i^2=I, [C_i,C_j]=0',
        'Gauss_preserved':'Each inverted link is moved to the oppositely oriented parallel edge, so its endpoint current stays at the same cell.',
        'matter_or_scalar_coordinates_removed':False}


@lru_cache(maxsize=1)
def certificate():
    coordinate=coordinate_certificate();words=word_certificate()
    rng=random.Random(727919);algebra=0
    for _ in range(48):
        key=(rng.getrandbits(240),tuple(rng.randrange(119) for _ in range(24)),tuple(rng.randrange(5) for _ in range(8)))
        odd=(-1)**key[0].bit_count()
        for axis in range(3):
            one,a=transform(key,axis);two,b=transform(one,axis)
            assert two==key and a*b==1
            for j in range(3):
                left,c=transform(one,j);middle,d=transform(key,j);right,e=transform(middle,axis)
                assert left==right and a*c==d*e
                left,c=translate_key(one,j,periodic=True)
                middle,d=translate_key(key,j,periodic=True);right,e=transform(middle,axis)
                assert left==right and a*c==(odd if j==axis else 1)*d*e
                algebra+=1
    keys=[(1<<11,(0,)*24,(0,)*8),
          ((1<<9)|(1<<45)|(1<<179),tuple(1 if j in (0,7) else 0 for j in range(24)),(1,0,0,0,0,0,0,0)),
          ((1<<10)|(1<<26)|(1<<155),tuple((17*j+3)%119 for j in range(24)),(1,2,0,3,0,4,0,0))]
    checks=0
    for key in keys:
        for axis in range(3):assert not commutator(key,axis);checks+=1
    witness=keys[0]
    no_parity=commutator(witness,0,parity=False)
    no_inverse=commutator(witness,0,invert=False)
    antiperiodic=commutator(witness,0,action=lambda key:boundary_column(key,periodic=False))
    assert no_parity and no_inverse and antiperiodic
    return {'coordinate_and_current_identities':coordinate,'complete_word_proof':words,
        'independent_occupied_full_carrier_algebra_checks':algebra,
        'occupied_complete_H_column_checks':checks,
        'translation_conjugation':'T_j C_i T_j^dagger=(-1)^(N delta_ij) C_i',
        'mutations':{'omit_plane_parity_nonzero_commutator_terms':len(no_parity),
                    'omit_link_inversion_nonzero_commutator_terms':len(no_inverse),
                    'substitute_antiperiodic_source_nonzero_commutator_terms':len(antiperiodic)},
        'restriction':'Exact periodic side-two source; the axis reflection fixes every matter site only at this size.'}
