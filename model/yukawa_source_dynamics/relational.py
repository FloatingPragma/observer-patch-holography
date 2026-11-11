"""Relational phase invariants with every quantum record retained.

No scalar branch phase or analyzer choice is supplied here. The cyclic
routing and the default weak-feedback word are consumed as specified by
their owning source. The weak word contains an inherited analyzer projector
and the inherited rational clock tangent1/64.
"""
from functools import lru_cache
import sympy as s
from model.standard_model_source_map import koide_source
from .complete_c3 import exact as cyclic, clean


@lru_cache(maxsize=1)
def exact():
    base=cyclic(); omega=-s.Rational(1,2)+s.I*s.sqrt(3)/2
    # Columns are normalized source output record vectors conditioned on the
    # three actual C3 characters. All128 fine records occur as rows.
    E=s.Matrix([[s.sqrt(w)*omega**(k*exponent) for k in range(3)]
                for raw,exponent,w,rotation in base['records']])
    gram=clean(E.H*E)
    wanted=s.eye(3)*s.Rational(5,24)+s.ones(3)*s.Rational(19,24)
    if gram!=wanted or any(gram[i,i]!=1 for i in range(3)):
        raise ArithmeticError('full quantum-record Gram differs from native weights')
    bargmann=s.simplify(gram[0,1]*gram[1,2]*gram[2,0])
    if bargmann!=s.Rational(19,24)**3:
        raise ArithmeticError('a retained relative cyclic phase was lost')
    # Evaluate the actual full weak record word. Its canonical basis is fixed
    # by X and the native C3 characters, not by a measured fermion spectrum.
    F=s.Matrix([[omega**(-i*j)/s.sqrt(3) for j in range(3)] for i in range(3)])
    Xframe=s.Matrix([[1,1],[1,-1]])/s.sqrt(2)
    basis=s.kronecker_product(Xframe,F)
    K=koide_source.weak_kraus()
    eigenrecords=[]
    for branch in K:
        diagonal=clean(basis.H*branch*basis)
        if not diagonal.is_diagonal():
            raise ArithmeticError('the exact source-character frame does not diagonalize the weak word')
        eigenrecords.append(list(diagonal.diagonal()))
    Eweak=s.Matrix(eigenrecords)
    weakgram=clean(Eweak.H*Eweak)
    if any(weakgram[i,i]!=1 for i in range(6)):
        raise ArithmeticError('weak relational records lost a quantum outcome')
    triples=[]
    for sign in range(2):
        ids=tuple(3*sign+i for i in range(3))
        invariant=s.simplify(weakgram[ids[0],ids[1]]*weakgram[ids[1],ids[2]]*weakgram[ids[2],ids[0]])
        if s.simplify(s.im(invariant))!=0:
            raise ArithmeticError('the inherited weak source has a nonzero three-family relative phase')
        triples.append(invariant)
    return dict(record_vectors=E,gram=gram,bargmann=bargmann,weak_records=Eweak,
                weak_gram=weakgram,weak_triples=tuple(triples))


def report():
    x=exact(); entries=lambda a:[[str(v) for v in row] for row in a.tolist()]
    return dict(retained_quantum_record_rows=128,conditional_family_characters=3,
                record_gram_exact=entries(x['gram']),
                oriented_three_record_invariant_exact=str(x['bargmann']),
                oriented_three_record_phase_exact='0',
                quantum_records_traced_or_postselected=False,
                record_Gram_spectrum_exact={'31/12':1,'5/24':2},
                uniform_family_input_Schmidt_probabilities_exact=['31/36','5/72','5/72'],
                common_record_unitaries_preserve_Gram_and_triple_phase=True,
                canonical_record_normalization='E*(E^dag E)^(-1/2) has Gram I3',
                weak_word_source='koide_source.weak_kraus(z=1,tangent=1/64)',
                weak_word_all_quantum_records=8,weak_word_input_modes=6,
                weak_word_record_gram_exact=entries(x['weak_gram']),
                weak_word_family_oriented_invariants_exact=list(map(str,x['weak_triples'])),
                weak_word_family_oriented_invariants_imaginary_part_exact=['0','0'],
                inherited_default_word_contains_analyzer_projector=True,
                scalar_phase_inserted=False,empirical_inputs=[],
                arbitrary_coherent_readouts_excluded=False,
                specified_source_relational_records_derive_three_nonzero_distinct_Yukawas=False)
