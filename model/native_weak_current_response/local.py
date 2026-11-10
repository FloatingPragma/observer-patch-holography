"""Exact Higgs-dressed weak-current moments on the complete local ground."""
from pathlib import Path
from functools import lru_cache
import json
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from model.native_weak_current_interface.model import dressed_weak_polynomial,dressed_weak_matrix
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source,polynomial_map
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import field_data,field_value,polynomial,encoded
from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import evaluate,outward

EIGENVALUES=(s.Integer(-1),s.Rational(-1,2),s.Integer(0),s.Rational(3,2),s.Integer(2))

def projectors(W):
    identity=s.eye(W.rows);out=[]
    annihilator=identity
    for value in EIGENVALUES:annihilator=annihilator*(W-value*identity)
    assert annihilator==s.zeros(W.rows)
    for value in EIGENVALUES:
        P=identity
        for other in EIGENVALUES:
            if other!=value:P=P*(W-other*identity)/(value-other)
        assert P*P==P
        out.append(s.ImmutableSparseMatrix(P))
    assert sum(out,s.zeros(W.rows))==identity
    for i,P in enumerate(out):
        for j,Q in enumerate(out):assert P*Q==(P if i==j else s.zeros(W.rows))
    return out

@lru_cache(None)
def calculate():
    W,D,S,Ka,NL=dressed_weak_polynomial()
    scalar_density=car.scale(car.multiply(car.monomial(matrix=S),NL),s.Rational(1,2))
    assert D==car.add(scalar_density,W)
    operators={'W':W,'W2':car.multiply(W,W),'D':D,'scalar_density':scalar_density}
    e,f,field=field_data();eroot=f.intervals(eps=s.Rational(1,10**100))[0][0]
    all_moments=[];source_controls=[]
    for magnetic in (-1,1):
        block=source(magnetic);H=block['H'];G=block['gram'];matrices={}
        for name,operator in operators.items():
            raw=polynomial_map(block['keys'],block['keys'],operator)
            image=raw*block['basis']
            reduced=block['left']*image.extract(block['coordinates'],range(block['basis'].cols))
            assert block['basis']*reduced==image
            assert G*reduced==reduced.H*G
            matrices[name]=s.ImmutableSparseMatrix(reduced)
        assert matrices['W2']==matrices['W']**2
        assert matrices['D']==matrices['W']+matrices['scalar_density']
        matrices.update({'probability_'+str(value):P for value,P in zip(EIGENVALUES,projectors(matrices['W']),strict=True)})
        x=s.Symbol('x');p=s.Poly(H.charpoly(x).as_expr(),x).sqf_part()
        powers=[s.eye(H.rows)]
        for _ in range(p.degree()-1):powers.append(H*powers[-1])
        coefficients=[]
        for i in range(p.degree()):
            c=s.Poly(sum(p.nth(j)*x**(j-i-1) for j in range(i+1,p.degree()+1)),x)
            coefficients.append(field_value(c)/field_value(p.diff()))
        def expectation(O):
            return sum((c*field.convert(s.trace(power*O)) for c,power in zip(coefficients,powers,strict=True)),field.zero)
        assert expectation(s.eye(H.rows))==field.one
        moments={name:expectation(O) for name,O in matrices.items()}
        assert sum((moments['probability_'+str(value)] for value in EIGENVALUES),field.zero)==field.one
        assert sum((field.convert(value)*moments['probability_'+str(value)] for value in EIGENVALUES),field.zero)==moments['W']
        assert sum((field.convert(value**2)*moments['probability_'+str(value)] for value in EIGENVALUES),field.zero)==moments['W2']
        all_moments.append(moments)
        source_controls.append({'twice_magnetic_spin':magnetic,'complete_source_dimension':H.rows,
          'full_polynomial_image_inclusion_checks':len(operators),'metric_Hermiticity_checked':True,
          'source_squarefree_polynomial':str(p.as_expr()),'ground_projector_rank':1})
    assert all_moments[0]==all_moments[1]
    moments=all_moments[0];moments['variance']=moments['W2']-moments['W']**2
    full=dressed_weak_matrix();Ps=projectors(full)
    assert full.H==full and full.shape==(80,80)
    ranks=[int(s.trace(P)) for P in Ps]
    assert sum(ranks)==80 and all(rank>0 for rank in ranks)
    result={'source_controls':source_controls,'full_active_dimension':80,
      'complete_full_eigenvalues':list(map(str,EIGENVALUES)),'complete_full_projector_ranks':ranks,
      'actual_local_ground_moments':{name:{'coefficients_in_ground_energy':encoded(value),
        'outward_interval':outward(evaluate(polynomial(value),eroot),100)} for name,value in moments.items()},
      'both_local_ground_magnetic_components_agree_exactly':True,
      'full_source_complement_discarded':False}
    return result

if __name__=='__main__':
    value=calculate()
    print(json.dumps({k:v['outward_interval'] for k,v in value['actual_local_ground_moments'].items()},indent=2))
