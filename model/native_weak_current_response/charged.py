"""Actual charged atomic target insertion of the original weak density."""
from pathlib import Path
from functools import lru_cache
import json
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from model.native_weak_current_interface.model import dressed_weak_polynomial
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.model import target,weight_basis
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import polynomial_map,validate_annihilator,reduced,source
from model.coherent_defect_wave.chiral_joint.local_reference import physical
from model.coherent_defect_wave.chiral_joint.local_reference.reduction import quark_module,lepton_operators,coordinate_kernel
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.local import positive_response
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import evaluate,outward,mul
from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
from .local import EIGENVALUES,projectors

@lru_cache(None)
def complete_target(magnetic):
    if magnetic in (0,2):return target('ec',magnetic)
    assert magnetic==-2
    keys=weight_basis(4,magnetic,0,-2);raised=weight_basis(4,magnetic,2,-2)
    q=quark_module();l=lepton_operators();scalar=higgs_memory()[1]
    weak=physical.tensor_map(keys,raised,quark=q['currents'][0]+s.I*q['currents'][1],
      lepton=l['currents'][0]+s.I*l['currents'][1],scalar=scalar[8]+s.I*scalar[9])
    basis,_=coordinate_kernel(weak);assert basis.cols==len(keys)-len(raised)
    return reduced(keys,basis)

@lru_cache(None)
def calculate():
    W,D,S,Ka,NL=dressed_weak_polynomial();local=positive_response('ec');q=local['q'];x=q.gen
    selected=s.CRootOf(q.as_expr(),0);field=s.QQ.algebraic_field(selected)
    minimal=s.Poly(field.ext.minpoly.as_expr().subs(field.ext.minpoly.gen,x),x)
    assert minimal.degree()==4 and q.rem(minimal).is_zero
    interval=minimal.intervals(eps=s.Rational(1,10**100))[0][0]
    broad=local['rows'][0]['energy_interval'];assert broad[0]<=interval[0]<interval[1]<=broad[1]
    def at(poly):
        value=field.zero
        for c in poly.all_coeffs():value=value*field.unit+field.convert(c)
        return value
    def polynomial(value):return s.Poly.from_list([s.Rational(c.numerator,c.denominator) for c in value.to_list()],x)
    coefficients=[]
    for i in range(q.degree()):
        c=s.Poly(sum(q.nth(j)*x**(j-i-1) for j in range(i+1,q.degree()+1)),x)
        coefficients.append(at(c)/at(q.diff()))
    controls=[];moments=[]
    for magnetic in (-2,0,2):
        block=complete_target(magnetic);H=block['H'];G=block['gram'];matrices={}
        for name,operator in {'W':W,'W2':car.multiply(W,W)}.items():
            raw=polynomial_map(block['keys'],block['keys'],operator)
            image=raw*block['basis']
            reduced=block['left']*image.extract(block['coordinates'],range(block['basis'].cols))
            assert block['basis']*reduced==image
            assert G*reduced==reduced.H*G
            matrices[name]=s.ImmutableSparseMatrix(reduced)
        assert matrices['W2']==matrices['W']**2
        matrices.update({'probability_'+str(value):P for value,P in zip(EIGENVALUES,projectors(matrices['W']),strict=True)})
        validate_annihilator(H,q)
        powers=[s.eye(H.rows)]
        for _ in range(q.degree()-1):powers.append(H*powers[-1])
        def trace(O):return sum((c*field.convert(s.trace(power*O)) for c,power in zip(coefficients,powers,strict=True)),field.zero)
        assert trace(s.eye(H.rows))==field.one
        row={name:trace(O) for name,O in matrices.items()};moments.append(row)
        controls.append({'twice_magnetic_spin':magnetic,'complete_charged_target_dimension':H.rows,
          'selected_energy_projector_rank':1,'full_polynomial_images_and_metric_Hermiticity_checked':True})
    assert moments[0]==moments[1]==moments[2]
    for magnetic in (-1,1):
        sb=source(magnetic)
        for spin in (0,1):
            tb=complete_target(magnetic-(1-2*spin))
            raw=polynomial_map(sb['keys'],tb['keys'],car.monomial(((9+15*spin,False),)))
            image=raw*sb['basis']
            mapped=tb['left']*image.extract(tb['coordinates'],range(sb['basis'].cols))
            assert tb['basis']*mapped==image
    row=moments[0];row['variance']=row['W2']-row['W']**2
    assert sum((row['probability_'+str(value)] for value in EIGENVALUES),field.zero)==field.one
    values={name:evaluate(polynomial(value),interval) for name,value in row.items()}
    atomic_weight=local['rows'][0]['weight_interval']
    return {'complete_target_controls':controls,'selected_charged_energy_minimal_polynomial':str(minimal.as_expr()),
      'selected_charged_energy_interval':outward(interval,100),
      'charged_target_ground_current':{name:{'coefficients_in_charged_energy':[str(c) for c in polynomial(value).all_coeffs()],
          'outward_interval':outward(values[name],100)} for name,value in row.items()},
      'rank_one_compression_identity':'P_e W P_e=tr(P_e W)P_e in all three complete magnetic target blocks',
      'both_initial_ground_magnetic_states_and_both_native_removal_maps_checked':4,
      'atomic_pumped_bin_norm_squared_interval':outward(atomic_weight,100),
      'atomic_unnormalized_charged_current_insertion_interval':outward(mul(atomic_weight,values['W']),100),
      'atomic_unnormalized_charged_current_second_moment_interval':outward(mul(atomic_weight,values['W2']),100),
      'all_source_local_reference_and_target_complements_retained':True}

if __name__=='__main__':
    value=calculate()
    print(json.dumps({k:v['outward_interval'] for k,v in value['charged_target_ground_current'].items()},indent=2))
