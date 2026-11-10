"""Independent exact contractions of the conjugated full weak density."""
from functools import lru_cache
import sympy as s
from model.native_weak_current_interface.model import dressed_weak_polynomial
from model.native_weak_current_response.charged import complete_target
from model.native_weak_current_response.local import EIGENVALUES, projectors
from model.coherent_defect_wave.chiral_joint.neutral_local_spectrum.symmetries import conjugate_polynomial
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import polynomial_map, validate_annihilator
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.local import positive_response
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import evaluate, outward, mul
from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.coherent_defect_wave.chiral_particle_map.model import pauli


@lru_cache(None)
def polynomial():
    W, _, _, _, _ = dressed_weak_polynomial()
    conjugated = conjugate_polynomial(W)
    phi = higgs_memory()[0]; sigma = pauli(); terms = []
    for a in range(3):
        scalar = sum((s.conjugate(sigma[a][j,i])*phi[i].H*phi[j]
                      for i in range(2) for j in range(2)), s.zeros(5))
        one = s.zeros(30)
        for spin in range(2):
            for i in range(2):
                for j in range(2):
                    one[15*spin+10+i,15*spin+10+j] = sigma[a][i,j]/2
        terms.append(car.multiply(car.monomial(matrix=scalar), car.one_body(one)))
    formula = car.scale(car.add(*terms), -1)
    assert conjugated == formula and conjugated != W
    assert car.adjoint(conjugated) == conjugated
    assert conjugate_polynomial(conjugated) == W
    for a in range(12): assert not car.commutator(car.gauge_current(a), conjugated)
    for a in range(3): assert not car.commutator(car.spin_current(a), conjugated)
    return conjugated


@lru_cache(None)
def calculate():
    Wbar = polynomial(); atom = positive_response('ec'); q = atom['q']; x = q.gen
    selected = s.CRootOf(q.as_expr(), 0); field = s.QQ.algebraic_field(selected)
    minimal = s.Poly(field.ext.minpoly.as_expr().subs(field.ext.minpoly.gen, x), x)
    assert minimal.degree() == 4 and q.rem(minimal).is_zero
    interval = minimal.intervals(eps=s.Rational(1,10**100))[0][0]
    broad = atom['rows'][0]['energy_interval']
    assert broad[0] <= interval[0] < interval[1] <= broad[1]

    def at(poly):
        value = field.zero
        for c in poly.all_coeffs(): value = value*field.unit+field.convert(c)
        return value

    def as_poly(value):
        return s.Poly.from_list([s.Rational(c.numerator,c.denominator) for c in value.to_list()],x)

    # Full q(H)=0 Riesz projector, not a selected numerical eigenvector.
    coefficients = [at(s.Poly(sum(q.nth(j)*x**(j-i-1)
                         for j in range(i+1,q.degree()+1)),x))/at(q.diff())
                    for i in range(q.degree())]
    rows=[]; controls=[]
    for magnetic in (-2,0,2):
        block = complete_target(magnetic); H=block['H']; G=block['gram']; matrices={}
        for name,operator in {'W':Wbar,'W2':car.multiply(Wbar,Wbar)}.items():
            raw=polynomial_map(block['keys'],block['keys'],operator)
            image=raw*block['basis']
            reduced=block['left']*image.extract(block['coordinates'],range(block['basis'].cols))
            assert block['basis']*reduced == image
            assert G*reduced == reduced.H*G
            matrices[name]=s.ImmutableSparseMatrix(reduced)
        assert matrices['W2'] == matrices['W']**2
        matrices.update({'probability_'+str(value):P for value,P in
                         zip(EIGENVALUES,projectors(matrices['W']),strict=True)})
        validate_annihilator(H,q); powers=[s.eye(H.rows)]
        for _ in range(q.degree()-1): powers.append(H*powers[-1])

        def trace(operator):
            return sum((c*field.convert(s.trace(power*operator))
                        for c,power in zip(coefficients,powers,strict=True)),field.zero)

        assert trace(s.eye(H.rows)) == field.one
        rows.append({name:trace(O) for name,O in matrices.items()})
        controls.append({'twice_magnetic_spin':magnetic,
                         'complete_target_dimension':H.rows,'spectral_projector_rank':1,
                         'complete_polynomial_image_and_metric_Hermiticity_checked':True})
    assert rows[0] == rows[1] == rows[2]
    row=rows[0]; row['variance']=row['W2']-row['W']**2
    assert sum((row['probability_'+str(value)] for value in EIGENVALUES),field.zero)==field.one
    assert row['W2'] == -row['W']/field.convert(2)
    values={name:evaluate(as_poly(value),interval) for name,value in row.items()}
    weight=atom['rows'][0]['weight_interval']
    return {'operator':'Wbar=Theta W Theta^-1, distinct from W',
        'full_native_Gauss_commutators_checked':12,'full_native_spin_commutators_checked':3,
        'complete_target_controls':controls,
        'selected_charged_energy_minimal_polynomial':str(minimal.as_expr()),
        'selected_charged_energy_interval':outward(interval,100),
        'atomic_addition_branch_W':{name:{
            'coefficients_in_charged_energy':[str(c) for c in as_poly(value).all_coeffs()],
            'outward_interval':outward(values[name],100)} for name,value in row.items()},
        'atomic_pumped_bin_norm_squared_interval':outward(weight,100),
        'atomic_unnormalized_W_insertion_interval':outward(mul(weight,values['W']),100),
        'atomic_unnormalized_W2_insertion_interval':outward(mul(weight,values['W2']),100),
        'same_values_for_every_spin_magnetic_block':True,
        'opposite_charge_means_equal_for_same_untransformed_W':False}
