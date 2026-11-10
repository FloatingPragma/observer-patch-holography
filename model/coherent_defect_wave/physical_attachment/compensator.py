"""A finite charged-memory controller compensating the Yukawa determinant charge.

The small executed controls isolate the total-degree conservation law. The
full spin/internal tensor uses the invariant lines in chiral.py.
"""
import sympy as s


def raising():
    return s.Matrix([[0,0],[1,0]])


def fermion_creation(index):
    if type(index) is not int or index not in (0,1):
        raise ValueError('one of two fermionic modes required')
    out=s.zeros(4)
    for mask in range(4):
        if not mask>>index&1:
            out[mask|(1<<index),mask]=(-1)**((mask&((1<<index)-1)).bit_count())
    return out


def component_control(degrees,higgs_charge):
    if tuple(degrees) not in ((2,2),(2,4),(4,2)) or higgs_charge not in (-1,1):
        raise ValueError('a retained auxiliary pair and Higgs character required')
    if sum(degrees)+higgs_charge!=5:
        raise ValueError('the unbalanced vertex must have determinant charge five')
    c0,c1=fermion_creation(0),fermion_creation(1)
    pair=c0*c1
    matter_degree=s.diag(0,degrees[0],degrees[1],sum(degrees))
    hdegree=s.diag(0,higgs_charge);xdegree=s.diag(0,-5)
    unbalanced=s.kronecker_product(pair,raising())
    before_degree=s.kronecker_product(matter_degree,s.eye(2))+s.kronecker_product(s.eye(4),hdegree)
    assert before_degree*unbalanced-unbalanced*before_degree==5*unbalanced
    vertex=s.kronecker_product(unbalanced,raising())
    total=s.kronecker_product(before_degree,s.eye(2))+s.kronecker_product(s.eye(8),xdegree)
    generator=vertex+vertex.H
    parity=s.kronecker_product(s.diag(1,-1,-1,1),s.eye(4))
    assert generator==generator.H and generator**3==generator
    assert total*generator==generator*total
    assert parity*generator==generator*parity
    angle=s.pi/4
    evolution=s.eye(16)+(s.cos(angle)-1)*generator**2-s.I*s.sin(angle)*generator
    assert (evolution.H*evolution-s.eye(16)).applyfunc(s.simplify)==s.zeros(16)
    output=evolution*s.eye(16)[:,0]
    probabilities=[s.simplify(s.conjugate(x)*x) for x in output]
    assert sum(probabilities)==1 and probabilities[0]==probabilities[15]==s.Rational(1,2)
    assert all(total[i,i]==0 for i,p in enumerate(probabilities) if p!=0)
    return {'auxiliary_matter_degrees':list(degrees),'Higgs_component_charge':higgs_charge,
            'memory_charge_transition':[0,-5],
            'uncompensated_vertex_charge':5,
            'generator_commutes_with_total_degree_and_fermion_parity':True,
            'declared_dimensionless_control_angle':'pi/4',
            'complete_sixteen_detector_probabilities':list(map(str,probabilities))}


def certificate():
    x=raising();state=s.Matrix([1,1])/s.sqrt(2)
    assert (state.H*x*state)[0]==s.Rational(1,2)
    assert x*x==s.zeros(2)
    return {'memory_representation':'trivial plus inverse determinant of U(3) x U(2)',
            'condition':'source-to-field map retains the full lifted U(3) x U(2) phase action',
            'unique_or_universal_Higgs_requirement':False,
            'memory_transition_phase_charges':[-1,-1,-1,-1,-1],
            'memory_is_even_fermion_parity':True,
            'component_controls':[component_control((2,2),1),
                                  component_control((2,4),-1),
                                  component_control((4,2),-1)],
            'coherent_memory_transition_expectation':'1/2',
            'charge_eigenstate_transition_expectation':'0',
            'finite_transition_nilpotent':True,
            'declared_inputs':['charged even memory register','Hermitian coupling strength',
                               'Higgs register and state','memory preparation'],
            'full_gauge_tensor_condition':'multiply every complete invariant Yukawa vertex by the same inverse-determinant memory transition',
            'two_mode_component_control_is_full_gauge_invariant_vertex':False,
            'constant_Yukawa_or_persistent_condensate_derived':False,
            'GeV_mass_derived':False}
