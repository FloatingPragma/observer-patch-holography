"""An exact finite spectral left-Weyl clock with all cutoff modes retained."""
from functools import lru_cache
from itertools import product
import sympy as s

from ..physical_attachment.chiral import representation
from ..chiral_source.model import currents as source_currents, native_basis, source_descriptor
from ..mass_instrument.model import annihilator

PROGRAM_ID='coherent-native-spectral-weyl-v1'


def clean(matrix):
    return s.ImmutableMatrix(matrix.applyfunc(s.simplify))


def pauli():
    return (s.ImmutableMatrix([[0,1],[1,0]]),
            s.ImmutableMatrix([[0,-s.I],[s.I,0]]),
            s.ImmutableMatrix([[1,0],[0,-1]]))


def exact_real(value):
    value=s.sympify(value)
    if value.has(s.Float) or value.is_real is not True:
        raise ValueError('exact finite real data required')
    if value.has(s.oo,-s.oo,s.zoo,s.nan):
        raise ValueError('finite data required')
    return value


def spectral_momenta(side,length=2*s.pi):
    if type(side) is not int or side<2 or side%2:
        raise ValueError('positive even side at least two required')
    length=exact_real(length)
    if length.is_positive is not True:raise ValueError('positive physical length required')
    return tuple(s.simplify(2*s.pi*(s.Rational(2*n+1,2))/length)
                 for n in range(-side//2,side//2))


def symbol(momentum):
    if len(momentum)!=3:raise ValueError('three physical momenta required')
    p=tuple(exact_real(v) for v in momentum)
    return clean(-sum((q*matrix for q,matrix in zip(p,pauli(),strict=True)),s.zeros(2)))


def energy(momentum):
    p=tuple(exact_real(v) for v in momentum)
    if len(p)!=3:raise ValueError('three physical momenta required')
    return s.sqrt(sum(x*x for x in p))


def projectors(momentum):
    e=energy(momentum)
    if e.is_positive is not True:raise ValueError('nonzero momentum required for filled-band projector')
    h=symbol(momentum)
    return clean((s.eye(2)+h/e)/2),clean((s.eye(2)-h/e)/2)


def fourier_matrix(side):
    momenta=spectral_momenta(side,2*s.pi)
    return s.ImmutableMatrix(side,side,lambda x,n:
        s.expand_complex(s.exp(2*s.pi*s.I*momenta[n]*s.Rational(x,side)))/s.sqrt(side))


def momentum_matrix(side,length=2*s.pi):
    """Dense position derivative with antiperiodic spin structure."""
    spectral_momenta(side,length)
    return s.ImmutableMatrix(side,side,lambda x,y:0 if x==y else
        s.I*s.pi*s.Integer(-1)**(x-y+1)/(length*s.sin(s.pi*s.Rational(x-y,side))))


def cubic_hamiltonian(side,length=2*s.pi):
    """Position-major spin ordering; internal R15 is a spectator factor."""
    p=momentum_matrix(side,length);eye=s.eye(side)
    axes=(s.kronecker_product(p,eye,eye),s.kronecker_product(eye,p,eye),
          s.kronecker_product(eye,eye,p))
    return clean(-sum((s.kronecker_product(q,sigma) for q,sigma in zip(axes,pauli(),strict=True)),s.zeros(2*side**3)))


def spectral_certificate():
    p=s.symbols('p_x p_y p_z',real=True);h=symbol(p);r2=sum(v*v for v in p)
    assert clean(h*h-r2*s.eye(2))==s.zeros(2)
    _,_,currents,_=representation();y=currents[-1]
    trace5=s.simplify(s.trace(y**5))
    assert native_basis()==tuple(tuple(row) for row in representation()[0])
    h30=s.ImmutableSparseMatrix(s.kronecker_product(symbol((1,2,2)),s.eye(15)))
    for current in source_currents():
        assert h30*current-current*h30==s.zeros(30)
    assert trace5!=0 and s.trace((-y.T)**5)==-trace5
    helicity_checks=0
    for p in [(s.Rational(1,2),)*3,(s.Rational(-3,2),s.Rational(1,2),s.Rational(-1,2)),(1,2,2)]:
        e=energy(p);hp=sum((p[i]*pauli()[i] for i in range(3)),s.zeros(2))/(2*e)
        plus,minus=projectors(p)
        assert clean(hp*plus+plus/2)==s.zeros(2)
        assert clean(hp*minus-minus/2)==s.zeros(2)
        # A hole has physical momentum -p and conjugate spin generators -J^T.
        assert clean(hp.T*minus.T-minus.T/2)==s.zeros(2)
        helicity_checks+=3
    return {'program':PROGRAM_ID,'sitewise_source_program':source_descriptor()['program'],
            'actual_native_current_kinetic_commutators_checked':12,'internal_dimension':15,'Cauchy_components_per_cell':30,
            'one_body_symbol':'-sum_i sigma_i p_i tensor I15','squared_energy':'|p|^2',
            'particle_helicity':'-1/2','antiparticle_helicity':'1/2',
            'particle_current':'T_A','hole_current':'-T_A^T',
            'helicity_projector_identities':helicity_checks,
            'Tr_Y_power5':str(trace5),'Tr_conjugate_Y_power5':str(-trace5),
            'hole_is_independent_same_charge_right_Weyl_mirror':False}


def position_certificate():
    rows=[]
    for side in (2,4,6):
        f=fourier_matrix(side);p=momentum_matrix(side);values=spectral_momenta(side)
        assert clean(f.H*f)==s.eye(side)
        assert clean(p.H-p)==s.zeros(side)
        assert clean(p*f-f*s.diag(*values))==s.zeros(side)
        assert set(values)==set(-v for v in values)
        rows.append({'side':side,'one_axis_momenta':[str(v) for v in values],
                     'nonzero_derivative_entries':sum(v!=0 for v in p),
                     'all_N_cubed_momenta_retained':side**3,
                     'one_family_Cauchy_dimension':30*side**3,
                     'one_family_full_Fock_dimension_exponent':30*side**3})
    # A direct complete small cube, without relying solely on a symbol formula.
    side=2;h=cubic_hamiltonian(side);f=fourier_matrix(side)
    cubic=s.kronecker_product(f,f,f,s.eye(2));blocks=s.diag(*(symbol(p) for p in product(spectral_momenta(side),repeat=3)))
    assert clean(h*cubic-cubic*blocks)==s.zeros(16)
    # The non-neighbor matrix element is an explicit finite nonlocality witness.
    far=momentum_matrix(6)[0,3]
    assert far==-s.I/2 and far!=0
    return {'axis_executions':rows,'literal_side_two_cubic_columns_checked':16,
            'finite_nonlocality_witness':{'side':6,'length':'2*pi','from':3,'to':0,
                                          'momentum_matrix_element':str(far),'nearest_neighbor':False},
            'finite_strict_light_cone_proved':False,
            'arbitrary_range_derivative_is_declared_architecture':True}


def filled_reference_certificate():
    rows=[]
    for side in (2,4,6):
        momenta=list(product(spectral_momenta(side),repeat=3))
        squared=[sum(q*q for q in p) for p in momenta]
        assert min(squared)==s.Rational(3,4)
        # One negative and one positive spin eigenvalue for each internal row.
        assert all(q>0 for q in squared)
        rows.append({'side':side,'physical_length':'2*pi','retained_momenta':len(momenta),
                     'occupied_negative_modes':15*len(momenta),
                     'source_charge_relative_to_reference':'N-15*V',
                     'full_excitation_gap':'sqrt(3)/2','neutral_pair_threshold':'sqrt(3)'})
    return {'finite_ground':'unique pure Slater determinant filling all negative bands',
            'sea_energy':'-15*sum_p |p|','positive_excitation_energy':'sum_p |p|*(b_dagger*b+d_dagger*d)',
            'general_finite_gap':'sqrt(3)*pi/L','gap_vanishes_in_infinite_volume':True,
            'infinite_volume_rest_mass':'0','finite_volume_gap_is_rest_mass':False,
            'rows':rows}


def instrument_certificate(momentum=(s.Rational(1,2),)*3):
    """All outcomes of a two-mode block, including unsuccessful preparations."""
    h=symbol(momentum);e=energy(momentum);_,pm=projectors(momentum)
    c=[annihilator(i,2) for i in range(2)];dim=4
    H=clean(sum((h[i,j]*c[i].H*c[j] for i,j in product(range(2),repeat=2)),s.zeros(dim))+e*s.eye(dim))
    number=s.diag(0,1,1,2);rho=s.zeros(dim)
    for i,j in product(range(2),repeat=2):rho[1<<i,1<<j]=pm[i,j]
    rho=clean(rho)
    assert s.trace(rho)==1 and clean(rho*rho-rho)==s.zeros(dim) and clean(H*rho)==s.zeros(dim)
    effects=[]
    for j in range(3):
        ep=s.eye(dim)
        for k in range(3):
            if j!=k:ep=clean(ep*(H/e-k*s.eye(dim))/(j-k))
        for charge in (-1,0,1):
            cq=s.diag(*(int(m.bit_count()-1==charge) for m in range(dim)))
            effect=clean(ep*cq)
            if effect!=s.zeros(dim):effects.append((j,charge,effect))
    assert len(effects)==4 and clean(sum((p for _,_,p in effects),s.zeros(dim)))==s.eye(dim)
    n=c[0].H*c[0]
    prep=(('add',c[0].H/s.sqrt(2)),('add_fail',-n/s.sqrt(2)),
          ('remove',c[0]/s.sqrt(2)),('remove_fail',(s.eye(dim)-n)/s.sqrt(2)))
    assert clean(sum((k.H*k for _,k in prep),s.zeros(dim)))==s.eye(dim)
    records=[]
    for label,k in prep:
        out=clean(k*rho*k.H)
        for j,q,effect in effects:
            prob=s.simplify(s.trace(effect*out))
            assert prob.is_nonnegative is True
            records.append({'preparation':label,'energy':str(j*e),'relative_charge':q,'probability':str(prob)})
    assert sum(s.sympify(row['probability']) for row in records)==1
    return {'momentum':[str(v) for v in momentum],'energy':str(e),'records':records,
            'complete_records':len(records),'nonzero_records':sum(row['probability']!='0' for row in records),
            'failure_records_retained':True,'rest_mass_inferred_from_finite_momentum_energy':False}


def tail_bound(cutoff,sobolev_order,norm):
    cutoff=exact_real(cutoff);norm=exact_real(norm)
    if cutoff.is_positive is not True or norm.is_nonnegative is not True or type(sobolev_order)!=int or sobolev_order<0:
        raise ValueError('positive cutoff, nonnegative norm and integer Sobolev order required')
    return s.simplify(norm/(1+cutoff**2)**s.Rational(sobolev_order,2))
