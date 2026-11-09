"""A declared coframe coupling with exact primary flat action and strain records."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_geometric_response import geometry
from model.particle_local_wilson_density import density
from model.particle_primary_gravity_premises import model as modular
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS = (
    'model/particle_primary_open_coframe/__init__.py',
    'model/particle_primary_open_coframe/__main__.py',
    'model/particle_primary_open_coframe/model.py',
    'model/particle_primary_open_coframe/test_model.py',
    'model/particle_primary_open_coframe/PROOF.md',
    'model/particle_primary_open_coframe/README.md',
    'model/particle_common_contract/PROOF.md',
    'model/particle_geometric_response/PROOF.md',
    'model/particle_primary_gravity_premises/PROOF.md',
    'model/particle_primary_stress/PROOF.md',
    'model/matched_einstein_limit/PROOF.md',
    'Lean/QFT/CoarseEnergyWard.lean',
    'Lean/Thermodynamics/FirstLawIdentity.lean',
)


def shift_data(shifts):
    """Parent covariant partial shifts, before product link compression."""
    if len(shifts) != 3:
        raise ValueError('three open shifts are required')
    shifts = tuple(geometry.clean(t) for t in shifts)
    size = shifts[0].rows
    if not size or any(t.shape != (size, size) for t in shifts):
        raise ValueError('common nonempty shift carrier required')
    identity = s.eye(size)
    holes = tuple(geometry.clean(identity-t.H*t) for t in shifts)
    if any(p != p.H or p*p != p for p in holes):
        raise ValueError('open shifts must have projection initial defects')
    forward = tuple(identity-t for t in shifts)
    derivative = tuple(geometry.clean((t-t.H)/(2*s.I)) for t in shifts)
    laplace = tuple(tuple(geometry.clean(
        (forward[i].H*forward[j]+forward[j].H*forward[i])/2
        + (holes[i] if i == j else s.zeros(size)))
        for j in range(3)) for i in range(3))
    for i in range(3):
        if laplace[i][i] != geometry.clean(2*identity-shifts[i]-shifts[i].H):
            raise ArithmeticError('the open boundary completion failed')
    return {'shifts': shifts, 'holes': holes, 'D': derivative, 'L': laplace}


def matter(v, shifts):
    """All four Dirac rows on a finite coordinate control; bare spacing one."""
    frame = geometry.coframe_data(v); data = shift_data(shifts)
    gamma = density.gamma_data(); beta = gamma[3]
    alpha = tuple(geometry.clean(s.I*beta*g) for g in gamma[:3])
    size = shifts[0].rows; h = s.zeros(4*size)
    for A, i in product(range(3), repeat=2):
        h += frame['inverse_frame'][A, i]*s.kronecker_product(alpha[A], data['D'][i])
        h += frame['inverse_metric'][A, i]*s.kronecker_product(beta, data['L'][A][i])/2
    return geometry.clean(h)


def strain(q, shifts):
    """Exact symmetric spatial derivative of the boundary-completed matter term."""
    q = s.ImmutableMatrix(q)
    if q.shape != (3, 3) or q != q.T or q.has(s.Float) or any(x.is_real is not True for x in q):
        raise ValueError('an exact real symmetric strain is required')
    data = shift_data(shifts); gamma = density.gamma_data(); beta = gamma[3]
    alpha = tuple(geometry.clean(s.I*beta*g) for g in gamma[:3])
    h = s.zeros(4*shifts[0].rows)
    for A, i in product(range(3), repeat=2):
        h -= q[A, i]*s.kronecker_product(alpha[A], data['D'][i])
        h -= q[A, i]*s.kronecker_product(beta, data['L'][A][i])
    return geometry.clean(h)


def open_shifts(side=2, phases=None):
    side = ground._i(side, 'small open side', 1)
    if side > 3:
        raise ValueError('dense coordinate control is limited to side three')
    cells = tuple(product(range(side), repeat=3)); indices = {c:j for j,c in enumerate(cells)}
    phases = (s.Integer(1), s.Integer(1), s.Integer(1)) if phases is None else tuple(phases)
    if len(phases) != 3 or any(s.simplify(s.conjugate(z)*z-1) != 0 for z in phases):
        raise ValueError('three exact unit phases required')
    shifts = []
    for i in range(3):
        t = s.zeros(len(cells))
        for c in cells:
            if c[i]+1 < side:
                d = tuple(x+int(k==i) for k,x in enumerate(c))
                # T_(x,y) is a positively oriented covariant source link.
                t[indices[c], indices[d]] = phases[i]
        shifts.append(s.ImmutableSparseMatrix(t))
    return tuple(shifts)


def magnetic_gram(loops):
    """One common base point; None retains an absent open plaquette."""
    if len(loops) != 3:
        raise ValueError('three oriented plane slots required')
    present = [x for x in loops if x is not None]
    size = present[0].rows if present else 15
    if size != 15 or any(u.shape != (15,15) or geometry.clean(u.H*u) != s.eye(15) for u in present):
        raise ValueError('complete native15 parent loops required')
    curvatures = [s.zeros(15) if u is None else s.eye(15)-u for u in loops]
    return s.ImmutableMatrix(3,3,lambda i,j:
        s.simplify(s.re(s.trace(curvatures[i].H*curvatures[j]))/30
                   -int(i==j and loops[i] is not None)))


def strain_bounds(member, q):
    contract.verify_member(member); q = s.ImmutableMatrix(q)
    if q.shape != (3,3) or q != q.T or any(not x.is_Rational for x in q):
        raise ValueError('rational symmetric source strain required')
    source = member['source']; L=source['side']; V=L**3
    a=F(source['spacing']); ell=source['per_link_word_degree']
    e=F(source['electric_coefficient']); k=F(source['magnetic_coefficient'])
    total=sum((abs(F(x)) for x in q),F()); diagonal=sum(abs(F(q[i,i])) for i in range(3))
    c=2*q-s.trace(q)*s.eye(3)
    cdiag=sum(abs(F(c[i,i])) for i in range(3))
    coff=sum(abs(F(c[i,j])) for i,j in product(range(3),repeat=2) if i!=j)
    matter_cap=F(60*V)/a*(5*total+diagonal)
    electric_cap=F(V)*e/(2*a)*F(29,5)*ell**2*(cdiag+coff)
    magnetic_cap=F(V)*k/a*(cdiag+2*coff)
    return {'matter_full_Fock_norm_upper':matter_cap,
            'electric_norm_upper':electric_cap,'magnetic_norm_upper':magnetic_cap,
            'full_strain_norm_upper':matter_cap+electric_cap+magnetic_cap,
            'boundary_strain_norm_upper':F(30*L*L)/a*diagonal,
            'complete_record_values':[-(matter_cap+electric_cap+magnetic_cap),
                                      matter_cap+electric_cap+magnetic_cap],
            'effects':('(I plus_or_minus K_S/R_S)/2 on full physical support'
                       if total else 'I/2 for both zero-score records'),
            'instrument':'positive square roots of both effects; arbitrary idle reference retained',
            'root_service_declared':True,'postselection_normalization':False}


@lru_cache(maxsize=1)
def flat_control():
    shifts=open_shifts(2,(s.Rational(3,5)+s.I*s.Rational(4,5),s.I,-1))
    h=matter(s.eye(3),shifts); gamma=density.gamma_data(); beta=gamma[3]
    direct=3*s.kronecker_product(beta,s.eye(8))
    for i,t in enumerate(shifts):
        hop=(-beta+beta*gamma[i])/2
        direct+=s.kronecker_product(hop,t)+s.kronecker_product(hop.H,t.H)
    if h != geometry.clean(direct):
        raise ArithmeticError('completed flat action differs from the actual source')
    data=shift_data(shifts)
    wilson=sum((s.kronecker_product(beta,data['L'][i][i])/2 for i in range(3)),s.zeros(32))
    if strain(s.eye(3),shifts) != geometry.clean(-h-wilson):
        raise ArithmeticError('the source dilation identity failed')
    return {'open_side':2,'full_Dirac_spatial_columns':32,
            'flat_source_equality':True,'dilation_identity':True,
            'native_internal_carrier_proof_dimension':15,
            'all_native_rows_retained_in_operator_proof':True,
            'coordinate_control_is_prepared_quantum_source':False}


def selected_response(member):
    contract.verify_member(member); flat_control()
    rows=[{'strain':name,**strain_bounds(member,q)} for name,q in geometry.strain_basis()]
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'coupling':'boundary-completed constant spatial coframe, declared architecture',
        'flat_primary_action_identical':True,'selected_flat_state_reference_clock_unchanged':True,
        'parent_products_before_link_compression':True,'source_Gauss_and_CAR_parity_preserved':True,
        'all_six_spatial_strains':rows,
        'dilation_identity':'K_I=-H-H_W with completed open H_W',
        'coordinate_pressure_mean':'Tr(rho*(H+H_W))/(3*V*a^3)',
        'uncompleted_flat_boundary_error':modular.boundary_obstruction(member),
        'complete_mean_detector_records':2,'operator_state_disturbance_retained':True,
        'metric_coupling_unique_or_forced':False,'proper_metric_recovered':False,
        'symmetric_Lorentz_stress_derived':False,'diamond_modular_or_area_law':False,
        'Einstein_equation_derived':False,'physical_particle_mass_table_GeV':[]})


def bind_member(member):
    result=selected_response(member)
    return {'result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):
        raise ValueError('open coframe or source-bound strain data changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_open_coframe.v1','members':[bind_member(contract.member(n)) for n in (0,1,2)],
            'flat_control':flat_control(),'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:
        raise ValueError('coframe inputs changed during replay')
    return result
