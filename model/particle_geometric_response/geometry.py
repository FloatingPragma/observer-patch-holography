"""Six spatial strain operators, native covariance and complete work records.

The coframe is a declared action input. Its mechanical response is not
identified with emergent proper geometry or a completed Einstein premise.
"""
from functools import lru_cache
from itertools import combinations
import sympy as s

from model.particle_local_wilson_density import density
from model.particle_observer_outputs import field_route

SOURCE_INPUTS = (
    'model/particle_geometric_response/geometry.py',
    'model/particle_geometric_response/PROOF.md',
    'model/particle_field_state_response/state_response.py',
    'model/particle_field_state_response/PROOF.md',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_local_wilson_density/density.py',
    'model/particle_gravity_response/PROOF.md',
    'model/particle_physical_causality/PROOF.md',
    'model/particle_field_causality/PROOF.md',
    'model/matched_einstein_limit/PROOF.md',
    'Lean/QFT/CoarseEnergyWard.lean',
    'Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean',
)


def clean(x):
    # Controls have Gaussian-rational entries (and a cancelled sqrt2 spin
    # chart). Expansion/cancellation avoids generic simplifier searches.
    def normal(z):
        z=s.expand(z)
        if z.is_Rational: return z
        re,im=z.as_real_imag()
        if re.is_Rational and im.is_Rational: return re+s.I*im
        return s.cancel(z)
    return s.ImmutableSparseMatrix(x).applyfunc(normal)


def coframe_data(v):
    """Exact real, orientation-preserving constant coframes only."""
    v = s.ImmutableMatrix(v)
    if v.shape != (3, 3) or v.atoms(s.Float) or any(x.is_real is not True for x in v):
        raise ValueError('a real3-by-3 coframe is required')
    volume = s.factor(v.det())
    if volume.is_positive is not True:
        raise ValueError('a strictly positive coframe determinant is required')
    metric = clean(v.T*v)
    return {'V': v, 'volume': volume, 'metric': metric,
            'inverse_frame': clean(v.inv().T),
            'inverse_metric': clean(metric.inv()),
            'electric_magnetic_coefficient': clean(metric/volume)}


@lru_cache(maxsize=1)
def strain_basis():
    result = []
    for i in range(3):
        q = s.zeros(3); q[i, i] = 1
        result.append((str(i)+str(i), s.ImmutableMatrix(q)))
    for i, j in ((0, 1), (0, 2), (1, 2)):
        q = s.zeros(3); q[i, j] = q[j, i] = s.Rational(1, 2)
        result.append((str(i)+str(j), s.ImmutableMatrix(q)))
    return tuple(result)


@lru_cache(maxsize=1)
def native_links():
    return (density.native_weak_link(0, s.Rational(3, 5), s.Rational(4, 5)),
            density.native_weak_link(1, s.Rational(5, 13), s.Rational(12, 13)),
            density.native_weak_link(2, s.Rational(8, 17), s.Rational(15, 17)))


@lru_cache(maxsize=8)
def link_data(links):
    """The complete native15 rows at a one-cell coordinate configuration."""
    links = tuple(clean(x) for x in links)
    if len(links) != 3 or any(x.shape != (15, 15) or clean(x.H*x) != s.eye(15)
                              for x in links):
        raise ValueError('three native fifteen-row unitary links are required')
    difference = tuple(clean((u-u.H)/(2*s.I)) for u in links)
    forward = tuple(s.eye(15)-u for u in links)
    laplace = tuple(tuple(clean((forward[i].H*forward[j]
                                +forward[j].H*forward[i])/2)
                          for j in range(3)) for i in range(3))
    loops = tuple(clean(links[j]*links[k]*links[j].H*links[k].H)
                  for j, k in ((1, 2), (2, 0), (0, 1)))
    curvature = tuple(s.eye(15)-u for u in loops)
    magnetic = s.ImmutableMatrix(3, 3, lambda i, j:
        s.simplify(s.re(s.trace(curvature[i].H*curvature[j]))/30-int(i == j)))
    return {'links': links, 'D': difference, 'L': laplace,
            'loops': loops, 'magnetic_gram_minus_identity': magnetic}


def native_symbol(v, links=None, kappa=s.Rational(1, 7)):
    """At bare spacing1; all60 native spin/internal columns are retained."""
    kappa = s.sympify(kappa)
    if kappa.has(s.Float) or kappa.is_real is not True or kappa.is_nonnegative is not True:
        raise ValueError('a nonnegative exact magnetic coefficient is required')
    d = coframe_data(v); link = link_data(native_links() if links is None else links)
    gamma = density.gamma_data(); beta = gamma[3]
    alpha = tuple(clean(s.I*beta*g) for g in gamma[:3])
    h = s.zeros(60)
    for A in range(3):
        for i in range(3):
            if d['inverse_frame'][A,i]:
                h += d['inverse_frame'][A, i]*s.kronecker_product(alpha[A], link['D'][i])
    for i in range(3):
        for j in range(3):
            if d['inverse_metric'][i,j]:
                h += d['inverse_metric'][i, j]*s.kronecker_product(beta, link['L'][i][j])/2
    magnetic = s.simplify(kappa*sum(d['electric_magnetic_coefficient'][i, j]
                                  *link['magnetic_gram_minus_identity'][i, j]
                                  for i in range(3) for j in range(3)))
    return {'matter': clean(h), 'magnetic': magnetic,
            'total': clean(h+magnetic*s.eye(60)), 'coframe': d, 'link_data': link}


def native_strain(q, links=None, kappa=s.Rational(1, 7)):
    q = s.ImmutableMatrix(q)
    if q.shape != (3, 3) or q != q.T or any(x.is_real is not True for x in q):
        raise ValueError('strain must be a real symmetric3-by-3 matrix')
    link = link_data(native_links() if links is None else links)
    gamma = density.gamma_data(); beta = gamma[3]
    alpha = tuple(clean(s.I*beta*g) for g in gamma[:3])
    h = s.zeros(60)
    for A in range(3):
        for i in range(3):
            if q[A,i]:
                h -= q[A, i]*s.kronecker_product(alpha[A], link['D'][i])
    for i in range(3):
        for j in range(3):
            if q[i,j]:
                h -= q[i, j]*s.kronecker_product(beta, link['L'][i][j])
    coefficient = 2*q-s.trace(q)*s.eye(3)
    magnetic = s.simplify(kappa*sum(coefficient[i, j]
                                  *link['magnetic_gram_minus_identity'][i, j]
                                  for i in range(3) for j in range(3)))
    return {'matter': clean(h), 'magnetic': magnetic,
            'total': clean(h+magnetic*s.eye(60))}


@lru_cache(maxsize=1)
def native_control():
    links = native_links(); base = native_symbol(s.eye(3), links)
    gamma=density.gamma_data(); beta=gamma[3]
    # The original Wilson projectors provide an independent source expression.
    # Its square-decomposition theorem is not needed to test this identity.
    dw=3*s.eye(60)
    for i,u in enumerate(links):
        dw-=s.kronecker_product((s.eye(4)-gamma[i])/2,u)
        dw-=s.kronecker_product((s.eye(4)+gamma[i])/2,u.H)
    source=clean(s.kronecker_product(beta,s.eye(15))*dw)
    if base['matter'] != source:
        raise ArithmeticError('flat coframe is not the actual Wilson source')
    expected_mag = -s.Rational(1, 7)*sum(s.re(s.trace(u))/15
                                       for u in base['link_data']['loops'])
    if s.simplify(base['magnetic']-expected_mag) != 0:
        raise ArithmeticError('flat magnetic action changed')
    strains = [(name, native_strain(q, links)) for name, q in strain_basis()]
    if any(row['total'] != row['total'].H for _, row in strains):
        raise ArithmeticError('a complete strain operator is not Hermitian')
    # Independent gauge and physical-frame rotations on a nonorthogonal frame.
    v = s.Matrix([[1, s.Rational(1, 3), 0], [0, 1, s.Rational(1, 5)], [0, 0, 1]])
    value = native_symbol(v, links)
    gauge = links[2]
    ug = clean(s.kronecker_product(s.eye(4), gauge))
    changed = native_symbol(v, tuple(clean(gauge*u*gauge.H) for u in links))
    gauge_ok = clean(ug*value['total']*ug.H-changed['total']) == s.zeros(60)
    alpha = tuple(clean(s.I*beta*g) for g in density.gamma_data()[:3])
    spin_numerator = clean(s.eye(4)-alpha[0]*alpha[1])
    rotation = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    ur = clean(s.kronecker_product(spin_numerator, s.eye(15)))
    rotated = native_symbol(rotation*v, links)
    rotation_ok = clean(ur*value['total']*ur.H/2-rotated['total']) == s.zeros(60)
    if not gauge_ok or not rotation_ok:
        raise ArithmeticError('native gauge or frame covariance failed')
    return {'native_internal_rows': 15, 'Dirac_rows': 4,
            'complete_one_body_columns': 60,
            'all_six_strain_labels': [name for name, _ in strains],
            'all_six_strain_HS_squared': [s.simplify(s.trace(r['total']**2)) for _, r in strains],
            'flat_Wilson_source_equality': True,
            'flat_magnetic_source_equality': True,
            'native_gauge_covariance': gauge_ok, 'constant_frame_rotation_covariance': rotation_ok,
            'magnetic_scalar_is_second_quantized': False,
            'coordinate_symbol_is_a_normalizable_quantum_preparation': False,
            'full_Fock_lift_is_algebraic_not_dense_execution': True}


def rotor_strains():
    from model.particle_field_state_response.state_response import exact_control
    d = exact_control()
    he, hw, hd = d['HE'], d['HW'], d['HD']
    if clean(he+hw+hd-d['H']) != s.zeros(208):
        raise ArithmeticError('thermal/strain source Hamiltonians differ')
    return {'source': d, 'longitudinal': clean(he-hd-2*hw),
            'isotropic': clean(-he-hd-2*hw), 'detector_bound': s.Integer(28)}


def exact_mean(rho, operator):
    return s.factor(sum(value*operator[j, i] for (i, j), value in rho.todok().items()))


@lru_cache(maxsize=1)
def rotor_record_control():
    from model.particle_field_state_response.state_response import prepared_density
    data = rotor_strains(); source = data['source']; prepared = prepared_density()
    rho = prepared['rho']; R = data['detector_bound']
    eps = s.Rational(prepared['preparation_trace_error_upper'].numerator,
                     prepared['preparation_trace_error_upper'].denominator)
    rows=[]
    for name in ('longitudinal', 'isotropic'):
        op=data[name]
        if op != op.H or clean(op*source['Gauss']-source['Gauss']*op) != s.zeros(208):
            raise ArithmeticError('strain lost Hermiticity or the source gauge constraint')
        row_bound=max(sum(abs(op[i,j]) for j in range(op.cols)) for i in range(op.rows))
        if row_bound > R:
            raise ArithmeticError('complete mean-reading effects are not certified positive')
        mean=exact_mean(rho, op)
        probabilities=[s.factor((1+mean/R)/2),s.factor((1-mean/R)/2)]
        if any(p.is_nonnegative is not True for p in probabilities) or sum(probabilities) != 1:
            raise ArithmeticError('complete strain records are not a probability law')
        rows.append({'strain':name,'record_labels':[R,-R],
                     'original_probabilities':probabilities,
                     'exact_prepared_work_slope':mean,
                     'operator_row_sum_bound':row_bound,
                     'target_Gibbs_mean_error_upper':2*R*eps,
                     'Gauss_commutator_zero':True})
    energy=exact_mean(rho,source['H'])
    pressure=-rows[1]['exact_prepared_work_slope']/3
    energy_error=44*eps; pressure_error=2*R*eps/3
    ratio=s.factor(pressure/energy)
    if energy<=energy_error:
        raise ArithmeticError('the energy denominator lacks a positive target bound')
    ratio_error=s.factor(pressure_error/(energy-energy_error)
                         +abs(pressure)*energy_error/(energy*(energy-energy_error)))
    components={name:exact_mean(rho,source[name]) for name in ('HE','HW','HD')}
    if s.factor(sum(components.values())-energy) != 0:
        raise ArithmeticError('the thermal energy split changed the common source')
    return {'source':'unchanged field_route.exact_rotor(6)',
            'magnetic_coefficient':0,
            'complete_dimension':208,'inverse_temperature_target':s.Rational(1,64),
            'same_clock':'exp(-it H); no coframe redefinition of duration',
            'prepared_density':'positive order12 amplitude polynomial, normalized without postselection',
            'target_state_trace_error_upper':eps,
            'all_CAR_occupations':16,'all_electric_fluxes':list(range(-6,7)),
            'stress_instruments':rows,'energy_components_prepared':components,
            'energy_prepared':energy,'pressure_prepared_at_declared_volume1':pressure,
            'pressure_over_energy_density_prepared':ratio,
            'target_Gibbs_pressure_error_upper':pressure_error,
            'target_Gibbs_pressure_ratio_error_upper':ratio_error,
            'global_Gibbs_entropy_slope_is_regional_geometric_entropy':False,
            'pressure_is_a_physical_cosmological_equation_of_state':False,
            'rotor_is_a_normalizable_reducing_subspace_of_full3D_quantum_field':False}


def serial(x):
    if isinstance(x, s.Basic): return str(x)
    if isinstance(x, dict): return {k: serial(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [serial(v) for v in x]
    return x


def exact_rational(value, name, positive=False):
    value=s.sympify(value)
    if value.is_Rational is not True or (value<=0 if positive else value<0):
        raise ValueError(name+' must be an exact '+('positive' if positive else 'nonnegative')+' rational')
    return value


def coframe_causal_bound(v, magnetic_scale=s.Rational(1,100), weight=2):
    """Rational parent interaction rate for one constant coframe.

    All collapsed periodic cells are allowed. The bound is on the parent
    interaction decomposition; a global word cutoff adds its separate error.
    """
    d=coframe_data(v)
    if any(e.is_Rational is not True for e in d['metric']):
        raise ValueError('the rate certificate requires a rational metric')
    k=exact_rational(magnetic_scale,'magnetic scale')
    lam=exact_rational(weight,'spatial weight',True)
    if lam<=1: raise ValueError('spatial weight must exceed one')
    C=d['inverse_metric']; A=d['electric_magnetic_coefficient']
    if any(e.is_Rational is not True for e in A):
        raise ValueError('the rate certificate requires rational G/det(V)')
    diagonal=120*lam*sum(max(s.Integer(1),C[i,i]) for i in range(3))
    mixed=1080*lam**2*sum(abs(C[i,j]) for i,j in combinations(range(3),2))
    magnetic=k*lam**2*(16*s.trace(A)+64*sum(abs(A[i,j]) for i,j in combinations(range(3),2)))
    return {'nearest_bond_rate_upper':diagonal,'mixed_Wilson_rate_upper':mixed,
            'magnetic_rate_upper':magnetic,'weighted_rate_times_spacing_upper':diagonal+mixed+magnetic,
            'spatial_weight':lam,'magnetic_scale_upper':k,
            'electric_cross_terms_are_tail_cell_onsite':True,
            'self_loop_terms_absorbed_into_onsite_generator':True,
            'collapsed_periodic_neighbors_allowed':True}


def compact_coframe_family(min_singular=s.Rational(1,2), max_singular=2,
                         magnetic_scale=s.Rational(1,100), weight=2):
    """A rate uniform in constant coframes with m^2 I <= V^T V <= L^2 I."""
    m=exact_rational(min_singular,'minimum singular value',True)
    L=exact_rational(max_singular,'maximum singular value',True)
    k=exact_rational(magnetic_scale,'magnetic scale')
    lam=exact_rational(weight,'spatial weight',True)
    if L<m or lam<=1: raise ValueError('require L>=m>0 and spatial weight>1')
    W=360*lam*max(s.Integer(1),m**-2)+3240*lam**2/m**2+240*k*lam**2*L**2/m**3
    mu=(lam-1)/lam
    return {'min_singular':m,'max_singular':L,'magnetic_scale_upper':k,
            'spatial_weight':lam,'weighted_rate_times_spacing_upper':W,
            'log_weight_rational_lower':mu,'physical_speed_rational_upper':2*L*W/mu,
            'physical_distance':'Euclidean on open carriers; quotient Euclidean on periodic carriers',
            'coframe_cube_to_coordinate_halfwidth_factor':2/m,
            'constant_coframes_only':True,
            'proper_geometry_Lorentz_covariance_or_limiting_state_proved':False}


def principal_minors(a):
    return [s.factor(a.extract(indices,indices).det()) for size in (1,2,3)
            for indices in combinations(range(3),size)]


def coframe_family_member(v, family):
    """Validate every principal minor, not only the leading minors."""
    expected=compact_coframe_family(family['min_singular'],family['max_singular'],
                                   family['magnetic_scale_upper'],family['spatial_weight'])
    if family!=expected: raise ValueError('compact coframe family was changed')
    d=coframe_data(v); G=d['metric']; m=family['min_singular']; L=family['max_singular']
    low=principal_minors(G-m*m*s.eye(3)); high=principal_minors(L*L*s.eye(3)-G)
    if any(x.is_nonnegative is not True for x in low+high):
        raise ValueError('coframe lies outside the declared compact family')
    point=coframe_causal_bound(v,family['magnetic_scale_upper'],family['spatial_weight'])
    if point['weighted_rate_times_spacing_upper']>family['weighted_rate_times_spacing_upper']:
        raise ArithmeticError('pointwise coframe rate exceeds the compact majorant')
    return {'lower_metric_principal_minors':low,'upper_metric_principal_minors':high,
            'all_fourteen_principal_minors_nonnegative':True,'single_coframe_rate':point}


def coframe_causality_certificate():
    family=compact_coframe_family(magnetic_scale=s.Rational(1,7))
    V=s.Matrix([[1,s.Rational(1,3),0],[0,1,s.Rational(1,5)],[0,0,1]])
    flat=coframe_causal_bound(s.eye(3),magnetic_scale=s.Rational(1,7))
    if flat['weighted_rate_times_spacing_upper']!=s.Rational(5232,7):
        raise ArithmeticError('flat Wilson causal rate changed')
    return {'compact_family':family,'native_control_coframe':V.tolist(),
            'native_coordinate_control_spacing':1,
            'native_coordinate_control_magnetic_coefficient':s.Rational(1,7),
            'native_control_membership':coframe_family_member(V,family),
            'flat_source_rate_recovered':True,'flat_source_rate':flat,
            'shared_rotor_thermal_state_transferred_to_spatial_refinement':False}


def certificate():
    return serial({'schema': 'oph.spatial_coframe_response.v1',
                   'source_inputs': list(SOURCE_INPUTS),
                   'spatial_coframe_action': {
                       'metric': 'G=V^T V', 'volume': 'v=det(V)>0',
                       'Dirac_coefficient': 'V^(-T)', 'Wilson_coefficient': 'G^(-1)',
                       'electric_magnetic_coefficient': 'G/v',
                       'dilation_derivative': '-H_E-H_D-2H_W-H_B',
                       'metric_probe_is_declared': True},
                   'native_sixty_column_stress': native_control(),
                   'shared_quantum_rotor_work':rotor_record_control(),
                   'compact_coframe_causality':coframe_causality_certificate(),
                   'metric_identification_boundary': {
                       'same_flat_action_fixes_unique_metric_response': False,
                       'proper_geometric_stress_or_Einstein_limit_claimed': False,
                       'physical_particle_or_GeV_assignment': False}})
