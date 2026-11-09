"""Exact spatial Weyl hopping with the configured M15 quantum links."""
from functools import lru_cache
from itertools import combinations_with_replacement

import sympy as s

from model.common_source import native_matter_attachment as matter
from model.native_quantum_gauss import algebra as native


PAULI = (s.ImmutableMatrix([[0, 1], [1, 0]]),
         s.ImmutableMatrix([[0, -s.I], [s.I, 0]]),
         s.ImmutableMatrix.diag(1, -1))


def _axis(axis):
    if type(axis) is not int or not 0 <= axis < 3:
        raise ValueError('a spatial axis zero, one or two is required')
    return axis


@lru_cache(maxsize=1)
def _source():
    currents = tuple(s.ImmutableSparseMatrix(g) for g in native.native_currents())
    r = matter.representation()
    encoder = r['encoder']
    expected = tuple(encoder.T*matter.exterior_derivative(matter.spin_lift(g))*encoder
                     for g in matter.canonical_currents())
    if currents != expected or len(currents) != 12:
        raise ArithmeticError('the native M15 current attachment differs')
    blocks = tuple(indices for _, indices in native.native_blocks())
    labels = ((-1, -1),)+tuple((a, b) for block in blocks for a in block for b in block)
    if len(labels) != 60 or sorted(a for block in blocks for a in block) != list(range(15)):
        raise ArithmeticError('the full native link and matter blocks are required')
    return currents, blocks, labels


def currents():
    """Immutable matrices in the actual canonical port-current ordering."""
    return tuple(s.ImmutableSparseMatrix(g) for g in _source()[0])


def link_labels():
    return tuple(_source()[2])


def link_generator(generator, endpoint, *, orientation=1):
    if endpoint not in (0, 1) or type(endpoint) is not int:
        raise ValueError('one of the two endpoints is required')
    if type(orientation) is not int or orientation not in (-1, 1):
        raise ValueError('an exact orientation sign is required')
    generator = s.ImmutableSparseMatrix(generator)
    if generator.shape != (15, 15) or generator+generator.H != s.zeros(15):
        raise ValueError('an anti-Hermitian fifteen-state generator is required')
    labels = link_labels(); position = {label: i for i, label in enumerate(labels)}
    out = {}
    for column, (a, b) in enumerate(labels[1:], 1):
        for (i, j), value in generator.todok().items():
            if endpoint == 0 and i == a:
                target, coefficient = (j, b), -value
            elif endpoint == 1 and j == b:
                target, coefficient = (a, i), value
            else:
                continue
            if target not in position:
                raise ValueError('generator does not preserve the actual link blocks')
            out[position[target], column] = orientation*coefficient
    return s.ImmutableSparseMatrix(60, 60, out)


def hopping_coefficients(axis):
    """Normal-ordered c_p^dagger c_q tensor link coefficients.

    Each key is (creation mode, annihilation mode, link output, link input).
    There are thirty physical modes at each endpoint; no Fock cutoff occurs
    in this operator expansion. The Weyl coefficient is i sigma_axis/2.
    """
    spin = s.I*PAULI[_axis(axis)]/2
    out = {}
    for link, (a, b) in enumerate(link_labels()[1:], 1):
        for (u, v), coefficient in spin.todok().items():
            p, q = 15*u+a, 30+15*v+b
            out[p, q, link, 0] = coefficient
            out[q, p, 0, link] = s.conjugate(coefficient)
    return out


def coefficient_matrix(axis):
    """Faithful coefficient array, also the full one-particle/link operator."""
    return s.ImmutableSparseMatrix(3600, 3600,
        {(60*p+l, 60*q+k): value
         for (p, q, l, k), value in hopping_coefficients(axis).items()})


def gauss_coefficients(generator, endpoint, *, orientation=1):
    link = link_generator(generator, endpoint, orientation=orientation)
    entries = {}
    for spin in range(2):
        for (a, b), value in generator.todok().items():
            p, q = 30*endpoint+15*spin+a, 30*endpoint+15*spin+b
            for label in range(60):
                entries[60*p+label, 60*q+label] = value
    for (row, col), value in link.todok().items():
        for mode in range(60):
            key = 60*mode+row, 60*mode+col
            entries[key] = entries.get(key, 0)+value
    return s.ImmutableSparseMatrix(3600, 3600, entries)


def anomaly_certificate():
    """Every cubic Lie-algebra component, with one M15 per Weyl chirality."""
    gs = currents()
    linear = tuple(s.expand(s.trace(g)) for g in gs)
    entries = []
    for a, b, c in combinations_with_replacement(range(12), 3):
        value = s.expand(s.trace(gs[a]*(gs[b]*gs[c]+gs[c]*gs[b])))
        if value:
            raise ArithmeticError('the same-source cubic anomaly tensor is nonzero')
        entries.append((a, b, c))
    if any(linear):
        raise ArithmeticError('the same-source linear anomaly trace is nonzero')
    r = matter.representation(); ix = list(r['selected'])
    casimir = r['operators']['weak_C2'].extract(ix, ix)
    doublet_components = sum(value == s.Rational(3, 4) for value in casimir.diagonal())
    if doublet_components != 8:
        raise ArithmeticError('the source weak-doublet parity differs')
    return {'internal_dimension': 15, 'native_generator_count': 12,
            'linear_trace_components': len(linear),
            'symmetric_cubic_tensor_components': len(entries),
            'all_linear_and_cubic_components_zero': True,
            'weak_doublets_counting_color': doublet_components//2,
            'physical_spatial_spin_counted_as_an_extra_chiral_flavor': False,
            'regulated_chiral_quantum_measure_constructed': False}


def ward_certificate():
    hs = tuple(coefficient_matrix(j) for j in range(3))
    if any(h.H != h for h in hs):
        raise ArithmeticError('the configured Weyl/link hopping is not Hermitian')
    occupied = s.diag(0, *([1]*59))
    checked = 0; electric_checks = 0; wrong = 0
    for generator in currents():
        for endpoint in (0, 1):
            link = link_generator(generator, endpoint)
            if link*occupied != occupied*link:
                raise ArithmeticError('the electric link energy violates native Gauss')
            electric_checks += 1
            gauss = gauss_coefficients(generator, endpoint)
            for h in hs:
                if (gauss*h-h*gauss).is_zero_matrix is not True:
                    raise ArithmeticError('a full-CAR Weyl/link Ward coefficient is nonzero')
                checked += 1
            if generator == currents()[-1]:
                broken = gauss_coefficients(generator, endpoint, orientation=-1)
                wrong += int((broken*hs[0]-hs[0]*broken).is_zero_matrix is not True)
    if wrong != 2:
        raise ArithmeticError('wrong link orientation did not fail at both endpoints')
    return {'matter_components': 15, 'spatial_spin_components': 2,
            'physical_CAR_modes_per_vertex': 30, 'link_dimension': 60,
            'edge_Fock_dimension': str(2**60), 'coefficient_array_dimension': 3600,
            'full_Fock_Ward_identities': checked,
            'electric_projector_Ward_identities': electric_checks,
            'wrong_orientation_endpoint_failures': wrong,
            'proof_method': 'independent normal-ordered CAR bilinear coefficients',
            'all_link_states_retained': True,
            'physical_Fock_space_truncated_to_one_particle': False}


def _dressed_action(state, axis):
    """Literal owning link Hamiltonian with its scalar edge coefficient replaced."""
    spin = s.I*PAULI[_axis(axis)]/2
    out = {}
    for (spin_in, *raw), amplitude in state.items():
        key = tuple(raw)
        moved = native.hamiltonian_action({key: amplitude}, 15, link_energy=0,
                                          blocks=_source()[1])
        for target, value in moved.items():
            for spin_out in range(2):
                coefficient = (spin[spin_out, spin_in] if key[0] == 1
                               else spin.H[spin_out, spin_in])
                native.add(out, (spin_out, *target), coefficient*value)
        if key[-2] >= 0:
            native.add(out, (spin_in, *key), amplitude)
    return native.normalized(out)


def dressed_reduction_certificate():
    checks = 0
    for indices in _source()[1]:
        d = len(indices)
        a, b = native.reference_states(15, indices)
        basis = [{(spin, *key): value for key, value in state.items()}
                 for state in (a, b) for spin in range(2)]
        for vector in basis:
            for generator in currents():
                for endpoint in (0, 1):
                    by_spin = [{key[1:]: value for key, value in vector.items() if key[0] == spin}
                               for spin in range(2)]
                    if any(native.gauss_action(generator, part, endpoint) for part in by_spin):
                        raise ArithmeticError('a dressed spatial-spin state violates native Gauss')
                    checks += 1
        for axis in range(3):
            spin = s.I*PAULI[axis]/2
            h = s.zeros(4); h[:2, 2:] = d*spin.H; h[2:, :2] = spin
            h[2:, 2:] = s.eye(2)
            if h*h-h != d*s.eye(4)/4:
                raise ArithmeticError('the finite dressed energy polynomial differs')
            for column, vector in enumerate(basis):
                expected = native.sum_states(*(native.scale(v, h[row, column])
                                               for row, v in enumerate(basis)))
                if _dressed_action(vector, axis) != expected:
                    raise ArithmeticError('the literal native link restriction differs')
    return {'source_blocks': 5, 'literal_native_Gauss_state_checks': checks,
            'literal_Weyl_link_Hamiltonian_columns': 5*3*4,
            'restricted_dimension_per_block': 4,
            'restricted_Gram_diagonal': ['1', '1', 'd', 'd'],
            'restricted_polynomial_at_unit_parameters': 'H^2-H=d I/4',
            'level_separation_squared_at_unit_parameters': '1+d',
            'finite_transport_frequency_identified_with_particle_mass': False}


def dressed_history(axis, *, duration=s.Rational(1, 8)):
    """All block, position and spin records for the configured two-site clock.

    The initial matter has one physical fermion and a static dual reference,
    with source block weight dim(Vr)/15 and spatial spin |up>. This is an
    explicit finite preparation, with no interacting vacuum assertion.
    """
    axis = _axis(axis)
    duration = native.exact_real(duration, nonnegative=True)
    spin = s.I*PAULI[axis]/2
    records = []; checks = 0
    for indices in _source()[1]:
        d = len(indices)
        h = s.zeros(4); h[:2, 2:] = d*spin.H; h[2:, :2] = spin
        h[2:, 2:] = s.eye(2)
        gram = s.diag(1, 1, d, d)
        if h.H*gram != gram*h:
            raise ArithmeticError('the dressed Weyl restriction lost its Gram metric')
        u = ((s.eye(4)+s.I*duration*h/2).inv()
             *(s.eye(4)-s.I*duration*h/2)).applyfunc(s.simplify)
        if (u.H*gram*u-gram).applyfunc(s.simplify) != s.zeros(4):
            raise ArithmeticError('the same-clock Cayley step is not complete')
        output = u*s.Matrix([1, 0, 0, 0])
        for position in (0, 1):
            for spin_label in (0, 1):
                i = 2*position+spin_label
                probability = s.simplify(s.Rational(d, 15)*gram[i, i]*abs(output[i])**2)
                records.append({'block_indices': list(indices), 'position': 1-position,
                                'spin': spin_label, 'original_probability': str(probability)})
        checks += 1
    records.append({'block_indices': [], 'position': None, 'spin': None,
                    'original_probability': '0', 'outcome': 'other_occupation'})
    if sum(s.Rational(row['original_probability']) for row in records) != 1:
        raise ArithmeticError('complete original gauge-sector probabilities do not sum to one')
    return {'axis': axis, 'duration': str(duration), 'hopping_coefficient': '1',
            'link_energy': '1', 'blocks_retained': checks, 'records': records,
            'one_fermion_record_count': 20, 'other_occupation_record_count': 1,
            'measurement_complete_on_full_Fock': True,
            'zero_probability_records_retained': True,
            'initial_block_measurement_renormalized': False,
            'clock': 'Cayley step of the same configured Weyl/link Hamiltonian',
            'static_dual_reference_is_a_declared_input': True,
            'native_pulse_execution_claimed': False,
            'interacting_ground_or_particle_mass_inferred': False}


def certificate():
    return {'schema': 'oph.standard_model_source_map.gauge_attachment.v1',
            'full_source_anomaly_tensor': anomaly_certificate(),
            'full_CAR_local_gauge_Ward': ward_certificate(),
            'literal_dressed_source_reduction': dressed_reduction_certificate(),
            'complete_dressed_controls': [dressed_history(j) for j in range(3)],
            'scope': {
                'same_fifteen_source_components_and_twelve_currents': True,
                'executed_family_copies': 1, 'family_count_selected': False,
                'spatial_Weyl_matrices_commute_with_internal_gauge_action': True,
                'quantum_links_physical_CAR_and_spatial_hopping_are_declared_architecture': True,
                'quantum_links_identified_with_unitary_background_links': False,
                'source_repair_selects_this_gauge_action': False,
                'interacting_relativistic_chiral_continuum_proved': False,
                'physical_Standard_Model_particles_identified': False,
                'empirical_mass_or_coupling_inputs': []}}
