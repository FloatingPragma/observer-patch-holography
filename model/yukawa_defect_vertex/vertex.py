"""Complete local CAR/Higgs readouts, without a supplied family matrix.

The source is the existing 320-state unit-coupling lepton/Higgs law. This
module derives its finite spectral normalization and equation-of-motion
vertex; it does not identify that normalization with an on-shell coupling.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from ..coherent_defect_wave.chiral_higgs import model as source
from ..coherent_defect_wave.physical_attachment.chiral import representation

ROOT = Path(__file__).resolve().parents[2]
SOURCE_FILES = (
    'model/coherent_defect_wave/chiral_higgs/model.py',
    'model/coherent_defect_wave/physical_attachment/chiral.py',
    'model/common_source/chiral_candidate/action.py',
    'model/observer_particles/results/representation_identification_operators.json',
    'model/observer_particles/results/representation_identification.json',
    'model/observer_particles/results/unified_yukawa_operators.json',
    'model/observer_particles/results/unified_yukawa_action.json',
    'model/yukawa_defect_vertex/vertex.py',
    'model/yukawa_defect_vertex/PROOF.md',
)


def clean(matrix):
    matrix = s.ImmutableSparseMatrix(matrix)
    return s.ImmutableSparseMatrix(matrix.rows, matrix.cols,
        {key: value for key, entry in matrix.todok().items()
         if (value := s.simplify(entry)) != 0})


def zero(matrix):
    return not clean(matrix).todok()


def trace(matrix):
    return s.simplify(sum(value for (i, j), value in matrix.todok().items() if i == j))


@lru_cache(maxsize=16)
def eye(n):
    return s.ImmutableSparseMatrix(n, n, {(i, i): 1 for i in range(n)})


def kron(left, right):
    # The generic symbolic Kronecker implementation traverses zero entries.
    # Literal sparse products preserve exactly the same full operator.
    return s.ImmutableSparseMatrix(left.rows*right.rows, left.cols*right.cols,
        {(i*right.rows+k, j*right.cols+l): a*b
         for (i, j), a in left.todok().items()
         for (k, l), b in right.todok().items()})


def comm(a, b):
    return a*b-b*a


def trace_product(a, b):
    """Exact trace without forming a large intermediate matrix product."""
    return s.simplify(sum(value*b[j, i] for (i, j), value in a.todok().items()))


@lru_cache(maxsize=1)
def source_data():
    """Sparse replay of the literal existing source, on every occupation.

    There is no Y argument. The native current tensor, six actual CAR
    operators and five scalar levels fix every matrix entry.
    """
    scalar, scalar_currents, number, _ = source.higgs_memory()
    cs = tuple(source.annihilator(i) for i in range(6))
    native = representation()
    components = native[3]['charged_lepton']['components']
    internal = (9, 10, 11)
    pairs = []
    for component in components:
        component = component.extract(internal, internal)
        pair = s.zeros(64, 64)
        for i in range(3):
            for j in range(3):
                if component[i, j]:
                    pair += component[i, j]*cs[i]*cs[3+j]
        pairs.append(clean(pair))
    h = kron(eye(64), number)
    for pair, phi in zip(pairs, scalar, strict=True):
        h += kron(pair, phi.H)+kron(pair.H, phi)
    h = clean(h)
    if h != h.H:
        raise ArithmeticError('the literal full source is not Hermitian')
    currents = []
    for current, scalar_current in zip(native[2], scalar_currents, strict=True):
        current = current.extract(internal, internal)
        matter = s.zeros(64, 64)
        for spin in range(2):
            for i in range(3):
                for j in range(3):
                    if current[i, j]:
                        matter += current[i, j]*cs[3*spin+i].H*cs[3*spin+j]
        total = kron(clean(matter), eye(5))+kron(eye(64), scalar_current)
        if not zero(comm(total, h)):
            raise ArithmeticError('a native current does not preserve the full law')
        currents.append(clean(total))
    vacuum = list(range(0, 320, 5))
    occupied = [i for i in range(320) if i % 5]
    w = h.extract(occupied, vacuum)
    gram = clean(w.H*w)
    projectors = []
    for r in range(7):
        p = eye(64)
        for k in range(7):
            if k != r:
                p = p*(gram-k*eye(64))/(r-k)
        p = clean(p)
        if not zero(p*p-p) or not zero(gram*p-r*p):
            raise ArithmeticError('the exact complete source spectrum changed')
        projectors.append(p)
    if sum(projectors, s.zeros(64)) != eye(64):
        raise ArithmeticError('the complete source eigenspaces do not resolve identity')
    lift = s.MutableSparseMatrix(320, 64, {})
    for j, i in enumerate(vacuum):
        lift[i, j] = 1
    for (i, j), value in w.todok().items():
        lift[occupied[i], j] = -value/3
    ground = clean(s.Rational(3, 5)*lift*projectors[6]*lift.H)
    rho = ground/2
    if (trace(ground) != 2 or not zero(ground*ground-ground)
            or not zero((h+2*eye(320))*ground)
            or any(not zero(q*ground) for q in currents)):
        raise ArithmeticError('the actual rank-two singlet ground was not replayed')
    return {'H': h, 'rho': rho, 'ground': ground, 'currents': tuple(currents),
            'c': tuple(kron(c, eye(5)) for c in cs),
            'phi': tuple(kron(eye(64), p) for p in scalar),
            'pairs': tuple(kron(p, eye(5)) for p in pairs),
            'small_pairs': tuple(pairs), 'small_scalar': scalar,
            'W': w, 'spectral_projectors64': tuple(projectors),
            'vacuum_indices': vacuum, 'occupied_indices': occupied}


@lru_cache(maxsize=1)
def spectral_trace_forms():
    """Rational trace forms for the complete fourteen spectral projectors."""
    data = source_data()
    vacuum, occupied = data['vacuum_indices'], data['occupied_indices']
    w = data['W']

    def embed(matrix, rows, columns):
        return s.ImmutableSparseMatrix(320, 320,
            {(rows[i], columns[j]): x for (i, j), x in matrix.todok().items()})

    result = []
    for p in data['spectral_projectors64']:
        wp = w*p
        v = embed(p, vacuum, vacuum)
        x = embed(wp, occupied, vacuum)
        x = x+x.H
        z = embed(wp*w.H, occupied, occupied)
        result.append((v, x, z, int(trace(p))))
    scalar_occupied = s.ImmutableSparseMatrix(320, 320, {(i, i): 1 for i in occupied})
    return tuple(result), scalar_occupied


def spectral_weights(density):
    """All fourteen spectral outcomes, including zero probabilities.

    The argument may be an off-diagonal matrix unit. This gives the whole
    two-point matrix, including coherent reference entries, by linearity.
    """
    forms, occupied = spectral_trace_forms()
    values = [(trace_product(density, a), trace_product(density, b),
               trace_product(density, c), rank) for a, b, c, rank in forms]
    result = [(s.S.Zero, values[0][0], values[0][3])]
    dark = trace_product(density, occupied)-sum(values[r][2]/r for r in range(1, 7))
    result.append((s.S.One, s.simplify(dark), 198))
    for r in range(1, 7):
        v, x, z, rank = values[r]
        for sign in (-1, 1):
            energy = (1+sign*s.sqrt(1+4*r))/2
            d = energy-1
            value = s.simplify((v+x/d+z/d**2)/(1+r/d**2))
            result.append((energy, value, rank))
    if s.simplify(sum(row[1] for row in result)-trace(density)) != 0:
        raise ArithmeticError('a complete energy record lost original weight')
    return tuple(result)


@lru_cache(maxsize=1)
def two_point_data():
    """Exact matrix-valued addition/removal spectra of all active fields."""
    data = source_data()
    rho, cs, phis = data['rho'], data['c'], data['phi']
    sectors = {}
    for name, operators in (('fermion', cs), ('scalar', phis)):
        n = len(operators)
        rows = [{'energy': energy, 'gap': s.simplify(energy+2), 'rank': rank,
                 'addition': s.MutableSparseMatrix(n, n, {}),
                 'removal': s.MutableSparseMatrix(n, n, {})}
                for energy, _, rank in spectral_weights(rho)]
        for i in range(n):
            for j in range(n):
                # C+_ij=<O_i P_E O_j^dagger>, C-_ij=<O_j^dagger P_E O_i>.
                plus = spectral_weights(operators[j].H*rho*operators[i])
                minus = spectral_weights(operators[i]*rho*operators[j].H)
                for row, a, b in zip(rows, plus, minus, strict=True):
                    row['addition'][i, j], row['removal'][i, j] = a[1], b[1]
        for row in rows:
            row['addition'], row['removal'] = clean(row['addition']), clean(row['removal'])
        sectors[name] = rows
    fermion_sum = sum((r['addition']+r['removal'] for r in sectors['fermion']), s.zeros(6))
    scalar_zero = sum((r['addition']-r['removal'] for r in sectors['scalar']), s.zeros(2))
    scalar_first = clean(sum((r['gap']*(r['addition']+r['removal'])
                             for r in sectors['scalar']), s.zeros(2)))
    if not zero(fermion_sum-eye(6)) or not zero(scalar_zero):
        raise ArithmeticError('the retarded two-point normalization changed')
    direct_first = clean(s.Matrix(2, 2, lambda i, j:
        trace(rho*comm(phis[i], comm(data['H'], phis[j].H)))))
    if not zero(scalar_first-direct_first):
        raise ArithmeticError('spectral records and the double commutator disagree')
    if scalar_first != scalar_first[0, 0]*eye(2) or scalar_first[0, 0].is_positive is not True:
        raise ArithmeticError('positive complete doublet normalization required')
    sectors['fermion_CAR_sum'] = clean(fermion_sum)
    sectors['scalar_commutator_sum'] = clean(scalar_zero)
    sectors['scalar_first_moment'] = scalar_first
    return sectors


@lru_cache(maxsize=1)
def vertex_data():
    """Extract the full anomalous equation-of-motion vertex from three points.

    K_ij={ [c_i,H],c_j}=i d/dt {c_i(t),c_j}|0. Both CAR signs and
    scalar orientation are checked as 320-state operator identities.
    """
    data = source_data()
    rho, h, cs, phis = data['rho'], data['H'], data['c'], data['phi']
    covariance = clean(s.Matrix(2, 2, lambda a, b: trace(rho*phis[a]*phis[b].H)))
    if covariance.det().is_positive is not True:
        raise ArithmeticError('the measured scalar Gram is singular')
    inverse = covariance.inv()
    scalar_moment = two_point_data()['scalar_first_moment'][0, 0]
    rows = []
    for i in range(6):
        ci_h = comm(cs[i], h)
        for j in range(6):
            insertion = clean(ci_h*cs[j]+cs[j]*ci_h)
            three_point = s.Matrix([[trace(rho*insertion*phi.H) for phi in phis]])
            coefficients = clean(three_point*inverse)
            reconstructed = sum((coefficients[0, k]*phis[k] for k in range(2)), s.zeros(320))
            if not zero(insertion-reconstructed):
                raise ArithmeticError('the full equation-of-motion vertex did not reconstruct')
            if not zero(insertion+(comm(cs[j], h)*cs[i]+cs[i]*comm(cs[j], h))):
                raise ArithmeticError('fermion antisymmetry failed')
            if any(coefficients):
                rows.append({'left': i, 'right': j, 'mixed_three_point': three_point,
                             'source_coefficients': coefficients,
                             'spectrally_normalized_coefficients': clean(s.sqrt(scalar_moment)*coefficients)})
    return {'scalar_covariance': covariance, 'scalar_first_moment': scalar_moment,
            'normalized_unit_coefficient': s.sqrt(scalar_moment),
            'nonzero_components': rows, 'all_fermion_component_pairs_checked': 36}


def complete_probe_records():
    """Coherent charged probes, every energy and scalar no-operation outcome."""
    data, spectrum = source_data(), two_point_data()
    rows = []
    for name, divisor in (('fermion', 6), ('scalar', 4)):
        for branch in ('addition', 'removal'):
            for row in spectrum[name]:
                probability = s.simplify(trace(row[branch])/divisor)
                if probability.is_nonnegative is not True:
                    raise ArithmeticError('a spectral probability is negative')
                rows.append({'probe': name, 'branch': branch, 'energy': row['energy'],
                             'gap': row['gap'], 'probability': probability})
    # This complement is an actual instrument outcome, not dropped weight.
    no_op = kron(eye(64), s.diag(0, s.sqrt(3)/2, s.sqrt(3)/2, s.sqrt(3)/2, s.sqrt(3)/2))
    scalar_effect = sum((p.H*p+p*p.H for p in data['phi']), s.zeros(320))/4+no_op.H*no_op
    if not zero(scalar_effect-eye(320)):
        raise ArithmeticError('the scalar probe is incomplete')
    for energy, probability, _ in spectral_weights(no_op*data['rho']*no_op.H):
        rows.append({'probe': 'scalar', 'branch': 'no_operation', 'energy': energy,
                     'gap': s.simplify(energy+2), 'probability': probability})
    for name in ('fermion', 'scalar'):
        if s.simplify(sum(r['probability'] for r in rows if r['probe'] == name)-1) != 0:
            raise ArithmeticError('a complete probe failed to normalize')
    preparation = [{'energy': energy, 'probability': s.Rational(rank, 320),
                    'ground': energy == -2}
                   for energy, _, rank in spectral_weights(data['rho'])]
    complete = []
    for name in ('fermion', 'scalar'):
        complete.extend({'probe': name, 'stage': 'preparation_other_energy', **r}
                        for r in preparation if not r['ground'])
        complete.extend({**r, 'stage': 'ground_then_probe',
                         'probability': r['probability']/160} for r in rows if r['probe'] == name)
        if s.simplify(sum(r['probability'] for r in complete if r['probe'] == name)-1) != 0:
            raise ArithmeticError('the original preparation probability was discarded')
    return {'ground_conditional_records': rows, 'all_preparation_records': preparation,
            'unconditional_preparation_and_probe_records': complete,
            'fermion_reference': 'coherent direct sum of all six dual addition and six removal modes',
            'scalar_reference': 'coherent dual doublet in each orientation, plus no-operation branch'}


@lru_cache(maxsize=1)
def scalar_intervention():
    """A charged Higgs pulse and gauge-invariant pair pointer on its reference.

    The reference stays quantum. The two pointer effects retain both outcomes;
    their sums are identity, and all failed ground preparations are kept.
    """
    data = source_data()
    scalar = data['small_scalar']
    reference_raise = tuple(s.ImmutableSparseMatrix(3, 3, {(k+1, 0): 1}) for k in range(2))
    x_small = sum((kron(phi.H, r)+kron(phi, r.H)
                   for phi, r in zip(scalar, reference_raise, strict=True)), s.zeros(15))
    x_small = clean(x_small)
    x = kron(eye(64), x_small)
    b = sum((kron(pair, r.H) for pair, r in zip(data['pairs'], reference_raise, strict=True)),
            s.ImmutableSparseMatrix(960, 960, {}))
    # ||P_k||<=2 and ||R_k||=1 imply ||i(B-B^dagger)||<=8.
    pointer = clean(s.I*(b-b.H)/8)
    reference_vacuum = s.ImmutableSparseMatrix(3, 3, {(0, 0): 1})
    rho = kron(data['rho'], reference_vacuum)
    ward = 0
    scalar_currents = source.higgs_memory()[1]
    for q, qh in zip(data['currents'], scalar_currents, strict=True):
        t = qh[1:3, 1:3]
        qr = s.diag(s.zeros(1), -t.T)
        total = kron(q, eye(3))+kron(eye(320), qr)
        if not zero(comm(total, x)) or not zero(comm(total, pointer)):
            raise ArithmeticError('the scalar insertion or pair pointer changed gauge charge')
        ward += 2
    derivative = s.simplify(s.I*trace(rho*comm(x, pointer)))
    if derivative == 0:
        raise ArithmeticError('the actual charged-scalar response vanished')
    # Exact finite pulse from X^2 eigenvalues 0,1,2, preserving the null space.
    square = x_small*x_small
    p1, p2 = clean(2*square-square*square), clean((square*square-square)/2)
    if not zero(p1*p1-p1) or not zero(p2*p2-p2) or not zero(p1*p2):
        raise ArithmeticError('scalar pulse projectors failed')
    theta = s.Symbol('theta', real=True)
    u_small = (eye(15)+(s.cos(theta)-1)*p1+(s.cos(s.sqrt(2)*theta)-1)*p2
               -s.I*s.sin(theta)*x_small*p1-s.I*s.sin(s.sqrt(2)*theta)*x_small*p2/s.sqrt(2))
    # Only columns reached from the blank reference are needed for a trace;
    # these are exact columns of the full unitary, not a projected generator.
    source_columns = [3*j for j in range(5)]
    columns = u_small[:, source_columns]
    column_gram = clean(columns.H*columns)
    if not zero(column_gram-eye(5)):
        raise ArithmeticError('the scalar pulse lost a reference branch')
    operators = [kron(eye(64), columns[r::3, :]) for r in range(3)]
    # B maps reference k+1 to zero, so only these coherent blocks contribute.
    b_mean = sum(trace(pair*operators[k+1]*data['rho']*operators[0].H)
                 for k, pair in enumerate(data['pairs']))
    mean = s.simplify(s.I*(b_mean-s.conjugate(b_mean))/8)
    if s.simplify(s.diff(mean, theta).subs(theta, 0)-derivative) != 0:
        raise ArithmeticError('finite charged pulse and three-point response disagree')
    at = s.pi/(2*s.sqrt(2))
    expectation = s.simplify(mean.subs(theta, at))
    probabilities = {'minus': s.simplify((1-expectation)/2),
                     'plus': s.simplify((1+expectation)/2)}
    if sum(probabilities.values()) != 1:
        raise ArithmeticError('the complete pair pointer lost probability')
    all_records = [{'stage': 'preparation_other_energy', 'energy': energy,
                    'probability': s.Rational(rank, 320)}
                   for energy, _, rank in spectral_weights(data['rho']) if energy != -2]
    all_records.extend({'stage': 'ground_then_scalar_pulse', 'pointer': label,
                        'probability': probability/160}
                       for label, probability in probabilities.items())
    if s.simplify(sum(row['probability'] for row in all_records)-1) != 0:
        raise ArithmeticError('scalar intervention discarded failed preparations')
    return {'scalar_reference_dimension': 3, 'source_and_scalar_reference_dimension': 960,
            'total_Gauss_operator_identities': ward,
            'intervention': 'X=sum_k phi_k^dagger tensor |k dual><0| + h.c.',
            'pair_pointer': 'A=i(B-B^dagger)/8, B=sum_k P_k tensor |0><k dual|',
            'pointer_norm_bound': '1', 'pointer_effects': '(I +/- A)/2',
            'derivative_of_pair_expectation_at_zero_pulse': derivative,
            'full_finite_pulse_expectation': mean,
            'pulse_angle': at, 'pair_expectation': expectation,
            'conditional_pointer_probabilities': probabilities,
            'complete_preparation_and_pointer_records': all_records,
            'quantum_reference_coherences_retained': True,
            'scalar_cost_or_interaction_coefficient_differentiated': False}


def serial(value):
    if isinstance(value, s.MatrixBase):
        return [[str(x) for x in row] for row in value.tolist()]
    if isinstance(value, s.Basic):
        return str(value)
    if isinstance(value, dict):
        return {k: serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    return value


def certificate():
    data = source_data()
    return serial({'schema': 'oph.yukawa_defect_vertex.v1',
        'source_files_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_FILES},
        'construction': {'source': 'existing 320-state unit lepton-Higgs interaction',
            'dimension': data['H'].rows, 'CAR_modes': 6, 'scalar_levels': 5,
            'primitive_interaction_coefficient': '1', 'scalar_number_energy': '1',
            'supplied_family_matrix': False, 'measured_mass_inputs': [],
            'ground_energy': '-2', 'ground_rank': 2, 'ground_preparation_probability': '1/160'},
        'two_points': two_point_data(), 'vertex': vertex_data(),
        'complete_records': complete_probe_records(), 'actual_scalar_intervention': scalar_intervention(),
        'scope': {'actual_charged_Higgs_insertion': True,
            'complete_source_and_quantum_references_retained': True,
            'normalization_derived_from_full_addition_and_removal_records': True,
            'scalar_high_frequency_z_inverse_term_exactly_zero': True,
            'scalar_high_frequency_z_inverse_squared_normalization': True,
            'equation_of_motion_contact_vertex_not_full_1PI_vertex': True,
            'on_shell_pole_residue_or_relativistic_Higgs_kinetic_normalization': False,
            'primitive_unit_interaction_forced_by_native_idle_repair': False,
            'three_family_hierarchy_or_physical_particle_masses': False}})


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2, sort_keys=True))
