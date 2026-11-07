"""Exact complete transport instruments on the literal sixty-state quantum link.

The matter representation and link register are declared source extensions.
The calculation selects neither those extensions nor an evolution parameter.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from ..particle_gauge_band import quantum_link_completion as link


NAMES = {-4: 'u^c', -3: 'L', 1: 'Q', 2: 'd^c', 6: 'e^c'}
SOURCE_FILES = (
    'model/gauge_defect_transport/__init__.py',
    'model/gauge_defect_transport/transport.py',
    'model/gauge_defect_transport/test_transport.py',
    'model/gauge_defect_transport/PROOF.md',
    'model/particle_gauge_band/quantum_link_completion.py',
    'model/native_quantum_gauss/algebra.py',
    'model/common_source/native_matter_attachment.py',
    'model/common_source/native_sm_current.py',
    'model/defect_field_bridge/colored.py',
    'model/defect_field_bridge/full_domain.py',
)


def _block(a, b, c, d):
    """Sparse two-by-two block matrix without dense intermediate matrices."""
    if a.rows != b.rows or c.rows != d.rows or a.cols != c.cols or b.cols != d.cols:
        raise ValueError('compatible block sizes required')
    entries = {}
    for matrix, row, column in ((a, 0, 0), (b, 0, a.cols),
                                (c, a.rows, 0), (d, a.rows, a.cols)):
        entries.update({(i+row, j+column): value
                        for (i, j), value in matrix.todok().items()})
    return s.ImmutableSparseMatrix(a.rows+c.rows, a.cols+b.cols, entries)


def _exact_matrix(value):
    if not isinstance(value, s.MatrixBase) or value.rows != value.cols or value.rows < 1:
        raise ValueError('a nonempty square exact matrix is required')
    if any(x.has(s.Float) or x.free_symbols or x.is_finite is not True
           for x in value.todok().values()):
        raise ValueError('finite exact coefficients without free symbols are required')
    return s.ImmutableSparseMatrix(value)


def complete_partial_transport(value):
    """Two-position instrument and reversible completion of a partial isometry.

    Position zero is the target, position one the source. Each direction has
    a transported outcome and a held outcome, with no conditioning on success.
    """
    b = _exact_matrix(value)
    p = link.clean(b.H*b)
    r = link.clean(b*b.H)
    if not link.zero(p*p-p) or not link.zero(r*r-r):
        raise ValueError('the transport must be a partial isometry')
    z = s.SparseMatrix.zeros(b.rows)
    identity = s.ImmutableSparseMatrix(s.eye(2*b.rows))
    t = _block(z, b, z, z)
    source = _block(z, z, z, p)
    target = _block(r, z, z, z)
    hold_forward = identity-source
    hold_reverse = identity-target
    reversible = t+t.H+identity-source-target
    return {
        'forward': t, 'forward_hold': hold_forward,
        'reverse': t.H, 'reverse_hold': hold_reverse,
        'source': source, 'target': target, 'reversible': reversible,
    }


def coherent_transport(transport, cosine, sine):
    """Declared rotation; its dimensionless angle is an input."""
    c, a = s.sympify(cosine), s.sympify(sine)
    for value in (c, a):
        if (value.has(s.Float) or value.free_symbols or value.is_real is not True
                or value.is_finite is not True):
            raise ValueError('finite exact real rotation coefficients required')
    if s.simplify(c*c+a*a-1) != 0:
        raise ValueError('cosine squared plus sine squared must equal one')
    t, p, r = (transport[key] for key in ('forward', 'source', 'target'))
    return s.ImmutableSparseMatrix(s.eye(t.rows))+(c-1)*(p+r)-s.I*a*(t+t.H)


@lru_cache(maxsize=1)
def _literal_data():
    masks, blocks, currents = link.exterior_data()
    tensors, _ = link.clebsch_gordan(masks, blocks, currents)
    indices, labels, left, right = link.link_representation(blocks, currents)
    return blocks, currents, tensors, indices, labels, left, right


@lru_cache(maxsize=10)
def literal_transport(charge, completed=False):
    if type(charge) is not int or charge not in NAMES or type(completed) is not bool:
        raise ValueError('a retained primitive charge and a boolean completion flag are required')
    blocks, currents, tensors, indices, labels, left, right = _literal_data()
    selected = blocks[charge]
    d = len(selected)
    raising, coordinate = link.coordinate_operators(charge, blocks, indices, tensors)
    a = link.combined_coordinate(coordinate if completed else raising, d)
    b = link.clean(a/s.sqrt(d))
    maps = complete_partial_transport(b)
    source_projector = link.clean(b.H*b)
    target_projector = link.clean(b*b.H)
    checks = 0
    for current, lx, ly in zip(currents, left, right, strict=True):
        small = current.extract(selected, selected)
        matter = s.SparseMatrix(s.kronecker_product(small, s.eye(60)))
        gx = s.SparseMatrix(s.kronecker_product(s.eye(d), lx))
        gy = s.SparseMatrix(s.kronecker_product(s.eye(d), ly))
        # At the source the matter sits at y; at the target it sits at x.
        for incoming, outgoing in ((gx, gx+matter), (gy+matter, gy)):
            if not link.zero(outgoing*b-b*incoming):
                raise ArithmeticError('literal endpoint Gauss intertwining failed')
            if not link.zero(source_projector*incoming-incoming*source_projector):
                raise ArithmeticError('held source outcome violates Gauss covariance')
            if not link.zero(target_projector*outgoing-outgoing*target_projector):
                raise ArithmeticError('held target outcome violates Gauss covariance')
            checks += 3
    return {'charge': charge, 'dimension': d, 'completed': completed,
            'coordinate': s.ImmutableSparseMatrix(a), 'normalized': s.ImmutableSparseMatrix(b),
            'source_projector': s.ImmutableSparseMatrix(source_projector),
            'target_projector': s.ImmutableSparseMatrix(target_projector),
            'maps': maps, 'endpoint_checks': checks}


def car_hop_sign(mask, origin, target, flavors):
    """Independent bitstring CAR sign for one occupied origin and empty target.

    The two entries are (site, internal label); all spectator occupations
    remain in the full occupation mask.
    """
    if type(mask) is not int or mask < 0 or type(flavors) is not int or flavors < 1:
        raise ValueError('an exact occupation mask and positive internal dimension are required')
    for entry in (origin, target):
        if (type(entry) is not tuple or len(entry) != 2
                or any(type(x) is not int for x in entry)
                or entry[0] < 0 or not 0 <= entry[1] < flavors):
            raise ValueError('each endpoint requires a nonnegative site and valid internal label')
    i, j = (site*flavors+internal for site, internal in (origin, target))
    if i == j or not mask >> i & 1 or mask >> j & 1:
        raise ValueError('an occupied origin and a distinct empty destination are required')
    sign = (-1)**(mask & ((1 << i)-1)).bit_count()
    removed = mask ^ (1 << i)
    sign *= (-1)**(removed & ((1 << j)-1)).bit_count()
    return removed | (1 << j), sign


def neutralized_control(charge, completed=False):
    """Same prepared gauge singlet and detector distinguish two rotation inputs.

    A declared dual reference at y neutralizes the transported charge. State
    amplitudes are matrices with the reference label on the column axis.
    """
    value = literal_transport(charge, completed)
    blocks, currents, _, _, _, left, right = _literal_data()
    d, maps = value['dimension'], value['maps']
    n = 60*d
    initial = s.SparseMatrix(2*n, d, {(n+60*b, b): 1/s.sqrt(d) for b in range(d)})
    transported = maps['forward']*initial
    if s.trace(initial.H*initial) != 1 or s.trace(transported.H*transported) != 1:
        raise ArithmeticError('neutralized charged preparation lost normalization')
    checks = 0
    z = s.SparseMatrix.zeros(n)
    for current, lx, ly in zip(currents, left, right, strict=True):
        small = current.extract(blocks[charge], blocks[charge])
        matter = s.SparseMatrix(s.kronecker_product(small, s.eye(60)))
        gx = s.SparseMatrix(s.kronecker_product(s.eye(d), lx))
        gy = s.SparseMatrix(s.kronecker_product(s.eye(d), ly))
        endpoint_x = _block(gx+matter, z, z, gx)
        endpoint_y = _block(gy, z, z, gy+matter)
        for state in (initial, transported):
            # The dual reference action at y is -small.T. Vectorization
            # turns it into right multiplication by -small.
            if not link.zero(endpoint_x*state) or not link.zero(endpoint_y*state-state*small):
                raise ArithmeticError('complete preparation is not an endpoint Gauss singlet')
            checks += 2
    probabilities = []
    for cosine, sine in ((s.Rational(3, 5), s.Rational(4, 5)),
                         (s.Rational(5, 13), s.Rational(12, 13))):
        output = coherent_transport(maps, cosine, sine)*initial
        probability = s.simplify(s.trace(output.H*maps['target']*output))
        if probability != sine*sine:
            raise ArithmeticError('gauge-invariant endpoint detector contrast failed')
        probabilities.append(str(probability))
    return {'primitive_charge': charge, 'reference_dimension': d,
            'full_matter_position_link_reference_dimension': 120*d*d,
            'endpoint_singlet_checks': checks,
            'detector': 'target-support projector, with its full complementary outcome',
            'same_input_and_detector_rotation_probabilities': probabilities,
            'preparation_and_reference_are_declared': True}


def certificate():
    rows = []
    for completed in (False, True):
        for charge, name in NAMES.items():
            value = literal_transport(charge, completed)
            maps, d = value['maps'], value['dimension']
            identity = s.SparseMatrix(s.eye(120*d))
            for direction in ('forward', 'reverse'):
                k, held = maps[direction], maps[direction+'_hold']
                if not link.zero(k.H*k+held.H*held-identity):
                    raise ArithmeticError('complete transport effects do not sum to identity')
            if not link.zero(maps['reversible']**2-identity):
                raise ArithmeticError('reversible transport completion failed')
            rotation = coherent_transport(maps, s.Rational(3, 5), s.Rational(4, 5))
            if not link.zero(rotation.H*rotation-identity):
                raise ArithmeticError('coherent transport rotation is not unitary')
            rows.append({'matter_block': name, 'primitive_charge': charge,
                         'coordinate_completion': completed, 'matter_dimension': d,
                         'one_position_dimension': 60*d, 'two_position_dimension': 120*d,
                         'transport_rank': int(s.trace(value['source_projector'])),
                         'held_outcome_rank': 120*d-int(s.trace(value['source_projector'])),
                         'exact_endpoint_intertwiner_and_hold_checks': value['endpoint_checks']})
    root = Path(__file__).resolve().parents[2]
    # Preserve the owning current/matter producer's geometry and input-fixture
    # bindings as well as this proof, executable and independent test file.
    bindings = link.native.matter.source_bindings()
    bindings.update({p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in SOURCE_FILES})
    return {
        'schema': 'oph.gauge_defect_transport.v1',
        'declared_matter_dimension': 15, 'literal_link_dimension': 60,
        'complete_directional_outcomes': ['transported', 'held'],
        'proof_domain': 'One particle on a two-endpoint edge with its entire literal link register; arbitrary untouched references. A hard-core spectator contributes the independently verified CAR sign.',
        'exact_identity': 'B=A/sqrt(d), P=B^dagger B, R=B B^dagger; P^2=P, R^2=R.',
        'positioned_transport': 'T=|x><y| tensor B; instrument {T,I-T^dagger T}. Reverse instrument {T^dagger,I-T T^dagger}.',
        'reversible_completion': 'W=T+T^dagger+I-T^dagger T-T T^dagger; W^dagger=W and W^2=I.',
        'coherent_family': 'U(c,s)=I+(c-1)(T^dagger T+T T^dagger)-i s(T+T^dagger), c^2+s^2=1.',
        'rows': rows,
        'neutralized_parameter_controls': [neutralized_control(charge) for charge in NAMES],
        'scope': {
            'source_selected_representation_or_link_extension': False,
            'source_selected_action_or_clock': False,
            'full_Fock_collision_law_proved': False,
            'equivalence_to_original_C_E_repair_law_claimed': False,
            'site_local_Gauss_covariant_complete_instrument': True,
            'all_reference_extensions_and_finite_history_compositions': True,
            'physical_chiral_continuum_or_particle_identity': False,
        },
        'source_sha256': dict(sorted(bindings.items())),
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = certificate()
    path = Path(__file__).with_name('results.json')
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    elif args.check:
        if json.loads(path.read_text()) != result:
            raise SystemExit('gauge defect transport receipt drift')
        print('gauge defect transport receipt OK')
    else:
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
