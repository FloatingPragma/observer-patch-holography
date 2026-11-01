"""Exact current/basis join and a complete parity-controlled native instrument.

The parity-correlated code below is a new allowed memory design. It does
not equate global Fock parity with the exchange of two field particles.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from ..defect_identity import identity
from ..defect_protection import protection
from ..common_source.chiral_candidate import action

ROOT = Path(__file__).resolve().parents[2]


def _zero(matrix):
    return all(s.simplify(value) == 0 for value in s.SparseMatrix(matrix).todok().values())


def _q5(value):
    value = s.expand(value)
    b = value.coeff(s.sqrt(5))
    a = s.expand(value-b*s.sqrt(5))
    if not a.is_Rational or not b.is_Rational:
        raise ValueError('exact Q(sqrt(5)) source coefficient required')
    return identity.exterior.F5(str(a), str(b))


def _source_exterior(one_body):
    """Independent existing FockSpace.dgamma in B's actual subset order."""
    values = []
    for row in one_body.tolist():
        values.append([identity.exterior.C5(*map(_q5, s.expand(z).as_real_imag())) for z in row])
    fock = identity.exterior.FockSpace()
    full = fock.dgamma(values)
    indices = [fock.index[subset] for subset in identity.matter_basis()]
    return s.ImmutableSparseMatrix([[identity._scalar(full[i][j]) for j in indices] for i in indices])


@lru_cache(maxsize=1)
def representation():
    """Expose the actual retained operators, not a replacement label table."""
    subsets, g5, currents, operators, channels = action._representation()
    origin = identity.matter_basis()
    if len(set(subsets)) != 15 or set(subsets) != set(origin):
        raise ArithmeticError('the retained field is not the same exterior module')
    P = s.ImmutableSparseMatrix(15, 15, {(i, origin.index(row)): 1 for i, row in enumerate(subsets)})
    if P.T*P != s.eye(15):
        raise ArithmeticError('matter basis join is not unitary')
    source_currents = tuple(_source_exterior(g) for g in g5)
    for original, current in zip(source_currents, currents, strict=True):
        if not _zero(current*P-P*original):
            raise ArithmeticError('same-current exterior intertwiner failed')
    P90 = s.ImmutableSparseMatrix(s.kronecker_product(s.eye(2), P, s.eye(3)))
    W = s.ImmutableSparseMatrix.vstack(P90, -P90)
    if W.H*W != 2*s.eye(90):
        raise ArithmeticError('negative one-particle encoding lost norm')
    currents90 = tuple(s.ImmutableSparseMatrix(s.kronecker_product(s.eye(2), g, s.eye(3))) for g in currents)
    source90 = tuple(s.ImmutableSparseMatrix(s.kronecker_product(s.eye(2), g, s.eye(3))) for g in source_currents)
    for current, original in zip(currents90, source90, strict=True):
        if not _zero(s.kronecker_product(s.eye(2), current)*W-W*original):
            raise ArithmeticError('spin/family current join failed')
    # Bilinear mass kernels transform by P^T K P, not P^* K P.
    transported_channels = {}
    for name, channel in channels.items():
        transformed = tuple(s.ImmutableSparseMatrix(P.T*K*P) for K in channel['components'])
        if any(not _zero(P*K*P.T-target) for K, target in zip(transformed, channel['components'], strict=True)):
            raise ArithmeticError('retained Higgs bilinear join failed')
        transported_channels[name] = transformed
    return {'subsets': subsets, 'permutation15': P, 'permutation90': P90,
            'negative_one_particle_embedding': W,
            'source_one_body5': g5, 'source_currents15': source_currents,
            'currents15': currents, 'currents90': currents90,
            'identification15': operators, 'Higgs_channels': channels,
            'Higgs_channels_in_B_basis': transported_channels}


def parity_code(modes):
    """Small exact full-Fock code W=sqrt(2)E=[I;P]; no state truncation."""
    if type(modes) is not int or not 1 <= modes <= 4:
        raise ValueError('explicit full-Fock oracle is bounded to one through four modes')
    d = 1 << modes
    parity = s.diag(*[(-1)**mask.bit_count() for mask in range(d)])
    W = s.ImmutableSparseMatrix.vstack(s.ImmutableSparseMatrix(s.eye(d)), s.ImmutableSparseMatrix(parity))
    annihilators = []
    for mode in range(modes):
        entries = {}
        for mask in range(d):
            if mask >> mode & 1:
                entries[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
        annihilators.append(s.ImmutableSparseMatrix(d, d, entries))
    return W, s.ImmutableSparseMatrix(parity), tuple(annihilators)


def small_complete_instrument(law, seam=(0, 5), modes=2):
    """Actual native branches coupled to an exact even field control.

    W and all branch matrices omit their stated factors sqrt(2). The field
    control is (-1)^n0 on native coin one, identity on coin zero. It is a
    complete unitary field operation, not a second-quantized native Kraus map.
    """
    W, parity, annihilators = parity_code(modes)
    d = parity.rows
    code = s.ImmutableSparseMatrix(132*d, 66*d, {
        ((r*66+x)*d+f, x*d+f): (parity[f, f] if r else 1)
        for r in (0, 1) for x in range(66) for f in range(d)})
    number0 = annihilators[0].H*annihilators[0]
    field_phase = s.eye(d)-2*number0
    rows = [row for row in protection.branches(law) if row.seam == tuple(seam)]
    if not rows:
        raise ValueError('an actual native seam and C/E law are required')
    total_full, total_encoded = s.zeros(132*d, 1), s.zeros(66*d, 1)
    outputs = []
    for row in rows:
        U = field_phase if row.coin else s.eye(d)
        if U.H*U != s.eye(d) or U*parity != parity*U:
            raise ArithmeticError('field control is not an even unitary')
        full_entries, encoded_entries = {}, {}
        for i, j in protection.lifted_pairs(row):
            for f in range(d):
                full_entries[j*d+f, i*d+f] = U[f, f]
                total_full[i*d+f] += 1
        for i, j, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
            for f in range(d):
                encoded_entries[j*d+f, i*d+f] = U[f, f]*(parity[f, f] if sign == -1 else 1)
                total_encoded[i*d+f] += 1
        full = s.ImmutableSparseMatrix(132*d, 132*d, full_entries)
        encoded = s.ImmutableSparseMatrix(66*d, 66*d, encoded_entries)
        if full*code != code*encoded:
            raise ArithmeticError('complete native-field branch fails the code intertwiner')
        outputs.append((row, full, encoded))
    # Every map is a partial injection times a diagonal unitary. Its effect
    # is therefore exactly its domain projector, and the Kraus square is1/2.
    if set(total_full) != {2} or set(total_encoded) != {2}:
        raise ArithmeticError('joint native-field instrument is not trace preserving')
    return {'code': code, 'parity': parity, 'branches': outputs,
            'kraus_normalization': '1/sqrt(2)', 'code_normalization': '1/sqrt(2)',
            'full_dimension': 132*d, 'encoded_dimension': 66*d}


def statistics_certificate():
    W, parity, fields = parity_code(3)
    d = parity.rows
    X, Z = s.Matrix([[0, 1], [1, 0]]), s.diag(1, -1)
    if W.H*W != 2*s.eye(d) or s.kronecker_product(X, s.eye(d))*W != W*parity:
        raise ArithmeticError('parity-correlated isometry identity failed')
    for field in fields:
        if parity*field != -field*parity or s.kronecker_product(Z, field)*W != W*field:
            raise ArithmeticError('dressed odd field fails the same-code CAR identity')
    # Exact Fock lift of exchanging two mode labels: |11> acquires det=-1.
    swap = s.Matrix([[1,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,-1]])
    _, parity2, _ = parity_code(2)
    if swap[3, 3] != -1 or parity2[3, 3] != 1:
        raise ArithmeticError('exchange-versus-parity obstruction disappeared')
    return {'code': 'E=[I;P]/sqrt(2), P=(-1)^N; Q E=E P',
            'same_code_odd_field': '(Z_orientation tensor c_a) E = E c_a',
            'same_code_even_operator': '(I_orientation tensor O_even) E = E O_even',
            'native_branch_compression': 'sum_x |f_h(x)><x| tensor P^p_h(x), times unchanged1/sqrt(2)',
            'one_particle_restriction': 'the fixed negative code of defect_identity, tensored with Spin2 and family3',
            'vacuum': 'even parity and positive native Q in this new correlated code',
            'ninety_mode_particle_number_sectors': [{'N': n, 'parity': (-1)**n,
                                                     'dimension': str(s.binomial(90, n))} for n in range(91)],
            'two_occupied_modes_control': {'native_odd_loop_under_this_code': 1, 'actual_fermionic_mode_exchange': -1},
            'native_loop_identified_as_field_particle_exchange': False,
            'fixed_chi_minus_equals_global_Fock_parity': False,
            'warning': 'This is a designed global-parity coupling. The map from native two-token exchange to exchange of two field particles remains unconstructed.'}


def certificate():
    rep = representation()
    instruments = {law: small_complete_instrument(law) for law in ('C', 'E')}
    files = ['model/defect_field_bridge/bridge.py', 'model/defect_identity/identity.py',
             'model/defect_protection/protection.py',
             'model/common_source/chiral_candidate/action.py',
             'code/a5_closure/super_tannakian_matter_lift_certificate.py',
             'code/a5_closure/port_current_inner_certificate.py', *action.SOURCE_FILES]
    return {
        'schema': 'oph.defect_field_bridge.representation.v1',
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files},
        'matter_basis_permutation': [list(rep['permutation15'].row(i)).index(1) for i in range(15)],
        'same_source_current_intertwiners_checked': 12,
        'same_Higgs_channels_checked': sorted(rep['Higgs_channels']),
        'chosen_spin_dimension': 2, 'chosen_family_dimension': 3, 'one_particle_dimension': 90,
        'negative_one_particle_embedding_dimension': [180, 90],
        'statistics': statistics_certificate(),
        'complete_small_instrument': {law: {'seam': [0, 5], 'field_modes': 2,
             'raw_slots': len(packet['branches']), 'zero_slots': sum(not row.inputs for row, _, _ in packet['branches']),
             'abort_slots': sum(row.abort for row, _, _ in packet['branches']),
             'full_dimension': packet['full_dimension'], 'encoded_dimension': packet['encoded_dimension'],
             'all_branch_intertwiners_exact': True, 'trace_preserving_on_full_and_encoded_carriers': True}
             for law, packet in instruments.items()},
        'design': 'Same retained chiral matter representation, chosen Spin2 and family3; conditional even field unitaries on complete native-record branches.',
        'general_complete_instrument_proof': 'For unitary U_h, (K_h tensor U_h)^*(K_h tensor U_h)=K_h^*K_h tensor I. Summation retains every native slot and equals I. Even U_h preserve E.',
        'reference_and_history': 'Exact branch equality tensors with arbitrary references and composes under finite record-adaptive schedules; no successful-record renormalization.',
        'not_claimed': ['physical pair exchange from global parity', 'physical Lorentz or family identification from tensor dimensions',
                        'successful masses from the old trial Yukawa law', 'native Kraus second quantization is trace preserving'],
    }


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    payload = certificate()
    path = Path(__file__).with_name('representation_results.json')
    if args.check:
        if json.loads(path.read_text()) != payload:
            raise SystemExit('defect field representation receipt drift')
        print('defect field representation receipt OK')
    elif args.write:
        path.write_text(json.dumps(payload, indent=2)+'\n')
        print(f'wrote {path}')
    else:
        print(json.dumps(payload, indent=2))
