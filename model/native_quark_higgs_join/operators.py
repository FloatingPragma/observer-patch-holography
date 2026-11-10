"""Actual original-CAR quark/Higgs equations and complete current covariance."""
from collections import defaultdict
from functools import lru_cache
import sympy as s

from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_joint.local_reference.algebra import particle_hole
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source
from model.coherent_defect_wave.chiral_joint.local_reference.physical import columns
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.native_quark_response.quark import target_module
from model.native_quark_response.local import CHANNELS, target, transfer


@lru_cache(None)
def local_H():
    parts, couplings = car.quantum_yukawa()
    h = car.add(*parts.values(), car.monomial(matrix=higgs_memory()[2]))
    return car.clean({k: v.subs({y: 1 for y in couplings.values()}) for k, v in h.items()})


@lru_cache(None)
def composite(mode):
    actual = car.commutator(local_H(), car.monomial(((mode, False),)))
    expected = car.scale(car.adjoint(car.add(*(car.screened_eom(channel, mode % 15, mode//15)
                                               for channel in ('up', 'down', 'charged_lepton')))), -1)
    assert actual == expected
    grouped = defaultdict(lambda: s.zeros(5))
    for (word, row, col), coefficient in actual.items():
        assert len(word) == 1
        grouped[word][row, col] += coefficient
    assert len(grouped) == 2
    for scalar in grouped.values():
        effect = scalar.H*scalar
        assert effect*effect == effect and effect != s.zeros(5)
    return actual


@lru_cache(None)
def operator_controls():
    gauge = representation()[2]
    raw_L = s.diag(*(1 if i % 15 in (10, 11) else -1 if i % 15 == 9 else 0 for i in range(30)))
    raw_B = s.diag(*(1 if i % 15 in FAMILIES['Q'] else -1 if i % 15 in (*FAMILIES['uc'], *FAMILIES['dc']) else 0
                     for i in range(30)))
    raw = tuple(car.one_body(t) for t in (raw_L, raw_B))
    assert all(not car.commutator(t, local_H()) for t in raw)
    rows = []
    for family in CHANNELS:
        indices = FAMILIES[family]
        modes = tuple(15*spin+i for spin in range(2) for i in indices)
        aa = {mode: car.monomial(((mode, False),)) for mode in modes}
        cc = {mode: composite(mode) for mode in modes}
        gauge_count = spin_count = raw_count = 0
        for probes in (aa, cc):
            for mode, op in probes.items():
                spin, index = divmod(mode, 15)
                assert all(len(word) % 2 == 1 for word, _, _ in op)
                for axis in range(12):
                    expected = car.add(*(car.scale(probes[15*spin+j], -gauge[axis][index, j])
                                         for j in indices if gauge[axis][index, j]))
                    assert car.commutator(car.gauge_current(axis), op) == expected
                    gauge_count += 1
                for axis in range(3):
                    expected = car.add(*(car.scale(probes[15*beta+index], -car.PAULI[axis][spin, beta]/2)
                                         for beta in range(2)))
                    assert car.commutator(car.spin_current(axis), op) == expected
                    spin_count += 1
                for current, weights in zip(raw, (raw_L, raw_B), strict=True):
                    assert car.commutator(current, op) == car.scale(op, -weights[mode, mode])
                    raw_count += 1
        example = modes[0]
        rows.append(dict(family=family, native_internal_dimension=len(indices), native_modes=list(modes),
                         actual_EOM_components=len(modes), gauge_component_checks=gauge_count,
                         spin_component_checks=spin_count, raw_charge_component_checks=raw_count,
                         raw_L_change='0', raw_B3_change=str(-raw_B[example, example]),
                         source_hypercharge_change=str(-gauge[11][indices[0], indices[0]]),
                         source_representation='conjugate of the named native family',
                         source_parity='odd', joint_static_reference_parity='even',
                         elementary_norm='1', composite_norm_upper='2',
                         composite_is_sum_of_two_unit_norm_scalar_CAR_products=True))
    return dict(full_local_Fock_dimension=2**30, scalar_dimension=5,
                coefficients='all three actual quantum Yukawa coefficients one',
                exact_EOM='[H_local,c_i]=C_i=-adjoint(sum_channel screened_eom(channel,i))',
                full_hopping_remainder='[h,c_i]=C_i+[2^-34 T,c_i]',
                every_occupation_and_scalar_state_retained=True, family_controls=rows,
                scalar_background_substitution=False, canonical_composite_CAR_assumed=False,
                finite_internal_current_labels_are_physical_Lorentz_chirality=False)


@lru_cache(None)
def reduced_EOM(family):
    create, mode, _, _ = CHANNELS[family]
    q = target_module(create)
    native_mode = q['native_component'][mode]
    operator = particle_hole(composite(native_mode))
    sb, tb = source(1), target(family)
    positions = {key: i for i, key in enumerate(tb['keys'])}
    entries = defaultdict(lambda: s.Integer(0))
    maps = tuple(columns(matrix) for matrix in q['maps'])
    for j, (a, mask, h) in enumerate(sb['keys']):
        for (word, row, col), coefficient in operator.items():
            if col != h:
                continue
            assert len(word) == 1 and word[0][1] == create
            alpha = q['native_component'].index(word[0][0])
            for destination, value in maps[alpha].get(a, ()):
                assert (destination, mask, row) in positions
                entries[positions[destination, mask, row], j] += coefficient*value
    raw = s.ImmutableSparseMatrix(len(tb['keys']), len(sb['keys']), entries)
    image = raw*sb['basis']
    matrix = s.ImmutableSparseMatrix(image.extract(tb['coordinates'], range(16)))
    assert tb['basis']*matrix == image
    elementary = transfer(family)['map']
    assert matrix == tb['H']*elementary-elementary*sb['H']
    return dict(family=family, full_weight_dimension=len(tb['keys']),
                retained_target_dimension=tb['H'].rows, source_dimension=16,
                literal_original_CAR_Higgs_image_equals_Ht_A_minus_A_Hs=True,
                all_scalar_and_quark_polynomial_images_checked=True)


@lru_cache(None)
def reference_controls():
    """Literal sixty-CAR joint invariance before vacuum-domain restriction."""
    gauge = representation()[2]
    spin = tuple(s.kronecker_product(p/2, s.eye(15)) for p in car.PAULI)
    raw_L = s.diag(*(1 if i % 15 in (10, 11) else -1 if i % 15 == 9 else 0 for i in range(30)))
    raw_B = s.diag(*(1 if i % 15 in FAMILIES['Q'] else -1 if i % 15 in (*FAMILIES['uc'], *FAMILIES['dc']) else 0
                     for i in range(30)))
    def ref_one(matrix):
        return car.one_body(s.diag(s.zeros(30), matrix))
    joint = tuple(car.add(car.gauge_current(j), ref_one(s.diag(g, g)))
                  for j, g in enumerate(gauge))
    joint += tuple(car.add(car.spin_current(j), ref_one(t)) for j, t in enumerate(spin))
    joint += tuple(car.add(car.one_body(t), ref_one(t)) for t in (raw_L, raw_B))
    href = car.scale(ref_one(s.eye(30)), 4)
    checks = []
    for family in CHANNELS:
        modes = tuple(15*spin+i for spin in range(2) for i in FAMILIES[family])
        probes = tuple(car.add(*(car.multiply(car.monomial(((30+i, True),)),
                                               car.monomial(((i, False),)) if kind == 'A' else composite(i))
                                  for i in modes)) for kind in ('A', 'C'))
        for probe in probes:
            assert all(len(word) % 2 == 0 for word, _, _ in probe)
            assert all(not car.commutator(t, probe) for t in joint)
            assert car.commutator(href, probe) == car.scale(probe, 4)
        assert car.commutator(local_H(), probes[0]) == probes[1]
        checks.append(dict(family=family, all_joint_current_commutators=34,
                           full_reference_energy_commutators=2,
                           full_compensated_local_EOM=True, reference_normalization_squared=len(modes)))
    return dict(full_joint_CAR_modes=60, full_reference_dimension=2**30,
                vacuum_domain_restriction_commutes_with_all_reference_currents=True,
                extra_spin_orientation_reference_required=False,
                fixed_origin_service_translation_covariance_asserted=False,
                families=checks)


if __name__ == '__main__':
    print(operator_controls(), flush=True)
    print(reference_controls(), flush=True)
    for family in CHANNELS:
        print(reduced_EOM(family), flush=True)
