"""An exact quark particle-hole symmetry of the complete joined law."""
from functools import lru_cache
import sympy as s

from model.defect_mass_join.source import (
    car, compact_points, local_interactions, one_cell_polynomial,
    representation, spatial_operators,
)

PAIRS = ((0, 12), (1, 13), (4, 14), (2, 3), (5, 6), (7, 8))
MAP = {a: (b, 1) for a, b in PAIRS} | {b: (a, -1) for a, b in PAIRS}


def transform(operator):
    """Local D=1 map; each spatial word also acquires its literal D factors."""
    terms = []
    for (word, row, col), value in operator.items():
        mapped, phase = [], 1
        for mode, create in word:
            spin, internal = divmod(mode, 15)
            if internal in MAP:
                following, sign = MAP[internal]
                phase *= sign*(1 if spin == 0 else -1)
                mapped.append((15*(1-spin)+following, not create))
            else:
                mapped.append((mode, create))
        for ordered, sign in car.normal_order(tuple(mapped)):
            terms.append({(ordered, row, col): phase*sign*value})
    return car.add(*terms)


@lru_cache(maxsize=1)
def certificate():
    parts = local_interactions()
    assert transform(parts['up']) == parts['down']
    assert transform(parts['down']) == parts['up']
    assert transform(parts['charged_lepton']) == parts['charged_lepton']
    assert transform(one_cell_polynomial()) == one_cell_polynomial()
    baryon_weights = [s.Rational(-1, 3) if i in (0, 1, 4, 12, 13, 14)
                      else s.Rational(1, 3) if i in (2, 3, 5, 6, 7, 8) else 0
                      for i in range(15)]
    baryon = car.one_body(s.kronecker_product(s.eye(2), s.diag(*baryon_weights)))
    hypercharge = car.gauge_current(11)
    assert transform(baryon) == car.scale(baryon, -1)
    assert transform(hypercharge) == car.add(hypercharge, car.scale(baryon, -1))
    signs = []
    for index in range(11):
        current = car.gauge_current(index)
        changed = transform(current)
        options = [sign for sign in (1, -1) if changed == car.scale(current, sign)]
        assert len(options) == 1
        signs.append(options[0])
    assert not car.commutator(baryon, one_cell_polynomial())
    spatial = spatial_operators(compact_points())
    entries = {}
    for x, point in enumerate(spatial['points']):
        phase = 1 if sum(point) % 2 == 0 else -1
        entries[2*x, 2*x+1], entries[2*x+1, 2*x] = phase, -phase
    c = s.ImmutableSparseMatrix(686, 686, entries)
    assert c.T*c == s.eye(686)
    assert spatial['h']*c+c*spatial['h'].conjugate() == s.zeros(686)
    return {'full_spatial_h_particle_hole_identity_exact': True,
        'all_three_local_quantum_channels_mapped_exactly': True,
        'scalar_Omega_unchanged': True,
        'native_quark_internal_signed_pairs': [list(pair) for pair in PAIRS],
        'CAR_map': 'c_quark(x,s,i) -> D_x epsilon_st S_ij c_quark(x,t,j)^dagger; leptons and scalar unchanged',
        'native_current_signs_first_eleven': signs,
        'hypercharge_map': 'Y -> Y-B', 'baryon_map': 'B -> -B',
        'global_native_singlet_and_B_zero_domain_preserved': True,
        'whole_global_native_singlet_domain_preserved_without_B_zero': False,
        'complete_spatial_physical_ground_B_zero_proved': False,
        'local_physical_ground_B_zero': 'the full local_reference theorem has B=0',
        'source_spectral_consequence': 'on an invariant B=0 reference, up removal and down addition are related by the exact spatial particle-hole frame; no spatial-ground symmetry assumption is made',
        'positive_Y_squared_energy_breaks_this_symmetry_on_charged_targets': True}
