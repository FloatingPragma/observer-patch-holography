"""Microscopic particle-hole symmetries of the complete compact source."""
from functools import lru_cache
import sympy as s

from model.defect_mass_join.source import (
    car, compact_points, local_interactions, one_cell_polynomial, spatial_operators,
)
from model.quark_gauge_dressing.symmetry import certificate as quark_particle_hole
from model.coherent_defect_wave.chiral_joint.neutral_local_spectrum.symmetries import conjugate_polynomial


def rotate_spin(operator):
    """The exact determinant-one local fermion rotation i sigma_x."""
    terms = []
    for (word, row, col), value in operator.items():
        mapped, phase = [], s.S.One
        for mode, create in word:
            spin, internal = divmod(mode, 15)
            mapped.append((15*(1-spin)+internal, create))
            phase *= -s.I if create else s.I
        for ordered, sign in car.normal_order(tuple(mapped)):
            terms.append({(ordered, row, col): phase*sign*value})
    return car.add(*terms)


@lru_cache(maxsize=1)
def certificate():
    quark = quark_particle_hole()
    h_local = one_cell_polynomial()
    if conjugate_polynomial(h_local) != h_local:
        raise ArithmeticError('native antiunitary local identity failed')
    for part in local_interactions().values():
        if rotate_spin(part) != part:
            raise ArithmeticError('spin rotation changed an actual native channel')
    for index in range(12):
        q = car.gauge_current(index)
        if conjugate_polynomial(q) != car.scale(q, -1) or rotate_spin(q) != q:
            raise ArithmeticError('complete global current transformation failed')
    source = spatial_operators(compact_points())
    if s.trace(source['h']) != 0:
        raise ArithmeticError('the antiunitary spatial constant is nonzero')
    points = source['points']
    address = {p: i for i, p in enumerate(points)}
    reflection, spin_reflection = {}, {}
    for i, (x, y, z) in enumerate(points):
        j = address[-x, y, z]
        reflection[j, i] = 1
        spin_reflection[2*j+1, 2*i] = s.I
        spin_reflection[2*j, 2*i+1] = s.I
        if (-1)**(x+y+z) != (-1)**(-x+y+z):
            raise ArithmeticError('the native pair stagger changed under reflection')
    p = s.ImmutableSparseMatrix(343, 343, reflection)
    r = s.ImmutableSparseMatrix(686, 686, spin_reflection)
    if r.H*r != s.eye(686) or r.H*source['h']*r != -source['h']:
        raise ArithmeticError('complete compact spin-reflection chiral identity failed')
    if p.T*source['K']*p != source['K']:
        raise ArithmeticError('scalar stiffness changed under the same reflection')
    return {'quark_particle_hole': quark,
            'full_compact_rows_checked': 686,
            'all_scalar_sites_checked': 343,
            'spin_reflection': 'R=(x -> -x) tensor i sigma_x',
            'exact_spatial_identity': 'R^dagger h R=-h',
            'exact_scalar_identity': 'P^T K P=K; hence P^T sqrt(K) P=sqrt(K)',
            'all_native_local_interactions_spin_rotation_invariant': True,
            'full_antiunitary_A': 'A c_i A^-1=i c_i^dagger, with scalar conjugation',
            'A_one_body_identity': 'A dGamma(h) A^-1=Tr(h)-dGamma(h)=-dGamma(h)',
            'combined_spatial_antiunitary_C': 'C=R A preserves the full joined H and reverses every global native current',
            'composed_up_down_removal_exchange': 'C after quark particle-hole maps up removal to down removal, with reflected site and unitary spin frame',
            'composed_charge_map': 'Y -> -Y+B, B -> B',
            'exact_equality_domain': 'invariant stationary density on global native singlets with B=0, at reflection-fixed sites or paired reflected probes',
            'actual_local_ground_in_equality_domain': True,
            'full_spatial_B_zero_ground_projector_in_equality_domain_by_construction': True,
            'unrestricted_spatial_global_ground_B_zero_proved': False,
            'Q_squared_changes_under_exchange': '(Y-B)^2 differs from Y^2 on charged target B!=0',
            'polarization_coupling_changes_under_exchange': 'Y X_A -> (-Y+B) X_A; no single memory sign removes the B term'}
