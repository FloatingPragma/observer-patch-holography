"""Exact local premises for a full-CAR thermodynamic ground construction.

The finite checks bind the local action and its boundary stabilization.
Compactness and positive-energy GNS construction are proved in README.md;
no finite array is called an infinite-volume state.
"""
from fractions import Fraction as F
from functools import lru_cache

from ..coarse_energy_ward.ward import Clifford, annihilator, commutator, energies, total
from ..full_fock_field.field import FieldGeometry
from ..full_fock_field.bounds import local_boundary_bound, stability_bounds


def embed(value, small, large):
    """Exact CAR inclusion for nested lexicographically ordered site lists."""
    if not isinstance(value, Clifford) or not set(small.sites) <= set(large.sites):
        raise ValueError('a Clifford polynomial and nested source geometries required')
    index = {x: j for j, x in enumerate(large.sites)}
    mapping = [4*index[x]+j for x in small.sites for j in range(4)]
    if mapping != sorted(mapping):
        raise ValueError('this inclusion requires the inherited lexicographic order')
    terms = {}
    for mask, coefficient in value.terms.items():
        if mask.bit_length() > len(mapping):
            raise ValueError('polynomial exceeds its stated source algebra')
        target = sum(1 << mapping[j] for j in range(mask.bit_length()) if mask >> j & 1)
        terms[target] = coefficient
    return Clifford(terms)


def local_monomials(geometry, site):
    if site not in geometry.sites:
        raise ValueError('site outside geometry')
    first = 4*geometry.sites.index(site)
    return tuple(Clifford({mask << first: F(1)}) for mask in range(16))


def stationarity_error(radius, support_radius, support_sites, duration, speed=F(1, 8)):
    """For any extended finite-box ground, |omega_R(tau_t(A))-omega_R(A)|.

    This is the original local-dynamics boundary bound for ||A||<=1,
    because omega_R(tau_t^R(A))=omega_R(A) exactly. It is NOT a rate of
    convergence of the finite-box states or a selected ground algorithm.
    """
    return local_boundary_bound(radius, support_radius, support_sites, duration, speed)


def folner_translation_error(radius, shift, observable_norm=F(1)):
    """State-independent defect of cube-averaging translates of any state."""
    if (type(radius) is not int or radius < 0 or type(shift) is not tuple or
            len(shift) != 3 or any(type(x) is not int for x in shift) or
            type(observable_norm) not in (int, F) or observable_norm < 0):
        raise ValueError('nonnegative integer cube, integer shift and exact norm required')
    side = 2*radius+1
    overlap = 1
    for value in shift:
        overlap *= max(0, side-abs(value))
    return 2*F(observable_norm)*(1-F(overlap, side**3))


@lru_cache(maxsize=1)
def source_certificate():
    small, large = FieldGeometry((3, 3, 3)), FieldGeometry((4, 4, 4))
    hs, hl = total(energies(small)), total(energies(large))
    checked = 0
    for a in local_monomials(small, (1, 1, 1)):
        if not (embed(commutator(hs, a), small, large)-
                commutator(hl, embed(a, small, large))).is_zero():
            raise ArithmeticError('an exterior even interaction affected a disjoint CAR probe')
        checked += 1
    # A nonzero contact contribution on an even observable survives in every
    # containing algebra and every nonzero unital GNS representation.
    a0, a1, a2, a3 = (annihilator(j) for j in range(4))
    n = [a.adjoint()@a for a in (a0, a1, a2, a3)]
    contact = n[0]@n[1]+n[2]@n[3]
    probe = a0.adjoint()@a2
    quartic = probe@(n[1]-n[3])
    if quartic.is_zero() or not (commutator(contact, probe)-quartic).is_zero():
        raise ArithmeticError('the original contact derivation was lost')
    degree_four = [mask for mask in quartic.terms if mask.bit_count() == 4]
    if not degree_four:
        raise ArithmeticError('contact witness has no genuine quartic coefficient')
    # Deliberately omitting an incident hopping destroys stabilization.
    origin = small.sites.index((1, 1, 1))
    neighbor = small.sites.index((2, 1, 1))
    left, right = annihilator(2*origin), annihilator(2*neighbor+1)
    term = (left.adjoint()@right).scale((0, F(1, 16)))
    edge = term+term.adjoint()
    if commutator(edge, left).is_zero():
        raise ArithmeticError('missing-boundary-term control failed')
    stability = stability_bounds((3, 3, 3))
    return {'source': 'H=dGamma[-c sigma.sin k]+U sum_x n_xup n_xdown; c=1/8,U=1,mu=0',
            'small_box_sites': 27, 'larger_box_sites': 64,
            'full_Fock_dimensions_symbolic': ['2^54', '2^128'],
            'complete_one_site_CAR_basis_checked': checked,
            'odd_probes_retained': True,
            'local_derivation_stabilizes_after_one_neighbor_layer': True,
            'nonzero_contact_derivation': '[D,a0^dagger a2]=a0^dagger a2(n1-n3)',
            'quartic_majorana_terms': len(degree_four),
            'nonzero_in_every_nonzero_unital_GNS_representation': True,
            'reason_for_representation_claim': 'The witness belongs to a full local matrix algebra; every unital representation of that simple algebra is injective.',
            'missing_incident_bond_control_nonzero': True,
            'full_Fock_stability_27_sites': {k: str(v) if isinstance(v, F) else v for k,v in stability.items()},
            'fixed_spacing_stationarity_budgets': [
                {'radius': r, 'support_radius': 0, 'support_sites': 1,
                 'time': '1/8', 'error_upper': str(stationarity_error(r, 0, 1, F(1, 8)))}
                for r in (4, 8, 16, 32, 64)],
            'translation_average_defects': [
                {'radius': r, 'shift': [1, 0, 0],
                 'error_upper': str(folner_translation_error(r, (1, 0, 0)))}
                for r in (1, 4, 16, 64)]}


def existence_certificate():
    """Theorem ledger tied to the complete proof, not a computed limit state."""
    return {
        'domain': 'Z^3 CAR with two modes/site; fixed finite c>=0,U>=0; original mu=0 H; no particle-number restriction',
        'finite_reference': 'normalized spectral projector of the absolute full-Fock minimum of each open-box H',
        'outside_extension': 'even product with vacuum on the complementary CAR algebra; no conditioning',
        'limit': 'a weak-star cluster state of spatially averaged nested-box full-ground extensions; countable local matrix entries permit a subsequence',
        'local_ground_inequality': 'omega(A^*[H,A])>=0 for every local A, including odd A',
        'stationary_under_original_thermodynamic_dynamics': True,
        'gauge_invariant_state_exists': True,
        'translation_invariant_ground_state_exists_by_Folner_average': True,
        'GNS': 'U_t pi(A)Omega=pi(tau_t(A))Omega; U_t=exp(i t L), L>=0, L Omega=0',
        'spectral_proof': 'Finite full-ground correlation measures have support[0,infinity). Compact-time locality convergence and weak-star convergence pass every Schwartz Fourier test supported on negative energies. Density of local GNS vectors makes the negative spectral projection vanish.',
        'interacting_derivation_represented_faithfully': True,
        'Hamiltonian_scalar_shift_or_chemical_potential_added': False,
        'ground_state_computed_or_prepared': False,
        'chosen_subsequence_or_convergence_rate': False,
        'state_functional_local_faithfulness_asserted': False,
        'uniqueness_or_purity_asserted': False,
        'uniform_gap_or_clustering_asserted': False,
        'complete_repair_stationarity_asserted': False,
        'spacing_zero_or_Lorentz_continuum_asserted': False,
        'nontrivial_scattering_or_vacuum_connected_four_point_asserted': False,
        'physical_particle_identifications_added': 0}
