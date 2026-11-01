"""Exact finite energy continuity on the complete Weyl/contact CAR algebra.

The half-bond local-energy convention is explicit. Integer real/imaginary
arithmetic certifies all operator identities, including every number sector.
The reference state, model time and action coefficients remain supplied.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

from scipy import sparse

from ..field_stress.energy import GaussianInteger, commutator, sum_exact
from ..full_fock_field.field import FieldGeometry, FullFock


@dataclass(frozen=True)
class EnergyOperators:
    geometry: FieldGeometry
    scale: int
    local: tuple[GaussianInteger, ...]
    total: GaussianInteger
    currents: dict[tuple[int, int], GaussianInteger]

    @property
    def dimension(self):
        return 1 << self.geometry.modes

    def region(self, sites):
        """Return energy and outward current numerators for a site subset."""
        if (type(sites) not in (tuple, list, set, frozenset) or
                any(type(x) is not int or not 0 <= x < self.geometry.volume for x in sites)
                or len(set(sites)) != len(sites)):
            raise ValueError('distinct existing integer sites required')
        sites = set(sites)
        energy = sum_exact([self.local[x] for x in sites], self.dimension)
        flux = GaussianInteger.zero(self.dimension)
        for (x, y), current in self.currents.items():
            if (x in sites) != (y in sites):
                flux = flux+current if x in sites else flux-current
        return energy, flux


@dataclass(frozen=True)
class GrandCanonicalOperators:
    field: EnergyOperators
    grand: EnergyOperators
    chemical_potential: Fraction
    local_numbers: tuple[GaussianInteger, ...]
    total_number: GaussianInteger
    number_currents: dict[tuple[int, int], GaussianInteger]


def _integer_multiple(operator, factor):
    if type(factor) is not int:
        raise ValueError('integer rescaling required')
    # Preserve the same conservative entry bound as from_exact; perform no
    # unbounded int64 product or floating reduction of a large coefficient.
    if operator.maximum()*abs(factor) > 2**30:
        raise OverflowError('rescaled Gaussian-integer entry bound exceeded')
    return GaussianInteger.from_exact(factor*operator.matrix())


def local_energies(geometry=FieldGeometry((2, 1, 1))):
    """Build h_x=U n_x0 n_x1 + 1/2 sum_(b incident x) H_b.

    The mathematical construction is finite-volume general. The explicit
    certificate admits at most four sites, exact binary-rational parameters,
    denominator at most 2**16 and bounded Gaussian-integer numerators.
    """
    if not isinstance(geometry, FieldGeometry) or geometry.volume > 4:
        raise ValueError('one through four complete field sites required')
    scale = math.lcm((Fraction(geometry.speed)/4).denominator,
                     Fraction(geometry.interaction).denominator)
    if scale > 2**16:
        raise ValueError('exact coefficient denominator exceeds certificate domain')
    space = FullFock(geometry.modes)
    a = [space.annihilator(p) for p in range(geometry.modes)]
    local = [sparse.csr_matrix((space.dimension, space.dimension), dtype=complex)
             for _ in geometry.sites]
    for group in geometry.groups:
        for p, q, z in group:
            term = z*(a[p].conj().T@a[q])
            half = (term+term.conj().T)/2
            local[p//2] += half
            local[q//2] += half
    for site in range(geometry.volume):
        local[site] += geometry.interaction*(space.number(2*site)@space.number(2*site+1))
    local = tuple(GaussianInteger.from_exact(scale*h) for h in local)
    total = sum_exact(local, space.dimension)
    direct = GaussianInteger.from_exact(scale*space.hamiltonian(geometry))
    if not (total-direct).is_zero():
        raise ArithmeticError('local energies do not sum to the actual full-Fock action')
    currents = {}
    for x in range(geometry.volume):
        for y in range(x+1, geometry.volume):
            current = commutator(local[x], local[y]).times_i()
            if not current.is_zero():
                currents[x, y] = current
    return EnergyOperators(geometry, scale, local, total, currents)


def certify(geometry=FieldGeometry((2, 1, 1))):
    """Certify exact complete-algebra identities, with no reference restriction."""
    operators = local_energies(geometry)
    dim = operators.dimension
    space = FullFock(geometry.modes)
    number = GaussianInteger.from_exact(sum((space.number(p) for p in range(geometry.modes)),
                                           sparse.csr_matrix((dim, dim))))
    maximum_range = 0
    for x, energy in enumerate(operators.local):
        if not (energy-energy.adjoint()).is_zero():
            raise ArithmeticError('non-Hermitian local energy')
        if not commutator(number, energy).is_zero():
            raise ArithmeticError('local energy fails number conservation')
        _, divergence = operators.region((x,))
        if not (commutator(operators.total, energy).times_i()+divergence).is_zero():
            raise ArithmeticError('local energy continuity failed')
    for (x, y), current in operators.currents.items():
        if not (current-current.adjoint()).is_zero():
            raise ArithmeticError('non-Hermitian energy current')
        if not commutator(number, current).is_zero():
            raise ArithmeticError('energy current fails number conservation')
        distance = sum(abs(a-b) for a, b in zip(geometry.sites[x], geometry.sites[y], strict=True))
        if distance > 2:
            raise ArithmeticError('energy current exceeds overlapping-bond support')
        maximum_range = max(maximum_range, distance)
    # Exhaust all site subsets: cancellation uses oriented currents including
    # range-two terms induced by the half-bond local energy assignment.
    for mask in range(1 << geometry.volume):
        region = tuple(x for x in range(geometry.volume) if mask >> x & 1)
        energy, outward = operators.region(region)
        if not (commutator(operators.total, energy).times_i()+outward).is_zero():
            raise ArithmeticError('regional continuity failed')
    return {'geometry': list(geometry.sides), 'field_modes': geometry.modes,
            'full_Fock_dimension': dim,
            'number_sector_dimensions': [math.comb(geometry.modes, n)
                                         for n in range(geometry.modes+1)],
            'local_energy_denominator': operators.scale,
            'energy_current_denominator': operators.scale**2,
            'all_site_subsets_checked': 1 << geometry.volume,
            'all_unordered_site_pairs_checked': geometry.volume*(geometry.volume-1)//2,
            'nonzero_energy_current_pairs': len(operators.currents),
            'maximum_current_range': maximum_range,
            'complete_operator_identities_exact': True,
            'number_sector_projection_applied': False,
            'reference_state_required_for_continuity': False,
            'physical_time_energy_or_gravity_identified': False}, operators


def grand_canonical(geometry=FieldGeometry((2, 1, 1)), chemical_potential=Fraction(1, 2)):
    """Exact k_x=h_x-mu*n_x, number currents and their grand-energy current.

    Number-current numerators have field.scale as denominator; grand-current
    numerators have grand.scale**2 as denominator. The number-current formula
    uses the declared half-bond energy convention, not an arbitrary partition.
    """
    if type(chemical_potential) not in (int, Fraction):
        raise ValueError('exact rational chemical potential required')
    mu = Fraction(chemical_potential)
    energy = local_energies(geometry)
    scale = math.lcm(energy.scale, mu.denominator)
    if scale > 2**16:
        raise ValueError('grand-energy coefficient denominator exceeds certificate domain')
    space = FullFock(geometry.modes)
    numbers = tuple(GaussianInteger.from_exact(space.number(2*x)+space.number(2*x+1))
                    for x in range(geometry.volume))
    total_number = sum_exact(numbers, space.dimension)
    coefficient = int(mu*scale)
    grand_local = tuple(_integer_multiple(h, scale//energy.scale)-_integer_multiple(n, coefficient)
                        for h, n in zip(energy.local, numbers, strict=True))
    grand_total = sum_exact(grand_local, space.dimension)
    expected = _integer_multiple(energy.total, scale//energy.scale)-_integer_multiple(total_number, coefficient)
    if not (grand_total-expected).is_zero():
        raise ArithmeticError('grand-canonical partition identity failed')
    number_currents, grand_currents = {}, {}
    for x in range(geometry.volume):
        for y in range(x+1, geometry.volume):
            current_n = (commutator(energy.local[x], numbers[y])+
                         commutator(numbers[x], energy.local[y])).times_i()
            current_k = commutator(grand_local[x], grand_local[y]).times_i()
            current_h = energy.currents.get((x, y), GaussianInteger.zero(space.dimension))
            # i[k_x,k_y] = i[h_x,h_y] - mu*j_N(x,y), since the n_x commute.
            expected = (_integer_multiple(current_h, (scale//energy.scale)**2)-
                        _integer_multiple(current_n, coefficient*(scale//energy.scale)))
            if not (current_k-expected).is_zero():
                raise ArithmeticError('grand-current j_K=j_H-mu*j_N identity failed')
            if not current_n.is_zero():
                number_currents[x, y] = current_n
            if not current_k.is_zero():
                grand_currents[x, y] = current_k
    grand = EnergyOperators(geometry, scale, grand_local, grand_total, grand_currents)
    for x in range(geometry.volume):
        divergence_n = GaussianInteger.zero(space.dimension)
        for (a, b), current in number_currents.items():
            if a == x:
                divergence_n = divergence_n+current
            elif b == x:
                divergence_n = divergence_n-current
        if not (commutator(energy.total, numbers[x]).times_i()+divergence_n).is_zero():
            raise ArithmeticError('number continuity failed')
        if not commutator(total_number, energy.local[x]).is_zero():
            raise ArithmeticError('H-local energy does not conserve total number')
        _, divergence_k = grand.region((x,))
        if not (commutator(grand_total, grand_local[x]).times_i()+divergence_k).is_zero():
            raise ArithmeticError('grand-energy continuity failed')
    return GrandCanonicalOperators(energy, grand, mu, numbers, total_number, number_currents)


def charged_probe_certificate(bundle):
    """Check [K,O]=[H,O]-mu*q*O for every elementary charge q=+/-1 probe."""
    if not isinstance(bundle, GrandCanonicalOperators):
        raise ValueError('exact grand-canonical operators required')
    space = FullFock(bundle.field.geometry.modes)
    qscale = bundle.grand.scale
    for mode in range(space.modes):
        lower = GaussianInteger.from_exact(space.annihilator(mode))
        for probe, charge in ((lower, -1), (lower.adjoint(), 1)):
            if not (commutator(bundle.total_number, probe)-_integer_multiple(probe, charge)).is_zero():
                raise ArithmeticError('elementary CAR probe has wrong number charge')
            expected = (_integer_multiple(commutator(bundle.field.total, probe),
                                          qscale//bundle.field.scale)-
                        _integer_multiple(probe, int(bundle.chemical_potential*qscale)*charge))
            if not (commutator(bundle.grand.total, probe)-expected).is_zero():
                raise ArithmeticError('charged probe violates K/H energy identity')
    return {'chemical_potential': str(bundle.chemical_potential),
            'charge_changing_CAR_probes_checked': 2*space.modes,
            'number_continuity_and_grand_current_exact': True,
            'transition_identity': 'Delta K = Delta H - mu Delta N',
            'grand_current_identity': 'j_K = j_H - mu j_N',
            'field_and_grand_energy_identified_with_each_other': False}
