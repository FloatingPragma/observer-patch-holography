"""All-parameter charged-spectrum certificate, independent of eigenvalue fits."""
from __future__ import annotations

import json
import sympy as sp


def certificate():
    # Work in the onsite-untwisted Hubbard frame d_1=i sigma_axis c_1.
    # Each CAR matrix is independently assembled from the wedge sign.
    t = sp.symbols('t', positive=True)
    u = sp.symbols('U', nonnegative=True)
    z = sp.symbols('z')
    operators = []
    for mode in range(4):
        a = sp.zeros(16)
        for mask in range(16):
            if mask >> mode & 1:
                a[mask ^ (1 << mode), mask] = (-1)**((mask & ((1 << mode)-1)).bit_count())
        operators.append(a)
    h = sp.zeros(16)
    number = sp.zeros(16)
    densities = [a.T*a for a in operators]
    for n in densities:
        number += n
    for spin in range(2):
        hop = operators[spin].T*operators[spin+2]
        h += t*(hop+hop.T)
    h += u*(densities[0]*densities[1]+densities[2]*densities[3])
    expected = [z, (z*z-t*t)**2,
                z**3*(z-u)*(z*z-u*z-4*t*t),
                ((z-u)**2-t*t)**2, z-2*u]
    polynomials = []
    for sector in range(5):
        indices = [i for i in range(16) if i.bit_count() == sector]
        polynomial = h.extract(indices, indices).charpoly(z).as_expr()
        # charpoly normalizes the generator symbol's assumptions.
        polynomial = polynomial.subs(next(s for s in polynomial.free_symbols if s.name == 'z'), z)
        if sp.expand(polynomial-expected[sector]) != 0:
            raise ArithmeticError('charged-family characteristic polynomial mismatch')
        polynomials.append(str(sp.factor(expected[sector])))
    if h*number != number*h:
        raise ArithmeticError('number conservation failed')
    # Check the eigenvalue/root identities also at the noninteracting endpoint.
    # With t>0, d>U, so the N=2 low singlet is simple even when U=0.
    d = sp.sqrt(u*u+16*t*t)
    roots = ((0,), (-t, t), (0, u, (u-d)/2, (u+d)/2),
             (u-t, u+t), (2*u,))
    for polynomial, sector_roots in zip(expected, roots, strict=True):
        for root in sector_roots:
            if sp.simplify(polynomial.subs(z, root)) != 0:
                raise ArithmeticError('all-sector root identity failed')
            if sp.simplify(polynomial.subs(z, root).subs(u, 0)) != 0:
                raise ArithmeticError('noninteracting endpoint root identity failed')
    # Let d=sqrt(U²+16t²). Then d>=4t for U>=0; no sampling enters this proof.
    # E0(K)=-(U+d)/2, and the N=1,3 eigenvalues are -U/2 +/-t.
    # Thus g_-=d/2-t>=t>0 and 3g_- - g_+=d-4t>=0.
    return {
        'schema': 'oph.dimer_family.v1',
        'domain': 't>0, U>=0, mu=U/2, unique N=2 ground; complete 16-state field',
        'parameter_trials': 'symbolic universal certificate; no finite parameter scan',
        'sector_characteristic_polynomials': polynomials,
        'ground_energy_of_H_minus_mu_N': '-(U+sqrt(U^2+16t^2))/2',
        'half_filled_ground_separations': {
            'N0_and_N4': '(U+d)/2>0',
            'N1_and_N3_lowest': 'd/2-t>=t>0',
            'N2_triplet': '(d-U)/2>0 since d^2-U^2=16t^2>0',
            'N2_other_singlets': '(U+d)/2>0 and d>0',
        },
        'strict_N2_ground_interval': 'abs(mu-U/2)<d/2-t, d=sqrt(U^2+16t^2)',
        'ground_interval_proof': 'N1/N3 gaps are g_minus +/- (mu-U/2); N0/N4 are positive inside because d-U<=4t; the N2 internal gaps are independent of mu',
        'off_center_charge_one_gaps': {
            'addition': 'd/2 +/- t - (mu-U/2)',
            'removal': 'd/2 +/- t + (mu-U/2)',
        },
        'charge_one_gaps': ['sqrt(U^2+16t^2)/2-t', 'sqrt(U^2+16t^2)/2+t'],
        'multiplicity_per_gap_per_charge': 2,
        'gap_ratio_lower_strict': 1,
        'gap_ratio_upper_inclusive': 3,
        'proof': 'd^2-(4t)^2=U^2>=0; g_minus>=t>0; 3g_minus-g_plus=d-4t>=0',
        'higher_local_residue': '(1-4t/sqrt(U^2+16t^2))/4; zero at U=0',
        'maximum_distinct_energies_per_charge': 2,
        'three_same_charge_distinct_mass_families_possible': False,
        'assumptions': ['one common positive linear energy-to-mass conversion',
                        'mass assigned to a ground-relative charge-one gap',
                        'same half-filled reference and action for all species'],
        'boundary': 'Other references, nonlinear dictionaries, larger carriers and other actions are outside the ratio theorem. The two-level count also holds for arbitrary mu while this N=2 ground persists.',
        'physical_identification': False,
    }


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
