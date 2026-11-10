"""Recompute the information preserved and lost by the six-class register."""
from collections import Counter
from functools import lru_cache

import sympy as s

from model.standard_model_source_map.algebra import source_components, match_components


def unit(n, a, b):
    m = s.zeros(n)
    m[a, b] = 1
    return m


def exterior_rows():
    """All exterior occupations, before attaching particle or statistics names."""
    rows = []
    for mask in range(32):
        n3 = (mask & 7).bit_count()
        n2 = (mask >> 3).bit_count()
        q6 = 3 * n2 - 2 * n3
        t3 = s.Rational(((mask >> 3) & 1) - ((mask >> 4) & 1), 2)
        rows.append(dict(mask=mask, n3=n3, n2=n2, q6Y=q6,
                         Y=str(s.Rational(q6, 6)), T3=str(t3), Q=str(t3+s.Rational(q6, 6)),
                         centre_class=q6 % 6, colour_triality=n3 % 3,
                         weak_centre_parity=n2 % 2,
                         exterior_degree_parity=(n3+n2) % 2))
    return rows


def centre_information():
    rows = exterior_rows()
    groups = []
    for k in range(6):
        selected = [r for r in rows if r['centre_class'] == k]
        assert all(r['colour_triality'] == k % 3 and r['weak_centre_parity'] == k % 2
                   for r in selected)
        groups.append(dict(centre_class=k, colour_triality=k % 3,
                           weak_centre_parity=k % 2, exterior_dimension=len(selected),
                           q6Y_values=sorted({r['q6Y'] for r in selected}),
                           Q_values=sorted({r['Q'] for r in selected}),
                           exterior_degree_parities=sorted({r['exterior_degree_parity'] for r in selected}),
                           native_particle_identity=None))
    # This is the existing current representation, recomputed before comparison.
    actual = source_components()
    by_mask = {r['mask']: r for r in rows}
    assert len(actual) == 15
    for r in actual:
        x = by_mask[r['occupation_mask']]
        assert s.Rational(r['q6Y']) == x['q6Y']
        assert r['Y'] == x['Y'] and r['T3'] == x['T3'] and r['Q'] == x['Q']
        assert (x['n3']+x['n2']) in (2, 4)
    return dict(rows=rows, centre_classes=groups,
                actual_fifteen_current_match=match_components(),
                actual_fifteen_auxiliary_exterior_parity=0,
                auxiliary_exterior_parity_is_physical_CAR_parity=False)


@lru_cache(maxsize=1)
def adjoint_witness():
    """Full Lie-current Casimir distinguishes centre-neutral octet and singlet."""
    generators = []
    for a, b in ((0, 1), (0, 2), (1, 2)):
        generators.extend(((unit(3,a,b)+unit(3,b,a))/2,
                           -s.I*(unit(3,a,b)-unit(3,b,a))/2))
    generators.extend((s.diag(1,-1,0)/2, s.diag(1,1,-2)/(2*s.sqrt(3))))
    assert all(s.trace(a*b) == (s.Rational(1,2) if i == j else 0)
               for i,a in enumerate(generators) for j,b in enumerate(generators))
    adj = [s.kronecker_product(t,s.eye(3))-s.kronecker_product(s.eye(3),t.T)
           for t in generators]
    casimir = sum((a*a for a in adj),s.zeros(9))
    identity = s.eye(3).reshape(9,1)
    singlet = identity*identity.T/3
    assert casimir == 3*(s.eye(9)-singlet)
    assert all(a*identity == s.zeros(9,1) for a in adj)
    charged = unit(3,0,1).reshape(9,1)
    assert casimir*charged == 3*charged
    # A centre zI acts on End(C^3) by z X / z = X, on both summands.
    return dict(tensor_product='3 x 3bar = 1 + 8', common_centre_class=0,
                singlet_rank=singlet.rank(), non_singlet_rank=casimir.rank(),
                Casimir_singlet='0', Casimir_octet='3',
                centre_neutral_implies_full_Gauss_singlet=False)


def pair_witness():
    """An even diquark and an odd d-conjugate can have identical gauge currents."""
    # Antisymmetric colour pair: wedge^2(3)=3bar. Antisymmetric weak pair: wedge^2(2)=1.
    colour = [s.diag(1,-1,0),s.diag(0,1,-1)]
    pairs = [(0,1),(0,2),(1,2)]
    weights = {tuple(str(t[a,a]+t[b,b]) for t in colour) for a,b in pairs}
    assert weights == {('-1','0'),('1','-1'),('0','1')}
    weak = s.Matrix([0,1,-1,0])/s.sqrt(2)
    for p in (s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)):
        assert (s.kronecker_product(p,s.eye(2))+s.kronecker_product(s.eye(2),p))*weak == s.zeros(4,1)
    q6_pair = 1+1
    fermion_parity_pair = (1+1) % 2
    assert q6_pair == 2 and fermion_parity_pair != 1
    return dict(incoming='Q_L x Q_L', channel='(3bar,1,Y=1/3)',
                comparison='d^c internal gauge representation',
                total_q6Y=q6_pair, centre_class=q6_pair % 6,
                tensor_colour_weights=sorted([list(x) for x in weights]),
                pair_CAR_parity=0, single_d_conjugate_CAR_parity=1,
                pair_baryon_number='2/3', single_d_conjugate_baryon_number='-1/3',
                direct_pair_to_single_fermion_preserves_parity=False,
                same_gauge_representation_implies_same_particle=False)


def conjugation_witness():
    rows = exterior_rows()
    for r in rows:
        c = rows[31 ^ r['mask']]
        assert c['q6Y'] == -r['q6Y'] and s.Rational(c['Q']) == -s.Rational(r['Q'])
        assert c['centre_class'] == -r['centre_class'] % 6
    # Class three equals its inverse, although its full hypercharge is nonzero.
    return dict(complement_pairs=16, class_three_is_self_inverse=True,
                class_three_nonzero_q6Y=sorted({r['q6Y'] for r in rows if r['centre_class']==3}),
                centre_self_inverse_implies_Majorana_or_self_conjugate=False)


def result():
    return dict(schema='oph.particle-catalogue-correspondence.v1',
                information=centre_information(), adjoint=adjoint_witness(),
                pair=pair_witness(), conjugation=conjugation_witness(),
                source_derived_particle_masses_GeV={},
                required_transfer=[
                    'An operator map intertwining every full gauge current, with boundary charges retained.',
                    'The actual physical CAR grading and exchange instrument, distinct from auxiliary exterior degree.',
                    'One preparation, interacting Hamiltonian, reference and clock, with complete detector outcomes.',
                    'Momentum-dependent spectral support and controlled errors for the claimed particle or resonance.',
                    'A physical energy convention and input ancestry for comparison with measured masses.'
                ])
