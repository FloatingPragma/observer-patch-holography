"""Assemble the exact class identification, its named checks and the hashes of its inputs."""
from collections import Counter
from fractions import Fraction as F
import hashlib
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from . import flux_costs as fc
from . import model as m
from . import statistics as st
from .census import census

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
SCHEMA = 'oph.defect_particle_classes.v1'


def _s(x):
    return str(x) if isinstance(x, F) else x


def class_rule():
    rows = m.patterns()
    checks = dict(
        six_Y_equals_3n2_minus_2n3=all(r['six_Y'] == 3 * r['n2'] - 2 * r['n3'] for r in rows),
        class_formula=all(r['k'] == (3 * (r['n2'] % 2) - 2 * (r['n3'] % 3)) % 6 for r in rows),
        class_mod_3_is_triality=all(r['k'] % 3 == r['n3'] % 3 for r in rows),
        class_mod_2_is_weak_parity=all(r['k'] % 2 == r['n2'] % 2 for r in rows),
        electric_charge_is_T3_plus_Y=all(r['Q'] == r['T3'] + r['Y'] for r in rows))
    multiplicities = [sum(r['k'] == k for r in rows) for k in range(6)]
    table = [dict(mask=r['mask'], occupied=r['occupied'], n3=r['n3'], n2=r['n2'], six_Y=int(r['six_Y']),
                  Y=_s(r['Y']), T3=_s(r['T3']), Q=_s(r['Q']), k=r['k'], name=r['name']) for r in rows]
    return rows, table, multiplicities, checks


def kernel(rows, fields):
    on_patterns = all(m.kernel_phase_turns(r['n3'], r['n2'], r['Y']).denominator == 1 for r in rows)
    on_fields = all(m.kernel_phase_turns(m.SU3_TRIALITY[f['su3']], m.SU2_TWICE_J[f['su2']], f['Y']).denominator == 1
                    for f in fields)
    centre_characters = all(f['k'] % 3 == m.SU3_TRIALITY[f['su3']] and f['k'] % 2 == m.SU2_TWICE_J[f['su2']] % 2
                            for f in fields)
    return dict(identity_on_all_32_patterns=on_patterns, identity_on_every_field=on_fields,
                field_class_gives_triality_and_weak_parity=centre_characters)


def class_table(rows, fields):
    # patterns grouped by gauge labels must be exactly the fermion fields, component by component
    dims3 = {'1': 1, '3': 3, '3b': 3, '8': 8}
    dims2 = {'1': 1, '2': 2, '2b': 2, '3': 3}
    weak_dim = {'1': 1, '2': 2}
    by_label = Counter((r['colour'], weak_dim[r['weak']], r['Y']) for r in rows)
    fermions = [f for f in fields if f['lorentz'] in 'LR']
    field_labels = Counter()
    for f in fermions:
        field_labels[(f['su3'], dims2[f['su2']], f['Y'])] += dims3[f['su3']] * dims2[f['su2']]
    patterns_are_one_generation = by_label == field_labels
    multiplet_of_pattern = all(
        any(f['name'] == r['multiplet'] and f['k'] == r['k'] and f['Y'] == r['Y'] for f in fermions) for r in rows)
    out = {}
    for k in range(6):
        pats = [r for r in rows if r['k'] == k]
        out[str(k)] = dict(
            colour_triality=k % 3, weak_centre_parity=k % 2, conjugate_class=(-k) % 6, patterns=len(pats),
            pattern_fields=dict(sorted(Counter(r['name'] for r in pats).items())),
            fermion_multiplets=sorted({f['name'] for f in fermions if f['k'] == k}),
            left_handed_names=sorted({f['left_handed_name'] or f['name'] for f in fermions
                                      if f['k'] == k and f['lorentz'] == 'L'}),
            bosons=sorted(f['name'] for f in fields if f['k'] == k and f['lorentz'] in 'SV'))
    left_generation = {'Q': 1, 'u^c': 2, 'd^c': 2, 'L': 3, 'e^c': 0, 'nu^c': 0, 'H': 3}
    by_left = {(f['left_handed_name'] or f['name']): f['k'] for f in fields}
    expected = {'0': (['anti-e_R', 'anti-nu_R', 'e_R', 'nu_R'], ['B', 'G', 'W']), '1': (['Q'], []),
                '2': (['anti-d_R', 'anti-u_R'], []), '3': (['L', 'anti-L'], ['H', 'anti-H']),
                '4': (['d_R', 'u_R'], []), '5': (['anti-Q'], [])}
    checks = dict(
        patterns_are_one_generation_with_nu_R_and_conjugates=patterns_are_one_generation,
        every_pattern_in_its_multiplet_class=multiplet_of_pattern,
        left_handed_generation_classes=all(by_left[n] == k for n, k in left_generation.items()),
        class_table_as_stated=all(
            (out[k]['fermion_multiplets'], out[k]['bosons']) == v for k, v in expected.items()))
    return out, left_generation, checks


def flux_visibility():
    phase, seen = m.aharonov_bohm()
    checks = dict(colour_flux_seen_by_quark_classes_only=seen[2] == seen[4] == [1, 2, 4, 5],
                  weak_flux_seen_by_doublet_classes_only=seen[3] == [1, 3, 5],
                  generating_flux_seen_by_every_nonzero_class=seen[1] == seen[5] == [1, 2, 3, 4, 5])
    return dict(phase_turns={str(f): {str(k): str(v) for k, v in row.items()} for f, row in phase.items()},
                classes_seeing_flux={str(f): v for f, v in seen.items()}), checks


def reactions(fields):
    rows = {}
    for r in sorted(set(m.REACTION_TYPES) | set(m.PAIR_CHANNELS), key=lambda x: (x not in m.REACTION_TYPES, x)):
        a, b, c = m.parse(r)
        reals = m.realizations(r, fields)
        kinds = sorted({x['vertex'] for x in reals})
        rows[r] = dict(vertex_form=f'{a}+{b}->{c}', conserves_class=(a + b - c) % 6 == 0,
                       in_sixteen_types=r in m.REACTION_TYPES, conjugate=m.conjugate_reaction(r),
                       claimed_vertex=m.CLAIMED_VERTEX[r], reading=m.READING[r],
                       realization_vertex_types=kinds, realizations=reals)
    # conjugation closure of the 16 types and the pair channels gives the 27 channels of the grading system
    def norm(t):
        return (t[0],) + tuple(sorted(t[1:])) if len(t) == 3 else tuple(sorted(t))
    fusions, splits = set(), set()
    for r in m.REACTION_TYPES:
        a, b, c = m.parse(r)
        if c == 0:
            continue
        lhs = r.split('->')[0]
        for x, y, w in ((a, b, c), ((-a) % 6, (-b) % 6, (-c) % 6)):
            if '+' in lhs:
                fusions.add(norm((w, x, y)))
            else:
                splits.add(norm((w, x, y)))
    closure_ok = (fusions == {norm((c, a, b)) for a, b, c in st.FUSIONS}
                  and splits == {norm(t) for t in st.SPLITS})
    ycheck = m.yukawa_checks(fields)
    # every gauge-invariant triple among the listed fields conserves the class
    every = True
    for a in fields:
        for b in fields:
            for c in fields:
                n, y_ok = m.gauge_singlets(a, b, c)
                if n and y_ok and (a['k'] + b['k'] - c['k']) % 6:
                    every = False
    checks = dict(
        sixteen_reaction_types=len(m.REACTION_TYPES) == 16 and len(set(m.REACTION_TYPES)) == 16,
        all_reactions_conserve_class=all(x['conserves_class'] for x in rows.values()),
        every_reaction_has_a_gauge_invariant_field_realization=all(x['realizations'] for x in rows.values()),
        claimed_vertex_realized=all(x['claimed_vertex'] in x['realization_vertex_types'] for x in rows.values()),
        diquark_channels_are_colour_epsilon_three_fermion_singlets_only=all(
            x['realization_vertex_types'] == ['three-fermion gauge singlet']
            and all(y['colour_epsilon'] for y in x['realizations'])
            for x in rows.values() if x['claimed_vertex'] == 'three-fermion gauge singlet'),
        yukawa_operators_hypercharge_neutral_and_gauge_invariant=all(v['neutral'] for v in ycheck.values()),
        every_gauge_invariant_field_triple_conserves_class=every,
        conjugation_closure_is_12_fusions_12_splits=closure_ok and len(fusions) == len(splits) == 12)
    return rows, ycheck, checks


def carriers():
    stokes = m.carrier_stokes()
    neutral = m.neutral_multisets(4)
    listed = [tuple(sorted(c)) for c in m.CARRIER_CONTENTS]
    checks = dict(
        carrier_is_closed_icosahedron=(stokes['ports'], stokes['seams'], stokes['faces']) == (12, 30, 20),
        every_seam_on_two_faces_with_opposite_orientation=(stokes['faces_per_seam'] == [2]
                                                           and stokes['seam_column_sums'] == [0]),
        listed_contents_class_neutral=all(sum(c) % 6 == 0 for c in listed),
        listed_contents_among_neutral_multisets=all(c in neutral for c in listed),
        odd_classes_even_in_number_on_neutral_contents=all(sum(x % 2 for x in c) % 2 == 0 for c in neutral))
    return dict(stokes=stokes, listed_contents=[list(c) for c in listed],
                neutral_nonzero_multisets_up_to_four=len(neutral),
                neutral_nonzero_multisets_up_to_three=[list(c) for c in neutral if len(c) <= 3]), checks


def spin_statistics():
    cov = st.cover()
    zdiag = cov['rho_z_diagonal']
    pats = st.z_on_patterns(zdiag)
    grad = st.gradings()
    comp = st.composites()
    links = st.links_two_sided()
    flux = st.port_flux_parity()
    rot = grad['rotation_grading']
    checks = dict(
        cover_has_order_120=cov['order'] == 120 and cov['unit_quaternions'] and cov['closed'],
        rotations_are_special_orthogonal=cov['rotations_special_orthogonal'],
        cover_centre_is_plus_minus_one=cov['centre_size'] == 2 and cov['centre_is_plus_minus_one'],
        galois_conjugate_cover_is_a_group=(cov['galois_conjugation_is_homomorphism']
                                           and cov['galois_image_order'] == 120 and cov['galois_image_closed']),
        code_representation_faithful=cov['representation_kernel_size'] == 1,
        central_involution_is_diag_11100=zdiag == [1, 1, 1, -1, -1],
        central_involution_equals_exp_i_pi_6Y=zdiag == st.exp_i_pi_6Y_on_code(),
        unique_nontrivial_central_involution_of_G=st.central_involutions_of_G() == [[1, 1, 1, -1, -1], [1, 1, 1, 1, 1]],
        z_on_patterns_is_minus_one_to_n2=all(p['z'] == p['minus_one_to_n2'] for p in pats),
        z_on_patterns_is_minus_one_to_6Y=all(p['z'] == p['minus_one_to_6Y'] for p in pats),
        z_on_patterns_is_class_3_holonomy=all(p['z'] == p['class3_holonomy'] for p in pats),
        minus_one_to_degree_outside_G=(-1) ** 5 != 1,
        z_fixes_every_defect_placement=(cov['rotations_permute_faces'] and cov['z_fixes_every_face']
                                        and cov['elements_fixing_every_face'] == 2),
        z_two_sided_on_every_link_block=all(v['two_sided'] == [1] and len(v['z_values']) == 1
                                            for v in links.values()),
        port_flux_parity_reproduces_rho_z=flux['signs'] == zdiag and flux['pattern_sector_parity_equals_n2'],
        vertex_algebra_has_16_gradings=grad['consistent'] and grad['number_of_gradings'] == 16,
        four_gradings_on_physical_states=sorted(grad['physical_generators']) in ([(0, 1), (1, 0)], [(0, 1), (1, 1)], [(1, 0), (1, 1)], [(0, 1), (1, 0), (1, 1)]),
        defect_class_gradings_trivial_and_k_mod_2=grad['defect_class_gradings'] == [[0, 0, 0, 0, 0], [1, 0, 1, 0, 1]],
        rotation_grading_is_a_vertex_grading=grad['rotation_grading_is_a_grading'],
        rotation_grading_all_patterns_odd=set(rot['quanta_by_block'].values()) == {1} and rot['quanta_constant_on_blocks'],
        rotation_grading_register_defects_links_even=(rot['register'] == [0, 0] and set(rot['defects'].values()) == {0}
                                                      and rot['links'] == [0]),
        other_rules_differ_from_rotation=not any(grad['rules_equal_to_rotation'].values()),
        baryon_odd_meson_even_flux_contents_even=(comp['baryon_three_quanta'] == [1]
                                                  and comp['meson_quark_antiquark'] == [0]
                                                  and comp['pure_flux_contents'] == [0]))
    return dict(cover=cov, z_on_code=zdiag, z_on_patterns=pats, link_blocks=links, port_flux_parity=flux,
                gradings=grad, composites=comp), checks


def result():
    rows, table, multiplicities, c1 = class_rule()
    fields = m.fields()
    c2 = kernel(rows, fields)
    ctable, left, c3 = class_table(rows, fields)
    ab, c4 = flux_visibility()
    react, yuk, c5 = reactions(fields)
    carr, c6 = carriers()
    spin, c7 = spin_statistics()
    imported = census()
    c8 = dict(imported_contents_class_neutral=all(imported['carrier_contents_class_neutral'].values()))
    flux = fc.costs()
    c9 = fc.checks(flux)
    checks = {}
    for prefix, group in (('1_class_rule', c1), ('2_z6_kernel', c2), ('3_class_table', c3), ('4_aharonov_bohm', c4),
                          ('5_reactions', c5), ('6_closed_carriers', c6), ('7_spin_statistics', c7),
                          ('8_imported_census', c8), ('9_central_flux_costs', c9)):
        for name, ok in group.items():
            checks[f'{prefix}.{name}'] = bool(ok)
    field_rows = [dict(name=f['name'], su3=f['su3'], su2=f['su2'], Y=str(f['Y']), lorentz=f['lorentz'], k=f['k'],
                       left_handed_name=f['left_handed_name']) for f in fields]
    return dict(
        schema=SCHEMA,
        scope=('Centre classes k = 6Y mod 6 of the five-level code and of the Standard Model fields, flux '
               'visibility, class conservation of the reaction types, closed carriers, and the grading fixed by '
               'the 2 pi rotation. The class does not carry generation, T3 or rest mass, and colour singlet and '
               'octet share a class.'),
        patterns=table, class_multiplicities=multiplicities, fields=field_rows, class_table=ctable,
        left_handed_generation_classes=left, aharonov_bohm=ab, reactions=react, yukawa_operators=yuk,
        closed_carriers=carr, spin_statistics=spin, imported_census=imported, central_flux_costs=flux,
        premises=PREMISES, checks=checks, all_checks_pass=all(checks.values()))


PREMISES = {
    'P1_pointer': 'each quantum carries a spin doublet (pointer) on which the 2 pi rotation acts as -1',
    'P2_gauss_sector': ('on a Gauss-invariant state the code, register and link-end factors of z multiply to +1, so '
                        'the 2 pi rotation on physical states equals (-1)^(number of quanta)'),
    'P3_vertices': ('the local vertices are the twelve native currents, the register toggle, gauge-covariant '
                    'hopping inside a block, the four Yukawa vertices, CPT, and the 27 register channels'),
    'P4_lift': 'the occupations are lifted antisymmetrically; the tick alone admits both lifts',
    'P5_port_sectors': 'vector and spinor levels carry port loads 2 and -3, seen at flux 11 x load mod 20',
    'P6_closed_carrier': 'a carrier is a closed icosahedral screen and register classes are face fluxes of seam links',
}


def certificate():
    seeds = {str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py', '.md')}
    seeds.add('model/requirements.txt')
    value = result()
    value['input_sha256'] = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
                             for p in sorted(import_closure(seeds))}
    return value
