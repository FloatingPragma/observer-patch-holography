"""Assemble the register, its weak round, its clock unit and its read counts, with named checks and input hashes."""
from fractions import Fraction as F
import hashlib
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from . import carrier, classes, clock, hexagon, kick, register

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
SCHEMA = 'oph.lepton_register_round.v1'
SECTIONS = ('1_register', '2_double_loads', '3_weak_round', '4_clock_unit', '5_kick', '6_carrier', '7_classes')


def result():
    reg, c1 = register.result()
    dbl, c2, expected = hexagon.section_double_loads()
    weights = reg['object_weights']['C']
    c2['double_load_weights_equal_the_register_object_weights'] = (
        (F(weights['idle']), F(weights['agree']), F(weights['repair +']), F(weights['repair -']))
        == (expected['double'], expected['in-dimer'], expected['forward'], expected['backward']))
    rnd, c3, _ = hexagon.section_weak_round()
    clk, c4 = clock.result()
    kck, c5 = kick.result()
    car, c6 = carrier.result()
    cls, c7 = classes.result()
    checks = {}
    for prefix, group in zip(SECTIONS, (c1, c2, c3, c4, c5, c6, c7)):
        for name, ok in group.items():
            checks[f'{prefix}.{name}'] = bool(ok)
    return dict(
        schema=SCHEMA,
        scope=('The charged-lepton generation register on the coherent seam law: the 432-state tick and its '
               'compression (sqrt2/3) C(theta), the read count 2/9 of two independent objects, the weak round with '
               'the inverse tick, the unit of the declared clock, the kick as the commutator of the oriented toggle '
               'with the tick of the pair, the counts on the sixty-address carrier, and the selection of the pair on '
               'one port by a class-3 flux. The register, its preparation, its analyzer, the kick and the fresh '
               'rounds are declared constructions on the public law; the lepton scale and a source of the double '
               'load are outside this package.'),
        register=reg, double_loads=dbl, weak_round=rnd, clock_unit=clk, kick=kck, carrier=car, classes=cls,
        premises=PREMISES, checks=checks, all_checks_pass=all(checks.values()))


PREMISES = {
    'P1_register': ('two seam objects, C on the seam (2, 3) and E on the seam (2, 10), each with its own load, '
                    'record and two-orientation address; a neutral register F carrying the third-turn R of the '
                    'face (4, 5, 11); T applies R or R^-1 when C repairs at orientation + or -; the clock '
                    'K_theta = exp(i theta J) acts on the oriented current J = |1><1|_(r_C) sigma_z(o_C); V is the '
                    'oriented record, Z_r at the reversed address'),
    'P2_preparation_and_analyzer': ('each object starts with activity 1/3, its load at the screen tail of its '
                                    'seam, a blank record and the uniform address; the analyzer is the Fourier '
                                    'zero over the flags E agrees with C idle, C repairs + with E idle, and C '
                                    'repairs - with E idle'),
    'P3_double_load_objects': ('two independent hexagon objects, each with two loads, the four-valued record of '
                               'the public completion and its own twelve-address register reflected about its '
                               'uniform state; the standing state S = (D_2 + D_10)/sqrt2 has blank records'),
    'P4_weak_round': ('a round is one tick, the kick exp(-i h J (x) Qx) on the exclusive read Qx of fresh '
                      'records, the inverse tick and a return test projecting on the standing state; N rounds '
                      'per unit clock time with h = 1/N'),
    'P5_declared_clock': ('the clock H = sum_x (I - Re U_x) + sum_edges (pi/4) dGamma(I - S_xy) with evolution '
                          'exp(-i t H), as declared in model/coherent_defect_wave/energy'),
    'P6_oriented_toggle': ('the named interaction T(eta) = exp(i eta J (x) Z_F Z_I) between the clocked '
                           'register and the flags of the pair'),
    'P7_two_cells': ('two cells with independent sixty-address registers, each holding two loads moved by the '
                     'public completion at the addressed seam in the order of the address'),
    'P8_flux': ('a class-k flux through the hexagon is the link phase exp(i pi k/3) on the hexagon seam (2, 3), '
                'entering as the gauge conjugation of the completion at both addresses of that seam'),
}


def certificate():
    seeds = {str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py', '.md')}
    seeds.add('model/requirements.txt')
    value = result()
    value['input_sha256'] = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
                             for p in sorted(import_closure(seeds))}
    return value
