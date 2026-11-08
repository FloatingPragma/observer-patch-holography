"""Independent exact identities on both complete interacting sectors."""
from fractions import Fraction as F

from .source import Controls, Sector, inputs


BLOCKS = ((0, 3, F(4)), (2, 6, F(25, 8)), (3, 3, F(3)))


def scaled(vector, value):
    return {i:x*value for i, x in vector.items() if x*value}


def keep(vector, indices):
    return {i:x for i, x in vector.items() if i in indices}


def certificate():
    data = inputs()
    rows = []
    controls = (Controls(), Controls(F(3, 2), F(2, 3), F(5, 4), F(7, 5), F(3, 2)),
                Controls(F(5, 4), F(-3, 2), F(1, 3), F(2, 7), F(7, 3)))
    for q in (0, 1):
        make = lambda c: Sector(q, c, data)
        base = make(Controls())
        electric = make(Controls(1, 0, 0, 0, 0))
        loop = make(Controls(0, -1, 0, 0, 0))
        hopping = make(Controls(0, 0, 1, 0, 0))
        higgs = make(Controls(0, 0, 0, 1, 0))
        scalar = make(Controls(0, 0, 0, 0, 1))
        seed = {base.seed:F(1)}
        s_seed, t_seed = higgs.apply(seed), hopping.apply(seed)
        if base.inner(s_seed, s_seed) != 12-2*q or base.inner(t_seed, t_seed) != q:
            raise ArithmeticError('source Higgs/hopping seed norm')
        for j in range(base.dimension):
            e = {j:F(1)}
            if higgs.apply(loop.apply(e)) != loop.apply(higgs.apply(e)):
                raise ArithmeticError('full Higgs/loop commutator')
        colored = {0, 2, 3}
        signs = [-1 if x['link0_block'] in colored and x['link1_block'] in colored else 1
                 for x in base.labels]
        if signs[base.seed] != 1:
            raise ArithmeticError('loop sign involution must fix the seed')
        for op, expected in ((electric, 1), (hopping, 1), (higgs, 1), (scalar, 1), (loop, -1)):
            if any(signs[i]*signs[j] != expected for i, j in op.numerators):
                raise ArithmeticError('full loop-sign involution identity')
        moment_rows = []
        for r, d, c in BLOCKS:
            indices = base.block(r)
            ell = scaled(keep(loop.apply(seed), indices), F(1, d))
            if base.inner(ell, ell) != 1 or electric.apply(ell) != scaled(ell, 2*c):
                raise ArithmeticError('normalized native loop or two-link energy')
            if any(hopping.apply({j:F(1)}) for j in indices):
                raise ArithmeticError('hopping must vanish on the entire colored block')
            u = higgs.apply(ell)
            if base.inner(u, u) != 12-2*q:
                raise ArithmeticError('complete Higgs loop complement')
            for control in controls:
                sector = make(control)
                m, e, h, g = control.magnetic, control.electric, control.hopping, control.higgs
                powers = sector.powers(3)
                if keep(powers[1], indices) != scaled(ell, -m*d):
                    raise ArithmeticError('first full-H projected amplitude')
                expected = scaled(ell, -2*m*d*e*c)
                for i, value in scaled(u, -2*m*d*g).items():
                    expected[i] = expected.get(i, 0)+value
                expected = {i:x for i, x in expected.items() if x}
                if keep(powers[2], indices) != expected:
                    raise ArithmeticError('second full-H projected amplitude')
                third = -m*d*(4*e*e*c*c+54*m*m+3*g*g*(12-2*q)+q*h*h)
                if sector.inner(ell, powers[3]) != third:
                    raise ArithmeticError('third full-H return amplitude')
                moments = sector.moments(r)
                if moments[2] != d*d*m*m or moments[4] != -d*d*m*m*(e*e*c*c/3+18*m*m+q*h*h/3):
                    raise ArithmeticError('complete quartic coefficient identity')
                moment_rows.append({'block':r, 'controls':control.packet(),
                                    'moments':{str(k):str(v) for k, v in moments.items()}})
        full = {r:base.moments(r) for r, _, _ in BLOCKS}
        removed = make(Controls(higgs=F(0)))
        differences = {str(r):full[r][6]-removed.moments(r)[6] for r, _, _ in BLOCKS}
        expected_differences = {str(r):-d*d*F(27-5*q, 15) for r, d, _ in BLOCKS}
        if differences != expected_differences:
            raise ArithmeticError('interacting sixth-order color response')
        no_hop, neither = make(Controls(hopping=F(0))), make(Controls(hopping=F(0), higgs=F(0)))
        if any(no_hop.moments(r)[6] != neither.moments(r)[6] for r, _, _ in BLOCKS):
            raise ArithmeticError('Higgs dependence must vanish without hopping')
        missing = {str(r):full[r][4]-base.moments(r, nh=0)[4] for r, _, _ in BLOCKS}
        if missing != {str(r):F(d*d*(12-2*q)) for r, d, _ in BLOCKS}:
            raise ArithmeticError('omitted Higgs probability did not change quartic response')
        rows.append({'Qstar':q, 'dimension':base.dimension,
                     'full_Higgs_loop_commutator_columns':base.dimension,
                     'loop_sign_involution_checked_on_every_generator_entry':True,
                     'Higgs_seed_norm_squared':12-2*q, 'hopping_seed_norm_squared':q,
                     'complete_moment_checks':moment_rows,
                     'Higgs_induced_sixth_order_color_change':{r:str(v) for r, v in differences.items()},
                     'quartic_probability_lost_if_Higgs_outcomes_removed':{r:str(v) for r, v in missing.items()}})
    return rows
