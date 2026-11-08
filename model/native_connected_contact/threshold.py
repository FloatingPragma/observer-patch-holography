"""A primitive Rayleigh trial for the actual connected comparison threshold."""
from fractions import Fraction as F
from model.native_spatial_scaling import geometry as g
from model.native_contact_response import contact

ZERO = (F(), F())


def scale(x, a): return x[0]*a, x[1]*a
def sub(x, y): return x[0]-y[0], x[1]-y[1]
def square(x): return g.multiply(x, x)
def total(xs):
    result = ZERO
    for x in xs: result = g.add(result, x)
    return result


def upper(level):
    contact.integer(level)
    return F(550, 240*4**level)


def trial(level):
    contact.integer(level, 0, 2)
    graph = g.graph(level); rows, denominator = g.coordinates(level); v = graph.ports
    divisor = (F(5*denominator*denominator, 2), F(denominator*denominator, 2))
    raw = tuple((F(row[0]), F(row[1])) for row in rows)
    if total(raw) != ZERO:
        raise ArithmeticError('affine coordinate is not exactly centered')
    products = lambda a, b: g.divide(g.multiply(raw[a], raw[b]), divisor)
    squares = tuple(products(i, i) for i in range(v)); a = total(squares)
    fourth = total(square(x) for x in squares)
    # Audit the signed/cyclic symmetry consequences on literal coordinate data.
    axes = []
    for j in range(3):
        axis = tuple((F(row[2*j]), F(row[2*j+1])) for row in rows)
        if total(axis) != ZERO: raise ArithmeticError('coordinate centering changed')
        axes.append(total(g.divide(square(x), divisor) for x in axis))
    if axes != [a]*3 or not g.sign(sub(a, (F(v, 5), F()))) > 0 or not g.at_most(a, F(v, 3)):
        raise ArithmeticError('the source affine variance bounds changed')
    u = tuple(scale(g.add(scale(x, -1), scale(a, F(1, 2*(v-1)))), F(1, v-2)) for x in squares)
    removed = scale(sub(fourth, scale(square(a), F(1, 2*(v-1)))), F(1, v-2))
    norm = sub(scale(sub(square(a), fourth), F(1, 2)), removed)
    energy = ZERO; projected = ZERO; contact_energy = ZERO; coordinate_energy = ZERO
    for i, j in graph.seams:
        delta = g.divide(square(sub(raw[i], raw[j])), divisor)
        coordinate_energy = g.add(coordinate_energy, delta)
        energy = g.add(energy, g.multiply(delta, sub(sub(a, squares[i]), squares[j])))
        projected = g.add(projected, scale(square(sub(squares[i], squares[j])), F(1, v-2)))
        h = sub(sub(products(i, j), u[i]), u[j])
        contact_energy = g.add(contact_energy, scale(square(h), 2))
    action_energy = sub(energy, projected)
    numerator = g.add(action_energy, contact_energy); quotient = g.divide(numerator, norm)
    if (g.sign(norm) <= 0 or not g.at_most(coordinate_energy, 4)
            or g.sign(sub(norm, (F(v*v, 75), F()))) < 0
            or not g.at_most(quotient, upper(level))):
        raise ArithmeticError('the exact primitive connected threshold trial failed')
    return {'level': level, 'ports': v, 'coordinate_squared_norm': g.display(a),
            'coordinate_energy': g.display(coordinate_energy),
            'primitive_squared_norm': g.display(norm),
            'exchange_energy': g.display(action_energy), 'contact_energy': g.display(contact_energy),
            'comparison_Rayleigh_quotient': g.display(quotient),
            'compressed_threshold_upper': str(upper(level))}


def certificate():
    return {'literal_trials': [trial(level) for level in range(3)],
            'all_level_upper': '550/V = 55/(24*4^level)',
            'formal_level200_upper': str(upper(200)),
            'trial_is_in_global_primitive_degree_two': True,
            'vanishing_lower_bound_used_as_collapse_proof': False,
            'isolated_pole_above_connected_comparison_bottom': True,
            'all_bound_states_or_resonances_excluded': False,
            'all_level_native_floating_constructor_return_claimed': False}
