"""Independent outward observable of the recorded finite positive filters."""
from fractions import Fraction as F
from interval_math import *
from interval_reader import read
from counter_core import E, B, BD, lift, counter, cp


def rational(value):
    x = F(value)
    return real(x.numerator)/real(x.denominator)


def matrix(counts):
    return [[complex_(rational(x)) for x in row] for row in counts]


def compiled(record, spin=False):
    mean = matrix(record['counts'])
    coordinate = record.get('coordinate_system', 'energy_target_order' if not spin else None)
    require(coordinate in ('energy_target_order', 'spin_vertex_order'), 'undeclared source coordinates')
    basis = [E[1], E[2], E[0]] if coordinate == 'spin_vertex_order' else E
    t = mm(mm(adj(basis), mean), basis)
    if record['preboost']:
        t = add(scale(power(t, 3), 4), scale(t, -3))
    alpha = 2*(1-rational(record['phase_record']['real']))
    square = power(t, 2)
    return add(eye(2), scale(mm(add(eye(2), scale(square, -1)),
                               power(add(eye(2), scale(square, -alpha)), 2)), -1))


def preparation(record, gram):
    contrast = compiled(record['contrast'])
    cool = compiled(record['cool'])
    leaf = [[conjugate(x)] for x in E[1]]
    v = mm(mm(power(cool, 24), contrast), leaf)
    rho = mm(v, adj(v))
    preparation_success = hermitian_trace(rho)
    weights = [rho[i][i].real/preparation_success for i in range(2)]
    projected = cp(gram, lift(rho))
    position = [row[:2] for row in projected[:2]]
    require(hermitian_trace(position) > 0, 'counter/decode success')
    return position, preparation_success, weights


def scalar(record):
    counts = [[F(x) for x in row] for row in record['counts']]
    require(sum(counts[0]) == 1 and sum(counts[i][i] for i in range(3)) == 1,
            'control optical trace')
    u = (sum(counts[i][j]*counts[j][i] for i in range(3) for j in range(3))-1)/2
    if record['preboost']:
        u *= (4*u-3)**2
    alpha = 2*(1-F(record['phase_record']['real']))
    result = 1-(1-u)*(1-alpha*u)**2
    require(0 < result <= 1, 'finite scalar range')
    return result


def arms(catalog, m, h, turns, direction, stride, scalar_law):
    require(3*m % stride == 0, 'filter stride domain')
    if scalar_law:
        values = [scalar(rec) for rec in catalog]
        a, b = F(1), F(1)
        for step in range(1, 3*m*turns+1):
            if step % stride == 0:
                a *= values[(2*m+direction*step) % (3*m)]**h
                b *= values[2*m]**h
        return a, b
    matrices = [power(compiled(rec, spin=True), h) for rec in catalog]
    a, b = eye(2), eye(2)
    for step in range(1, 3*m*turns+1):
        if step % stride == 0:
            a = mm(matrices[(2*m+direction*step) % (3*m)], a)
            b = mm(matrices[2*m], b)
    return a, b


def returned(gram, position, arm_pair, turns, scalar_law):
    a, b = arm_pair
    sign = -1 if turns % 2 else 1
    if scalar_law:
        positive = cp(gram, lift(position))
        total = hermitian_trace(positive)*rational((a*a+b*b)/2)
        visibility = rational(2*a*b/(a*a+b*b))
        selected = (a+sign*b)**2/4
        tau = scale(positive, rational(selected))
    else:
        aa, bb = sandwich(a, position), sandwich(b, position)
        total = hermitian_trace(cp(gram, lift(scale(add(aa, bb), real(1)/2))))
        cross = trace(cp(gram, lift(sandwich(a, position, b)))).real
        visibility = cross/total
        selected = scale(add(a, scale(b, sign)), real(1)/2)
        tau = cp(gram, lift(sandwich(selected, position)))
    require(total > 0, 'returned analyzer denominator')
    success = hermitian_trace(tau)
    if scalar_law and selected == 0:
        return visibility, real(0), None, total
    require(success > 0, 'returned selected clock absent')
    spectral = scale(mm(mm(BD, tau), B), 1/success)
    weights = [spectral[i][i].real for i in range(4)]
    require(all(x > 0 for x in weights), 'signed clock weights')
    require(1 in sum(weights, real(0)), 'signed clock normalization')
    return visibility, success, weights, total


def pair_records(left, right, z, k):
    weights = [[a*b for b in right] for a in left]
    off = 1-weights[0][0]-weights[1][1]
    result = {}
    for stride in (1, 8):
        steps = stride*2**k
        moment = sum((left[i]*z[i]**steps for i in range(4)), complex_())
        moment *= sum((right[i]*z[i]**steps for i in range(4)), complex_())
        result[str(stride)] = {'x_plus': (1+moment.real)/2, 'x_minus': (1-moment.real)/2,
                               'y_plus': (1+moment.imag)/2, 'y_minus': (1-moment.imag)/2}
    return weights, off, result


