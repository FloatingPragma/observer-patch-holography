"""Native-record block encoding and complete amplified choice instruments.

Only the coefficient s comes from the native census. No ray, rotation angle,
spin matrix, desired phase or visibility enters the compiler.
"""
import math
from fractions import Fraction
import numpy as np
from native_source import native

OMEGA = np.ones(3) / math.sqrt(3)
OPTICAL = np.array([[0., -2/math.sqrt(6)],
                    [1/math.sqrt(2), 1/math.sqrt(6)],
                    [-1/math.sqrt(2), 1/math.sqrt(6)]])
TRIANGLE = ((0, 1), (1, 2), (0, 2))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def routing(slot, m):
    vertex, step = divmod(slot % (3*m), m)
    other = (vertex+1) % 3
    return [e for e in TRIANGLE if vertex in e]*(m-step) + [e for e in TRIANGLE if other in e]*step


def census(edges, law):
    """Known basis/choice inputs; observer reconstructs output only from currents."""
    width = 2*len(edges)
    require(width > 0 and width & (width-1) == 0, 'Choice count must be dyadic')
    permutation = np.empty(3*width, dtype=int)
    mean = np.zeros((3,3))
    records = []
    for initial in range(3):
        for ticket, edge in enumerate(edges):
            for coin in (0,1):
                q = tuple(int(p == initial) for p in range(3))
                output_audit, event = native(q, edge, coin, law)
                a,b,c,current = event
                reconstructed = list(q)
                reconstructed[a] -= current
                reconstructed[b] += current
                require(tuple(reconstructed) == output_audit, 'Current reconstruction failed')
                final = reconstructed.index(1)
                memory_coin = c ^ int(current != 0)
                permutation[initial*width+2*ticket+c] = final*width+2*ticket+memory_coin
                mean[final,initial] += 1/width
                records.append([initial,ticket,a,b,c,current,final,memory_coin])
    require(np.array_equal(np.sort(permutation), np.arange(3*width)), 'Native reversal is not a permutation')
    require(np.array_equal(permutation[permutation], np.arange(3*width)), 'Native reversal is not involutive')
    require(np.linalg.norm(mean-mean.T) < 1e-13, 'Mean is not symmetric')
    require(np.linalg.norm(mean@OMEGA-OMEGA) < 1e-13, 'Mean does not preserve constant sector')
    return mean, permutation, records


def hadamard_choices(columns, width):
    """Tensor product of fair native blank-bit transforms on the choice register."""
    values = columns.reshape(3,width,-1).copy()
    step = 1
    while step < width:
        for start in range(0,width,2*step):
            a = values[:,start:start+step,:].copy()
            b = values[:,start+step:start+2*step,:].copy()
            values[:,start:start+step,:] = (a+b)/math.sqrt(2)
            values[:,start+step:start+2*step,:] = (a-b)/math.sqrt(2)
        step *= 2
    return values.reshape(3*width,-1)


def compile_instrument(edges, law, fixed_calibration=None):
    mean, permutation, records = census(edges,law)
    width = 2*len(edges)
    counts = [[Fraction(float(mean[i,j])) for j in range(3)] for i in range(3)]
    trace = sum(counts[i][i] for i in range(3))-1
    square_trace = sum(counts[i][j]*counts[j][i] for i in range(3) for j in range(3))-1
    discriminant = 2*square_trace-trace*trace
    require(discriminant >= 0,'Invalid optical characteristic polynomial')
    norm_record = {'optical_trace':str(trace),'optical_square_trace':str(square_trace),
                   'discriminant':str(discriminant)}
    s = (abs(float(trace))+math.sqrt(float(discriminant)))/2
    spectrum = np.array([(float(trace)-math.sqrt(float(discriminant)))/2,
                         (float(trace)+math.sqrt(float(discriminant)))/2])
    if fixed_calibration is not None:
        s = fixed_calibration['s']
        norm_record = {'fixed_native_calibration':fixed_calibration['norm_record']}
    require(.249999999999 < s <= 1.000000000001, 'Outside admitted norm range')
    # The same census-derived rule is used on both laws. No law branch here.
    preboost = s < .5-1e-12
    boosted_norm = abs(4*s**3-3*s) if preboost else s
    require(boosted_norm >= .5-1e-12, 'One preliminary boost is insufficient')
    cosine = 1-1/(2*boosted_norm**2)
    sine2 = 1-cosine*cosine
    require(sine2 >= -1e-12, 'Invalid matching phase')
    # At the exact endpoint boosted_norm=1/2, sin(phi)=0 analytically.
    sine = math.sqrt(sine2) if sine2 > 0 else 0.
    phase = complex(cosine,sine)
    initial = np.zeros((3*width,3),complex)
    initial[np.arange(3)*width,np.arange(3)] = 1

    def primitive(x):
        v = hadamard_choices(x,width)
        out = np.empty_like(v)
        out[permutation] = v
        return hadamard_choices(out,width)

    def mark(x,z):
        out = x.copy()
        out[np.arange(3)*width] *= z
        return out

    def signal(x):
        if not preboost:
            return primitive(x)
        return primitive(mark(primitive(mark(primitive(x),-1.)), -1.))

    def amplified(z):
        return signal(mark(signal(mark(signal(initial),z)),z))

    positive = amplified(phase).reshape(3,width,3).transpose(1,0,2)
    inverse = amplified(phase.conjugate()).reshape(3,width,3).transpose(1,0,2)
    for instrument in (positive,inverse):
        completeness = sum(op.conj().T@op for op in instrument)
        require(np.linalg.norm(completeness-np.eye(3)) < 1e-11, 'Incomplete amplified source instrument')
    require(np.linalg.norm(inverse[0]-positive[0].conj().T) < 1e-11, 'Inverse successful map mismatch')
    block = inverse[0]@positive[0]
    require(np.linalg.norm(block-block.conj().T) < 1e-11, 'Success filter is not positive Hermitian')
    effective = 4*mean@mean@mean-3*mean if preboost else mean
    predicted = np.eye(3)-(np.eye(3)-effective@effective)@(np.eye(3)-effective@effective/boosted_norm**2)@(np.eye(3)-effective@effective/boosted_norm**2)
    require(np.linalg.norm(block-predicted) < 1e-10, 'Source-derived polynomial disagrees with full circuit')
    return {'P':mean,'permutation':permutation,'records':records,'s':s,
            'phase':[phase.real,phase.imag],'preboost':preboost,
            'native_reversals_per_filter':18 if preboost else 6,
            'choice_width':width,'forward':positive,'inverse':inverse,
            'K':block,'polynomial':predicted,'spectrum':spectrum,'norm_record':norm_record}
