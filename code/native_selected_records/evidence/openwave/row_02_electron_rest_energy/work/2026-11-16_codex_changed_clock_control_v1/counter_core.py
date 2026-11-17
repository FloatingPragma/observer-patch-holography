"""Four-mode counter map with measured rational source means as explicit inputs.
No producer imports; actual finite preparations and arms are supplied downstream.
"""
from interval_math import *

def basis():
    s2, s6 = iv.sqrt(real(2)), iv.sqrt(real(6))
    e = [[complex_(), complex_(-2/s6)],
         [complex_(1/s2), complex_(1/s6)],
         [complex_(-1/s2), complex_(1/s6)]]
    b = zeros(4)
    for i in range(2):
        b[i][2*i] = b[i][2*i+1] = complex_(1/s2)
        b[2+i][2*i] = complex_(0, -1/s2)
        b[2+i][2*i+1] = complex_(0, 1/s2)
    return e, b

E, B = basis()
BD = adj(B)

def lift(position):
    rho = zeros(4)
    for i in range(2):
        for j in range(2):
            rho[i][j] = position[i][j]
    return rho

def counter(k, means):
    """Exact geometric finite word, enclosed at every arithmetic operation."""
    q = 256*2**k
    eps = real(1)/(2**k)
    z = []
    for mean in means:
        p = 1-eps*eps*(1-real(mean.numerator)/mean.denominator)
        t = iv.sqrt(1-p*p)
        z.extend([complex_(p, t), complex_(p, -t)])
    numer = [(1-a**q)/q for a in z]
    root = root_of_unity(q)
    phase = complex_(1)
    gram = zeros(4)
    for r in range(1, q//2):
        phase *= root
        vals = [n/(1-a*phase) for n, a in zip(numer, z)]
        for i in range(4):
            for j in range(i, 4):
                gram[i][j] += vals[i]*conjugate(vals[j])
    for i in range(4):
        for j in range(i):
            gram[i][j] = conjugate(gram[j][i])
    return gram, z

def cp(gram, rho):
    spectral = mm(mm(BD, rho), B)
    out = [[spectral[i][j]*gram[i][j] for j in range(4)] for i in range(4)]
    return mm(mm(B, out), BD)
