"""Source-blind two-stride reader with outward inverse-branch enclosures."""
from interval_math import iv, real, require, endpoints, acos_point, pi_bound


def principal(record, eta, delta):
    r = 2*record['x_plus']-1
    lo = max(endpoints((r-2*delta-eta)/(1-eta))[0], -1)
    hi = min(endpoints((r+2*delta+eta)/(1-eta))[1], 1)
    require(lo <= hi, 'inconsistent cosine record')
    left = endpoints(acos_point(real(hi)))[0]
    right = endpoints(acos_point(real(lo)))[1]
    return left, right


def read(first, later, eta='0.01', probability_error='0.00002', stride=8):
    eta, delta = real(eta), real(probability_error)
    require(eta >= 0 and eta < 1 and delta >= 0 and stride >= 1, 'reader parameters')
    lo, hi = principal(first, eta, delta)
    a, b = principal(later, eta, delta)
    branches = []
    for j in range(-stride, stride+1):
        for sign in (1, -1):
            band = (2*pi_bound()*j+sign*iv.mpf([a, b]))/stride
            x, y = endpoints(band)
            x, y = max(lo, x), min(hi, y)
            if x <= y:
                branches.append((x, y))
    branches.sort()
    merged = []
    for x, y in branches:
        if merged and x <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(y, merged[-1][1]))
        else:
            merged.append((x, y))
    return [iv.mpf([x, y]) for x, y in merged]
