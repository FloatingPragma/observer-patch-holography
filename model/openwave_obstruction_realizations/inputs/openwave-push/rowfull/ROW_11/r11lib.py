"""Row 11, shared definitions (own code; imports nothing from the theory repository or from a frozen producer).

Circulant on the three characters of a third-turn:
    lambda_k = 1 + 2 beta mod cos(theta + 2 pi k / 3),  k = 0, 1, 2,
beta^2 = (|b|/a)^2 the balance, theta the phase in radians, mod a modulus factor (1 in the limit of many rounds).
Reading: masses m_k = lambda_k^2, sorted ascending.
"""
import hashlib
import json
from fractions import Fraction as F

import mpmath as mp

DPS = 45                      # working digits; 40 are reported
mp.mp.dps = DPS

# Measured values, as cited by reverse-engineering-reality/model/lepton_cycle/physical_comparison.json:
# CODATA 2022 (muon/electron ratio, electron energy) and PDG 2024 (tau mass, Phys. Rev. D 110, 030001).
RMU_MEAS = F('206.7682827')
RMU_SIG = F('0.0000046')
ME_MEV = F('0.51099895069')
ME_SIG = F('0.00000000016')
MTAU_MEV = F('1776.93')
MTAU_SIG = F('0.09')
RTAU_MEAS = MTAU_MEV / ME_MEV
RTAU_SIG = MTAU_SIG / ME_MEV

BAL = F(1, 2)                 # the register's balance (|b|/a)^2
PHASE = F(2, 9)               # the register's phase


def mpf(x):
    if isinstance(x, F):
        return mp.mpf(x.numerator) / mp.mpf(x.denominator)
    return mp.mpf(x)


def levels(beta2, theta, mod=1):
    b = mp.sqrt(mpf(beta2)) * mpf(mod)
    th = mpf(theta)
    return [1 + 2 * b * mp.cos(th + 2 * mp.pi * k / 3) for k in range(3)]


def reading(beta2, theta, mod=1):
    """Sorted squared levels and the derived numbers of one alternative."""
    lam = levels(beta2, theta, mod)
    m = sorted(x * x for x in lam)
    out = dict(lam=lam, min_level=min(lam), positive=all(x > 0 for x in lam))
    out['koide'] = sum(m) / sum(mp.sqrt(x) for x in m) ** 2
    if m[0] == 0:
        out['r_mu'] = out['r_tau'] = None
    else:
        out['r_mu'] = m[1] / m[0]
        out['r_tau'] = m[2] / m[0]
    return out


def distance(r_mu, r_tau):
    """d = max |ln(r/r_meas)| over the two ratios (preregistered ranking distance)."""
    if r_mu is None:
        return mp.inf
    return max(abs(mp.log(r_mu / mpf(RMU_MEAS))), abs(mp.log(r_tau / mpf(RTAU_MEAS))))


def hits(r_mu, r_tau):
    """The preregistered hit levels H1..H4 and H4tau."""
    if r_mu is None:
        return dict(H1=False, H2=False, H3=False, H4=False, H4tau=False)
    e_mu = abs(r_mu / mpf(RMU_MEAS) - 1)
    e_tau = abs(r_tau / mpf(RTAU_MEAS) - 1)
    s_mu = abs(r_mu - mpf(RMU_MEAS)) / mpf(RMU_SIG)
    s_tau = abs(r_tau - mpf(RTAU_MEAS)) / mpf(RTAU_SIG)
    return dict(H1=bool(e_mu <= mp.mpf('1e-2') and e_tau <= mp.mpf('1e-2')),
                H2=bool(e_mu <= mp.mpf('1e-3') and e_tau <= mp.mpf('1e-3')),
                H3=bool(e_mu <= mp.mpf('1e-4') and e_tau <= mp.mpf('1e-4')),
                H4=bool(s_mu <= 3 and s_tau <= 3),
                H4tau=bool(s_tau <= 3 and e_mu <= mp.mpf('1e-4')))


def literal_gate(r_mu, r_tau):
    """G2a: |r_mu - 207| <= 1/2 and |r_tau - 3477| <= 1/2."""
    if r_mu is None:
        return False
    return bool(abs(r_mu - 207) <= mp.mpf('0.5') and abs(r_tau - 3477) <= mp.mpf('0.5'))


def s(x, n=40):
    """Decimal string with n significant digits."""
    return mp.nstr(mpf(x), n)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def dump(path, obj):
    with open(path, 'w') as fh:
        json.dump(obj, fh, indent=1, ensure_ascii=False, sort_keys=False)
        fh.write('\n')
