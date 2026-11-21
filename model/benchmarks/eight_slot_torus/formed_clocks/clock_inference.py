"""Read elapsed time or frequency from local four-port counts.

The paired records compare independent terminal cohorts of the same
enrollment at one fixed receiver. At a mature even lapse all protected
spatial amplitudes return, while the rail coherence acquires a common
phase. Its initial direction and magnitude can differ at each receiver.
The inputs start with +X,+Y,-X,-Y counts; subsequent cells are retained
misses. Absolute appointments, source amplitudes and hidden lapses do
not enter the elapsed-time decoder. Rate inference uses independently
supplied local elapsed intervals and no source frequency.

For local incidence I and coherence Z, write u+iv=2Z/I. The four effects
have probabilities I(1+u,1+v,1-u,1-v)/4. Conditional on its two X
counts' total the positive X count is binomial with probability (1+u)/2;
the Y pair reads v in the same way. Retaining the original-start misses
also permits an input-law error allowance at arbitrarily small incidence.
A confidence rectangle containing zero supplies no phase constraint.
All enrolled words receive simultaneous coverage before selecting any
phase, lag or candidate time from their records.
"""
from __future__ import annotations

import math
from fractions import Fraction
from functools import lru_cache
from numbers import Integral

import numpy as np
from mpmath.ctx_iv import MPIntervalContext
from scipy.special import betaincinv

PI = math.pi
TURN = 2*PI
PAD = 3e-11
_iv = MPIntervalContext()
_iv.dps = 52


@lru_cache(maxsize=8192)
def _log_choose(n, k):
    return _iv.loggamma(n+1)-_iv.loggamma(k+1)-_iv.loggamma(n-k+1)


@lru_cache(maxsize=32768)
def _lower(k, n, tail):
    """Conservatively enclose a one-sided Clopper-Pearson lower endpoint.

    At a trial p the binomial right tail from k is at most mass(k)/(1-r),
    where r=mass(k+1)/mass(k)<1 bounds all its successive mass ratios.
    Outward interval evaluation verifies that upper bound against the
    requested tail area. The inverse-beta seed affects only tightness.
    """
    if k == 0:
        return 0.
    area = _iv.mpf(tail.numerator)/tail.denominator
    log_area = _iv.log(area)
    norm = _log_choose(n,k)
    def safe(p):
        if p == 0:
            return True
        if p >= 1:
            return False
        z = _iv.mpf(p)
        ratio = (n-k)*z/((k+1)*(1-z))
        if not ratio < 1:
            return False
        log_mass = norm+k*_iv.log(z)+(n-k)*_iv.log(1-z)
        return bool(log_mass-_iv.log(1-ratio) <= log_area)

    high = float(betaincinv(k, n-k+1, float(tail)))
    if not math.isfinite(high):
        high = k/n
    if safe(high):
        return high
    width = max(math.ulp(high)*2, .005*math.sqrt(high*(1-high)/n))
    low = max(0., high-width)
    while not safe(low):
        width *= 2
        low = max(0., high-width)
    for _ in range(8):
        mid = (low+high)/2
        if safe(mid):
            low = mid
        else:
            high = mid
    return low


def _settings(observed, family_size, alpha):
    size = observed if family_size is None else family_size
    if isinstance(size, bool) or not isinstance(size, Integral) or size < observed or not 0 < alpha < 1:
        raise ValueError("all enrolled count words and a failure probability in (0,1) required")
    # Each word has four binomial port tails and two lower axis-incidence tails.
    return Fraction(str(alpha))/(6*int(size))


@lru_cache(maxsize=32768)
def _hoeffding_lower(k, n, tail):
    """A one-sided, outward-rounded binomial endpoint for a fixed shot budget."""
    if not k:
        return 0.
    area = _iv.mpf(tail.numerator)/tail.denominator
    bound = _iv.mpf(k)/n - _iv.sqrt(-_iv.log(area)/(2*n))
    return max(0., math.nextafter(float(bound.a), -math.inf))


def _words(data, ndim):
    a = np.asarray(data)
    if (a.ndim != ndim or a.shape[-1] < 4 or 0 in a.shape
            or not np.issubdtype(a.dtype, np.integer) or np.any(a < 0)):
        raise ValueError("nonnegative integer words with four ordered ports required")
    return a


def _arc(word, tail, cell_error, lower_endpoint=_lower):
    """An unwrapped covering arc; None represents the entire circle."""
    total = sum(int(v) for v in word)
    if cell_error and len(word) == 4:
        raise ValueError("source-error correction requires the original-start miss count")
    coordinates = []
    for i in (0,1):
        plus, minus = int(word[i]), int(word[i+2])
        n = plus+minus
        if not n:
            coordinates.append((-1.,1.))
            continue
        lower = 2*lower_endpoint(plus,n,tail)-1
        upper = 1-2*lower_endpoint(minus,n,tail)
        correction = 0.
        if cell_error:
            incidence = lower_endpoint(n,total,tail)
            if incidence <= 2*cell_error:
                coordinates.append((-1.,1.))
                continue
            # |u_source-u_sample| <= 4*cell_error/p_sample(axis).
            correction = 4*cell_error/incidence
        coordinates.append((max(-1.,lower-correction-PAD), min(1.,upper+correction+PAD)))
    (xl,xh),(yl,yh) = coordinates
    near_x, near_y = max(0.,xl,-xh), max(0.,yl,-yh)
    if near_x**2+near_y**2 > 1+PAD:
        return ()
    if xl <= 0 <= xh and yl <= 0 <= yh:
        return None
    # The angle of a rectangle disjoint from the origin has width < pi.
    # Each extreme is attained at a corner. Using the whole rectangle is
    # conservative even where a corner lies outside the physical unit disk.
    center = math.atan2(yl+yh,xl+xh)
    angles = [center+math.remainder(math.atan2(y,x)-center,TURN)
              for x in (xl,xh) for y in (yl,yh)]
    return min(angles)-PAD, max(angles)+PAD


def _merge(intervals):
    ordered = []
    for left,right in sorted(intervals):
        if ordered and left <= ordered[-1][1]+PAD:
            ordered[-1][1] = max(ordered[-1][1],right)
        else:
            ordered.append([left,right])
    return ordered


def _restrict(domain, difference, lag):
    """Clip frequency intervals by a modular lag*rate constraint."""
    low, high = difference
    if high-low >= TURN:
        return domain
    result = []
    for left,right in domain:
        first = math.ceil((lag*left-high)/TURN-PAD)
        last = math.floor((lag*right-low)/TURN+PAD)
        for k in range(first,last+1):
            lo = max(left,(low+TURN*k)/lag-PAD)
            hi = min(right,(high+TURN*k)/lag+PAD)
            if lo <= hi:
                result.append((lo,hi))
    return _merge(result)


def _difference(a,b):
    return b[0]-a[1]-PAD, b[1]-a[0]+PAD


def _error_budget(cell_error):
    if not math.isfinite(cell_error) or not 0 <= cell_error < 1:
        raise ValueError("finite absolute original-start cell error in [0,1) required")


def decode_elapsed(early, late, *, rate, cap=8191, even=True, alpha=1e-6,
                   family_size=None, rate_radius=0., cell_error=0., interval="cp"):
    """Decode candidate elapsed *integer services* using a calibrated signed rate.

    Each row is the same local receiver in the two cohorts. All count
    words share one simultaneous envelope, also when a receiver has no
    detected hand. A calibrated rate interval can be supplied. Its own
    calibration failure probability must be accounted for separately.
    Without a booked same-parity comparison every odd lapse is untested.
    """
    a, b = _words(early,2), _words(late,2)
    _error_budget(cell_error)
    if a.shape != b.shape or not math.isfinite(rate) or not math.isfinite(rate_radius) or rate_radius < 0:
        raise ValueError("matched station words and a finite signed calibration required")
    if type(cap) is not int or not 0 <= cap <= 8191 or type(even) is not bool:
        raise ValueError("integer pre-wrap cap in [0,8191] and a parity booking required")
    tail = _settings(2*len(a), family_size, alpha)
    if interval not in ("cp", "hoeffding"):
        raise ValueError("unknown binomial interval")
    lower_endpoint = _lower if interval == "cp" else _hoeffding_lower
    rate = math.remainder(rate, TURN)
    phases = [[-PI,PI]]
    informative = 0
    invalid = 0
    for left,right in zip(a,b):
        before = _arc(left,tail,cell_error,lower_endpoint)
        after = _arc(right,tail,cell_error,lower_endpoint)
        if before == () or after == ():
            phases = []
            invalid += 1
            continue
        if before is None or after is None:
            continue
        informative += 1
        phases = _restrict(phases,_difference(before,after),1)
    candidates = []
    for n in range(0,cap+1,2):
        center = math.remainder(n*rate,TURN)
        half = n*rate_radius+PAD
        if any(center-half <= hi+shift and center+half >= lo+shift
               for lo,hi in phases for shift in (-TURN,0.,TURN)):
            candidates.append(n)
    runs = []
    for n in candidates:
        if runs and n == runs[-1][1]+2:
            runs[-1][1] = n
        else:
            runs.append([n,n])
    return {"even_runs_step_two":runs, "untested_odd_run_step_two":
            ([1,cap if cap%2 else cap-1] if not even and cap else []),
            "phase_increment_intervals":phases, "informative_stations":informative,
            "stations":len(a), "station_incompatible_with_unit_disk":invalid,
            "eligible_even_count":len(candidates),
            "original_start_count_total":sum(int(x) for x in a.flat)+sum(int(x) for x in b.flat)}


def infer_rate(words, ticks, *, alpha=1e-6, family_size=None, cell_error=0.):
    """Cover every signed per-service rate consistent with booked local spacings.

    The common rate is constrained by every same-parity difference. The
    station and parity may each have their own unknown phase offset. The
    returned union is an outer region, without imposing a source-rate prior.
    At an even-only cadence rates separated by pi are indistinguishable.
    """
    data = _words(words,3)
    _error_budget(cell_error)
    ts = np.asarray(ticks)
    if (ts.ndim != 1 or ts.size != data.shape[1] or not np.issubdtype(ts.dtype,np.integer)
            or any(int(y) <= int(x) for x,y in zip(ts,ts[1:]))
            or int(ts[-1])-int(ts[0]) >= 8192):
        raise ValueError("increasing local integer timestamps in a pre-wrap window required")
    tail = _settings(data.shape[0]*data.shape[1], family_size, alpha)
    by_lag = {}
    invalid = 0
    for station in data:
        arcs = [_arc(row,tail,cell_error) for row in station]
        invalid += sum(arc == () for arc in arcs)
        for j, after in enumerate(arcs):
            if not after:
                continue
            for i in range(j):
                lag = int(ts[j])-int(ts[i])
                if lag % 2 == 0 and arcs[i]:
                    by_lag.setdefault(lag,[]).append(_difference(arcs[i],after))
    bands = [] if invalid else [[-PI,PI]]
    for lag in sorted(by_lag):
        # The preimage of all station constraints for a lag can be lifted
        # once, even if thousands of full turns are compatible initially.
        phase = [[-PI,PI]]
        for difference in by_lag[lag]:
            phase = _restrict(phase,difference,1)
        bands = _merge([v for region in phase for v in _restrict(bands,region,lag)])
    return {"per_service_intervals":bands, "used_even_lags":sorted(by_lag),
            "words_incompatible_with_unit_disk":invalid,
            "station_count":data.shape[0], "local_appointments":len(ts),
            "original_start_count_total":sum(int(x) for x in data.flat)}
