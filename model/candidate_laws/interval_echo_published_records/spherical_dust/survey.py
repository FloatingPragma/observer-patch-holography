"""Proper-ball areas from five rounded clock readings at each surveyed endpoint."""
from dataclasses import dataclass
from decimal import Decimal as D, getcontext, localcontext
from functools import lru_cache

from model.candidate_laws.interval_echo_partition_repair.published_strings.finite_rods import RADIAL


ZERO = D(0)
ONE = D(1)
BEACON = D(1) / 8
AXES = tuple(tuple(D(i == j) for i in range(3)) for j in range(3))
INNER, OUTER = D(1)/4, D(1)/2
SCALE_FIELDS = ("clock", "inner_radius", "outer_radius")


def cube_root(x):
    """Positive local solution, including the exact flat calibration."""
    if x <= 0:
        raise ValueError("clock interval reaches the focusing time")
    root = ONE
    for _ in range(getcontext().prec.bit_length()+6):
        better = (2*root + x/(root*root))/3
        if better == root:
            break
        root = better
    return root


def conformal(t, h):
    y = cube_root(1-3*h*t/2)
    return 3*t/(1+y+y*y)


def proper(eta, h):
    return eta*(1-h*eta/2+(h*eta)**2/12)


def scale(t, h):
    return cube_root(1-3*h*t/2)**2


@lru_cache(maxsize=256)
def _series(w, precision):
    """Clock and spatial ray coefficients in h*s from the spacelike launch."""
    # p=f**a satisfies f*p'=a*f'*p. Here f=1-3*y/2, y=h*t.
    def power_term(n, exponent, coefficients, clock):
        if n == 0:
            return ONE
        return sum(((exponent+1)*k-n)*(-3*clock[k]/2)*coefficients[n-k]
                   for k in range(1, n+1))/n

    clock = [ZERO, w]
    flight = [ZERO]
    clock_power, flight_power = [], []
    for n in range(14):
        clock_power.append(power_term(n, -D(7)/3, clock_power, clock))
        flight_power.append(power_term(n, -D(4)/3, flight_power, clock))
        clock.append((1+w*w)*clock_power[n]/((n+1)*(n+2)))
        flight.append(flight_power[n]/(n+1))
    return tuple(clock[1:]), tuple(flight[1:])


def _polynomial(coefficients, x):
    out = ZERO
    for coefficient in reversed(coefficients):
        out = coefficient + x*out
    return out


def comoving_distance(s, h, clock_slope):
    """Distance multiplier of the conserved spatial launch momentum."""
    _, flight = _series(clock_slope, getcontext().prec)
    return s*_polynomial(flight, h*s)


def ray(s, h, gamma, boost, direction):
    """Plan an endpoint appointment from the independently measured scale history."""
    w = boost*direction[2]
    clock, _ = _series(w, getcontext().prec)
    z = h*s
    t = s*_polynomial(clock, z)
    distance = comoving_distance(s, h, w)
    return t, (distance*direction[0], distance*direction[1], distance*gamma*direction[2])


def timestamp_event(t, point, h, tick):
    """Coincident local dust clock plus four received, retarded beacon stamps."""
    origin2 = sum(x*x for x in point)
    squared = [origin2] + [origin2 + BEACON*BEACON - 2*BEACON*x for x in point]
    eta = conformal(t, h)
    emissions = [proper(eta-d.sqrt(), h) for d in squared]
    return tuple(int((v/tick).to_integral_value()) for v in (t, *emissions))


def locate(stamps, h, tick):
    if len(stamps) != 5 or any(not isinstance(s, int) for s in stamps):
        raise ValueError("one local and four emission timestamps are required")
    t = stamps[0]*tick
    eta = conformal(t, h)
    ranges = [eta-conformal(k*tick, h) for k in stamps[1:]]
    if min(ranges) < 0:
        raise ValueError("a beacon timestamp follows reception")
    point = tuple((ranges[0]**2-r*r+BEACON**2)/(2*BEACON) for r in ranges[1:])
    return t, point


def shell_radii(t, h, tick):
    """Two comoving circumferences in 2π length units, on the dust clocks' common slice."""
    a = scale(t, h)
    return tuple(int((q*a/tick).to_integral_value()) for q in (INNER, OUTER))


@dataclass(frozen=True)
class DustEmitter:
    """The prepared source exposes emitted clock and circumference records."""
    h: D

    def endpoint(self, t, point, tick):
        """The transponder on this dust line returns the actual five clock stamps."""
        return timestamp_event(t, point, self.h, tick)

    def circumferences(self, t, tick):
        return shell_radii(t, self.h, tick)

    def calibration(self, denominator):
        units = denominator**10
        with localcontext() as ctx:
            ctx.prec = max(ctx.prec, 36+2*len(str(units)))
            tick = ONE/units
            epochs = []
            for clock in (-(units//8), -(units//16), 0):
                inner, outer = self.circumferences(clock*tick, tick)
                epochs.append({"clock": clock, "inner_radius": inner, "outer_radius": outer})
        return {"ticks_per_unit": units, "epochs": epochs}


def decode_scale(record):
    """Read h and a conservative absolute error from shell and clock ticks alone.

    The two tagged radii have ratio two. Their separation divided by its
    value at the last epoch measures a. Then a(t)^(3/2) = 1 - 3*h*t/2.
    The nearer past epoch independently checks the constant inferred from
    the longer baseline. All three readings precede or coincide with t=0.
    """
    if not isinstance(record, dict) or set(record) != {"ticks_per_unit", "epochs"}:
        raise ValueError("a scale clock and three raw circumference records are required")
    units = record["ticks_per_unit"]
    if type(units) is not int or units < 32**10:
        raise ValueError("invalid scale clock units")
    epochs = record["epochs"]
    if (not isinstance(epochs, list) or len(epochs) != 3 or
            any(not isinstance(epoch, dict) or set(epoch) != set(SCALE_FIELDS) for epoch in epochs)):
        raise ValueError("three clock and circumference records are required")
    clocks = [epoch["clock"] for epoch in epochs]
    if any(type(clock) is not int for clock in clocks) or clocks != [-(units//8), -(units//16), 0]:
        raise ValueError("three registered scale epochs are required")
    widths = []
    for epoch in epochs:
        inner, outer = epoch["inner_radius"], epoch["outer_radius"]
        if (type(inner) is not int or type(outer) is not int or not 0 < inner < outer
                or abs(outer-2*inner) > 1):
            raise ValueError("the circumference stamps disagree on the two tagged shells")
        widths.append(outer-inner)
    if abs(4*epochs[-1]["inner_radius"]-units) > 2 or abs(2*epochs[-1]["outer_radius"]-units) > 1:
        raise ValueError("shell radii disagree with their registered lengths at the event")
    with localcontext() as ctx:
        ctx.prec = max(ctx.prec, 36+2*len(str(units)))
        tick = ONE/units
        estimates = []
        for clock, width in zip(clocks[:-1], widths[:-1], strict=True):
            a = D(width)/widths[-1]
            if not ONE/2 < a < D(3)/2:
                raise ValueError("shell scale is outside the calibrated range")
            estimates.append(2*(a*a.sqrt()-1)/(-3*clock*tick))
        # Each width has error <= 1 tick. On the calibrated source h in [0,1],
        # a <= (1+3/16)^(2/3) < 1.13, the ratio error is < 9*tick and
        # x^(3/2) is 2-Lipschitz here. Even the shorter baseline is > 1/17.
        error = 256*tick
        if abs(estimates[0]-estimates[1]) > 2*error:
            raise ValueError("the two shell clocks infer incompatible scale histories")
        if not -error <= estimates[0] <= 1+error:
            raise ValueError("measured collapse rate is outside the supported branch")
        return max(ZERO, estimates[0]), error


def mass_tally(mass, denominator):
    """Enclosed material rest-action tickets on two tagged shells at the preparation."""
    count_unit = 33*64*denominator**6
    inner = (count_unit*mass)//64
    outer = (count_unit*mass)//8
    return count_unit, inner, outer


def frame_clocks(frame, h, ell):
    """Measure six optical phases between two ticks of a moving receiver at the event."""
    gamma, *p = (D(x.numerator)/x.denominator for x in frame)
    tick, duration = ell**10, ell**2
    arms = []
    for component in p:
        for sign in (1, -1):
            stamps = []
            for end in (-1, 1):
                eta = conformal(end*gamma*duration, h)
                phase = eta*(1-sign*component/gamma)
                stamps.append(int((phase/tick).to_integral_value()))
            arms.append(stamps)
    rates = [(end-start)*tick/(2*duration) for start, end in arms]
    reconstructed = [sum(rates)/6, *((rates[i+1]-rates[i])/2 for i in (0, 2, 4))]
    return arms, reconstructed


def _angular(emitter, measured_h, gamma, s, ell):
    """Six positive/negative axes, with actual neighboring launches on each."""
    if not s:
        return ONE, ONE
    boost = (gamma*gamma-1).sqrt()
    tick = ell**10
    d = ell*ell
    cosine, sine = (1-d*d)/(1+d*d), 2*d/(1+d*d)

    def observed(direction):
        t, point = ray(s, measured_h, gamma, boost, direction)
        stamps = emitter.endpoint(t, point, tick)
        return locate(stamps, measured_h, tick)

    density, mass = ZERO, ZERO
    for i, axis in enumerate(AXES):
        for sign in (-1, 1):
            n = tuple(sign*x for x in axis)
            t, _ = observed(n)
            tangents = []
            for j, tangent in enumerate(AXES):
                if i == j:
                    continue
                endpoints = []
                for side in (-1, 1):
                    direction = tuple(cosine*x+side*sine*y for x, y in zip(n, tangent, strict=True))
                    time, position = observed(direction)
                    endpoints.append((time, *position))
                tangents.append(tuple((v-u)/(2*s*sine) for u, v in zip(*endpoints, strict=True)))
            inner, outer = emitter.circumferences(t, tick)
            a = (outer-inner)*tick/(OUTER-INNER)
            def metric(v, w):
                return -v[0]*w[0] + a*a*sum(v[k]*w[k] for k in range(1, 4))
            v, w = tangents
            jacobian = (metric(v, v)*metric(w, w)-metric(v, w)**2).sqrt()
            density += jacobian/6
            # Rest mass per proper volume evolves with the local volume of the dust.
            volume_ratio = (OUTER**3-INNER**3)/((outer*tick)**3-(inner*tick)**3)
            # Parallel transport along the orthogonal ray keeps a(t)*gamma(t) fixed.
            mass += jacobian*volume_ratio/(6*a**2)
    return density, mass


@dataclass(frozen=True)
class Ball:
    radius: D
    area_over_solid_angle: D
    volume_over_solid_angle: D
    mass_weight_over_solid_angle: D  # Before the measured boost contraction.
    sample_timestamps: tuple
    sample_tick: D
    measured_h: D
    scale_error_bound: D
    _scale_units: int
    _scale_epochs: tuple

    @property
    def scale_record(self):
        """A fresh JSON receipt, leaving the cached survey's integer tapes intact."""
        return {"ticks_per_unit": self._scale_units,
                "epochs": [dict(zip(SCALE_FIELDS, row, strict=True)) for row in self._scale_epochs]}


@lru_cache(maxsize=128)
def read_ball(geometric_mass, time_component, denominator):
    """Prepare a source, emit its scale record, then perform the observer survey."""
    if not isinstance(denominator, int) or denominator < 32:
        raise ValueError("the local clock survey requires a small proper ball")
    with localcontext() as ctx:
        ctx.prec = 36 + 12*len(str(denominator))
        m = D(geometric_mass.numerator)/geometric_mass.denominator
        if not 0 <= m <= 1:
            raise ValueError("compact mass outside the registered stock")
        emitter = DustEmitter((m/2).sqrt())
        return survey_ball(emitter, time_component, denominator, emitter.calibration(denominator))


def survey_ball(emitter, time_component, denominator, scale_record, *, target_radius=None):
    """Equal-volume stop and area from received records and an independent scale tape."""
    if type(denominator) is not int or denominator < 32 or scale_record["ticks_per_unit"] != denominator**10:
        raise ValueError("scale and spatial clock resolutions must agree")
    with localcontext() as ctx:
        ctx.prec = 36 + 12*len(str(denominator))
        measured_h, error = decode_scale(scale_record)
        gamma = D(time_component.numerator)/time_component.denominator
        if not 1 <= gamma <= 3:
            raise ValueError("timelike reader exceeds the calibrated velocity range")
        ell = ONE/denominator
        target = ell if target_radius is None else D(target_radius)
        if not ell <= target <= ONE/32:
            raise ValueError("choose an equal-volume target within the local clock resolution")
        quadrature = [(D(x.numerator)/x.denominator, D(w.numerator)/w.denominator) for x, w in RADIAL]

        def integrals(r):
            volume, work = ZERO, ZERO
            for x, weight in quadrature:
                jac, projected = _angular(emitter, measured_h, gamma, r*x, ell)
                volume += weight*x*x*jac
                work += weight*x*x*(1-x*x)*projected/2
            return volume, work

        r = target
        for _ in range(2):
            mean_volume, _ = integrals(r)
            r = target/cube_root(3*mean_volume)
        mean_volume, weighted = integrals(r)
        jac, _ = _angular(emitter, measured_h, gamma, r, ell)
        boost = (gamma*gamma-1).sqrt()
        t, point = ray(r, measured_h, gamma, boost, AXES[2])
        tick = ell**10
        return Ball(r, r*r*jac, r**3*mean_volume, r**4*weighted,
                    emitter.endpoint(t, point, tick), tick,
                    measured_h, error, scale_record["ticks_per_unit"],
                    tuple(tuple(epoch[k] for k in SCALE_FIELDS) for epoch in scale_record["epochs"]))
