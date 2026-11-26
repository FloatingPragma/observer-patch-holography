"""A falling sphere, its freely released probes and two stationary clocks."""
from dataclasses import dataclass
from decimal import Decimal as D, localcontext
from math import cos, log1p, sin, sqrt

from model.candidate_laws.interval_echo_partition_repair.published_strings.archive import source_densities


G = 1 / 4


def release_calibration(mass, denominator):
    """Infer G only from ordinary release radii and probe proper time on each refinement."""
    if mass <= 0 or denominator < 32:
        raise ValueError("the reference release needs positive mass and a short appointment")
    with localcontext() as ctx:
        ctx.prec = 32 + 12*len(str(denominator))
        m = D(mass.numerator)/mass.denominator
        dt = D(1)/denominator
        tick, start = dt**10, D(2)
        phase = dt*(m/(2*start**3)).sqrt()

        def sine(x):
            term = total = x
            for k in range(1, 7):
                term *= -x*x/((2*k)*(2*k+1))
                total += term
            return total

        angle = phase/2
        for _ in range(4):
            angle -= (angle+sine(2*angle)/2-phase)/(2*(1-sine(angle)**2))
        rest = int(start/tick)
        loss = int((start*sine(angle)**2/tick).to_integral_value())
        gain = 2*loss*tick*start**2/(m*dt**2)
        return {"duration": str(dt), "tick": str(tick), "start": rest, "released": rest-loss,
                "gain": float(gain), "eta": float(1/(4*gain))}


@dataclass(frozen=True)
class Body:
    ordinary: float
    collar: float
    collar_weight: float = 1.0

    def __post_init__(self):
        if min(self.ordinary, self.collar, self.collar_weight) < 0:
            raise ValueError("mass counts and their response must be nonnegative")

    @property
    def gravitating_count(self):
        return self.ordinary + self.collar_weight * self.collar

    def shell(self, label, time):
        """Areal radius of a marginally bound shell, initially at its positive label."""
        if label <= 0:
            raise ValueError("a shell has a positive initial radius")
        enclosed = self.gravitating_count * min(label**3, 1.0)
        clock = 1 - 1.5 * time * sqrt(2 * G * enclosed / label**3)
        if clock <= 0:
            raise ValueError("the survey ends before shell focusing")
        return label * clock**(2 / 3)

    def lapse(self, radius):
        if radius <= max(1, 2 * G * self.gravitating_count):
            raise ValueError("enroll stationary clocks outside the initial surface and horizon")
        return sqrt(1 - 2 * G * self.gravitating_count / radius)

    def outward_flight(self, inner, outer):
        """Time between two stationary spheres for an outgoing light signal."""
        if inner >= outer:
            raise ValueError("the receiver must be outside the emitter")
        self.lapse(inner)
        b = sqrt(2 * G * self.gravitating_count)
        difference = sqrt(outer) - sqrt(inner)
        return outer - inner + 2*b*difference + 2*b*b*log1p(difference/(sqrt(inner)-b))

    def radius_after_release(self, start, proper_time):
        """Radius of an exterior probe initially at rest in the stationary apparatus."""
        self.lapse(start)
        if self.gravitating_count == 0:
            return start
        phase = abs(proper_time) * sqrt(2 * G * self.gravitating_count / start**3)
        if phase >= 0.5:
            raise ValueError("release interval exceeds the enrolled short arc")
        angle = phase / 2
        for _ in range(6):
            angle -= (angle + sin(angle)*cos(angle) - phase)/(2*cos(angle)**2)
        radius = start * cos(angle)**2
        if radius <= 1:
            raise ValueError("the probe has reached the material surface")
        return radius

    def exterior_records(self, inner=2.0, outer=4.0):
        """Integer emission, reception and areal-radius records for the two experiments."""
        tick = 2.0**-44
        # Emit one inner-clock unit apart. Both outgoing rays solve the same null rule.
        emit_times = (0.0, 1 / self.lapse(inner))
        flight = self.outward_flight(inner, outer)
        arrivals = [round((t + flight) * self.lapse(outer) / tick) for t in emit_times]
        emissions = [round(t * self.lapse(inner) / tick) for t in emit_times]
        durations = (1/16, 1/32)
        radii = [round(self.radius_after_release(inner, dt) / tick) for dt in durations]
        start = round(inner / tick)
        return {
            "tick": tick, "emitted": emissions, "received": arrivals,
            "received_frequency_per_emitted": (emissions[1]-emissions[0])/(arrivals[1]-arrivals[0]),
            "release_proper_times": list(durations), "initial_radius_ticks": start,
            "released_radius_ticks": radii,
            "signed_radius_acceleration": [2*(r-start)*tick/dt**2 for r, dt in zip(radii, durations, strict=True)],
        }


def archive_bodies():
    return {name: Body(*densities) for name, densities in source_densities().items()}
