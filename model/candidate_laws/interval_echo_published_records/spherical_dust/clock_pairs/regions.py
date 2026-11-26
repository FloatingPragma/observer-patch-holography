"""Regional clock-pair entropy on a finite lattice covering the pressureless interior."""
import argparse
from dataclasses import dataclass
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from .. import survey
from ..body import archive_bodies, release_calibration
from ..experiment import FRAME, NAMES, measured_frame
from .channel import ClockPair, LagGrid


HERE = Path(__file__).resolve().parent
COEFFICIENT = 4*math.pi/15
SCALE_MAX = 2
CENTERS = ((0., 0., 0.), (1/8, -1/16, 1/32))
DIAGONAL = (F(2), F(1), F(1), F(1))


@dataclass(frozen=True)
class Band:
    """All settings depend only on the permanent separation band."""
    level: int

    @property
    def width(self):
        return 2.**-self.level

    @property
    def mesh(self):
        return self.width**2/8

    @property
    def radial_count(self):
        return 4 << self.level

    @property
    def angular_count(self):
        # Equal increments in cos(theta) give polar diameter O(1/sqrt(nz)).
        return 8 << (2*self.level)

    def source(self, mass):
        direction_count = 2*self.angular_count**2
        width = F(1, 1 << self.level)
        step = width/(2*self.radial_count)
        weight = F(2, 5)*F(mass)*(width**2/8)**3*step/direction_count
        bits = self.level+4
        return ClockPair(LagGrid(bits, (1 << bits)-1), weight, SCALE_MAX)


@lru_cache(maxsize=128)
def label_samples(level, power):
    """Fixed quadrature sample of origins, directions and separation labels."""
    return qmc.Sobol(6, scramble=True, seed=1802+level).random_base2(power)


def material_clock(points, center, velocity, a0, h, cutoff):
    """Invert the spatial ray map to schedule local clock and proper-distance reads."""
    gamma = velocity[0]
    spatial = velocity[1:]
    norm = np.linalg.norm(spatial)
    axis = spatial/norm if norm else np.zeros(3)
    x = (points-center)*a0
    y = x + (x @ axis)[:, None]*(1/gamma-1)*axis
    z = np.linalg.norm(y, axis=1)
    w = np.divide(y @ spatial, z, out=np.zeros_like(z), where=z>0)

    def time(distance):
        hz = h*distance
        return distance*(w+hz*(0.5+hz*(-w/4+hz*((w*w-1)/16+hz*w*(9-w*w)/160))))

    eta = time(z)
    a = (1-h*eta/2)**2
    midpoint = 1-h*time(z/2)/2
    # ds/dz=a². Simpson integration and the fifth-degree emitter have
    # uniform absolute errors O(z⁵) and O(z⁶) on the surveyed compact set.
    distance = z*(1+4*midpoint**4+a*a)/6
    return distance <= cutoff, eta/a0, a0*a


def _band_read(band, center, velocity, a0, h, epoch, radius, mass, power, timed):
    gamma, spatial = velocity[0], velocity[1:]
    norm = np.linalg.norm(spatial)
    axis = spatial/norm if norm else np.zeros(3)
    # This box contains the projected ball and the entire origin voxel of
    # every admitted address on the finite fixtures (verify the stated slack).
    half = radius*1.2/a0 + 2*band.mesh
    if band.width/2 > 2*gamma*radius*1.2/a0:
        return 0., 0
    samples = label_samples(band.level, power)
    y = (2*samples[:, :3]-1)*half
    origins = center+y+(y @ axis)[:, None]*(gamma-1)*axis
    origins = np.rint(origins/band.mesh)*band.mesh
    nz = band.angular_count
    mu = (np.floor(samples[:, 3]*nz)+0.5)*2/nz-1
    phi = (np.floor(samples[:, 4]*2*nz)+0.5)*math.pi/nz
    transverse = np.sqrt(1-mu*mu)
    directions = np.column_stack((transverse*np.cos(phi), transverse*np.sin(phi), mu))
    d = band.width/2*(1+(np.floor(samples[:, 5]*band.radial_count)+0.5)/band.radial_count)
    ia, eta_a, scale_a = material_clock(origins, center, velocity, a0, h, radius)
    ib, eta_b, scale_b = material_clock(origins+d[:, None]*directions, center, velocity, a0, h, radius)
    enrolled = ia & ib
    count = int(np.count_nonzero(enrolled))
    if not count or not mass:
        return 0., count
    pair = band.source(mass)
    eta_a, eta_b, d = eta_a[enrolled]+epoch, eta_b[enrolled]+epoch, d[enrolled]
    units = pair.grid.bins*pair.grid.subphases
    tick_a = np.floor(eta_a*units/d).astype(np.int64)
    tick_b = np.floor(eta_b*units/d).astype(np.int64)
    scale_tick = band.width/128
    scale_a = np.rint(scale_a[enrolled]/scale_tick)*scale_tick
    scale_b = np.rint(scale_b[enrolled]/scale_tick)*scale_tick
    lag = tick_b-tick_a if timed else np.zeros_like(tick_a)
    changes = pair.entropy_from_tick_differences(lag, scale_a, scale_b)
    scaled_read = float(changes.sum()/len(samples)/float(pair.weight))
    # Convert normalized label sampling into the entropy sum of its population.
    return (2*mass/5)*scaled_read*(2*half)**3*gamma*(band.width/2), count


@lru_cache(maxsize=128)
def geometry(mass, frame, denominator, epoch_time):
    with localcontext() as ctx:
        ctx.prec = 65
        ell, t0 = D(1)/denominator, D(str(epoch_time))
        geometric = D(mass.numerator)/mass.denominator
        source_h = (geometric/2).sqrt()
        emitter = survey.DustEmitter(source_h)
        # Locate the time origin and the change in two circumferences independently.
        scale_tape = emitter.calibration(denominator)
        read_h, _ = survey.decode_scale(scale_tape)
        unit = ell**10
        inner, outer = emitter.circumferences(t0, unit)
        a0 = (outer-inner)*unit/(survey.OUTER-survey.INNER)
        local_source = survey.DustEmitter(local_source_h(source_h, t0))
        local_tape = local_source.calibration(denominator)
        h, _ = survey.decode_scale(local_tape)
        ball = survey.survey_ball(local_source, frame[0], denominator, local_tape)
        flat = survey.read_ball(F(0), frame[0], denominator)
        phase_stamps, _ = survey.frame_clocks(frame, local_source_h(source_h, t0), ell)
        observed = np.array(list(map(float, measured_frame(phase_stamps, denominator))))
        area = 4*math.pi*float((ball.area_over_solid_angle-flat.area_over_solid_angle)/ell**4)
        return dict(a0=float(a0), h=float(h), epoch=float(survey.conformal(t0, read_h)),
                    radius=float(ball.radius), area=area, velocity=observed,
                    phase_stamps=phase_stamps, shell_width_ticks=int(outer-inner))


def local_source_h(h, t):
    return h/(1-3*h*t/2)


def read_region(stock, denominator, frame, label, center, *, epoch_time=0., power=12,
                timed=True, material_collar=True, gravity_collar=F(1)):
    body = archive_bodies()[stock]
    ordinary, collar = (F(x).limit_denominator(33) for x in (body.ordinary, body.collar))
    tickets = [survey.mass_tally(value, denominator) for value in (ordinary, collar)]
    tallies = [F(64*(hi-lo), 7*n) for n, lo, hi in tickets]
    source = sum(tallies) if material_collar else tallies[0]
    meter_mass = float(sum(tallies))
    g = geometry(ordinary+gravity_collar*collar, tuple(frame), denominator, epoch_time)
    ell = 1/denominator
    total, count = 0., 0
    # All longer bands are present too. Their joint overlap vanishes.
    depth = math.ceil(math.log2(8*denominator**2))
    for j in range(1, depth+1):
        s, n = _band_read(Band(j), np.asarray(center), g["velocity"], g["a0"], g["h"],
                          g["epoch"], g["radius"], source, power, timed)
        total += s
        count += n
    eta = release_calibration(F(1, 3), denominator)["eta"]
    stress = 3*meter_mass*g["velocity"][0]**2/(4*math.pi*g["a0"]**3)
    entropy = total/ell**4
    return dict(stock=stock, inverse_radius=denominator, frame=label, center=list(center),
                proper_time=epoch_time, mounted_radius=g["radius"], measured_scale=g["a0"],
                phase_stamps=g["phase_stamps"], shell_width_ticks=g["shell_width_ticks"],
                mass_tickets=tickets, finite_bands=depth, label_sample_power=power,
                sampled_joint_pairs=count, entropy_over_r4=entropy, area_over_r4=g["area"],
                stress=stress, bulk_over_r4=2*math.pi*COEFFICIENT*stress,
                balance_over_r4=entropy+eta*g["area"])


def measure(power=14):
    rows = [read_region("held", 128, frame, label, CENTERS[i%2], power=power)
            for i, (frame, label) in enumerate(zip(FRAME, NAMES, strict=True))]
    trials = {}
    # Same nonsymmetric center, unconstrained radius, and velocity withheld
    # from the ten-direction reconstruction. All modes use the same clocks.
    for arm in ("held", "ordinary", "blank", "frozen_pair_clock", "blind_collar_work", "half_collar_force"):
        trials[arm] = read_region(arm if arm in ("blank", "ordinary") else "held", 96,
                                  DIAGONAL, "body diagonal", CENTERS[1], epoch_time=1/32,
                                  power=power, timed=arm!="frozen_pair_clock",
                                  material_collar=arm!="blind_collar_work",
                                  gravity_collar=F(1, 2) if arm=="half_collar_force" else F(1))
    refinement = [read_region("held", den, DIAGONAL, "body diagonal", CENTERS[1],
                              epoch_time=1/32, power=power) for den in (64, 192)]
    return dict(schema="oph.dust_clock_pair_regions.v1", sampled_population=True,
                ten_frames=rows, held_out=trials, held_out_refinement=refinement)


def verify():
    data = measure()
    saved = json.loads((HERE/"regions.json").read_text())
    if data.keys()!=saved.keys() or data["held_out"].keys()!=saved["held_out"].keys():
        raise ValueError("changed clock-pair aperture menu")
    def records(packet):
        return packet["ten_frames"]+list(packet["held_out"].values())+packet["held_out_refinement"]
    for actual, recorded in zip(records(data), records(saved), strict=True):
        for key, value in actual.items():
            old = recorded[key]
            if isinstance(value, float):
                if not math.isclose(value, old, rel_tol=1e-7, abs_tol=1e-8):
                    raise ValueError(f"changed regional reading at {actual['frame']}/{key}")
            elif key == "mass_tickets":
                if [list(t) for t in value] != old:
                    raise ValueError("changed regional mass counter")
            elif value != old:
                raise ValueError(f"changed finite regional record {key}")
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    if parser.parse_args().write:
        data = measure()
        (HERE/"regions.json").write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    else:
        verify()
