"""Experiments on the law family: the scalar continuum and the two-clock journey."""
from __future__ import annotations

import time

import numpy as np

from .law import LAWS, Law
from .sector import continuum_field, field, plane_modes, prepare_field, sector_step
from .spectral import Packet, beat_phase, track

# The free reference speed measured on the zero-private-momentum wave of one constituent,
# in radar units per home tick. The two-clock forecasts use it.
REFERENCE_SPEED = 1 / 2.002604166666667

# The registered field menu of the continuum run: amplitude and chart wave vector K
# (cycles per chart side). It mixes axial, oblique and diagonal patterns.
CONTINUUM_MODES = (
    (0.55, (1, 0, 0)),
    (0.45, (0, 1, 1)),
    (0.40, (1, -1, 1)),
    (0.30, (2, 1, 0)),
    (0.25, (0, 0, 2)),
)
CONTINUUM_PRIVATE_ORDER = 1


def observer_stations(size: int):
    """Nine stations in chart coordinates: the observer at the centre, eight at the cube corners of
    a half-side box around it. Positions are integer addresses on a carrier of ``size``."""
    centre = size // 2
    half = size // 4
    stations = [(centre, centre, centre)]
    for dx in (-half, half):
        for dy in (-half, half):
            for dz in (-half, half):
                stations.append((centre + dx, centre + dy, centre + dz))
    return np.array(stations, int)


def continuum_run(size: int = 32, chart_time: float = 2.0, law_key: str = "torus",
                  modes=CONTINUUM_MODES, private_order: int = CONTINUUM_PRIVATE_ORDER,
                  keep_frames: bool = True):
    """Evolve the released field under the law and compare it with the Klein-Gordon limit.

    The carrier has ``size`` positions per spatial axis (chart spacing h = 1/size) and a
    private register whose sector angle is h Q, Q = 2 pi * private_order. One update is
    h chart time units. Returns the law's field at every update, the continuum field at
    the same chart times, and the largest deviation between them.
    """
    law = LAWS[law_key]
    h = 1.0 / size
    updates = int(round(chart_time * size))
    private_angle = 2 * np.pi * private_order * h
    initial = plane_modes(size, modes)
    chi = prepare_field(initial).astype(complex)
    x = np.arange(size) / size
    points = np.stack(np.meshgrid(x, x, x, indexing="ij"), axis=-1)
    stations = observer_stations(size)
    frames, station_law, errors = [], [], []
    started = time.time()
    norm = float(np.max(np.abs(initial)))
    for n in range(updates + 1):
        phi = field(chi)
        if n % max(1, size // 32) == 0 or n == updates:
            exact = continuum_field(modes, private_order, [n * h], points)[0]
            errors.append((n * h, float(np.max(np.abs(phi.real - exact))) / norm,
                           float(np.max(np.abs(phi.imag))) / norm))
        if keep_frames:
            frames.append(phi.real.astype(np.float32))
        station_law.append([float(phi.real[tuple(s)]) for s in stations])
        if n < updates:
            chi = sector_step(chi, law, private_angle)
    station_continuum = continuum_field(modes, private_order, np.arange(updates + 1) * h,
                                        stations / size)
    return {
        "law": law, "size": size, "updates": updates, "chart_time": chart_time,
        "private_order": private_order, "modes": modes, "norm": norm,
        "frames": frames, "stations": stations,
        "station_law": np.array(station_law), "station_continuum": station_continuum,
        "errors": errors, "max_error": max(e for _, e, _ in errors),
        "max_imaginary": max(i for _, _, i in errors),
        "seconds": time.time() - started,
    }


def clock_journey(route=(1, 0), size: int = 512, private_size: int = 4096,
                  reference_speed: float = REFERENCE_SPEED, origin=(0, 0, 0)):
    """The two-clock journey: beat and optical clocks on one route, read at home at reunion.

    Speed comes from the beat clock's own envelope track, the reunion time from that speed,
    and the forecast from the flat metric with the archived reference speed. Each clock's
    moving advance is divided by its own stationary twin's advance over the same home ticks.
    """
    kwargs = dict(size=size, private_size=private_size, origin=origin)
    moving = Packet("clock", route, **kwargs)
    times = np.arange(0, 400, 20)
    velocity = [np.polyfit(times, track(moving, axis, times), 1)[0] for axis in (0, 1)]
    speed = float(np.hypot(*velocity))
    lead = max(abs(velocity[0]), abs(velocity[1]))
    reunion = int(round(size / lead))
    forecast = float(np.sqrt(1 - (speed / reference_speed) ** 2))

    def beat_advance(packet):
        start = beat_phase(packet.ramsey(0, origin))
        end = beat_phase(packet.ramsey(reunion, origin))
        return float(np.angle(np.exp(1j * (end - start))))

    rest = Packet("clock", (0, 0), **kwargs)
    beat = beat_advance(moving) / beat_advance(rest)

    def optical_turns(packet):
        ticks = np.arange(0, reunion + 1, 4)
        if ticks[-1] != reunion:
            ticks = np.append(ticks, reunion)
        z = track(packet, 2, ticks)
        return float((z[-1] - z[0]) / size)

    optical_moving = Packet("optical", route, **kwargs)
    optical_rest = Packet("optical", (0, 0), **kwargs)
    optical = optical_turns(optical_moving) / optical_turns(optical_rest)
    return {"route": list(route), "speed": speed, "reunion": reunion, "forecast": forecast,
            "beat_ratio": beat, "optical_ratio": optical}


# ----------------------------------------------------------------------------
# Runs for the replay app: the successful experiments of steps 1 to 4.

def _xy_density(packet, t: float, stride: int) -> np.ndarray:
    """Probability over (x, y) sampled every ``stride`` addresses, traced over z, private and slots."""
    from .spectral import TAU
    n = packet.size
    xs = np.arange(0, n, stride) - packet.origin[0]
    ys = np.arange(0, n, stride) - packet.origin[1]
    ex = np.exp(1j * TAU * np.outer(xs, packet.orders[0]) / n)
    ey = np.exp(1j * TAU * np.outer(ys, packet.orders[1]) / n)
    amplitude = np.einsum("ai,bj,ijkps->abkps", ex, ey, packet.state(t))
    return np.sum(np.abs(amplitude) ** 2, axis=(2, 3, 4))


def _velocity(packet, survey=np.arange(0, 400, 20)):
    return [float(np.polyfit(survey, track(packet, axis, survey), 1)[0]) for axis in (0, 1)]


def _clock_packets(clock, law_key, route, size):
    from .spectral import PairPacket
    if clock == "rails":
        law = LAWS[law_key]
        return PairPacket(law, route, size=size), PairPacket(law, (0, 0), size=size)
    kind = "clock" if clock == "beat" else "optical"
    return (Packet(kind, route, size=size, private_size=8 * size),
            Packet(kind, (0, 0), size=size, private_size=8 * size))


def _reading(clock, packet, t, site):
    from .spectral import rail_phase
    if clock == "beat":
        return beat_phase(packet.ramsey(t, site))
    if clock == "rails":
        return rail_phase(packet.rails(t, site))
    raise ValueError(clock)


def journey_run(clocks=("beat",), law_key="torus", route=(1, 1), size=512,
                frames=180, stride=4):
    """A traveller and its stationary twin, every clock's advance, and the forecast.

    The density frames and the traveller's running reading are audit views: the
    traveller's hand is read at the address of its envelope's centre. The reunion
    ratios use the home instrument only, at departure and at reunion.
    """
    started = time.time()
    traveller, twin = _clock_packets(clocks[0], law_key, route, size)
    velocity = _velocity(traveller)
    speed = float(np.hypot(*velocity))
    reunion = int(round(size / max(abs(velocity[0]), abs(velocity[1]))))
    forecast = float(np.sqrt(1 - (speed / REFERENCE_SPEED) ** 2))
    ticks = np.unique(np.round(np.linspace(0, reunion, frames)).astype(int))
    density = np.stack([[_xy_density(p, t, stride) for p in (traveller, twin)] for t in ticks])
    centre = np.stack([track(traveller, axis, ticks) for axis in (0, 1)], axis=-1)
    home = (0, 0, 0)
    series = []
    for clock in clocks:
        moving, rest = (traveller, twin) if clock == clocks[0] else _clock_packets(clock, law_key, route, size)
        if clock == "optical":
            dense = np.arange(0, reunion + 1, 2)
            turns = lambda p: np.interp(ticks, dense, (track(p, 2, dense) - track(p, 2, [0])[0]) / size)
            rest_series, moving_series = turns(rest), turns(moving)
            ratio = float(moving_series[-1] / rest_series[-1])
            unit = "turns"
        else:
            def unwrap(packet, sites):
                raw = np.array([_reading(clock, packet, t, s) for t, s in zip(ticks, sites)])
                return np.unwrap(raw) - raw[0]
            moving_sites = [(int(round(x)) % size, int(round(y)) % size, 0) for x, y in centre]
            rest_series = unwrap(rest, [home] * len(ticks))
            moving_series = unwrap(moving, moving_sites)
            advance = lambda p: np.angle(np.exp(1j * (_reading(clock, p, reunion, home) - _reading(clock, p, 0, home))))
            ratio = float(advance(moving) / advance(rest))
            unit = "radians"
        series.append({"clock": clock, "unit": unit, "rest": rest_series.tolist(),
                       "moving": moving_series.tolist(), "ratio": ratio})
    return {"law": LAWS[law_key], "size": size, "route": list(route), "stride": stride,
            "ticks": ticks.tolist(), "density": density, "centre": (centre % size).tolist(),
            "speed": speed, "reunion": reunion, "forecast": forecast, "clocks": series,
            "seconds": time.time() - started}


def race_run(size: int = 1024, source: int = 256, distance: int = 128, stride: int = 2):
    """Field-only signalling: a compact free-wave launch at A read by a detector at B.

    The carrier is a strip of ``size`` addresses (the other axes have one position).
    The launch is one address in the equal-slot state. The law's strict front moves one
    address per update; the smooth cone moves at c = 1/2. B's detector probability is
    computed exactly at every update, precursors ahead of the smooth cone included.
    """
    from .law import step
    started = time.time()
    law = LAWS["torus"]
    psi = np.zeros((size, 1, 1, 1, 8), complex)
    psi[source, 0, 0, 0] = 1 / np.sqrt(8)
    late = 4 * distance
    receiver = source + distance
    detector, strip = [], []
    for n in range(late + 1):
        density = np.sum(np.abs(psi[:, 0, 0, 0]) ** 2, axis=-1)
        detector.append(float(density[receiver]))
        if n % stride == 0:
            strip.append(density)
        psi = step(psi, law)
    detector = np.array(detector)
    early = int(1.5 * distance)
    first = int(np.argmax(detector > 0)) if np.any(detector > 0) else None
    return {"law": law, "size": size, "source": source, "receiver": receiver, "distance": distance,
            "stride": stride, "strip": np.array(strip), "detector": detector,
            "deadlines": {"strict_front": distance, "early": early, "smooth": 2 * distance, "late": late},
            "first_nonzero": first,
            "largest_before_early": float(detector[: early + 1].max()),
            "largest_before_smooth": float(detector[: 2 * distance].max()),
            "received_by_late": float(detector[: late + 1].sum()),
            "seconds": time.time() - started}


def dispersion_scan(orders=(4, 8, 12, 16, 20, 24, 32, 40, 48), size=512):
    """The moving/rest beat ratio against v/c, from long to short wavelengths.

    The private register has 16 positions per spatial address and the doublet sits at
    orders P/32 and P/32 + 1, the private wave number of the two-clock clock with half its
    beat frequency, so the rest clock turns less than pi before the slowest reunion.
    """
    started = time.time()
    points = []
    private = 16 * size
    for order in orders:
        kwargs = dict(size=size, private_size=private, private_order=private // 32, carrier_order=order)
        moving = Packet("clock", (1, 0), **kwargs)
        rest = Packet("clock", (0, 0), **kwargs)
        vx = _velocity(moving)[0]
        reunion = int(round(size / abs(vx)))
        advance = lambda p: np.angle(np.exp(1j * (beat_phase(p.ramsey(reunion, (0, 0, 0)))
                                                   - beat_phase(p.ramsey(0, (0, 0, 0))))))
        points.append({"order": order, "v_over_c": abs(vx) / REFERENCE_SPEED,
                       "ratio": float(advance(moving) / advance(rest)),
                       "forecast": float(np.sqrt(max(0.0, 1 - (vx / REFERENCE_SPEED) ** 2)))})
    return {"points": points, "seconds": time.time() - started}


# Moving/rest rates the observers recorded on narrow-band diagonal flights of the bound pair,
# with the central orders of their packets.
STEP4_RECORDED = {12: 0.686818, 16: 0.577615, 22: 0.455911}


def pair_correction_scan(orders=(12, 16, 22), size=512):
    """Smooth forecast, finite law at the central wavelength, and the full packet."""
    from .spectral import PairPacket, pair_band
    started = time.time()
    law = LAWS["torus"]
    phi0, phi1 = law.rail_contact_phases
    points = []
    for order in orders:
        k = 2 * np.pi * order / size
        angles = np.array([k, k, 0.0, 0.0])
        e0, _ = pair_band(angles, phi0)
        e1, _ = pair_band(angles, phi1)
        central = float((e1 - e0) / ((phi1 - phi0) / 2))
        moving = PairPacket(law, (1, 1), size=size, carrier_order=order)
        rest = PairPacket(law, (0, 0), size=size, carrier_order=order)
        velocity = _velocity(moving)
        smooth = float(np.sqrt(1 - (np.hypot(*velocity) / REFERENCE_SPEED) ** 2))
        reunion = int(round(size / max(abs(velocity[0]), abs(velocity[1]))))
        from .spectral import rail_phase
        advance = lambda p: np.angle(np.exp(1j * (rail_phase(p.rails(reunion, (0, 0, 0)))
                                                   - rail_phase(p.rails(0, (0, 0, 0))))))
        points.append({"order": order, "smooth": smooth, "central": central,
                       "packet": float(advance(moving) / advance(rest)),
                       "recorded": STEP4_RECORDED.get(order)})
    return {"points": points, "seconds": time.time() - started}
