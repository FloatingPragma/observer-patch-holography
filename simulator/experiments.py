"""Experiments on the law family: the scalar continuum and the two-clock journey."""
from __future__ import annotations

import time

import numpy as np

from .law import LAWS, Law
from .sector import continuum_field, field, plane_modes, prepare_field, sector_step
from .spectral import Packet, beat_phase, track

# The free reference speed measured on the zero-private-momentum wave in step 1,
# in radar units per home tick. The step-3 forecasts use it.
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


def continuum_run(size: int = 32, chart_time: float = 2.0, law_key: str = "reflected-routing",
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
    """The step-3 journey: beat and optical clocks on one route, read at home at reunion.

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
