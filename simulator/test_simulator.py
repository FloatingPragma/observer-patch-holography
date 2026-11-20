"""Exactness and regression checks of the simulator.

Each check returns a short summary and raises on failure. ``python -m simulator
check`` prints them; pytest runs them as tests.
"""
from __future__ import annotations

import numpy as np

from .experiments import clock_journey, continuum_run
from .law import LAWS, PHI0, dispersion, step, step_bound, step_contact, step_pair, walk_matrix
from .sector import sector_step
from .spectral import Packet, PairPacket, rail_phase, track

RNG_SEED = 20261120


def _random_state(shape, rng):
    psi = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    return psi / np.linalg.norm(psi)


def check_unitarity():
    rng = np.random.default_rng(RNG_SEED)
    worst = 0.0
    for law in LAWS.values():
        psi = _random_state((5, 4, 3, 6, 8), rng)
        for _ in range(7):
            psi = step(psi, law)
        worst = max(worst, abs(np.linalg.norm(psi) - 1))
    assert worst < 1e-12, worst
    return f"norm preserved by all {len(LAWS)} members, worst drift {worst:.1e}"


def check_sector_reduction():
    rng = np.random.default_rng(RNG_SEED)
    law = LAWS["oph"]
    size, private, order = 6, 8, 3
    chi = _random_state((size, size, size, 8), rng)
    phase = np.exp(2j * np.pi * order * np.arange(private) / private)[None, None, None, :, None]
    full = chi[:, :, :, None, :] * phase
    for _ in range(9):
        chi = sector_step(chi, law, 2 * np.pi * order / private)
        full = step(full, law)
    error = float(np.max(np.abs(full - chi[:, :, :, None, :] * phase)))
    assert error < 1e-12, error
    return f"private-sector update equals the literal update, max deviation {error:.1e}"


def check_field_identity():
    rng = np.random.default_rng(RNG_SEED)
    equal = np.ones(8) / np.sqrt(8)
    worst = 0.0
    for _ in range(20):
        angles = rng.uniform(-np.pi, np.pi, 4)
        energy, walk, state = dispersion(angles), walk_matrix(angles), equal.astype(complex)
        for n in range(1, 16):
            state = walk @ state
            worst = max(worst, abs(np.vdot(equal, state) - np.cos(n * energy)))
    assert worst < 1e-12, worst
    return f"<e|W^n e> = cos(nE) on 20 random wave vectors, max deviation {worst:.1e}"


def check_encounter():
    rng = np.random.default_rng(RNG_SEED)
    carrier = (4, 3, 1, 2)
    a = _random_state(carrier + (8,), rng)
    b = _random_state(carrier + (8,), rng)
    pair = np.einsum("ijklm,nopqr->ijklmnopqr", a, b)
    free = LAWS["oph"]
    met = LAWS["control-encounter"]
    product, interacting = pair.copy(), pair.copy()
    for _ in range(5):
        a, b = step(a, free), step(b, free)
        product = step_pair(product, free)
        interacting = step_pair(interacting, met)
    independent = float(np.max(np.abs(product - np.einsum("ijklm,nopqr->ijklmnopqr", a, b))))
    differs = float(np.max(np.abs(interacting - product)))
    norm = abs(np.linalg.norm(interacting) - 1)
    assert independent < 1e-12 and differs > 1e-3 and norm < 1e-12
    return (f"zero encounter phase evolves two excitations independently ({independent:.1e}); "
            f"phase pi/2 changes the pair state by up to {differs:.3f} and keeps its norm")


def check_spectral_engine():
    """Spectral readings against the literal update on a small carrier."""
    law = LAWS["oph"]
    size, origin, sites = 8, (2, 4, 3), [(2, 4, 3), (5, 6, 7)]
    worst = 0.0
    for kind in ("optical", "clock", "reference"):
        packet = Packet(kind, (1, -1), size=size, private_size=size, carrier_order=1,
                        private_order=1, width=0.8, support=1, transverse_sign=-1, origin=origin)
        fourier = np.zeros((size,) * 4 + (8,), complex)
        fourier[np.ix_(*(o % size for o in packet.orders))] = packet.spinors
        psi = np.roll(np.fft.ifftn(fourier, axes=(0, 1, 2, 3), norm="ortho"), origin, axis=(0, 1, 2))
        for t in range(6):
            density = np.sum(np.abs(psi) ** 2, axis=(3, 4))
            for axis in range(3):
                projected = density.sum(axis=tuple(a for a in range(3) if a != axis))
                worst = max(worst, float(np.max(np.abs(projected - packet.marginal(t, axis)))))
            if kind == "clock":
                for site in sites:
                    local = np.sum(np.abs(psi[site]) ** 2, axis=-1)[:: size // 4]
                    worst = max(worst, float(np.max(np.abs(local - packet.ramsey(t, site)))))
            psi = step(psi, law)
    assert worst < 1e-12, worst
    return f"spectral marginals and four-aperture readings equal the literal update, max deviation {worst:.1e}"


def check_continuum():
    coarse = continuum_run(16, 2.0, keep_frames=False)["max_error"]
    fine = continuum_run(32, 2.0, keep_frames=False)["max_error"]
    ratio = coarse / fine
    assert 3.6 < ratio < 4.4 and fine < 0.012, (coarse, fine)
    return (f"largest deviation from Klein-Gordon {coarse:.4f} at N=16 and {fine:.4f} at N=32, "
            f"ratio {ratio:.2f} (second order)")


def check_axial_journey():
    j = clock_journey((1, 0))
    assert j["reunion"] == 1459 and abs(j["forecast"] - 0.71128) < 2e-4
    assert abs(j["beat_ratio"] - j["forecast"]) < 0.035 and abs(j["optical_ratio"] - j["forecast"]) < 0.035
    return (f"axial two-clock journey: reunion {j['reunion']}, forecast {j['forecast']:.5f}, "
            f"beat {j['beat_ratio']:.5f}, optical {j['optical_ratio']:.5f}")


def _embed_bound(chi):
    carrier = chi.shape[:4]
    volume = int(np.prod(carrier))
    psi = np.zeros((volume, 8, volume, 8), complex)
    flat = chi.reshape(volume, 8)
    for v in range(volume):
        psi[v, np.arange(8), v, np.arange(8)] = flat[v]
    return psi.reshape(carrier + (8,) + carrier + (8,))


def check_bound_pair():
    """The contact law: exact bound-sector reduction and no formation from outside it."""
    rng = np.random.default_rng(RNG_SEED)
    carrier = (3, 2, 1, 2)
    chi = _random_state(carrier + (8,), rng)
    psi, reduced = _embed_bound(chi), chi.copy()
    for _ in range(6):
        psi, reduced = step_contact(psi, PHI0), step_bound(reduced, PHI0)
    reduction = float(np.max(np.abs(psi - _embed_bound(reduced))))
    volume = int(np.prod(carrier))
    outside = _random_state(carrier + (8,) + carrier + (8,), rng).reshape(volume, 8, volume, 8)
    for v in range(volume):
        outside[v, np.arange(8), v, np.arange(8)] = 0
    outside = (outside / np.linalg.norm(outside)).reshape(carrier + (8,) + carrier + (8,))
    for _ in range(20):
        outside = step_contact(outside, PHI0)
    flat = outside.reshape(volume, 8, volume, 8)
    captured = max(float(np.max(np.abs(flat[v, np.arange(8), v, np.arange(8)]))) for v in range(volume))
    assert reduction < 1e-12 and captured == 0.0
    return (f"bound sector equals the literal two-constituent update ({reduction:.1e}); "
            f"a pair started outside it has bound amplitude {captured} after 20 updates")


def check_bound_clock():
    """The bound pair on the forward diagonal: reunion 1787 and the local home reading 0.577424."""
    law = LAWS["oph"]
    moving, rest = PairPacket(law, (1, 1)), PairPacket(law, (0, 0))
    survey = np.arange(0, 400, 20)
    velocity = [np.polyfit(survey, track(moving, a, survey), 1)[0] for a in (0, 1)]
    reunion = int(round(512 / max(abs(velocity[0]), abs(velocity[1]))))
    advance = lambda p: np.angle(np.exp(1j * (rail_phase(p.rails(reunion, (0, 0, 0)))
                                               - rail_phase(p.rails(0, (0, 0, 0))))))
    ratio = advance(moving) / advance(rest)
    spectator = LAWS["control-spectator"]
    flat = advance(PairPacket(spectator, (1, 1))) / advance(PairPacket(spectator, (0, 0)))
    assert reunion == 1787 and abs(ratio - 0.577424) < 2e-6 and abs(flat - 1) < 1e-9
    return f"bound pair on the diagonal: reunion {reunion}, rail clock {ratio:.6f}; spectator rails {flat:.6f}"


def check_formation():
    """Formation: the branch evolution against a literal update with explicit marker wires."""
    from .formation import ANGLE, COLLISION, _strip_packet, formation_run
    from .law import EQUAL
    size, wire, updates = 8, 16, 10
    phi = LAWS["oph"].rail_contact_phases[0]
    a = _strip_packet(size, size // 4, 1, support=2)
    b = _strip_packet(size, 3 * size // 4, -1, support=2)
    state = {0: np.einsum("abcde,fghij->abcdefghij", a, b)}
    c, s = np.cos(ANGLE), np.sin(ANGLE)
    blank = np.zeros_like(state[0])
    literal = []
    for n in range(updates + 1):
        literal.append(sum(float(np.linalg.norm(v) ** 2) for m, v in state.items() if m))
        if n == updates:
            break
        new = {m: v.copy() for m, v in state.items()}
        for z in range(size):
            bit = 1 << (z * wire)
            for low in {m & ~bit for m in state}:
                high = low | bit
                lo, hi = state.get(low, blank), state.get(high, blank)
                x = lo[0, 0, z, 0, COLLISION[0], 0, 0, z, 0, COLLISION[1]]
                y = np.sum(hi[0, 0, z, 0, np.arange(8), 0, 0, z, 0, np.arange(8)] * EQUAL.conj())
                if x == 0 and y == 0:
                    continue
                new.setdefault(low, blank.copy())[0, 0, z, 0, COLLISION[0], 0, 0, z, 0, COLLISION[1]] = c * x - s * y
                target = new.setdefault(high, blank.copy())
                target[0, 0, z, 0, np.arange(8), 0, 0, z, 0, np.arange(8)] += (s * x + (c - 1) * y) * EQUAL
        def delay(m):
            out = 0
            for z in range(size):
                bits = (m >> (z * wire)) & ((1 << wire) - 1)
                out |= (((bits << 1) & ((1 << wire) - 1)) | (bits >> (wire - 1))) << (z * wire)
            return out
        state = {delay(m): step_contact(v, phi) for m, v in new.items()}
    branch = formation_run(size, updates, orders=(1, -1), support=2)
    error = max(abs(x - y) for x, y in zip(literal, branch["captured"]))
    norm = sum(float(np.linalg.norm(v) ** 2) for v in state.values())
    off = formation_run(size, updates, orders=(1, -1), angle=0.0, support=2)["captured"][-1]
    assert error < 1e-12 and abs(norm - 1) < 1e-12 and literal[-1] > 0.01 and off == 0.0
    return (f"capture probability {literal[-1]:.4f} equals the literal marker-wire update to {error:.1e}; "
            f"with the vertex off, capture {off}")


def check_row30():
    """OpenWave Row 30: native frequencies, the closed form, and the observer's kinetic test."""
    from .row30 import row30_run
    r = row30_run()
    law_rows = [x for x in r["rows"] if x["law"] == "oph"]
    changed = [x for x in r["rows"] if x["law"] == "changed"]
    assert r["identity_error"] < 1e-12 and abs(r["coefficient"] - 0.24919629) < 1e-8
    assert abs(r["rest"] - np.pi / 32) < 1e-12
    assert all(x["passes"] for x in law_rows) and not any(x["passes"] for x in changed)
    worst = max(abs(x["kinetic"] - 1) for x in law_rows)
    deficit = [1 - x["kinetic"] for x in changed]
    return (f"42 of 42 probes of the law within {worst:.4f} of the smooth relation; changed law "
            f"deficit {min(deficit):.4f}-{max(deficit):.4f}; curvature {r['coefficient']:.8f}")


def check_row04():
    """OpenWave Row 4: the circulation persists at every age; independent transport loses 3/4."""
    from .row04 import row04_run
    r = row04_run(size=6, updates=6, ages=(0, 127, 128))
    assert r["profile_change"] < 1e-12 and r["capture_projection"] < 1e-12 and r["eigen_error"] < 1e-12
    assert min(r["contact_matched"]) > 1 - 1e-12
    assert abs(r["control_matched"][1] - 0.25) < 1e-12
    return ("circulation: eigenvalue -exp(i phi/2) at every age, capture projection 0, matched weight 1; "
            f"independent transport keeps {r['control_matched'][1]:.6f} after one update")


def check_app_pages():
    """The replay app gives every experiment and every OpenWave row its own page."""
    import gzip
    import json
    import tempfile
    from pathlib import Path

    from . import runs, webapp
    with tempfile.TemporaryDirectory() as tmp:
        run = Path(tmp) / "03-example.run.json.gz"
        doc = {"format": runs.FORMAT, "kind": "scan", "step": "Clocks", "title": "Example", "summary": "One run. More."}
        run.write_bytes(gzip.compress(json.dumps(doc).encode()))
        row = dict(doc, openwave_row=30, title="Row example")
        run30 = Path(tmp) / "08-openwave-row-30.run.json.gz"
        run30.write_bytes(gzip.compress(json.dumps(row).encode()))
        out = Path(tmp) / "site"
        pages = webapp.build([run, run30], out)
        rows = json.loads(webapp.ROWS.read_text())["rows"]
        expected = {"", "experiments/", "experiments/example/", "openwave/", "openwave/row-30/replay/"}
        expected |= {f"openwave/row-{r['number']}/" for r in rows}
        assert set(pages) == expected, sorted(set(pages) ^ expected)
        for page in pages:
            text = (out / page / "index.html").read_text()
            assert "{{" not in text, page
        deep = (out / "openwave/row-30/replay/index.html").read_text()
        assert 'src="../../../data/openwave-row-30.js"' in deep and '"root":"../../../"' in deep
        assert (out / "data/example.js").read_text().startswith("window.OPH_RUNS.push(")
        assert '"run":"example"' in (out / "index.html").read_text() and (out / "favicon.svg").is_file()
    return f"{len(pages)} pages"


CHECKS = [
    ("unitarity", check_unitarity),
    ("sector reduction", check_sector_reduction),
    ("field identity", check_field_identity),
    ("encounter phase", check_encounter),
    ("spectral engine", check_spectral_engine),
    ("continuum limit", check_continuum),
    ("two-clock journey", check_axial_journey),
    ("bound pair", check_bound_pair),
    ("bound clock", check_bound_clock),
    ("formation", check_formation),
    ("OpenWave row 30", check_row30),
    ("OpenWave row 4", check_row04),
    ("app pages", check_app_pages),
]


def test_unitarity():
    check_unitarity()


def test_sector_reduction():
    check_sector_reduction()


def test_field_identity():
    check_field_identity()


def test_encounter():
    check_encounter()


def test_spectral_engine():
    check_spectral_engine()


def test_continuum():
    check_continuum()


def test_axial_journey():
    check_axial_journey()


def test_bound_pair():
    check_bound_pair()


def test_bound_clock():
    check_bound_clock()


def test_formation():
    check_formation()


def test_row30():
    check_row30()


def test_row04():
    check_row04()


def test_app_pages():
    check_app_pages()
