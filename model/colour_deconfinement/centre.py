"""Exact statements on the centre: the Z6 classes on the five-level code, the colour factor of the plaquette
action, and the action of a centre transformation on links, plaquettes and Polyakov loops."""
from fractions import Fraction

import numpy as np

from . import su3

# hypercharge on V = C^3 + C^2 (three colour levels, two weak levels)
Y_COLOUR = Fraction(-1, 3)
Y_WEAK = Fraction(1, 2)


def central_element(f):
    """Phases, in turns, of exp(2 pi i f Y) on the colour levels and on the weak levels."""
    return (f * Y_COLOUR) % 1, (f * Y_WEAK) % 1


def flux_table():
    """For each flux class f: its action on colour and weak levels and the classes k whose loop it rotates."""
    rows = {}
    for f in range(6):
        pc, pw = central_element(f)
        rows[str(f)] = dict(colour_phase_turns=str(pc), weak_phase_turns=str(pw),
                            acts_on_colour=pc != 0, acts_on_weak=pw != 0,
                            classes_rotated=[k for k in range(6) if Fraction(f * k, 6) % 1 != 0],
                            loop_phase_turns={str(k): str(Fraction(f * k, 6) % 1) for k in range(6)})
    return rows


def five_trace_cost(V3=None, V2=None, phase=0.0):
    """Seam cost 1 - Re tr_5(U)/5 of the holonomy U = exp(i phase Y) (V3 + V2) on the five-level code."""
    V3 = np.eye(3) if V3 is None else V3
    V2 = np.eye(2) if V2 is None else V2
    t = np.exp(1j * phase * float(Y_COLOUR)) * np.trace(V3) + np.exp(1j * phase * float(Y_WEAK)) * np.trace(V2)
    return 1.0 - t.real / 5.0


def colour_factor_identity(samples=200, seed=1):
    """Largest deviation of 1 - Re tr_5 diag(V, 1)/5 from (3/5)(1 - Re tr V / 3) over random V in SU(3)."""
    rng = np.random.default_rng(seed)
    V = su3.random_su3(samples, rng)
    lhs = np.array([five_trace_cost(V3=v) for v in V])
    rhs = 0.6 * (1.0 - np.trace(V, axis1=1, axis2=2).real / 3.0)
    return float(np.abs(lhs - rhs).max())


def thin_flux_cost():
    """Seam cost of a plaquette whose holonomy is the central element of class f, as exact fractions."""
    out = {}
    for f in range(6):
        pc, pw = central_element(f)
        cos = {Fraction(0): Fraction(1), Fraction(1, 2): Fraction(-1), Fraction(1, 3): Fraction(-1, 2),
               Fraction(2, 3): Fraction(-1, 2)}
        out[str(f)] = str(1 - (3 * cos[pc] + 2 * cos[pw]) / 5)
    return out


def lattice_checks(seed=7):
    """Centre transformation on a thermalised 4^4 lattice. Returns the deviations of each exact statement."""
    lat = su3.Lattice((4, 4, 4, 4), 5.7, seed=seed, start='hot')
    lat.compound_sweeps(2, 1)
    p0, l0, f0 = lat.plaquette(), lat.polyakov(), lat.polyakov_field().copy()
    out = {}
    g = su3.random_su3(lat.volume, lat.rng).reshape(lat.dims + (3, 3))
    lat.gauge_transform(g)
    out['gauge_invariance'] = max(abs(lat.plaquette() - p0), abs(lat.polyakov() - l0))
    lat.centre_transform(t_slice=2, k=1)
    out['global_plaquette_change'] = abs(lat.plaquette() - p0)
    out['global_loop_mismatch'] = abs(lat.polyakov() - su3.OMEGA * l0)
    lat.centre_transform(t_slice=2, k=2)                       # back
    # spatial links of one slice: the temporal loop is blind to it
    lat.U[0][1] *= su3.OMEGA
    out['spatial_transformation_loop_change'] = abs(lat.polyakov() - l0)
    out['spatial_transformation_plaquette_change'] = abs(lat.plaquette() - p0)
    lat.U[0][1] *= su3.OMEGA ** 2
    # region: loops inside rotate, loops outside stay
    region = np.zeros(lat.dims[:3], bool)
    region[1:3, 1:3, 1:3] = True
    lat.centre_transform(t_slice=0, k=1, region=region)
    f1 = lat.polyakov_field()
    out['region_inside_mismatch'] = float(np.abs(f1[region] - su3.OMEGA * f0[region]).max())
    out['region_outside_change'] = float(np.abs(f1[~region] - f0[~region]).max())
    return {k: float(v) for k, v in out.items()}


def vortex_checks():
    """A region transformation on the cold lattice: the flux it creates and its closure."""
    lat = su3.Lattice((6, 6, 6, 4), 1.0, start='cold')
    region = np.zeros(lat.dims[:3], bool)
    region[1:4, 2:4, 1:5] = True                               # a 3 x 2 x 4 block
    lat.centre_transform(t_slice=1, k=1, region=region)
    flux = {(mu, nu): lat.centre_flux(mu, nu) for mu in range(4) for nu in range(mu + 1, 4)}
    all_centre = all((f >= 0).all() for f in flux.values())
    spatial = sum(int((flux[p] != 0).sum()) for p in ((0, 1), (0, 2), (1, 2)))
    temporal = {i: flux[(i, 3)] for i in range(3)}
    count = sum(int((f != 0).sum()) for f in temporal.values())
    only_slice = all(int((f[:, :, :, [0, 2, 3]] != 0).sum()) == 0 for f in temporal.values())
    surface = 0
    for i in range(3):
        surface += int((region != np.roll(region, -1, axis=i)).sum())
    # closure: on every cube the six plaquette classes, with orientation, add to zero modulo 3
    def n(mu, nu):
        return flux[(mu, nu)] if mu < nu else (-flux[(nu, mu)]) % 3
    worst = 0
    for mu, nu, rho in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
        s = (np.roll(n(mu, nu), -1, axis=rho) - n(mu, nu) + np.roll(n(nu, rho), -1, axis=mu) - n(nu, rho)
             + np.roll(n(rho, mu), -1, axis=nu) - n(rho, mu)) % 3
        worst = max(worst, int(s.max()))
    loops = lat.polyakov_field()
    return dict(all_plaquettes_central=bool(all_centre), flux_plaquettes=count, region_surface_links=surface,
                spatial_flux_plaquettes=spatial, flux_only_on_the_slice=bool(only_slice), cube_closure_violation=worst,
                loops_inside=[float(np.abs(loops[region] - su3.OMEGA).max())],
                loops_outside=[float(np.abs(loops[~region] - 1.0).max())])
