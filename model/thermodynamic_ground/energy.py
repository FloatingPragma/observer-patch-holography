"""Ground energy density with explicit finite-box surface brackets."""
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path


def _side(value):
    if type(value) is not int or value<2 or value%2:
        raise ValueError('positive even cube side required')
    return value


def intertile_bonds(large_side,tile_side):
    large_side,tile_side=_side(large_side),_side(tile_side)
    if large_side%tile_side:raise ValueError('large side must be a multiple of tile side')
    return 3*(large_side//tile_side-1)*large_side**2


def finite_energy_bracket(side,energy_lower,energy_upper,speed=F(1,8)):
    """Conditional on supplied certified FULL-ground energy interval, not a gap."""
    side=_side(side)
    if (any(type(x) not in (int,F) for x in (energy_lower,energy_upper,speed)) or
            energy_lower>energy_upper or speed<0):
        raise ValueError('ordered exact ground-energy interval and nonnegative speed required')
    return F(energy_lower)/side**3-3*F(speed)/side,F(energy_upper)/side**3


def remainder_norm_bound(large_side,tile_side,speed=F(1,8),interaction=F(1)):
    """Conservative norm of all terms touching untiled boundary sites."""
    large_side,tile_side=_side(large_side),_side(tile_side)
    if large_side<tile_side or any(type(x) not in (int,F) or x<0 for x in (speed,interaction)):
        raise ValueError('ordered cube sides and exact nonnegative coefficients required')
    filled=(large_side//tile_side)*tile_side
    remainder=large_side**3-filled**3
    # At most six incident spatial bonds per remainder site, including
    # boundaries to the tiled region; each contact norm is U.
    return (6*F(speed)+F(interaction))*remainder


def certificate():
    root=Path(__file__).resolve().parents[2]
    path='model/interacting_vacuum/results.json'
    prior=json.loads((root/path).read_text())
    if (prior['all_pass'] is not True or prior['global_ground']['global_ground_energy']!='-sqrt(3)/4' or
            prior['global_ground']['full_Fock_dimension']!=65536):
        raise ArithmeticError('finite energy example lost its full-ground certificate')
    for name,digest in prior['source_bindings'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:
            raise ArithmeticError('finite energy input source binding is stale: '+name)
    scale=1<<80
    k=isqrt(3*scale*scale)
    sqrt_lo,sqrt_hi=F(k,scale),F(k+1,scale)
    if not sqrt_lo**2<=3<=sqrt_hi**2:raise ArithmeticError('exact square-root enclosure failed')
    bracket=finite_energy_bracket(2,-sqrt_hi/4,-sqrt_lo/4)
    return {'domain':'fixedc>0,U>0 originalmu0; full-Fock grounds on even-sided open cubes',
            'energy_density_exists':True,
            'definition':'e_inf=lim_(even m->infinity) E_m/m^3',
            'finite_box_bracket':'E_ell/ell^3-3c/ell<=e_inf<=E_ell/ell^3 for every evenell',
            'exact_intertile_count':'3(m/ell-1)m^2 for m divisible byell',
            'untiled_remainder_bound':'||(terms touching remainder)||<=(6c+U)[m^3-(ell floor(m/ell))^3]',
            'existence_argument':'For each fixed evenell the surface remainder vanishes and limsup-liminf<=3c/ell. Takingell->infinity forces equality.',
            'constructed_ground_local_energy':'omega(h_0)=e_inf, h_0=UD_0+one_half sum of incident hopping bonds',
            'finite_average_identity':'nu_L(h_0)=E_(2L)/(2L)^3 exactly because exterior vacuum has zero cross-boundary hopping expectation',
            'global_bounds':'-sqrt3*c<=e_inf<=-sqrt3*c/4',
            'finite_example':{'side':2,'c':'1/8','U':'1','full_ground_energy':'-sqrt3/4',
                              'source_receipt':path,'all_source_hashes_verified':True,
                              'resulting_lower_exact':str(bracket[0]),'resulting_upper_exact':str(bracket[1]),
                              'sharper_universal_lower_exact':str(-sqrt_hi/8)},
            'finite_ground_solver_run_for_growing_boxes':False,
            'bracket_requires_certified_finite_ground_input':True,
            'state_convergence_rate_or_particle_rest_mass':False}
