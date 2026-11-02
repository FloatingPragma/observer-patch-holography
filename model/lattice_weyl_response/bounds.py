"""Exact rational bounds; all lengths are relative to shell radius L."""
from fractions import Fraction as F
from functools import lru_cache
from .algebra import sobolev
from model.weyl_gravity_response.bounds import GAMMA_PI3, BALL_ENTRY_ERROR, INVERSE_NORM


def rational(value):
    if type(value) not in (int,F):
        raise ValueError('an exact rational is required')
    return F(value)


@lru_cache(maxsize=1)
def constants():
    k1,k2,k3,k4=sobolev()['H1_H2_H3_H4_integer_majorants']
    k0=k3+7
    kd=F(k4+7*k1)+F(k3,6)
    A=16*k1*k0+24*kd+48*k2+144*k1**3
    lattice=60/GAMMA_PI3*(A/2+16*k1)
    return {'K1':k1,'K2':k2,'K3':k3,'K4':k4,'K0':k0,'KD':kd,
            'stress_A':A,'lattice_entry':lattice}


def admissible_spacing(radius_ratio):
    rho=rational(radius_ratio)
    if not 0<rho<=F(1,2**28):
        raise ValueError('0<R/L<=2^-28 required for the certified tomography family')
    c=constants()
    return min(F(1,16),F(1,6*c['K1']),rho**7/(1+c['lattice_entry']))


def errors(spacing_ratio,radius_ratio):
    q,rho=rational(spacing_ratio),rational(radius_ratio);c=constants()
    if not 0<rho<=F(1,8) or not 0<q<=min(F(1,16),F(1,6*c['K1'])):
        raise ValueError('positive radius and spacing within the local normalized domain required')
    projection=q*c['K0'];derivative=q*c['KD']
    lattice=c['lattice_entry']*q/rho**4
    entry=lattice+BALL_ENTRY_ERROR*rho**2
    return {'projection_norm_upper':projection,'L_times_derivative_norm_upper':derivative,
            'normalization_defect_upper':9*q*q*c['K1']**2,
            'normalized_lattice_entry_upper':lattice,
            'normalized_total_entry_upper':entry,
            'Neumann_budget':9*INVERSE_NORM*entry}


def certificate():
    c=constants();rho=F(1,2**28);q=admissible_spacing(rho);error=errors(q,rho)
    if error['Neumann_budget']>=F(1,4):
        raise ArithmeticError('the joint lattice/radius family is not invertible')
    return {'constants':{k:str(v) for k,v in c.items()},
            'maximum_radius_ratio':str(rho),'example_spacing_ratio':str(q),
            'example_error_bounds':{k:str(v) for k,v in error.items()},
            'joint_family':'rho_m=2^-28/(m+1); q_m=min(1/16,1/(6K1),rho_m^7/(1+C_lattice))',
            'inverse_infinity_norm_upper':'272/15',
            'spacing_constraint_is_sufficient_not_necessary':True,
            'numerical_lattice_extrapolation_used_for_proof':False}
