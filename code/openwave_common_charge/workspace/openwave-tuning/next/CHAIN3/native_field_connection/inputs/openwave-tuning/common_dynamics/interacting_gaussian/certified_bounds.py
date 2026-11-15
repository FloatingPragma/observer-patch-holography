"""Rational bounds for every branch, with no numerical spectral oracle.

The frozen integer spectral certificate supplies the positive gap. Route
norms are recomputed from exact integer currents. These conservative bounds
hold along exact free-field intervals, not arbitrary finite-N split states.
"""
from fractions import Fraction as F
from itertools import combinations,product
from functools import lru_cache
import json
import math
import branches as b


def sqrt_upper(x,bits=40):
    x=F(x)
    if x<0:raise ValueError('nonnegative radicand required')
    den=1<<bits;n=math.isqrt(x.numerator*den*den//x.denominator)
    return F(n if F(n*n,den*den)==x else n+1,den)


def sqrt_lower(x,bits=40):
    x=F(x)
    if x<0:raise ValueError('nonnegative radicand required')
    den=1<<bits
    return F(math.isqrt(x.numerator*den*den//x.denominator),den)


@lru_cache(None)
def route_bounds(particles):
    if particles not in (1,2):raise ValueError('one/two disjoint root charges required')
    currents={}
    for cell in range(8):
        currents[cell]=[]
        for signs in product((-1,1),repeat=3):
            current=b.np.zeros(288,dtype=b.np.int64);here=cell
            for axis,sign in enumerate(signs):
                here,path,J=b.matter.route(here,axis,sign);current+=J
            currents[cell].append(current)
    cells=list(combinations(range(8),particles))
    j2=0;particular2=0
    P=b.matter.field.current_walk.geometry()['P']
    for occupied in cells:
        particular=sum((P[:,12*cell] for cell in occupied),start=b.np.zeros(288,dtype=b.np.int64))
        particular2=max(particular2,int(particular@particular))
        for choice in product(*(currents[cell] for cell in occupied)):
            J=sum(choice,start=b.np.zeros(288,dtype=b.np.int64))
            j2=max(j2,int(J@J))
    return j2,particular2


def contract(macros,particles=1,initial_cells=1):
    if type(macros)is not int or macros<0 or type(particles)is not int or particles not in (1,2):
        raise ValueError('nonnegative macro count and one/two particles required')
    if type(initial_cells)is not int or initial_cells not in (1,8) or particles==2 and initial_cells!=1:
        raise ValueError('one localized tuple or one-particle8-cell plane required')
    spec=json.loads((b.HERE.parent/'quadratic_field/FOCK.json').read_text())['spectral_certificate']
    lower=F(spec['minimum_eigenvalue'])-F(spec['actual_curvature_operator_error_upper'])
    assert lower>0
    wmin=sqrt_lower(F(3,2)*lower);wmax=sqrt_upper(F(3,2)*F(1504,17))
    J2,P2=route_bounds(particles)
    occupation=F(macros*macros*J2,4)/wmin
    branches=initial_cells*8**(particles*macros)
    weight=sqrt_upper(branches)
    polynomial=256*occupation**4+1792*occupation**3+3360*occupation**2+1680*occupation+105
    # Each face's variance squared <=188/51 by the actual K majorant.
    each_face=weight*F(188,51)*sqrt_upper(polynomial)
    vE=wmax
    meanE=sqrt_upper(P2)+2*sqrt_upper(vE*occupation)+macros*sqrt_upper(J2)
    vQco=F(1,4)/wmin
    vQ=max(vQco,F(1,2)+F((macros+1)**2,8))
    meanQ=2*sqrt_upper(vQco*occupation)+F(macros*(macros+1),2)*sqrt_upper(J2)
    return {'macros':macros,'particles':particles,'initial_cell_terms':initial_cells,
            'all_branches':branches,'branch_l1_upper_exact':str(weight),
            'actual_macro_current_norm_squared_max':J2,'integer_particular_norm_squared_max':P2,
            'coexact_frequency_lower_exact':str(wmin),'coexact_frequency_upper_exact':str(wmax),
            'coexact_occupation_upper_exact':str(occupation),
            'each_face_fourth_power_norm_upper_exact':str(each_face),
            'all208_face_fourth_power_norm_sum_upper_exact':str(208*each_face),
            'electric_variance_upper_exact':str(vE),'electric_mean_absolute_upper_exact':str(meanE),
            'theta_variance_upper_exact':str(vQ),'theta_mean_absolute_upper_exact':str(meanQ),
            'domain':'Every branch of the specified routed source through n macros and the following free interval of length1. Pure reference, real field, all190coexact and3harmonic modes; arbitrary spectator references by the same triangle/operator bounds.',
            'not_transferred':'Finite-N electric/magnetic splitting, compact periodization, native finite-register words and interacting spatial refinement are separate maps.'}


def photon_tail_squared(contract,cap):
    if type(cap)is not int or cap<0:raise ValueError('nonnegative integer cap required')
    lam=F(contract['coexact_occupation_upper_exact']);B=contract['all_branches']
    if lam==0:return F(0)
    if cap+2<=lam:return F(1)
    # exp(-lambda)<=1; the remaining positive tail is bounded geometrically.
    tail=lam**(cap+1)/math.factorial(cap+1)*F(cap+2)/(F(cap+2)-lam)
    return min(F(1),B*tail)


def cosine_error_coefficient(contract,substeps):
    if type(substeps)is not int or substeps<1:raise ValueError('positive splitting count required')
    return F(contract['all208_face_fourth_power_norm_sum_upper_exact'])/(8*substeps)
