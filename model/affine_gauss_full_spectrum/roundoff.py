"""Residual enclosure on the complete one-charge carrier.

The candidate eigenvector is untrusted. Outward source coefficient intervals,
the exact finite Peierls polynomial and explicit binary64 evaluation budgets
bound its actual-H residual; the independent full-complement E2 bound isolates
the charged ground. This is not a physical rest-mass certificate.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math
import numpy as np
import mpmath as mp
from .operator import HERE,ROOT,ActualOperator,algebraic_frame
from model.affine_gauss_spectrum import model as coarse

U=F(1,2**53)


def gamma(n):
    return n*U/(1-n*U)


def q(x):return F.from_float(float(x))


def up(x):
    return np.nextafter(float(x),np.inf)


def matrix_absolute_norm(M):
    # Row sums of exact binary64 entries, rounded outward once. The matrices
    # here have symmetric absolute values, so this bounds the induced2-norm.
    maximum=F(0)
    for i in range(M.shape[0]):
        row=M.data[M.indptr[i]:M.indptr[i+1]]
        maximum=max(maximum,sum((abs(q(x)) for x in row),F(0)))
    return maximum


def product_entry_error(value,error,factors):
    iv=mp.iv;maximum=F(0)
    for a,ea in zip(np.asarray(value).ravel(),np.asarray(error).ravel()):
        exact=iv.mpf([str(q(a)-q(ea)),str(q(a)+q(ea))])
        for n in factors:
            stored=float(a)*math.sqrt(n)
            target=exact*iv.sqrt(n)
            maximum=max(maximum,abs(q(stored)-coarse.endpoint(target)),
                        abs(q(stored)-coarse.endpoint(target,True)))
    return maximum


def evaluation_budget(frame,operator):
    """Uniform norm error per unit input for this actual floating implementation."""
    mp.iv.dps=80
    field_entry=product_entry_error(frame['field'],frame['field_error'],(1,2,3))
    # Each creation row has at most3 entries, each column at most30;
    # combined creation−annihilation has at most33 per row and column.
    D_error=33*field_entry
    omega_entry=product_entry_error(frame['omega'],frame['omega_error'],(1,2,3,4,9))
    # At most3 occupied annihilation modes times30 creation targets. The
    # assembled diagonal sums at most3 terms, adding its CSR summation error.
    omega_sum=3*max(abs(q(x)) for x in frame['omega'].ravel())
    photon_matrix_error=90*omega_entry+gamma(6)*omega_sum
    photon_abs=matrix_absolute_norm(operator.photons)
    photon_degree=max(np.diff(operator.photons.indptr))
    photon_eval_error=gamma(2*int(photon_degree))*photon_abs
    # The four-level collective q norm is sigma sqrt(3+sqrt6)<sqrt(1/2)*2.335<5/3.
    exact_D_norm=F(5,3);M=exact_D_norm+D_error
    assert M<2
    D_abs=max(matrix_absolute_norm(D) for D in operator.Ds)
    D_degree=max(int(max(np.diff(D.indptr))) for D in operator.Ds)
    D_round=gamma(2*D_degree)*D_abs
    polynomial_errors=[];polynomial_norms=[]
    for c,ce in zip(frame['peierls'],frame['peierls_error']):
        c=list(map(q,c));ce=list(map(q,ce))
        assert max(map(abs,c))<2
        coefficient_error=sum((ce[k]*exact_D_norm**k for k in range(9)),F(0))
        argument_error=D_error*sum((abs(c[k])*k*M**(k-1) for k in range(1,9)),F(0))
        # Horner source: multiply highest coefficient, then eight sparse
        # products, coefficient multiplications and additions. Signs exact.
        beta=(1+U)*abs(c[8]);error=U*abs(c[8])
        for k in range(7,-1,-1):
            multiplication=D_round*beta
            addition=U*((M+D_round)*beta+(1+U)*abs(c[k]))
            error=M*error+multiplication+U*abs(c[k])+addition
            beta=(1+U)*((M+D_round)*beta+(1+U)*abs(c[k]))
        total=coefficient_error+argument_error+error
        polynomial_errors.append(total);polynomial_norms.append(1+total)
    W_error=max(polynomial_errors);W_eval_norm=max(polynomial_norms)
    # Two initial summands: Hph x and the static diagonal. Every vertex
    # then accumulates at most12 hopping terms, each divided by24.
    static_error=max(map(q,frame['static_energy_error']))
    static_norm=max(abs(q(x)) for x in frame['static_energy'])
    initial_error=photon_matrix_error+photon_eval_error+static_error+U*static_norm
    initial_error+=U*(photon_abs+photon_eval_error+(1+U)*static_norm)
    hopping_error=F(1,2)*(W_error+U*W_eval_norm)
    final_sum_error=gamma(12)*(photon_abs+static_norm+initial_error+F(1,2)*(1+U)*W_eval_norm)
    # On the checked input-norm range[1/2,2], fewer than2^30 scalar
    # operations and less than2^10 subsequent Horner amplification admit
    # this loose per-unit-input underflow budget. It is not a homogeneous
    # bound on arbitrarily tiny input vectors.
    underflow=F(2**46,2**1022)
    total=initial_error+hopping_error+final_sum_error+underflow
    return {'H_apply_norm_error_per_unit_input':str(total),
        'exact_D_norm_upper':str(exact_D_norm),'rounded_D_error':str(D_error),
        'D_absolute_row_norm':str(D_abs),'D_maximum_row_entries':D_degree,
        'photon_absolute_row_norm':str(photon_abs),'photon_maximum_row_entries':int(photon_degree),
        'photon_matrix_error':str(photon_matrix_error),
        'Peierls_error_per_unit_input':list(map(str,polynomial_errors)),
        'arithmetic':'IEEE754 binary64 round-to-nearest; finite normal coefficients; explicit gamma_n and underflow bounds'}


def squared_norm_exact(v):
    # Exact rational sum of squares of the binary trial and residual entries.
    return sum((q(x)**2 for x in np.asarray(v).ravel()),F(0))


def sqrt_bounds(value):
    iv=mp.iv.mpf(value.numerator)/value.denominator
    return coarse.endpoint(mp.iv.sqrt(iv)),coarse.endpoint(mp.iv.sqrt(iv),True)
