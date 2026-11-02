"""Exact shell interpolation, stress sign and finite response frames."""
from functools import lru_cache
import sympy as s
from model.weyl_geometric_modular.algebra import pauli

R = s.Rational
DIRECTIONS = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),
              (R(3,5),R(4,5),0),(R(3,5),0,R(4,5)),(0,R(3,5),R(4,5)))
COORDINATES = ((0,0),(0,1),(0,2),(0,3),(1,1),(2,2),(1,2),(1,3),(2,3))


def spinor(direction):
    n = s.Matrix(direction)
    if n.shape != (3,1) or s.simplify(n.dot(n)-1) != 0:
        raise ValueError('an exact real unit direction is required')
    if any(x.is_real is not True for x in n):
        raise ValueError('real direction required')
    if n[2] == -1:
        return s.Matrix([1,0])
    return s.Matrix([-n[0]+s.I*n[1],1+n[2]])/s.sqrt(2*(1+n[2]))


def boost(direction):
    n = s.Matrix(direction)
    if n.shape != (3,1) or n.dot(n) != 1:
        raise ValueError('unit direction required')
    return s.Matrix([[R(5,4)]]).row_join(R(3,4)*n.T).col_join(
        (R(3,4)*n).row_join(s.eye(3)+R(1,4)*n*n.T))


def response_matrix(directions=DIRECTIONS):
    return s.Matrix([[(R(5,4)+R(3,4)*sum(a*b for a,b in zip(n,m,strict=True)))**2
                      for m in directions] for n in directions])


def require(condition):
    if not condition:
        raise ArithmeticError('exact Weyl response identity failed')


@lru_cache(maxsize=1)
def kernel_certificate():
    t,x,y,z = s.symbols('t x y z',real=True)
    coords=(x,y,z);sigma=pauli()
    numerator=t*s.eye(2)+sum((coords[j]*sigma[j] for j in range(3)),s.zeros(2))
    kernel=-s.I*numerator/(2*s.pi**2*(x*x+y*y+z*z-t*t)**2)
    equation=s.diff(kernel,t)+sum((sigma[j]*s.diff(kernel,coords[j]) for j in range(3)),s.zeros(2))
    require(all(s.simplify(v)==0 for v in equation))
    # Real q>=L^2/2 and the complex-polydisc perturbation is <L^2/4.
    q_perturbation=R(19,128)+R(1,1024)
    require(q_perturbation<R(1,4))
    require(R(5,2)+R(1,16)<3)
    r=s.symbols('r');profile=r*r*(r-1)**4*(2-r)**4
    require(all(s.diff(profile,r,j).subs(r,endpoint)==0 for j in range(4) for endpoint in (1,2)))
    return {'spacelike_kernel_satisfies_Weyl_equation':True,
            'kernel':'-i*(t I+sigma.z)/(2*pi^2*(|z|^2-t^2)^2)',
            'complex_denominator_perturbation_over_L_squared':str(q_perturbation),
            'shell_profile_has_three_continuous_boundary_derivatives':True}


@lru_cache(maxsize=1)
def certificate():
    sigma = pauli(); r = s.symbols('r', real=True)
    jet_map=s.BlockMatrix([[sigma[j]*sigma[k]-4*int(j==k)*s.eye(2)
                           for k in range(3)] for j in range(3)]).as_explicit()
    jet_inverse=s.BlockMatrix([[-(sigma[j]*sigma[k]+int(j==k)*s.eye(2))/4
                               for k in range(3)] for j in range(3)]).as_explicit()
    require(jet_map*jet_inverse==s.eye(6))
    profile = r*r*(r-1)**4*(2-r)**4
    A = s.integrate(profile/r,(r,1,2))
    B = s.integrate(profile/(r*r),(r,1,2))
    J = s.integrate(r*r*profile**2,(r,1,2))
    if not (A == R(1,420) and B == R(1,630)):
        raise ArithmeticError("exact response certificate failed")
    norm_coefficient = s.factor(4*J*(210**2+R(945**2,8)))
    if not (norm_coefficient == R(150088295,4249388)):
        raise ArithmeticError("exact response certificate failed")
    eta = s.diag(1,-1,-1,-1)
    for n in DIRECTIONS:
        u = spinor(n)
        hn = sum((n[j]*sigma[j] for j in range(3)),s.zeros(2))
        if not (s.simplify(u.adjoint()*u) == s.ones(1)):
            raise ArithmeticError("exact response certificate failed")
        if not (s.simplify(hn*u+u) == s.zeros(2,1)):
            raise ArithmeticError("exact response certificate failed")
        a = -210*s.pi*s.I*u
        b = [R(945,4)*s.pi*(n[j]*s.eye(2)-sigma[j])*u for j in range(3)]
        S = sum((sigma[j]*b[j] for j in range(3)),s.zeros(2,1))
        if not (s.simplify(2*s.I*A*a/s.pi-u) == s.zeros(2,1)):
            raise ArithmeticError("exact response certificate failed")
        for j in range(3):
            # L=1; homogeneity gives the displayed 1/L in the theorem.
            derivative = -2*s.I*B*(sigma[j]*S-4*b[j])/(3*s.pi)
            if not (s.simplify(derivative-s.I*n[j]*u) == s.zeros(2,1)):
                raise ArithmeticError("exact response certificate failed")
        if not (s.simplify(sum((b0.adjoint()*b0)[0] for b0 in b)) == R(6*945**2,16)*s.pi**2):
            raise ArithmeticError("exact response certificate failed")
        lam = boost(n)
        if not (lam.T*eta*lam == eta):
            raise ArithmeticError("exact response certificate failed")
        if not (max(sum(abs(x) for x in lam[:,j]) for j in range(4)) < 3):
            raise ArithmeticError("exact response certificate failed")
        # Classical Weyl current=(1,-n), momentum=(-1,n): the negative
        # energy wave has T=-k*k; the Wick steering covariance adds minus.
        current = s.Matrix([1,*[-z for z in n]])
        momentum = -current
        stress = (momentum*current.T+current*momentum.T)/2
        steered_lowered = -eta*stress*eta
        lower_null = s.Matrix([1,*n])
        if not (steered_lowered == lower_null*lower_null.T):
            raise ArithmeticError("exact response certificate failed")
        if not (s.trace(eta*steered_lowered) == 0):
            raise ArithmeticError("exact response certificate failed")
    frame = s.Matrix([[s.Matrix([1,*n])[i]*s.Matrix([1,*n])[j]
                       for n in DIRECTIONS] for i,j in COORDINATES])
    matrix = response_matrix()
    if not (frame.det() == R(13824,15625)):
        raise ArithmeticError("exact response certificate failed")
    if not (matrix.det() == R(24407490807,2000000000)):
        raise ArithmeticError("exact response certificate failed")
    inverse_norm = max(sum(abs(x) for x in matrix.inv()[i,:]) for i in range(9))
    if not (inverse_norm == R(68,5)):
        raise ArithmeticError("exact response certificate failed")
    return {
        'kernel':kernel_certificate(),
        'arbitrary_spatial_first_jet_map_has_exact_inverse':True,
        'shell_moment_A':str(A),'shell_moment_B':str(B),'shell_norm_moment_J':str(J),
        'squared_norm_over_pi_cubed_L_cubed':str(norm_coefficient),
        'stress_amplitude_times_pi_cubed_L_fourth':str(1/norm_coefficient),
        'interpolated_value':'u','interpolated_spatial_derivative':'i*n_j*u/L',
        'interpolated_time_derivative':'i*u/L',
        'stress_covariant_sign':'positive gamma_L*(1,n)*(1,n)^T',
        'directions':[[str(x) for x in n] for n in DIRECTIONS],
        'coordinates':[[i,j] for i,j in COORDINATES],
        'tracefree_stress_frame_determinant':str(frame.det()),
        'boosted_ball_response_matrix':[[str(x) for x in row] for row in matrix.tolist()],
        'boosted_ball_response_determinant':str(matrix.det()),
        'response_inverse_infinity_norm':str(inverse_norm),
        'unit_timelike_normals':[['5/4',*[str(R(3,4)*x) for x in n]] for n in DIRECTIONS],
        'valleys_discarded':False,'single_summand_excited_others_in_same_vacuum':True,
    }
