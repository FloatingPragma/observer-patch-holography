# An interior marginal modulus at positive selection field

## 1. The actual source and the interpolation

Use the physical spin frame of `native_ordered_ground`: on a finite
bipartite cubic subgraph, each edge has

    h_xy=J(S_x.S_y-I/4),  ||h_xy||=J,
    H_B=sum_edges h_xy-B sum_x epsilon_x S_x^1,  J>0, B>0.

Choose an interior box Lambda and split a larger finite graph into
Lambda and its complement. Let D be the crossing edges, M=|D|, and

    H(s)=H_cut+s V,  V=sum_(e in D) h_e,  0<=s<=1.           (1)

Every field is retained at its original value. Only the crossing bonds
vary. In the Marshall frame followed by global Hadamard, each nonnegative
edge strength J_e has Lee-Yang coefficients

    w_z=w_x=J_e/4,  w_y=-J_e/4,  w_xy=w_yx=0,  h=B/2.

The condition in Definition1 of Rayudu--Takahashi,
[arXiv:2607.10765v1](https://arxiv.org/html/2607.10765v1), is an equality.
The scalar bond shift has no effect on ground projections. Their Theorem11
therefore gives a unique ground projection P_s and uniform gap

    g=B/8                                                     (2)

throughout (1). Zero-strength crossing edges are admitted by that theorem.
At s=0 the ground is the tensor product of the two cut-region grounds.
This is an explicit source substitution into an imported finite theorem;
it is not an independent proof of the Lee-Yang theorem.

## 2. Gaussian approximate spectral flow with its error

Write tau_t^s(V)=exp(itH(s)) V exp(-itH(s)). For alpha>0 set

    W_alpha(t)=-sign(t) erfc(sqrt(alpha)|t|)/2,
    D_alpha(s)=integral_R W_alpha(t) tau_t^s(V) dt.             (3)

The integral converges in norm, D_alpha is Hermitian, and

    integral_R W_alpha(t) exp(i omega t) dt
       =-i [1-exp(-omega^2/(4alpha))]/omega,                  (4)
    integral_R |W_alpha(t)| dt=1/sqrt(pi alpha).

The value at omega=0 in (4) is its continuous value zero. Integration by
parts of the half-line sine integral proves (4); its remaining cosine
integral is the Fourier transform of a Gaussian.

Choose the differentiable normalized ground vector in parallel phase,
<Omega_s,Omega_s'>=0. With Q_s=I-P_s and
R_s=Q_s(H(s)-E_s)^(-1)Q_s, differentiation of the eigenvalue equation gives
Omega_s'=-R_s V Omega_s. Equations(2)-(4) then give

    ||Omega_s'+iD_alpha(s)Omega_s||
       <= ||V|| exp(-g^2/(4alpha))/g.                        (5)

The sign matters: the approximate ground-vector equation is
Omega_s'=-iD_alpha(s)Omega_s. In particular every fixed observable A obeys

    |d_s <A>_s-i<[D_alpha(s),A]>_s|
       <=2||A|| ||V|| exp(-g^2/(4alpha))/g.                  (6)

This estimates the derivative of the actual ground state, not of a Gibbs
or trial reference. No gap above a regional density matrix is used.

## 3. Applying the source locality bound

Let A be any operator supported in X inside Lambda, m=|X|, with graph
distance at least d>=1 from every crossing-edge support. The source
commutator-chain bound from `spin_canonical_field` gives

    ||[tau_t^s(A),h_e]||
       <= ||A|| J m 2^(-d) [exp(48J|t|)-1].                (7)

It applies uniformly in s: the interpolated edge strengths are between
zero and J. Passing to the onsite-field interaction picture preserves
supports and bond norms, so (7) is independent of B and of exterior size.
The trivial bound is 2J||A||. Put

    v=72J,  T=d/(2v)=d/(144J),  alpha=g/(2T).                (8)

For |t|<=T, log(2)>=2/3 makes the right side of (7) at most
||A|| J m exp(-d/3). The integral of |W_alpha| over this interval is
at most T. On the two tails,

    integral_(|t|>T) |W_alpha(t)| dt
      = exp(-alpha T^2)/sqrt(pi alpha)-T erfc(sqrt(alpha)T)
      <= exp(-alpha T^2)/sqrt(pi alpha).                    (9)

Summing the M crossing bonds in (3), and then using (6), yields

    |d_s <A>_s|/||A||
      <= JM [m T exp(-d/3)
            +(2/sqrt(pi alpha)+2/g) exp(-gT/2)].            (10)

Here alpha T^2=g^2/(4alpha)=gT/2. Neither a trace over the exterior nor
a Hilbert-space dimension enters the estimate. The number M is only the
number of bonds crossing the chosen interior boundary.

## 4. The marginal bound and the infinite-volume state at fixed B

Integrate (10) from zero to one. It holds for every Hermitian A on X,
with its operator norm, so trace-norm duality gives directly

    (1/2)||rho_X(full)-rho_X(cut)||_1 <= min(1,F),
    F=JM[(mT/2) exp(-d/3)
          +(sqrt(2T/(pi g))+1/g) exp(-gT/2)].              (11)

No Pauli-tomography dimension factor is required. The cut marginal on X
is exactly that of the smaller box's unique ground. For an interior box
[-R,R]^3, an observation support inside [-r,r]^3 has d>=R-r and at most
M=6(2R+1)^2 crossing bonds. At fixed J,B,r,m, the bound tends to zero
as R tends to infinity, independently of the exterior volume or shape.
The factors are polynomial in R and inverse B, and the second exponent
is exp(-Bd/(2304J)).

Consequently the finite open-graph ground marginals are Cauchy on every
fixed local region. They define a common positive-field infinite-volume
state, independent of the exhausting cubic subgraphs, and (11) bounds
the distance of the smaller-box marginal to this limit. The argument
uses only exhaustions retaining the same local source couplings and
positive field, with no added boundary terms. It does not assert a
bound for arbitrary boundary fields outside this family.

The executable planner uses a rational upper certificate. Set
eta=min(1/3,g/(288J)). Since pi>=3, sqrt(x)<=1+x for x>=0,
and e>=2, equation(11) is bounded by

    JM[mT/2+1+2T/(3g)+1/g] 2^(-floor(eta d)).                (12)

Doubling the source radius until (12) is below the requested rational
accuracy terminates: a fixed-degree polynomial times this decaying
exponential tends to zero. The search evaluates integers and rational
bounds, without constructing a source state or Hilbert-space matrix.
It provides a sufficient radius, not an optimal or practically small one.

For a complete instrument supported on X, including reference systems
initially independent of the source, trace-distance contraction carries
(11) to its entire output record and retained X/reference state. The source
exterior is not included in that output-distance claim. The bound is not an
all-input channel bound and does not cover an arbitrarily correlated
source/exterior preparation. Source gate errors and any clock-truncation
error add separately.

## 5. What remains at zero field

This construction supplies a computable volume modulus for each B>0.
It can be combined with the finite positive-field adiabatic preparation
of the same source. Its constants diverge when B tends to zero. The
convergence rate from the positive-field infinite-volume state to the
maximal-order zero-field state is a different premise and is not supplied
by (11). No zero-field gap, sharp dispersion, particle mass, Lorentz boost
or geometric stationarity follows from this positive-field modulus.
