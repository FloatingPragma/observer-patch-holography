# Original spectral residue and source memory

The centered full occupation source admits a canonical response measure at
every magnetic momentum. Its nodal zero-energy residue is an exact
function of the omitted channel's zero-energy weight. Finite-broadening
operator inequalities retain noncommuting real and imaginary memory. The
actual default source also has positive original spectral weight below
its atomic contact energy on a fixed nodal neighbourhood. These results
supply no nonzero weight in a shrinking ballistic energy window or
physical particle assignment.

## Source measure and Bloch charts

Use the proper-cubic average of the invariant centered full ground in
`MOMENTUM_RESPONSE.md`, with the same action, positive physical excitation
clock and signed retarded generator L in `COMPOSITE_CAR.md`. Write
`alpha=U/2>0`, and retain both original addition and removal weights. Their
combined signed spectral measure has mass one for every normalized original
coefficient packet. Separate smeared charge weights need not equal one half.

The coefficient momentum is that of the commuting squares of the actual
magnetic translations. In the physical-cell chart,

    h(k)=c sum_j cos(k_j) Gamma_j,
    h(k)^2=c^2 S(k) I,  S(k)=sum_j cos(k_j)^2,
    C(k)=kappa h(k),  |kappa|<=3.                       (1)

The original Weyl sin-symbol labels `pi b` require their specified onsite
and Bloch transport before insertion into (1). In a periodic magnetic
Bloch convention, finite-range kernels have trigonometric-polynomial
matrices. The physical-cell cos convention has unitary gluing at the
Brillouin-zone boundaries. Continuity below means continuity in a local
chart, with the same gluing applied to the probe. A boundary node is not
a physical discontinuity.

For every integer n the kernel of `V0*L^n V0` has range at most n.
Its entries are expectations of the anticommutator of an original local
annihilator with the n-fold source derivative at the other site, with
the retarded sign convention. The derivative has radius n, and disjoint
odd operators anticommute. Thus each moment matrix `H_n(k)` is continuous
in its Bloch chart. This argument holds before a ground expectation is
taken. The executable checks it in the complete 64-dimensional source
algebra, including a nonzero range-two coefficient and the actual Bloch
gluing matrices.

The translation spectral decomposition supplies a positive matrix measure
almost everywhere in k. The all-packet even-moment bounds of `UV_TAIL.md`
are uniform in momentum. Along momenta from that full-measure set tending
to any k, tightness gives a weak subsequential measure. Uniform higher
moments imply uniform integrability of each fixed power, so its moments
are precisely `H_n(k)`. The exponential-moment determinacy in `UV_TAIL.md`
makes this limit unique. Apply this construction to all coefficient
vectors and polarize; positivity and total matrix mass I are retained.
This proves a canonical matrix measure at every k, weakly continuous in
each local chart and covariant under chart gluing.

For fixed nonreal z, the original resolvent `G00(z,k)` is consequently
continuous. The exact identity `V1=(LV0-V0h)/alpha` expresses the other
three blocks through the same measure and polynomial moments. Uniform
higher moments control the weighted integrals. Moreover

    L W=W A+[0,R],  W=[V0,V1],
    A=[[h,alpha I],[alpha I,C]],  ||R||^2<=18c^2.        (2)

A uniform second-moment bound for W is
`M2_W=(3 sqrt(3)c+alpha)^2+18c^2`. Scalar Jensen gives
`-Im G_W(i eta,k)>=eta/(eta^2+M2_W) I`. Hence its inverse and the exact
Schur memory extend continuously at every fixed `eta>0`. This defines a
nodal value without choosing an arbitrary representative of an a.e.
fiber. At such a momentum, realize L as multiplication by energy in the
canonical matrix-measure cyclic space, V0 as constant vectors, and
`V1=(LV0-V0h)/alpha`. The moment identities make these two isometries
orthogonal and define bounded R; the same selfadjoint-complement argument
then defines D intrinsically. The construction supplies no joint
`eta->0, k->k_*` modulus and does not make atom weights continuous.

## Exact nodal residue

On the nodal response fiber `h=C=0`, the actual source decomposition is

    L0=[[0,alpha I,0],[alpha I,0,R*],[0,R,D]],
    D=L_perp,
    Sigma(z)=R*(z-D)^(-1)R.                            (3)

D is selfadjoint on its natural domain by the bounded off-diagonal
subtraction argument in `COMPOSITE_CAR.md`; R is bounded. Away from a node,
L0 denotes only the bounded comparison that sets h,C to zero and keeps
the actual R,D. It is not an alternative selected source.

Let `P0=1_{0}(D)` and `M0=R*P0R`. A zero vector of L0 has second component
zero, third component z in ker D, and first component `-R*z/alpha`.
Projection onto that graph gives the exact original zero-energy residue

    Z0=V0*1_{0}(L0)V0=M0(alpha^2 I+M0)^(-1).            (4)

Hidden kernel vectors orthogonal to R contribute no original weight.
The identity does not replace M0 by the full Gram R*R. In particular,
strictly positive total omitted weight does not establish a zero-energy
residue.

## Broadening bounds with noncommuting memory

For `eta>0` set

    M_eta=R* eta^2(D^2+eta^2)^(-1)R,
    Q_eta=-eta Im[V0*(i eta-L0)^(-1)V0],
    beta=alpha^2+eta^2,  r2=||R||^2.                   (5)

Both Q and M are positive operators. M_eta converges strongly to M0 as
eta decreases to zero. Write `H=beta I+M_eta` and
`J=eta Re Sigma(i eta)`. Block inversion gives
`Q_eta=I-alpha^2 Re(H+iJ)^(-1)`. For the selfadjoint
`T=H^(-1/2)JH^(-1/2)`,

    Re(H+iJ)^(-1)=H^(-1/2)(I+T^2)^(-1)H^(-1/2)<=H^(-1).

Therefore

    Q_eta >= (eta^2 I+M_eta)(beta I+M_eta)^(-1).         (6)

No commutation of J and M is used. There is a complementary upper bound.
For every f, spectral calculus gives

    ||Jf||^2
      <=r2 ||eta D(D^2+eta^2)^(-1)Rf||^2
      <=r2 <f,M_eta f>,

so `J^2<=r2 M_eta`. Since `I-(I+T^2)^(-1)<=T^2`, the correction to
`I-alpha^2 H^(-1)` is at most
`alpha^2 H^(-1)JH^(-1)JH^(-1)`. Use `H>=beta I` and the commutation of H
with M to bound it by `alpha^2 r2 M_eta/beta^3`. For a unit f, write
`m_eta=<f,M_eta f>`, `q_eta=<f,Q_eta f>`. The result is

    eta^2/beta + alpha^2 m_eta/[beta(beta+r2)]
      <=q_eta
      <=eta^2/beta + alpha^2(beta+r2)m_eta/beta^3.       (7)

The lower linear bound follows from (6) and `0<=M_eta<=r2 I`.
The code checks (6)-(7) by exact Hermitian principal minors on a complex
matrix example whose real and imaginary memory do not commute. It also
checks (4) with a matrix-valued nontrivial hidden kernel.

On a reducing momentum region with `||h||<=delta` and `C=kappa h`, the
actual source differs from L0 by `diag(h,C,0)`. Resolvent identity implies

    |q_eta(actual)-q_eta(L0)|
         <=max(||h||,||C||)/eta<=3delta/eta.            (8)

For unit hopping, `k=k_*+a p`, `|p|<=K` gives `delta<=aK`. At
`eta=a tau`, (8) is at most `3K/tau`; it is generally O(1), not a vanishing
regulator error. `score_bounds` retains this entire allowance.

## An actual source lower-energy bound

Use unit hopping `c=1,U=8`, so `alpha=4`; homogeneity restores other c.
On `S<=3/250`, `NODAL_MEMORY.md` proves

    R*R>=18I/125,  ||D R||<=5175,  R*R<=18I.            (9)

These are bounds on the specified full ground, with no dimer ground
substitution. Let mu be the original signed measure of a unit coefficient
packet supported there, and put `X=E^2`. The exact source moments give

    E_mu X=16+<S>=:a,
    E_mu X^2=<(16+S)^2+16[(h+C)^2+R*R]>,
    Var_mu(X)>=16(18/125)=288/125.                     (10)

The extra variance of S is nonnegative. A third moment of X supplies a
quantitative lower-energy conclusion. From (2),

    L^3 V0=
       V0[h^3+2alpha^2 h+alpha^2 C]
      +V1 alpha[h^2+alpha^2 I+C(h+C)+R*R]
      +alpha[R(h+C)+D R].                            (11)

The three terms occupy orthogonal subspaces. On `S<=1`, their norms are
at most 81,188,20780 respectively, using `sqrt(18)<5`, (1) and (9).
Consequently the exact admitted bounds are

    E X<=17,  E X^2<=833,
    E X^3<=81^2+188^2+20780^2=431850305.

For `Y=X-a`,

    E|Y|^3<=E(X+a)^3<=431912440=:T0,
    E(Y_-)=E|Y|/2 >= Var(X)^2/(2T0)=:b0.              (12)

The last inequality is Cauchy--Schwarz applied to
`|Y|^(1/2)|Y|^(3/2)`. Its hypotheses use only the finite source moments.
Exactly,

    b0=5184/843578984375,
    delta0=b0/4=1296/843578984375.

Restrict the packet to the smaller positive-measure region `S<=delta0`.
Then `a<=16+delta0`. Outside the event `X<=16-delta0`, the negative
part `Y_-` is at most `2delta0=b0/2`; everywhere it is at most `a<=17`.
Equation (12) proves the actual original probability bound

    mu(|E|<=sqrt(16-delta0))
          >=b0/34=2592/14340842734375>0.              (13)

Thus the full default source has response below the atomic contact
energy 4 on this nodal region. The construction retains both charge
components and never divides by their success probability. It supplies
no positive lower edge for the event in (13); a zero-energy component
is allowed. Restoring units multiplies its energy bound by c. Under
ballistic `c=v/a_spacing`, this is a diverging physical window, not the
shrinking unit-hopping window required for a finite physical branch.

## Low-frequency boundary and exact implication controls

The same source proof supplies, for every eta>0,

    M_eta >= eta^2 r0^2/(eta^2 r0+M^2) I,
    r0=18/125,  M=5175.                              (14)

It follows by Jensen from (9). This lower bound is O(eta^2), so (7)
does not establish a nonzero low-frequency limit. The all-spectrum
force-preparation probability in `NODAL_MEMORY.md` also supplies no
lower bound on M_eta in a shrinking window.

A finite exact control shows why energy reflection and the displayed
moments do not decide this limit. Real zero-diagonal Jacobi matrices
with successive couplings `(4,1,1)` and `(4,1,1,sqrt(3))` both
anticommute with the alternating diagonal grading. They share alpha=4,
`R*R=1`, `R*D^2R=1`, and original moments through degree six

    (1,0,16,0,272,0,4640).

Their characteristic polynomials are

    lambda^4-18lambda^2+16,
    lambda^5-21lambda^3+67lambda.

The first original measure has no zero atom. The second has `M0=3/4`
and original zero weight `3/67`, as (4) requires. Their degree-eight
moments differ. These are implication controls for the response identities
and chiral symmetry; they are not realizations of the full Hubbard source.
Tensoring the finite controls with internal symmetry multiplicities
preserves the comparison. No independent source modes are appended.

The extensive finite-graph charged separation in `strong_coupling_refinement`
is positive under `U>2c|edges|`; its canonical-isometry and dynamics bound
requires `U>=4c|edges|`. At U/c=8 these permit fewer than four edges and
at most two edges respectively, excluding growing cubic sources.
The complete default cube has its separately certified
finite charged gap; it supplies no volume-uniform gap. The positive-field
spin gap concerns a different selected action and state. None of these
premises discharges a centered full-CAR low-frequency estimate.

A concrete sufficient source estimate is a uniform original-packet floor
`<f_a,M_(a tau0)f_a>>=m_*>0` on `|p|<=K`. Monotonicity in eta and (7)-(8)
give `liminf q_(a tau)(actual)>=m_*/34-3K/tau`. Choosing
`tau>204K/m_*`, `tau>=tau0`, gives a positive floor. For an original
unit-mass measure,

    q_eta<=mu(|E|<=N eta)+1/(1+N^2).                  (15)

With `1/(1+N^2)<m_*/136`, the original physical-window weight is at least
`m_*/136` eventually along that family. No effective onset is asserted.
`window_from_score` keeps the sharper exact original-probability inequality.
The source floor in this paragraph is unproved. A dispersion graph and
relativistic field attachment also require the scaled response
`a G00(a zeta,k_*+a p)` and its retained spectral weights. Equations
(4),(7),(13) do not identify that limit or a Standard Model particle.
