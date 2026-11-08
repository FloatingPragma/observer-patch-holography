# Energy polarization on the actual ordered source

## 1. Source and cutoff

Use the same selected zero-field ground omega and spin source as
`spin_canonical_field`, `spin_canonical_energy` and `spin_current_spectrum`:

    H=J sum_<xy>(Sx.Sy-1/4),  e_x=(J/2)sum_(y~x)(Sx.Sy-1/4),
    L_p=K/(Jp),  Q_p(F)=p² sum_x epsilon_x F(px) Sx³,
    P_p(G)=-p sum_x G(px) Sx²,  X=px,  epsilon_x=(-1)^x.

Each unoriented source bond occurs once in H. The reference is the actual
maximal field-selected ground; a finite selected cube is not substituted
for it. Its common diagonal bond covariance has coefficient
`kappa=-(c_parallel+c_transverse)>0`. All formulas below retain that actual
coefficient. The stronger variational bounds can be used without assigning
its numerical value.

Let chi_R equal one on a growing spatial plateau, have compact support and
lie between zero and one. Define the finite source operator

    B_(i,p,R)=(1/J) sum_x x_i chi_R(px)e_x.

This is the candidate weighted energy in the same length and clock units:
`x_i e_x/J = X_i e_x/(Jp)`. Its coefficient on a bond is the average of
its two site weights. It is not assumed to generate Lorentz boosts.

If a local spin-rotation preparation is supported in the max-norm box
of radius r, only bonds touching this support change under conjugation.
For plateau radius at least r+p, its **energy-polarization increment**
`U* B_(i,p,R) U-B_(i,p,R)` is exactly independent of R. The same is true
of its commutator with any operator on that preparation support. No
large-cutoff norm limit of B itself is needed for these local identities.
The half-bond density, source energy Ward identity and its block transport
are the ones already verified in `native_spin_clock_join/energy.py` and
`QFT/CoarseEnergyWard.lean`.

## 2. Exact weighted work and controlled limit

For a transverse single-axis rotation `U_theta=exp(i sum theta_x Sx^a)`,
the actual source bond covariance gives the exact expectation

    Delta B_i = kappa sum_<xy> ((x_i+y_i)/2)
                            [1-cos(theta_x-theta_y)].              (1)

The ground expectation of the sine-current term vanishes by the inherited
real-state and spin-rotation symmetries. Multiplying each bond identity by
its actual weight proves (1); no product-state approximation is made.
The fixed ground energy scalar cancels from this difference.

For compact C3 tests, theta=p² epsilon F(px) for Q and theta=-pG(px) for P.
Taylor bounds and midpoint Riemann sums yield, respectively,

    Delta B_i(Q_p(F)) -> 6 kappa integral X_i F²,
    Delta B_i(P_p(G)) -> (kappa/2) integral X_i |grad G|².             (2)

These are limits of the actual microscopic expectation differences. They
are not asserted to equal matrix elements of an already reconstructed
continuum boost. In particular `spin_canonical_energy` proves only an
upper bound for continuum energy when first energy can escape the limiting
representation. Weighted source expectation convergence does not remove
that issue.

Here are explicit mesh constants. Suppose every derivative of order j is
bounded in absolute value by M_j, the support lies in [-r,r]^3 and p<=1.
Put V=(2r+2p)^3 and

    H_Q=4M0 M1+2(r+p)(M1²+M0 M2),
    H_P=4M1 M2+2(r+p)(M2²+M1 M3).

The absolute error in (2), divided by kappa, is at most

    Q: (27/4)p² H_Q V +(3/2)(r+p)p² M0 M2 V
                                      +2(r+p)V p^4 M0^4,
    P: (27/16)p² H_P V +(r+p)V p² M1 M3/8
                                      +(r+p)V p^4 M1^4/8.         (3)

For Q, the endpoint average differs from the midpoint value by at most
p²M2/8; the squared endpoint sum differs from four midpoint squares by
at most p²M0M2. For P, the edge difference quotient differs from the
midpoint derivative by at most p²M3/24. Its squared error is at most
p²M1M3/12. For each midpoint cell, Taylor's formula bounds the quadrature
error by `(9/8)p²` times a uniform Hessian bound times cell volume.
The Hessian bounds are H_Q and H_P. All affected cells fit in the stated
volume V. Finally `|1-cos u-u²/2|<=u^4/24`, with relative angles at most
2p²M0 and p²M1. These observations give every term of (3).

`samples.py` additionally computes the quadratic lattice sums exactly for
`F(X)=product_j (1-(X_j-a_j)²)_+^4`, with a lattice-aligned center. It
retains an exact fourth-angle remainder. Its values are per kappa, not
numerically selected energies.

The actual symmetric Weyl preparation W_p(F,G) also has the sum of the
two limits in (2). Reuse the source SO(3) linearization in
`spin_canonical_energy`: its covariance-weighted mixed quadratic term
is exactly zero, and its unweighted nonlinear energy error is at most
`3J t sum_edges[(Delta a)²+(Delta b)²]`, where a=-pG, b=p²epsilon F
and t<=p M0_G+p² M0_F<=1. On the preparation support the weighted bond
factor is at most (r+p)/p. Thus the absolute polarization remainder is
at most

    9 t (r+p) V (M1_G²+4 M0_F²),                            (3a)

which tends to zero with p. Add the two mesh errors (3), multiplied by
the actual kappa. This proves the combined source functional
`D_i(F,G)=6kappa integral X_iF²+(kappa/2)integral X_i|grad G|²`.
It does not infer a quadratic continuum Hamiltonian or Gaussian dynamics.

There is an exact source translation identity before any limit. Translate
the preparation by a=n p, with n an allowed even source translation, and
place the cutoff plateau around both supports. Changing the bond indices
in (1) gives

    Delta B_i(U translated by a)
        =Delta B_i(U)+a_i Delta H(U)/(Jp).                         (4)

The same law holds at the operator-increment level under source translation.
The limit is the corresponding affine law for the two integrals in (2).
This relates the energy first moment to the **actual spatial translations**
and the common source energy; it does not identify the finite energy-flux
operator with the translation generator.

## 3. The exact remaining scalar-field action

Write j_xy^a=J(Sx cross Sy)^a on oriented positive-axis bonds. The existing
spin-current identity gives, for arbitrary real g and f,

    [B(g),A(f)] = [H/J,A(gf)]
             + (i/2J)sum_<xy>(g_y-g_x)(f_x+f_y)j_xy^a,           (5)

where B(g)=sum_x g_x e_x/J and A(f)=sum_x f_x Sx^a. The elementary edge
identity is checked as exact source matrices for all three spin axes.
For g=x_i on the preparation support, (5) becomes

    [B_i,Q_p(F)] = [L_p,Q_p(X_i F)] + i R_(i,p)(F),
    R_(i,p)(F)=(p²/2J)sum_x epsilon_x
                     [F(px)-F(p(x+e_i))] j_(x,i)^3.              (6)

Thus the candidate scalar-field boost action has a specific extra source
observable. Its vacuum norm is precisely the quadratic four-spin response

    ||R_(i,p)(F)Omega||²
      =p^4/(4J²) sum_(x,y) epsilon_x epsilon_y Delta_iF_x Delta_iF_y
                                      omega(j_(x,i)^3 j_(y,i)^3). (7)

The current mean is zero, so (7) is also a connected-current expression.
The elementary current norm1/2 gives only
`||R|| <= M1 (2r+2p)^3/4`. This estimate does not vanish with p.
Replacing one spin factor in the current by the ordered mean would change
(7) and is not part of this construction. The inherited longitudinal
current sum rule does not fix this individual directional covariance.

On any actual propagated source domain vector, (6) is exact. Its residual
norm there is obtained by replacing omega in (7) with that vector's
matrix functional. Consequently a vanishing bound on these stated
current responses, together with the appropriate generator-domain and
state-covariance limits, transfers the scalar-field action. Positive
spectral filtering may be retained explicitly when second energy moments
are unavailable. No unbounded product domain is inferred from a first
energy moment alone.

This reduces a physical boost premise to actual current matrix elements
and an actual vacuum/history covariance condition. The finite cube defect
in `spin_boost_stress` neither supplies nor rules out their scaling limits.
The finite-rotation/one-operational-boost cone theorem can be used after
these operational conditions, rather than from the weighted symbol alone.

## 4. Schwartz preparations and complete charge records

A compact Fourier annulus gives a Schwartz spatial test, not compact
spatial support. Let zeta_L be one for |X|<=L, zero for |X|>=2L, between
zero and one, and |grad zeta_L|<=b/L. Use F_L=zeta_L F and G_L=zeta_L G.
At each fixed L the preceding cutoff construction applies. If |F| and,
for G, |grad G| are bounded by C(1+|X|)^(-N), N>2, then the absolute
differences of the respective continuum polarizations are bounded by

    Q: 24 pi kappa C² L^(4-2N)/(2N-4),
    P: 4 pi kappa C²(1+b²/L²) L^(4-2N)/(2N-4),                 (8)

for L>=1. This follows by weighting the radial tail by |X|. For P use
`||grad G|²-|grad(zeta G)|²| <=2|grad G|²+2|G grad zeta|²`
on the tail. The exact finite-source comparison, without replacing the
angles by their limits, is

    |Delta B(theta)-Delta B(theta_L)|
      <=(kappa/2)sum_<xy>|(x_i+y_i)/2|
           |Delta theta-Delta theta_L|(|Delta theta|+|Delta theta_L|). (9)

It uses the scalar cosine identity. Schwartz decay makes the weighted
sum absolutely convergent and gives dominated lattice Riemann sums.
Removing the source weight cutoff is controlled by the same absolute
tail sum, since its omitted coefficient is bounded by the uncut weight.
Taking p to zero at fixed L and then L to infinity is therefore justified;
no plateau is claimed to cover an untruncated annular test.

The **complete charged preparation** also has a controlled truncation,
although its coherent energy matrix elements are a separate question.
Use the actual commuting configuration fields x=s Q3(F), y=s Q2(F), and
the bounded charged creators of `spin_charge_dispersion`:

    B_+=(1/2pi) integral exp(-i theta) f(x cos theta+y sin theta)dtheta,
    f(t)=t/(1+t²),  B_-=B_+*,  B_0=sqrt(I-B_+*B_+-B_-*B_-).

The same configuration fields with F_L commute with these fields. Since
|f'|<=1, the angular Cauchy-Schwarz inequality and source rotational
symmetry give

    ||(B_+(F)-B_+(F_L))Omega||² <= s² V(F-F_L),
    V(h)=3 integral |hhat(k)|²/|k| d³k/(2pi)^3.                 (10)

This is the actual source configuration variance bound, not a fitted
Gaussian vacuum. Both creators have norm at most1/2, so B_0>=I/sqrt(2).
Pointwise in the common configuration spectral calculus, the complement
difference is at most sqrt(2) times the creator difference. Hence the
**full three-record isometry** changes on Omega by at most
`2|s|sqrt(V(F-F_L))`. The zero record is included.

For an explicit spatial bound, split Fourier momentum at1 and use

    V(h)<=3[||h||_1²/(4pi²)+||h||_2²].                       (11)

If |F|<=C(1+|X|)^(-N), N>3, then
`||h||_1<=4pi C L^(3-N)/(N-3)` and
`||h||_2²<=4pi C² L^(3-2N)/(2N-3)`.
The producer keeps the rational pi bounds outward. This proves complete
vacuum-preparation convergence, not operator-norm convergence for arbitrary
inputs. Later complete channels contract its error.

Equations(8)-(9) concern individual Weyl kick energy increments.
Equations(10)-(11) concern coherent charged preparation and all its records.
They do not replace a coherent superposition's energy-polarization matrix
elements by a classical mixture of those kicks.

## 5. Finite observation and remaining physical attachment

On a finite source, B is a declared weighted sum `b0 I+sum_P b_P P` of the
actual exchange Pauli terms. If there are no nonconstant terms, every
record can receive the constant b0 and the increment is exactly zero.
Otherwise a uniform choice among M terms, followed by
its complete signed Pauli readout, has score `b0+M b_P sign`. All failure
outcomes may be assigned zero within the same bounded score and retained.
A complete-output TV error epsilon changes its mean by at most
`2 epsilon max(|b0+M b_P sign|,0)`. Comparing two preparations adds their
mean-error bounds; the scalar reference cancels. This is the same complete
source readout principle as `source_energy_effort`, with actual spatial
weights. The construction needs finite preparation/control precision,
which remains explicit; no growing native experiment is claimed executed.

The source supplies weighted work, exact spatial affine covariance,
cutoff and mesh control, and the explicit current residual. Physical boost
covariance additionally needs convergence of the resulting operation on
the actual field/history domain, the retained source state, and the same
clock. Lorentz stress, a particle mass shell and Einstein geometry are not
consequences of the first moment alone.
