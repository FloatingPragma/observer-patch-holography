# A resolved native target and finite coefficient certificates

This is a source attachment design, not a new execution. It avoids using
a one-cell target whose capacity shrinks with the light-field grid spacing.
The scalar clock, slow-body and Maxwell readouts use the same target and
native hit probability in every member.

## Resolved target and scale order

Choose a fixed physical target region of radius R0 and a light impact
b much larger than R0. On a cubic native chart of spacing a, the target
contains every selected cell in that region; escape is a disclosed large
outer boundary. The target is therefore resolved by many cells as a->0.
Native current records reconstruct the probe location. Conditional on a
move its six directions are uniform, so its target-before-escape probability
is the discrete harmonic Dirichlet solution, independently of native holds.
The controller still uses completed hit letters, never a numerical potential
table.

The allowed hierarchy is

    a << wavelength << R0 << b << escape radius,

with packet widths and flight distance inserted into the explicit error
bounds of PROOF.md. A single fixed lattice cell instead has physical
capacity proportional to a; using its shrinking gravitational signal
without a separate relative-error argument could invalidate a proposed
light comparison. A resolved target avoids that particular obstruction.

For a spherical continuum target and spherical escape boundary R_out,
the limiting harmonic function is

 h(r)=(1/r-1/R_out)/(1/R0-1/R_out).

Its coefficient C_out=(1/R0-1/R_out)^(-1) and offset-C_out/R_out are fixed
by the target/escape apparatus. The offset, its clock normalization and
finite-boundary effects are retained until R_out tends to infinity. A
cubic or otherwise nonspherical target has its own independently measured
capacity and multipoles; it is not assigned the spherical formula.

The standard finite-range random-walk tools used for this mathematical
refinement are the Green asymptotic, first-hit representation and interior
first/second difference bounds in Lawler–Limic, Theorem4.3.1 and
Theorem6.3.8, [Random Walk: A Modern Introduction](https://www.math.uchicago.edu/~lawler/srwbook.pdf).
These results concern the derived cubic probe law. They supply no Maxwell
field or physical gravitational normalization. In a finite certificate,
source coefficient and difference errors are computed explicitly as below.

## Exact finite source certificates, without trusting a floating solve

For a box with integer coordinates max|z_i|<=B, the native one-step lazy
spatial walk has interior transition

    P(x,x)=1/2,       P(x,x+e)=1/12 for six neighbors.

This is the branch law of the pinned endpoint and six-address observer.
Let Q be its transient restriction after target and escape absorption,
and let b_hit be its one-step target-hit vector. The exact source law is

    (I-Q)h=b_hit.                                         (C1)

A rational approximate vector h_tilde, including exact boundary values,
can be certified by its exact rational residual r=(I-Q)h_tilde-b_hit.
The stopped squared-position martingale has drift1/2, so the expected
absorption time is at most6B^2. Therefore

    ||h_tilde-h||_infinity <=6B^2 ||r||_infinity.           (C2)

This follows from (I-Q)^(-1)r and its row sums. A floating linear solve
may propose h_tilde, but its binary floats are converted to exact rationals
before computing r. The residual, rather than an unexamined solver flag,
provides the certificate. No native experiment is replaced by that solve:
(C1) is the exact population of the authenticated branch kernel and the
actual controller continues to use native hit events.

For any finite edge/face stencil S, apply its linear coefficients to the
interval(C2). A source derivative certificate then has error at most
||S||_1 times that interval radius, with its actual powers of1/a. For
example, a first difference costs at most2 delta_h/a. It is not enough
to show only small potential error when a derivative enters a Maxwell
residual. An execution can tighten the rational solve to compensate and
compare the actual discrete differences against a smooth target profile.
Those finite checks give the required coefficient/residual bounds directly.

The cap tail is also constructive. Conditional on any interior history,
expected remaining absorption time is at most6B^2. Markov's inequality and
the strong Markov property give

    P(tau>12 k B^2)<=2^-k.                                (C3)

A fixed T=12B^2 ceil(log2(1/delta)) therefore bounds timeout probability
by delta. Native relocation between any two box positions has at most6B
Manhattan steps; R retries per step gives failure probability<=6B 2^-R.
Every failed move, padding slot and timeout remains in the archive. These
bounds determine real resource costs before a requested field error is
claimed.

## Same-source Maxwell coefficients and the reciprocal inner walk

The ideal electric weights M_xy=u_xu_y, u=1+h, are compiled through the
reciprocal weighted hit walk and its q flags. The latter has symmetric
edge transition probabilities at least1/48 on the connected box (before
absorption). Its finite Dirichlet law has the same exact residual method.
For a crude explicit mean-time bound, extend the comparison kernel to
the entire closed box, reflecting missing outer edges as holds, and use
the electrical commute-time identity for its symmetric stochastic
conductances. Its total vertex measure is N=(2B+1)^3. The target contains
the origin, and a coordinate path from any site to that fixed target
vertex has at most3B edges and resistance at most48 per edge. The hitting time is
bounded by the corresponding commute time, hence by144 N B. Treat loops
in the vertex measure as part of the stochastic conductance network.
Absorption at the whole target/escape set can only occur sooner.

This is a loose bound for the ideal symmetric kernel. Actual persistent
inner reset errors are coupled to that ideal kernel until their first
difference and charged by the same expected-time bound. One cannot apply
the symmetric commute identity directly to an arbitrary history-dependent
failed-reset process. Complete-word errors are tightened further before
rare compatibility conditioning, as required by the conditional-flux theorem.

All source probability intervals are propagated through the finite
products defining u, reciprocal q, electric M and magnetic beta. Since
1<=u<=2, these maps have explicit bounded derivatives on the entire admitted
range. The spatial Maxwell residual is evaluated with those intervals,
not with point estimates silently treated as exact coefficients.

## What still must be frozen for a source light run

A finite target and escape chart; its complete or locally authenticated
native branch dictionary; the rational residual and derivative intervals;
the actual nested probe and reset limits; electric/magnetic sampling points;
a source-derived packet preparation; the weighted discrete propagator and
its OS correction; source clock bins and detector windows; the WKB residual
and source-to-continuum error budget; and a changed controller/source control.
All of these precede source execution and belong in one preregistration.
The relative delay/deflection test can use the same certified capacity
without calling it G times target inertial mass.

## Explicit spherical boundary/continuum bound

For this subsection let the native grid be a Z^3 with a<=R0/2. Absorb
on first entry into |x|<=R0 or first exit to |x|>=Rout, with Rout>R0.
Let C=(1/R0-1/Rout)^(-1) and H(x)=C(1/|x|-1/Rout) in the continuum annulus.
A first inner landing has radius between R0-a and R0; a first outer
landing has radius between Rout and Rout+a. Thus the boundary mismatch
between H and the native hit values1,0 is at most2 C a/R0^2.

For one actual native proposal, P holds with probability1/2 and goes to
each nearest neighbor with probability1/12. Therefore

 (P-I)H=(1/12) sum_i[H(x+ae_i)+H(x-ae_i)-2H(x)].

The second-order contribution is a^2 Delta H/12=0. Each paired fourth
remainder is at most a^4 sup|partial_i^4 H|/12. The exact derivative
identity partial_i^4(1/r)=24 r^-5 P4(x_i/r) and |P4|<=1 give

 |(P-I)H| <= (a^4/144)*3*24 C/(R0-a)^5
           = C a^4/[2(R0-a)^5].                          (C4)

This retains the lazy-native normalization. The stopped process
|X_(t wedge tau)|^2-(a^2/2)(t wedge tau) is a martingale before absorption.
Its terminal radius is at most Rout+a, so monotone stopping gives

 E_x tau<=2[(Rout+a)^2-|x|^2]/a^2
         <=2(Rout+a)^2/a^2.                              (C5)

Finiteness follows already from a finite path-to-absorption bound; the
martingale then supplies the sharper polynomial estimate. Applying the
discrete stopped Dynkin identity to H, its boundary mismatch and(C4)-(C5)
gives the fully explicit source bound

 ||h_a-H||_infinity
   <=2 C a/R0^2 + C a^2(Rout+a)^2/(R0-a)^5 =:delta_h.     (C6)

This is a bound on the exterior grid and its first-landing boundary layer;
it does not extend the singular1/r formula through the target interior.
On a fixed exterior compact set, subsequently sending Rout to infinity
adds the exact finite-escape difference of the displayed radial functions.
For example, with r>=R0,

 |H(r)-R0/r|=R0(r-R0)/[r(Rout-R0)]<=R0/(Rout-R0).         (C7)

Choose a with a Rout->0 and then Rout->infinity (or an explicit joint
diagonal) to make(C6)-(C7) vanish. Light/clock error budgets use their
actual numerical values, not just this existence statement. The same C
is determined by the target/escape boundary apparatus before light data.
The bound does not infer source inertial mass from target size.

## Uniform coefficient error is enough in the first-order Maxwell form

A derivative of a noisy h_a estimate is unnecessary for the main Maxwell
wave comparison. Use a cubic Yee cochain pair, with electric mass
M_E,xy=u_a(x)u_a(y), and magnetic mass M_H,f=u_f^2. Define u_f as the
average of u_a at the four corners of face f. This average is itself
native-compilable: choose one of the four vertices by two fair bits, run
its plus-hit instrument, and repeat independently for the squared gate.
It uses the same source h_a and requires no off-grid potential oracle.
Its additional vertex-choice/reset costs must be included in a new task.
The resulting magnetic precision is exactly beta_f=M_H,f^-1.

For smooth continuum u in[1,2], let M1 bound its first derivatives and M2
its Hessian operator norm on the packet neighborhood. Comparing electric
edge midpoint and magnetic face-center samples gives

 |M_E-u(midpoint)^2| <=4 delta_h+a^2(2M2+M1^2)/4,
 |M_H-u(facecenter)^2| <=4 delta_h+a^2 M2.                 (C8)

The product bound follows by differentiating u(x-se)u(x+se), whose first
derivative vanishes at s=0. The four-corner average cancels the linear
Taylor term and has error<=a^2 M2/4 before squaring. Native instrument,
cap and word-limit errors are added separately to(C8).

Write C_a for curl incidence divided by a. The semidiscrete equations are

 M_E E_dot=C_a^T H,       M_H H_dot=-C_a E.                (C9)

The electric and magnetic mass matrices are positive diagonal. The
weighted discrete energy is conserved because the two curl cross terms
cancel exactly. For sampled smooth fields, the residuals are

 R_E=M_E E_app_dot-C_a^T H_app,
 R_H=M_H H_app_dot+C_a E_app.                             (C10)

Thus coefficient error multiplies the smooth time derivatives. It is
NOT differentiated and does not incur a delta_h/a factor. Centered Yee
curl differences act only on the sampled smooth fields; each component's
two-derivative truncation error is bounded by a^2 M3/12 when M3 bounds the
relevant third partial derivatives. Boundary cells require the declared
boundary discretization, or a packet supported away from them.

For r=(R_E,R_H), define its dual energy norm by

 ||r||_*^2=a^3 [R_E^T M_E^-1 R_E+R_H^T M_H^-1 R_H].

The exact semidiscrete solution differs from the approximate one by at
most its initial energy-norm mismatch plus integral||r||_*dt. This follows
by differentiating the error energy, Cauchy-Schwarz and the same curl
cancellation. It is the discrete counterpart of PROOF equation(8).
Equations(C6),(C8),(C10) therefore supply an explicit C0 native-source to
Maxwell wave residual, without importing an unproved pointwise derivative
convergence rate. The scalar operator has constant spatial coefficient
and variable time/mass multipliers, so the analogous coefficient comparison
also needs only uniform source error.

The actual OS propagation is then compared to this semidiscrete Maxwell
flow using the independent frequency and plaquette-CCR terms in PROOF(14).
A source finite-spectrum tail and the detector margin are retained. For
a cap-induced ternary-probability error, the electric precision derivative
is bounded by128 on p_E>=1/32, while M_H=8p_B has derivative8. These explicit
bounds propagate the ideal instrument error into(C8); history-dependent
reset discrepancies are still handled by whole-word coupling before
compatibility conditioning rather than renamed an exact IID field law.
