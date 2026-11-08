# Vanishing of the scalar energy-polarization remainder

The same selected ordered source and clock as `PROOF.md` are used.
For a smooth compact test F, let the plateau of the finite weighted-energy
operator B_i contain the test support and every incident bond. The exact
source identity from `spin_energy_polarization` is

\[
 [B_i,Q_p(F)]-[L_p,Q_p(X_iF)]=iR_{i,p}(F),
 \quad
 R_{i,p}(F)={p^2\over2}\sum_x\epsilon_x
 [F(px)-F(px+pe_i)]C_{x,i}^3,                               \tag{1}
\]

where `C_(x,i)^3=S_x^1S_(x+ei)^2-S_x^2S_(x+ei)^1` has norm 1/2.
Its coefficient depends on the actual source current and on the same
spin field; no replacement of one spin by its order parameter is made.

## The test removes the entire central current component

Let gamma be the internal pi rotation about axis two and
`beta_x=tau_x gamma^(x1+x2+x3 mod2)`. Gamma reverses C^3, since it
reverses S^1 and preserves S^2. Hence

\[
 \beta_x(C_{0,i}^3)=\epsilon_x C_{x,i}^3,\qquad
 R_{i,p}(F)=-{p^3\over2}\sum_xD_i^pF(px)\beta_x(C_{0,i}^3),
 \quad D_i^pF(X)={F(X+pe_i)-F(X)\over p}.                  \tag{2}
\]

For a bounded local observable C in the beta-invariant selected state,
the weighted mean-ergodic argument gives

\[
 p^3\sum_x h(px)\beta_x(C)\Omega
 \longrightarrow \left(\int h\right)P_{\rm inv}C\Omega    \tag{3}
\]

for smooth compact h. To prove this form of the argument, first use
indicator boxes. The spectral theorem for the three commuting beta
unitaries makes their normalized sums converge strongly to the joint
invariant projection. Translated and anisotropic boxes have the same
limit because their sides all diverge. Step functions approximate h
in l1 Riemann sums; the operator norm of the remainder is at most
`||C|| p^3 sum|h-h_step|`. This gives (3) by two successive limits.
Neither state purity nor ergodicity is required.

Replace the difference quotient in (2) by `partial_i F`. If F is
supported in `[-R,R]^3`, `R>=1`, `p<=1`, and its second derivative in
direction i is bounded by M2, at most `(2R+1)^3 p^-3` grid vertices
belong to the support box. Active positive-axis bond differences lie
in this set or its one-site translate, so their count is at most twice
that number. Taylor's formula and `||C||=1/2` therefore give
operator error at most `(2R+1)^3 M2 p/4`. But
`integral partial_i F=0`. Equation (3) therefore proves

\[
                    R_{i,p}(F)\Omega\longrightarrow0.     \tag{4}
\]

The possible central current variable is multiplied by zero. Its mean,
variance and scalarity do not enter this conclusion. The discrete
difference coefficients also sum exactly to zero at every finite mesh;
the code checks this on all three directions of the same compact source
test used by the weighted-energy calculation.

## The same vanishing holds on finite time histories

For `|partial_iF|<=M1`, the real bond coefficients of R obey

\[
 \|z\|_\infty\le M_1p^3/2,\qquad
 \sum|z|\le M_1(2R+1)^3,
 \qquad\|R_{i,p}(F)\|\le M_1(2R+1)^3/2.                 \tag{5}
\]

The overlap proof in `PROOF.md` applies to any selfadjoint two-site
operator of norm at most 1/2, including this current. Thus
`||[R,[H_J,R]]||<=198J ||z||_infinity sum|z|=O(Jp^3)`.
Its phase `exp(isR)` has the same energy-kick and time-transport estimate
as the bond coefficients, with `a=M1/2` and `V0=2(2R+1)^3`.
The squared transport error through a source time `T/(Jp)` on energy
at most `EJp` is bounded by

\[
 4|Ts|a\sqrt{99V_0E}\,p+990|T|s^2a^2V_0p^2.             \tag{6}
\]

Its equal-time commutators with existing Q and P Weyls are O(p^2)
and O(p). Move the phase through every factor of any fixed finite
bounded canonical time-history word W_p and then apply (4). This yields
`||(exp(isR)-I)W_p Omega||->0` for each fixed real s. Every word has
the needed source energy bound by the inherited energy-kick theorem.

The norm bound (5) is uniform. Choose a fixed nonzero s small enough
that `|s| ||R||<=1`. The elementary spectral bound
`|exp(isx)-1|>=|sx|/2` in that interval gives

\[
 \|R_{i,p}(F)W_p\Omega\|
 \le {2\over|s|}\|(e^{isR_{i,p}(F)}-I)W_p\Omega\|
 \longrightarrow0.                                     \tag{7}
\]

The zero operator case needs no division. Linear combinations of fixed
finite histories follow by the triangle inequality. The explicit
transport terms in (6) are supplemented by the qualitative microscopic
mean-ergodic convergence in (4). No computable rate for (4) is asserted.

Combining (1) and (7), the difference of the two actual source
commutators annihilates every such history in the scaling limit. This
closes that specific scalar-field-action remainder on the same source,
state and clock. It does not prove separate limits of the unbounded
commutators, a limiting boost generator, vacuum boost invariance, the
corresponding P-field action, Lorentz covariance or a sharp mass shell.

## An actual source control excludes operator-norm vanishing

The state and energy scope is essential. On the same spin source, take
the product input with unit spin directions
`n_x=(epsilon_x,sign(x_i),0)/sqrt(2)`, where `sign(0)=1`.
It has
`epsilon_x <C_(x,i)^3>=[sign(x_i)+sign(x_i+1)]/8`.
For `F(X)=product_j(1-X_j^2)_+^4`, write `b(t)=(1-t^2)_+^4`.
At mesh p=1/n the exact expectation in (1) is

\[
 {p^2\over8}[b(0)+b(-p)]\left(\sum_j b(pj)\right)^2
 \longrightarrow {1\over4}\left(\int_{-1}^1b(t)dt\right)^2
 ={16384\over99225}>0.                                   \tag{8}
\]

The spin-product input is a literal finite-source state for each
supported observable; the fixed vectors used in the selected-ground
proof are not replaced by it. The code checks its exact rational
expectations for four meshes. Equation (8) rules out a state-independent
operator-norm version of (7) on this same source. It does not conflict
with the selected low-energy history limit.
