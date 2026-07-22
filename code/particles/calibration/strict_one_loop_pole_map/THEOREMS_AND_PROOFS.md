# Strict-one-loop W/Z complex-pole map: theorems and proofs

Prepared: 2026-07-20  
Status: proved algebraically under the stated finite-order input contract; **not** an OPH-native physical prediction

## 1. Scope and frozen conventions

Let `eps` be a formal loop-counting parameter. The implementation suppresses
`eps` and stores each `Delta^(1)` with its one-loop factor already included.
For proofs it is useful to restore it.

The renormalized electroweak input at scale `Q` is

\[
\theta(Q)=(g,g',v_F,\ldots),
\qquad
w=\frac{g^2v_F^2}{4},
\qquad
z=\frac{(g^2+g'^2)v_F^2}{4},
\]

with `g>0`, `g'>0`, and `v_F>0`. The Pro package uses

\[
\Gamma^T(s)=s-m_0^2-\Delta^T(s).
\]

The earlier `survival-proof-3` convention was

\[
\Gamma^T(s)=s-m_0^2+\Pi^T(s).
\]

They are identical only after the explicit translation

\[
\boxed{\Delta^T=-\Pi^T.}
\]

This identity is emitted in every hardened receipt. No sign is inferred from
the numerical width.

For a decaying state the energy-pole coordinate is frozen as

\[
\sqrt{s_V}=M_V-\frac{i}{2}\Gamma_V,
\qquad M_V>0,\quad \Gamma_V\geq0.
\]

The strict perturbative coefficient and the exact square root of a truncated
pole are different mathematical objects and are reported separately.

## 2. Tree-mass lemma

**Lemma 2.1 (tree masses).** For the canonically normalized one-doublet
electroweak kinetic term `(D_mu H)^dagger(D^mu H)` and the background
`H=(0,v_F/sqrt(2))^T`, the transverse tree inverse propagators have nonzero
roots

\[
w=\frac{g^2v_F^2}{4},
\qquad
z=\frac{(g^2+g'^2)v_F^2}{4}.
\]

**Proof.** Expanding the covariant derivative around the background gives the
charged quadratic term `g^2 v_F^2 W^+ W^-/4`. In the neutral `(W^3,B)` basis,
the mass matrix is

\[
\frac{v_F^2}{4}
\begin{pmatrix}
g^2&-gg'\\
-gg'&g'^2
\end{pmatrix}.
\]

Its determinant is zero and its trace is `(g^2+g'^2)v_F^2/4`; hence its
eigenvalues are `0` and `z`. The orthogonal massive eigenvector is the `Z`
direction and the zero eigenvector is the photon direction. QED.

This lemma is conditional on the action and field normalization. It does not
prove that an OPH transmutation coordinate is the renormalized `v_F` in that
action.

## 3. Charged strict-one-loop pole theorem

Write the charged transverse inverse propagator as

\[
\Gamma_W^T(s,\epsilon)
=s-w-\epsilon\,\delta_{WW}^{(1)}(s)+O(\epsilon^2),
\]

where `delta_WW^(1)` is holomorphic in a neighborhood of `w` on the frozen
analytic sheet.

**Theorem 3.1 (charged pole coefficient).** If the tree root `w` is simple,
the formal pole through strict one loop is

\[
\boxed{s_W^{[1]}=w+\epsilon\,\delta_{WW}^{(1)}(w).}
\]

After suppressing `eps`, this is the implemented rule

\[
s_W^{[1]}=w+\Delta_{WW}^{(1)}(w).
\]

**Proof.** Put `s_W=w+eps sigma_W+O(eps^2)`. Taylor expansion gives

\[
\begin{aligned}
0&=\Gamma_W^T(s_W,\epsilon)\\
 &=\epsilon\sigma_W
   -\epsilon\left[\delta_{WW}^{(1)}(w)
   +\epsilon\sigma_W\partial_s\delta_{WW}^{(1)}(w)
   +O(\epsilon^2)\right]
   +O(\epsilon^2)\\
 &=\epsilon\left(\sigma_W-\delta_{WW}^{(1)}(w)\right)
   +O(\epsilon^2).
\end{aligned}
\]

The coefficient of `eps` vanishes iff
`sigma_W=delta_WW^(1)(w)`. The derivative/iteration term first contributes at
two loops. QED.

## 4. Neutral strict-one-loop pole theorem

In the tree-level photon/Z basis, restore loop counting in the transverse
neutral inverse matrix:

\[
\Gamma_N^T(s,\epsilon)=
\begin{pmatrix}
s-\epsilon\delta_{AA}^{(1)}(s)&-\epsilon\delta_{AZ}^{(1)}(s)\\
-\epsilon\delta_{ZA}^{(1)}(s)&s-z-\epsilon\delta_{ZZ}^{(1)}(s)
\end{pmatrix}
+O(\epsilon^2).
\]

**Theorem 4.1 (neutral massive pole coefficient).** Assume `z>0`, all four
one-loop entries are holomorphic near `z`, and the tree `Z` root is simple.
Then

\[
\boxed{s_Z^{[1]}=z+\epsilon\,\delta_{ZZ}^{(1)}(z).}
\]

The off-diagonal product contributes no strict-one-loop term.

**Proof.** The determinant through the displayed orders is

\[
D_N(s,\epsilon)=
\left(s-\epsilon\delta_{AA}^{(1)}\right)
\left(s-z-\epsilon\delta_{ZZ}^{(1)}\right)
-\epsilon^2\delta_{AZ}^{(1)}\delta_{ZA}^{(1)}+O(\epsilon^2),
\]

where the final remainder also contains genuine two-loop matrix entries. Put
`s=z+eps sigma_Z+O(eps^2)`. The first factor is `z+O(eps)` and the second is
`eps(sigma_Z-delta_ZZ^(1)(z))+O(eps^2)`. Therefore

\[
D_N(s,\epsilon)=
z\epsilon\left(\sigma_Z-\delta_{ZZ}^{(1)}(z)\right)
+O(\epsilon^2).
\]

Since `z` is nonzero, the coefficient vanishes iff
`sigma_Z=delta_ZZ^(1)(z)`. The product
`delta_AZ^(1) delta_ZA^(1)` carries `eps^2` and cannot be inserted into a
strict-one-loop root. QED.

**Corollary 4.2 (first omitted mixing term).** Eliminating the photon block by
a Schur complement gives the effective massive inverse entry

\[
s-z-\Delta_{ZZ}(s)
-\frac{\Delta_{ZA}(s)\Delta_{AZ}(s)}{s-\Delta_{AA}(s)}.
\]

Its leading off-diagonal term is

\[
-\frac{\Delta_{ZA}^{(1)}(s)\Delta_{AZ}^{(1)}(s)}{s},
\]

which is order two. A future two-loop map must include this term together with
genuine two-loop self-energies and pole-iteration derivatives; adding it alone
would be an inconsistent order mask.

The current numerical fixture lacks `AA`, `AZ`, and `ZA` values. Consequently
it validates the strict `ZZ` root formula but not a complete neutral-matrix
receipt, the photon Ward identities, or a neutral BRST identity.

## 5. Strict mass/width corollary

**Corollary 5.1 (linearized energy-pole coordinates).** Let

\[
s_V=m_{V,0}^2+\epsilon\delta_V+O(\epsilon^2),
\qquad m_{V,0}>0,
\]

and choose the branch continuously connected to `+m_V,0`. Then

\[
\sqrt{s_V}=m_{V,0}+\epsilon\frac{\delta_V}{2m_{V,0}}
+O(\epsilon^2).
\]

Writing the root as `M_V-i Gamma_V/2` gives

\[
\boxed{\delta M_V^{(1)}=
\frac{\operatorname{Re}\delta_V}{2m_{V,0}},}
\qquad
\boxed{\Gamma_V^{(1)}=
-\frac{\operatorname{Im}\delta_V}{m_{V,0}}.}
\]

**Proof.** Apply the Taylor series
`sqrt(m0^2+eps delta)=m0+eps delta/(2m0)+O(eps^2)` and equate real and
imaginary parts. QED.

The decaying-sheet condition at this order is `Im(delta_V)<=0`. The code
rejects a positive imaginary part instead of silently conjugating it.

## 6. Exact coordinate theorem for a truncated pole

**Theorem 6.1 (unique physical square-root coordinates).** Let `s=x+iy` with
`y<=0`, `s` not on the nonpositive real branch cut, and `s!=0`. There is a
unique pair `M>0`, `Gamma>=0` satisfying

\[
s=(M-i\Gamma/2)^2.
\]

It is

\[
M=\sqrt{\frac{|s|+x}{2}},
\qquad
\Gamma=\sqrt{2(|s|-x)}.
\]

**Proof.** Squaring gives

\[
x=M^2-\Gamma^2/4,
\qquad y=-M\Gamma.
\]

Also `|s|=M^2+Gamma^2/4`. Adding and subtracting the first identity gives the
displayed nonnegative square roots. The sign of `y` fixes the lower-half-plane
root. Positivity of `M` removes the remaining sign ambiguity. QED.

**Finite-order warning.** Applying this nonlinear formula exactly to
`s^[1]=m0^2+Delta^(1)` resums powers such as `(Delta^(1))^2`. Those are useful
coordinates of the truncated complex number, but they are not computed
two-loop physics. Therefore a finite-order gauge/Nielsen comparison must use
the coefficients in Corollary 5.1, not numerical equality of exact square-root
readouts.

## 7. Conditional gauge-independence theorem

Let `xi` denote any gauge-fixing parameter and let `D(s,xi,eps)` be the
relevant charged inverse entry or neutral determinant.

**Theorem 7.1 (Nielsen factorization implies pole independence).** Suppose in
a neighborhood of a simple pole that

\[
\partial_\xi D(s,\xi,\epsilon)
=C(s,\xi,\epsilon)D(s,\xi,\epsilon),
\]

with `C` regular at the pole. Then the exact pole `s_star(xi)` is independent
of `xi`.

**Proof.** Differentiate `D(s_star(xi),xi,eps)=0`:

\[
0=\partial_\xi D+(\partial_sD)\frac{ds_\star}{d\xi}.
\]

At the pole, factorization makes the first term zero. Simplicity gives
`partial_s D != 0`, so `ds_star/dxi=0`. QED.

**Corollary 7.2 (finite-order version).** If the factorization identity holds
as a formal series through order `eps^n`, the regularity assumptions hold
coefficientwise, and the pole equation is solved and re-expanded consistently
through that same order, then every pole coefficient through `eps^n` is gauge
independent.

This theorem explains what a gauge/BRST receipt must establish; it is not the
receipt itself. The supplied package has no symbolic extended-BRST identity,
no certified general-gauge comparison, and no complete neutral matrix. Its
`finite_order_BRST_ST_Ward_Nielsen_receipt` gate therefore remains false.
General all-order theory cannot certify that a particular finite diagram,
counterterm, tadpole, and analytic-sheet implementation is complete.

## 8. Conditional FJ/tadpole conversion theorem

Let `p_L` and `p_F` be two complete renormalized parameter coordinates related
at one loop by

\[
p_{L,a}=p_{F,a}+\epsilon\,\delta p_a^{(1)}+O(\epsilon^2).
\]

This equation fixes the orientation of `delta p`.  Write the same physical
observable in the two charts as

\[
O_X(p)=O^{(0)}(p)+\epsilon O_X^{(1)}(p)+O(\epsilon^2),
\qquad X\in\{F,L\}.
\]

**Theorem 8.1 (one-loop reparameterization).** For any differentiable tree
quantity `O^(0)(p)` with a consistently transformed one-loop coefficient,

\[
\boxed{O_F^{(1)}(p_F)=O_L^{(1)}(p_F)
+\sum_a\delta p_a^{(1)}
\frac{\partial O^{(0)}}{\partial p_a}(p_F).}
\]

**Proof.** Coordinate independence means
`O_F(p_F)=O_L(p_L)+O(epsilon^2)`.  Substitute
`p_L=p_F+epsilon delta p` and Taylor-expand.  The tree term contributes
`epsilon delta p . grad O^(0)(p_F)`, whereas shifting the argument of the
one-loop coefficient contributes only at order two.  Equating the coefficients
of `epsilon` gives the displayed identity. QED.

**Corollary 8.2 (scheme-equivalent strict pole).** If the finite transformation
is applied to every tree parameter, internal mass, counterterm,
wave-function term, and tadpole contribution before re-expansion, the direct
explicit-tadpole FJ route and the completely converted tadpole-free route
give the same strict-one-loop pole coefficient up to `O(eps^2)`.

This is a conditional equivalence theorem. The archive does not contain the
finite shift vector, its derivation, or two independently implemented routes.
It also documents that the current fixture uses a Landau-gauge
effective-potential VEV, not an OPH-selected FJ `v_F`. Hence it cannot be used
as the requested OPH-to-renormalized-VEV/FJ conversion receipt.

## 9. Certified pole-isolation theorem

**Theorem 9.1 (Rouche isolation certificate).** Let `D` and `D_hat` be
holomorphic inside and on a closed contour `C`. If

\[
|D(s)-\widehat D(s)|<|\widehat D(s)|
\quad\text{for every }s\in C,
\]

then `D` and `D_hat` have the same number of zeros inside `C`, counted with
multiplicity. In particular, if `D_hat` has exactly one simple zero there,
the contour certifies one pole of `D`.

**Proof.** This is the direct application of Rouche's theorem to
`D_hat+(D-D_hat)`. QED.

For a machine receipt, both sides must be enclosed on the entire contour by
certified complex balls or intervals. A sampled grid of small residuals is not
the hypothesis. The current strict-one-loop archive evaluates point
coefficients and does not claim a Rouche certificate.

## 10. Source-covariance transport theorem

Let `zeta` be a target-independent source variable with a licensed joint law,
mean `zeta_0`, and covariance `C_zeta`. Let the source, matching, and pole maps
be differentiable at the frozen point, with Jacobians `J_source`, `J_match`,
and `J_pole`. Represent complex poles by their ordered real and imaginary
parts.

**Theorem 10.1 (first-order source covariance).** The linearized source-induced
pole covariance is

\[
\boxed{C_{\mathrm{pole,source}}=
J_{\mathrm{pole}}J_{\mathrm{match}}J_{\mathrm{source}}
C_\zeta
J_{\mathrm{source}}^T J_{\mathrm{match}}^T
J_{\mathrm{pole}}^T.}
\]

**Proof.** Linearization gives
`delta y=J_pole J_match J_source delta zeta`. Taking
`E[delta y delta y^T]` gives the displayed congruence transformation. QED.

For a scalar pole equation `D(s,zeta)=0` with a simple root, the implicit
Jacobian is

\[
\frac{\partial s_\star}{\partial\zeta_i}
=-\frac{\partial_{\zeta_i}D}{\partial_sD}\bigg|_{s=s_\star}.
\]

This covariance does not include branch ambiguity, matching truncation, pole
truncation, scale variation, or numerical enclosure radius. Those must remain
separate frozen uncertainty components unless an independently justified
joint probabilistic model exists. If and only if the source is proved unique
and exact may one set `Law(zeta)=delta(zeta-zeta_0)` and `C_zeta=0`; this makes
source variance zero, not total theoretical uncertainty zero.

The present fixture explicitly reports source covariance as `not_licensed`.

## 11. Evidence non-self-attestation theorem

**Theorem 11.1 (untrusted booleans cannot certify external evidence).** Let a
producer receive an untrusted subject packet containing a boolean claim `b_E`
for an external property `E`. If a verifier checks only syntax,
self-consistency, and `b_E=true`, but resolves no independent witness whose
content is bound to the exact subject, then verifier acceptance does not imply
`E`.

**Proof.** Choose any subject for which `E` is false. Set its untrusted field
`b_E=true` and preserve every relation the verifier actually recomputes. Since
the verifier never evaluates a predicate that distinguishes the false subject
from a true one, it follows the same accepting path. Thus there exists an
accepted false instance, so acceptance cannot imply `E`. QED.

**Concrete regression.** The supplied upstream v1 checker bound a copied
fixture hash but did not compare the receipt's parameters or self-energies to
the fixture, and it trusted evidence booleans. A self-consistent receipt for a
different pole could therefore pass while claiming physical promotion. The
hardened integration now:

1. makes `physical_promotion_allowed=false` a schema constant for v1;
2. records caller booleans only as `unverified_evidence_claims`;
3. recomputes all v1 evidence gates, with every external gate false;
4. binds IDs, source, parameters, contribution mask, covariance, tolerances,
   and self-energy coefficients to the exact fixture;
5. caps tolerances; and
6. tests unrelated-receipt, tolerance-inflation, and promotion-flip attacks.

A future physical promotion must be performed by a separate aggregate
verifier that resolves hash-bound FJ, source-matching, covariance,
gauge/BRST, clock, and ancestry receipts and binds all of them to one exact
runtime subject.

## 12. Evaluated conditional fixture

The fixture reconstructs an archived SMDR v1.3 order-1 row at `Q=160 GeV`.
It is a regression fixture derived from differences of backend outputs, not an
independent self-energy calculation.

The stored one-loop coefficients are

\[
\Delta_{WW}^{(1)}(w)=
(64.72322786968198-160.53275277304510i)\ {\rm GeV}^2,
\]

\[
\Delta_{ZZ}^{(1)}(z)=
(-33.37572257933425-218.29276180643907i)\ {\rm GeV}^2.
\]

The resulting strict-one-loop truncated poles are

\[
s_W^{[1]}=
(6459.842027569383-160.532752773045i)\ {\rm GeV}^2,
\]

\[
s_Z^{[1]}=
(8222.835212344102-218.292761806439i)\ {\rm GeV}^2.
\]

| Quantity | W | Z |
|---|---:|---:|
| tree mass (GeV) | 79.969486678981 | 90.863694261919 |
| strict `M0+delta M` (GeV) | 80.374161202712 | 90.680036075608 |
| strict width (GeV) | 2.007425074735 | 2.402420059845 |
| exact `M` coordinate of truncated `s` (GeV) | 80.379345721647 | 90.687836666971 |
| exact `Gamma` coordinate of truncated `s` (GeV) | 1.997189095430 | 2.407078720028 |

These numbers certify the implementation of Theorems 3.1 through 6.1 for the
fixture. They do **not** certify target independence, an OPH/FJ VEV, threshold
matching, a source law, an independent loop engine, general-gauge BRST/Nielsen
identities, or a physical OPH clock.

## 13. Exact completion status

The strict finite-order pole-map kernel is now explicit, tested, and
fail-closed. Its current scientific classification remains

```text
CONDITIONAL_STRICT_1L_POLE_MAP_NOT_OPH_NATIVE_PHYSICAL
```

The following implication is proved:

```text
complete renormalized strict-1L input packet
  + frozen sign/order/sheet conventions
  -> strict W/Z complex poles and separated M/Gamma readouts.
```

The antecedent has not been produced from OPH. In particular, this package is
not the missing OPH-to-FJ conversion, source-matching certificate, source
covariance, or independent gauge/BRST receipt. It supplies the exact interface
and proof obligations those future artifacts must satisfy.

## 14. Primary references used to check the physics framing

- P. Gambino and P. A. Grassi, *The Nielsen Identities of the SM and the
  definition of mass*, [arXiv:hep-ph/9907254](https://arxiv.org/abs/hep-ph/9907254).
- D. Dūdėnas and M. Löschner, *Vacuum expectation value renormalization in the
  Standard Model and beyond*, [arXiv:2010.15076](https://arxiv.org/abs/2010.15076).
- S. Dittmaier and H. Rzehak, *Electroweak renormalization based on
  gauge-invariant vacuum expectation values of non-linear Higgs
  representations: 1. Standard Model*,
  [arXiv:2203.07236](https://arxiv.org/abs/2203.07236).
- S. P. Martin and D. G. Robertson, *Standard Model parameters in the
  tadpole-free pure MS-bar scheme*,
  [arXiv:1907.02500](https://arxiv.org/abs/1907.02500).
- S. Willenbrock, *Mass and width of an unstable particle*,
  [arXiv:2203.11056](https://arxiv.org/abs/2203.11056).
