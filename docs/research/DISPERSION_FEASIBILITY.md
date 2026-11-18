# Can the fixed edge dispersion support the proposed photon test?

This note contributes to the decision in [#1025](https://github.com/FloatingPragma/observer-patch-holography/issues/1025).
Its calculation addresses the kinematic uncertainty identified in the
[independent feasibility audit](https://github.com/FloatingPragma/observer-patch-holography/issues/1026#issuecomment-5977983054).
The executable check is [dispersion_feasibility](../../code/dispersion_feasibility/README.md).
It runs no flux simulation, Lean expansion or held-out-data comparison.

## Decision contribution

The proposed shared photon/electron/positron dispersion has an asymmetric
pair-production threshold. Its leading threshold differs from the exact,
direction-dependent, three-dimensional kinematic threshold by less than
`1e-13 eV` over the stated ultra-high-energy domain. The photon-only threshold
has the same error bound. The difference between the two variants therefore
survives the exact cosine law and arbitrary outgoing directions. Numerical
witnesses independently check conservation and stationarity; the error bound
comes from the analytic argument in the appendix, not an optimizer's success.

At hard momentum `1e19 eV`, the minimum soft-photon energy is about
`0.0001035098791 eV` for shared photon/electron/positron dispersion and
`0.1367736429 eV` for photon-only dispersion. Equal sharing overestimates the
shared threshold by a factor exceeding 990. Neither cancellation nor a
three-quarter rescaling transfers the photon-only observational bound.

There is also a concrete detector-consistency difference. In the photon-only
variant an electron of momentum `3e17 eV` can radiate a photon with momentum
magnitude `2e17 eV` in vacuum. The exact cosine law permits this channel in
every incident direction. In the shared-shell variant, the chord triangle
identity forbids both vacuum Cherenkov emission and photon decay into a
massive pair. These are kinematic statements; no emission rate is inferred.

**The published photon-only Auger coefficient bound does not decide the
shared-shell candidate.** Nor does a threshold shift establish a flux
exclusion for the photon-only candidate. A source/response model and an
interaction law are necessary for that inference. This contribution supplies
no positive start decision for #1026, no OPH exclusion and no empirical
support. It resolves a specific approximation question and identifies a
physical process that a proposed response calculation must include.

## Fixed hypotheses and input ancestry

Use natural units, additive energy and three-momentum, and the declared
cosmic-background rest frame. For the complete thirty-row unit edge orbit
`W`, set

\[
 \Lambda_a(p)=\frac1{5a^2}\sum_{w\in W}[1-\cos(a w\cdot p)],
 \qquad \Omega(p)=\sqrt{\Lambda_a(p)},\qquad d=a^2/20.
\]

The physical-frequency interpretation of this spatial symbol is a hypothesis.
A declared positive kinetic term gives two transverse oscillators with
frequency `Omega`; the quadratic Hamiltonian is nonnegative. This reuses the
[free-photon construction](../../code/a5_fingerprint/runtime/fz12_free_photon_hamiltonian_receipt.json).
Calling its parameter physical time and its modes photons is additional
physical content. The law contains no extra quadratic derivative terms.

The two controls differ only in the massive shell:

| Input | Photon-only control | Shared photon/electron/positron control |
| --- | --- | --- |
| Photon shell, both transverse polarizations | `Omega(p)` | `Omega(p)` |
| Electron and positron shell | `sqrt(m^2+|p|^2)` | `sqrt(m^2+Lambda_a(p))` |
| Meaning of “shared” | Not applicable | These three species only, not hadrons or all fields |

The massive shared shell is a specified hypothesis, not a derived electron
action or an implication of a photon principal symbol. Its full mass dependence
matters to the error proof. All three shells use positive energies. No electron
spin action, pair-production vertex, background population or quantum
transition probability is supplied by the free shells.

The nominal inputs are `P=1.6309682094039593`,
`l_P=1.616255e-35 m`, `hbar*c=1.973269804e-7 eV m`,
the CODATA electron rest energy `m` and the additional scale hypothesis `a=sqrt(P) l_P`.
They give `a≈1.04603e-28 eV^-1` and `d≈5.47094e-58 eV^-2`.
P inherits the exposed alpha calibration; the length conversion follows the
[existing diagnostic](../../code/a5_fingerprint/runtime/fz12_auger_threshold_diagnostic_receipt.json).
The error certificate
controls the approximation at these inputs, not their measurement uncertainty.

Hard momentum has magnitude `K` in `[1e17,1e20] eV`; a soft photon has momentum
`s u`, with `0<=s<=200 eV` and unit u. All final momenta obey `a|p|<=1`.
This explicit validity domain excludes unspecified ultraviolet branches.
The threshold minimizes over the soft direction and the full three-dimensional
final-pair domain. A head-on pair supplies an upper bound; the lower bound
allows every direction. At the minimizing energy all momenta are much
smaller than the cutoff. The continuum momentum and preferred frame are
physical hypotheses, not identified OPH record coordinates.

## Threshold and exact controls

The [joint-threshold receipt](../../code/a5_fingerprint/runtime/fz12_joint_threshold_receipt.json)
contains the leading relation. Write `y=x(1-x)`. For shared coefficients,
minimization reduces to the strictly convex function

\[
 \epsilon_0(K)=\min_{0<y\le1/4}
       \frac{m^2/y+3dK^4 y}{4K},\qquad
 y_* =\min\left(\frac14,\frac{m}{\sqrt{3d}K^2}\right).
\]

Thus, for `K_asym=(16m^2/(3d))^(1/4)≈2.24618e17 eV`,

\[
\epsilon_0(K)=
\begin{cases}
m^2/K+3dK^3/16,&K\le K_{\rm asym},\\
m\sqrt{3d}K/2,&K\ge K_{\rm asym}.
\end{cases}
\]

The smaller share is `x=2y*/(1+sqrt(1-4y*))`. The photon-only control instead
has `epsilon_0=m^2/K+dK^3/4` and equal sharing. Asymmetric thresholds are
established Lorentz-violating kinematics; see
[Mattingly, Jacobson and Liberati](https://arxiv.org/abs/hep-ph/0211466).
The contribution here is the uniform exact-symbol enclosure for this fixed
candidate and the implications for the proposed comparison.

| Hard momentum (eV) | Photon-only soft threshold (eV), rounded | Shared soft threshold (eV), rounded |
| ---: | ---: | ---: |
| `1e17` | `2.7479729e-6` | `2.7137795e-6` |
| `1e18` | `1.3703474e-4` | `1.0350988e-5` |
| `1e19` | `1.3677364e-1` | `1.0350988e-4` |
| `1e20` | `1.3677362e2` | `1.0350988e-3` |

The unrounded leading values enclose each exact threshold to `1e-13 eV`;
rounding in this table is not that enclosure. The retained report contains
30 exact-shell collinear witnesses at five momenta and three directions.
An independent checker uses a different edge-orbit construction. Repeating
representative solves at 70 and 100 decimal digits checks numerical stability.
It does not replace the global proof.

For vacuum channels define complex features
`Phi_w(p)=(exp(i a w.p)-1)/(sqrt(10)a)`. They satisfy
`||Phi(p)||=Omega(p)` and
`Phi(p+q)=Phi(p)+U(p)Phi(q)`, with diagonal unitary `U`.
Hence `Omega(p+q)<=Omega(p)+Omega(q)`. Adding the mass as an extra orthogonal
coordinate proves

\[
 E_m(p+q)<E_m(p)+\Omega(q) \quad (\Omega(q)>0),\qquad
 \Omega(p+q)<E_m(p)+E_m(q).
\]

These exclude shared-shell electron radiation and photon pair decay under
additive conservation. Within `a|q|<=1`, nonzero q has positive frequency.
For Lorentz-invariant leptons, photon decay is excluded by `Omega(k)<=|k|`,
as recorded in the joint-threshold receipt. Vacuum electron radiation is
different. At leading order its available energy for emitted share `z` is

\[
 \frac{zK}{2}\left[d z^2K^2-\frac{m^2}{K^2(1-z)}\right].
\]

The leading onset occurs at `z=2/3`,
`K=(27m^2/(4d))^(1/4)≈2.38243e17 eV`.
The exact existence claim uses `K=3e17`, emitted magnitude `2e17` and
recoil `1e17`, rather than asserting that approximate onset is exact.
The collinear surplus exceeds `7e-7 eV` in every direction, using the
same cosine and square-root remainder bounds. Opposite emission has negative
surplus. Continuity in emission angle supplies an energy-conserving channel.
This is a concrete reason to check shower response. Rates require a vertex.

## Consequence for the proposed rejection observation

[Auger sections 5.2 and 5.3](https://arxiv.org/html/2112.06773#S5.SS2)
give different outcomes for two cosmic-ray source scenarios. The reference
scenario supplies no electromagnetic bound; the alternative proton-containing
scenario supplies the quoted coefficient constraint. The direct
electromagnetic scan states no coefficient confidence level. Its result cannot
inherit the five-sigma confidence assigned to a separate hadronic fit.

Read physically, the published bound and the fixed-scale branch disagree.
With `a = sqrt(P) l_P`, the branch has `delta_gamma2 = -a^2/20 =
-5.47e-58 eV^-2`. The alternative scenario, whose sources include a
subdominant proton component extending to `1e20 eV`, gives
`delta_gamma2 > -1e-58 eV^-2` for photon-only Lorentz violation with
Lorentz-invariant electrons and positrons, ordinary additive conservation,
Breit-Wheeler opacity through the shifted threshold and the cosmic-background
rest frame. Under those published premises the photon-only fixed-scale branch
is excluded by a factor 5.47 in the coefficient; equivalently the bound
requires `a < 0.546 l_P`, against `a = 1.277 l_P`. The register reading is
separate: FZ-12 is unarmed and Auger-era data precede its freeze, so the
[frozen-prediction ladder](../FROZEN_PREDICTION_LADDER.md) records no verdict.

If both source scenarios are admitted, that comparison has no rejection
uniform over the nuisance set. This is an explicit surviving scenario, not an
uncertainty invented by the threshold calculation. It does not preclude a
different dataset or a justified narrower source model.

Free dispersion also fails to fix a pair-production rate. For example, in a
photon-only effective theory with a supplied Dirac electron, the Hermitian,
gauge-invariant interaction `g (bar(psi) psi) F_mu_nu F^mu_nu / M^3` changes
the tree-level pair-production amplitude while leaving all quadratic vacuum
shells unchanged. The restriction against extra quadratic derivative terms
does not exclude this interaction. Setting such coefficients to zero and
choosing a particular electromagnetic action are additional model hypotheses.
This tree-level counterexample is not a loop-renormalization or ultraviolet
completion claim. The need to specify interaction vertices and photon-shower
response is also explicit in
[Rubtsov, Satunin and Sibiryakov](https://arxiv.org/abs/1204.5782).

The finite decision for this contribution is therefore: **kinematics is
controlled; a rejection test is not established by the proposed transfer of
the Auger coefficient bound.** Before a positive #1025 selection, the decision
note must name an interaction/response prescription and an allowed source set
that actually permits an exclusion. If it retains both published scenarios
without further distinguishing evidence, this particular upper-limit test is
inconclusive. That is a stopping result for the comparison, not a reason to
create a prerequisite theorem queue or a verdict on all possible OPH models.

Concise corrections prepared for any decision note or summary making these
claims are:

| Claim to correct | Replacement |
| --- | --- |
| The pixel-scale coefficient is excluded by a factor 5.5 | Under the published photon-only premises of the proton-containing scenario, the fixed-scale branch is excluded by a factor 5.47; excluding the shared-shell variant or OPH requires the applicable source, interaction and response assumptions. |
| Universal propagation cancels the effect or is automatically unconstrained | Shared photon/electron/positron shells have a controlled asymmetric threshold; their flux constraints require their own calculation. |
| Two extra physical hypotheses arm FZ-12 | This fixed-scale extension is separately hypothesized. The immutable FZ-12 physical-attachment and comparison rules are unchanged. |
| Fixing the free coefficient fixes the measurable flux | The free shells constrain kinematics. Interaction vertices, source populations and detector response determine rates and recorded flux. |

These are prepared corrections, not edits to papers, book or deployed summaries.
All cited Auger information is exposed data. No prospective score, support
claim or change to FZ-13/FZ-15 or another frozen registration is made.

## Appendix: direction-uniform threshold enclosure

The support has unit vectors and exact moments
`sum(w.p)^2=10|p|^2`, `sum(w.p)^4=6|p|^4`.
Taylor inequalities for the cosine imply, for `ar<=1`,

\[
 r^2-dr^4\le\Lambda_a(rn)
 \le r^2-dr^4+a^4r^6/600.
\tag{1}
\]

The sixth moment is bounded by the fourth because `|w.n|<=1`. The same
support gives `Omega<=r`, unit Lipschitz continuity, and
`1-2dr^2<=partial_r Omega(rn)<=1` on this domain.
The lower derivative follows from `sin t>=t-t^3/6` for `0<=t<=1`, with odd
extension for the signed projections, and division by `2Omega<=2r`.

Let `Q=|K n+s u|` and let `M(Q,v)` be the minimum of the two exact massive energies
over all final vectors summing to `Qv=K n+s u` in the stated validity domain. A feasible
pair has energy at most `Q+2m`. By (1), each minimizing momentum is less than
`(Q+2m)/sqrt(19/20)<2Q`; it lies inside the cutoff. Compactness and continuity
therefore supply a minimum and a continuous energy function in this range.
One common compact domain is `|p|<=3K`; its other leg has norm below `5K`,
which is inside the cutoff. The energy bound puts every minimum inside it.

For shared shells put `f(r)=sqrt(m^2+r^2-dr^4)`.
It is increasing for `ar<=1`. Orthogonal projection followed by clamping onto
the segment `[0,Qv]` reduces both final norms. Thus

\[
 \min_{0\le p\le Q}[f(p)+f(Q-p)]\le M(Q,v).
\tag{2}
\]

An exact collinear pair supplies the upper bound. This sandwich accounts for
noncollinear and direction-dependent minimizers without declaring them
collinear.

For the uniform arithmetic take `Kmax=1e20+200`, `Qmin=9e16`, `b=1e13`,
`S=200`, `5e5<m<6e5`, and `5e-58<d<6e-58`, all in eV units as appropriate.
Define outward bounds

\[
 H=\frac{(20d)^2K_{\max}^5}{600},\quad
 T=\frac12\left(\frac{m^4}{b^3}+2m^2dK_{\max}+d^2K_{\max}^5\right),
 \quad R_\gamma=H+d^2K_{\max}^5/2.
\tag{3}
\]

For `0<=r<=Kmax`, (1) gives `0<=E_m(rn)-f(r)<=H`.
For `b<=r<=Kmax`, writing `z=m^2/r^2-dr^2`,

\[
 \left|f(r)-r-\frac{m^2}{2r}+\frac{dr^3}{2}\right|
 =\frac{r z^2}{2(\sqrt{1+z}+1)^2}\le T.
\tag{4}
\]

The analogous photon error is at most `R_gamma`.
No small-leg expansion is assumed: if a leg has `r<b`, the chord inequality
gives a massive pair gap above `Omega(Qv)` of at least
`m^2/(2b+m)`. Replacing both exact energies by f loses at most `2H`.
With outward constants this lower gap exceeds `0.01249 eV`. A collinear pair
at the leading optimum has gap at most
`2m^2/Qmin + m sqrt(3d) Kmax + 2(H+T)+R_gamma < 0.002589 eV`.
The leading optimum itself has both legs above b: on its asymmetric branch
`p>=m/(sqrt(3d)Q)>1e14`, and the symmetric branch has `p=Q/2`.
Consequently the minimum in (2) has neither leg below b, and (4) applies to it.

Write `D=M-Omega(Qv)` and `D0=2 epsilon_0`. Equations (2)--(4) prove

\[
 |D(Q,v)-D_0(Q)|\le R_U=2(H+T)+R_\gamma.
\tag{5}
\]

For photon-only leptons the exact minimum is `sqrt(Q^2+4m^2)`; the same
argument gives `R_L=2m^4/Qmin^3+R_gamma`.
Valid Lipschitz bounds for D0 are

\[
 L_U=2m^2/Q_{\min}^2+9dK_{\max}^2/8+m\sqrt{3d},\qquad
 L_L=2m^2/Q_{\min}^2+3dK_{\max}^2/2.
\]

For head-on incidence the exact energy residual is
`F(s)=M(K-s,n)-Omega(Kn)-Omega(sn)`. Using (5) and the radial derivative bound,

\[
 |F(s)-(D_0(K)-2s)|\le
 B_j=R_j+S(L_j+2dK_{\max}^2+dS^2).
\tag{6}
\]

For any soft direction u, the unit Lipschitz bound gives
`Omega(Kn+su)-Omega(Kn)>=-s`, and `|Q-K|<=s`.
Thus its residual is at least `D0(K)-2s-Rj-Lj*s`. Below `(D0-Bj)/2`
the energy is insufficient in every soft direction. At `(D0+Bj)/2` the head-on minimum
final energy is no greater than the incoming energy. Higher-energy final
pairs in the interior of the validity domain and continuity supply energy
conservation. Hence the first kinematic threshold lies between these two
values even after optimizing the soft direction; the conclusion does not need
global monotonicity of F. For the higher-energy endpoint one can use final
momenta `2Kn` and `-(K+s)n`; their energy exceeds the incoming energy by (1).
The endpoints are positive and below S throughout the band. Converting s
to soft energy costs at most `dS^3`.

Exact rational arithmetic with outward `sqrt(3d)<4.3e-29` gives threshold
errors below `8.240e-15 eV` for shared shells and `4.201e-15 eV` for
photon-only shells. Replacing K by the exact hard energy in the leading
formula costs less than `2.701e-15 eV`, since
`0<=K-Omega(Kn)<=dK^3` and `|epsilon_0'|<=L_j/2`.
Their sums are below the declared `1e-13 eV` bound. The bound is uniform in
incident direction and all admissible final momenta; it is not a propagated
physical-parameter uncertainty or a flux likelihood.
