# Gauge and Higgs inputs versus computed source responses

Several finite interaction coefficients are computable and operationally
readable in the declared source. They are retained here alongside the
conditional gauge-ray solution and the exact electroweak mass form.
Their identification with physical renormalized couplings and the Higgs
phase requires a response map on the same spatial source, state and clock.
Architectural coefficients may be declared; their uniqueness is not a
prerequisite for calculating that source's outputs.

## Finite couplings supplied by the configured sources

The [native color-link source](../source_strong_coupling/PROOF.md) sets
\(J=1\) and electric energy \(\epsilon=1\). Observers allowed to read
its specifications therefore have the exact ratio \(\epsilon/|J|=1\).
Its complete endpoint detector also determines the coefficients. On every
retained matter block of dimension \(d_r\), the actual Gauss-singlet
source states give

\[
 H_r=\begin{pmatrix}0&J\sqrt{d_r}\\J\sqrt{d_r}&\epsilon\end{pmatrix},
 \qquad p_r(t)=\frac{4d_rJ^2}{\epsilon^2+4d_rJ^2}
 \sin^2\!\left(\frac{\sqrt{\epsilon^2+4d_rJ^2}\,t}{2}\right).
\]

The complete probabilities determine
\(p_r''(0)=2d_rJ^2\) and
\(p_r''''(0)=-2d_rJ^2(\epsilon^2+4d_rJ^2)\).
The inherited `readout.witness` encloses \(J^2\), \(\epsilon^2\)
and their ratio using rational finite-time errors on the original two
endpoint outcomes. This package reuses that witness and checks the actual
source-function defaults. The sign of \(J\) is an endpoint phase
convention for these records. The assumed \(\epsilon\ge0\) fixes the
electric sign.

This operational coefficient is not missing. The source's occupied-link
energy is uniform across its five representations, however, and cannot
equal any expression \(aC_3+bC_2+cY^2\) on all of them. Thus replacing
it by a heat-kernel Casimir coefficient would change the action. The
source normalization and detector do not by themselves identify
\(\epsilon/|J|\) with a running QCD coupling.

The [common scalar source](../common_source/scalar_coupling.py) specifies
the neutral scalar-charge coupling \(\gamma=1/4\) and executes its
interaction through the twelve-port pulse compiler. For its dimension-two
scalar control,

\[
 H_{\rm int}=\tfrac14 I_2\otimes X\otimes
 \operatorname{diag}(1,-1,1,-1).
\]

This exact finite matrix does not commute with the same register's hopping
term. The owning all-column execution tests check its source words and
complete joint evolution. The present certificate checks the parameter
and independent dimension-two matrix form; its float residual is reported
as a diagnostic, and it does not repeat the source tape. This neutral
oscillator coupling is not an electroweak gauge coupling or a Higgs quartic
merely because it acts on a scalar.

The [quantum-Higgs source](HIGGS_RECORDS.md) supplies another exact
operational coupling. Its complete 307-state response has

\[
 g_{\rm pair}^2=4\operatorname{Tr}(Y^\dagger Y),\qquad
 \Delta^2=d^2+4g_{\rm pair}^2,
 \qquad T(t)=\frac{4g_{\rm pair}^2}{\Delta^2}
                 \sin^2(\Delta t/2).
\]

After the complete count inverse determines \(\Delta\), preparation
of the empty state gives

\[
 g_{\rm pair}^2=
 \frac{\Delta^2 T(t)}{4\sin^2(\Delta t/2)},\qquad
 d^2=\Delta^2-4g_{\rm pair}^2.
\]

For the exact fixture \(Y=\operatorname{diag}(1,2,3)\), \(d=1\),
the records at \(t=\pi/45\) give \(T=56/225\), \(\Delta=15\),
\(g_{\rm pair}^2=56\) and \(d^2=1\). These population records do
not determine the sign of a general \(d\); the fixture's nonnegative
cost is declared. A reproducible mixed preparation with a measured
\(\alpha=\operatorname{Tr}(\rho P_0)>0\) uses
\(T=J_{01}/\alpha\), with complete count and preparation errors
propagated through the same expression. No accepted-only history law is
substituted for the original counts.

For the configured coherent analyzer with pass amplitude one and unit
phase \(z=1\), meaning \(\theta=0\), the actual three-family coefficient gives

\[
 g_{\rm pair}^2=\frac{425+100\sqrt2}{1944},\qquad
 \Delta^2=\frac{911+100\sqrt2}{486}\quad(d=1).
\]

The certificate calculates this norm from the source's exact family
matrix and reuses the full-sector response theorem for that matrix. It
does not introduce a numerical coupling fit. These results specify a
local scalar-pair interaction and its spectral detector; they do not
identify \(g_{\rm pair}\) with an electroweak gauge coupling or infer
a continuum Higgs background from scalar occupation.

## Repair response and a controlled interacting spin limit

The [full repair-complement calculation](../common_source/full_repair_response/coupling.py)
eliminates all complementary ports of its specified one-proposal
Dirichlet operator. It gives an exact rank-one family kernel
\(S_0=\operatorname{diag}(0,0,\kappa_0)\) and

\[
 Z=I+(U_{\rm harm}-B)^{\mathsf T}(U_{\rm harm}-B)/G\ge I.
\]

Positive canonical normalization preserves the two zero family
coefficients. The unprotected uniform-repair control instead gives
\((5-\sqrt5)I_3/60\), with equal family coefficients. These are
computed source responses. Assigning them as Lorentzian Yukawa matrices
is a separate coupling prescription; the assignment does not yield
three nonzero unequal charged-lepton masses. The receipt reuse authenticates
its scientific source bindings. Research-note bindings are excluded from
evidence under the repository's code-evidence policy; the elimination
is not replayed by this package.

The [strong-interaction refinement](../strong_coupling_refinement/README.md)
also derives a finite spin exchange, rather than assuming one:

\[
 H_{\rm spin}=-PT^2P/U,
 \qquad J_{xy}=4t_{xy}^2/U=c^2/U.
\]

The owning bounds control its isometry and complete real-time dynamics
on each fixed graph. For example, \(c=128\), \(U=128^2\) gives
\(J_{\rm spin}=1\) with the displayed finite errors. This is an
interacting spin limit of its source. The word “strong” does not identify
that exchange with QCD's gauge coupling or create a chiral relativistic
continuum.

## Reuse of the conditional gauge-ray solution

The [coefficient-selection producer](../observer_particles/coefficient_selection.py)
has an isolated positive solution for the adopted gauge-ray and entropy
constraints. The exact ray is

\[
 x=c(10/3,2,2)+t(41/6,-19/6,-7),\quad
 x=(\alpha_Y^{-1},\alpha_2^{-1},\alpha_3^{-1}).
\]

With \(A=1/\alpha_{\rm Chain}\), same-scale tree electromagnetic
matching gives \(x_Y+x_2=A\), hence

\[
 c=(3A-11t)/16,\qquad
 x=(5A/8+109t/24,\;3A/8-109t/24,\;3A/8-67t/8).
\]

The positive domain is \(0\le t<3A/67\), and the adopted entropy
equation is \(\ell_2(1/x_2)+\ell_3(1/x_3)=P/4\). The source uses
\(\alpha_{\rm Chain}=(79+48\sqrt2+42\sqrt3+3\sqrt6)/31104\)
and \(P=(1+\sqrt5)/2+\sqrt\pi\alpha_{\rm Chain}\). Its representation
weights are \(d_r e^{-4\pi^2\alpha C_r}\), with one power of
representation dimension. This trace convention is part of the
conditional construction.

The owning numerical solution is approximately
\((\alpha_Y,\alpha_2,\alpha_3)=(0.0092278213,0.0348820290,0.1053610898)\),
with same-scale tree \(\sin^2\theta_W=0.2092009203\).
Its interval monotonicity and tail analysis establish the stated positive
branch; this package authenticates the producer and all eleven inherited
source bindings, checks the exact ray elimination and reuses that result.
It does not repeat the entropy-tail proof or label \(0.1053610898\) as
\(\alpha_s(M_Z)\).

The full hypotheses matter: the matter-trace ray must be the physical
kinetic boundary; the beta column assumes three fermion families and one
elementary scalar doublet; electromagnetic matching is at a common scale
at tree level without thresholds; and the representation-entropy law is
adopted. Applying this branch to the native quantum-link source also needs
its state, clock, response and scale map. A physical comparison needs a
stated renormalized coupling convention and threshold/radiative matching.
The one-loop coefficients alone do not supply these attachments. This
conditional solution is useful and is not replaced by a generic claim
that every coupling is unspecified.

## What symmetry fixes and what an action must supply

The [finite electroweak theorem](../../Lean/Screen/ElectroweakBreakingComposition.lean)
and the [owning Hessian calculation](../particle_atlas_horizontal/FIELD_NATURE.md)
give, for a chosen canonically normalized broken-phase action,

\[
 m_W^2=g^2v^2/4,\qquad m_Z^2=(g^2+g'^2)v^2/4,
 \qquad m_h^2=2\lambda v^2,
 \qquad m_\gamma^2=0.
\]

Consequently

\[
 \frac{m_h^2}{m_W^2}=\frac{8\lambda}{g^2},\qquad
 \frac{m_Z^2}{m_W^2}=1+\frac{g'^2}{g^2}.
\]

Two exact positive choices illustrate the dependency. The actions
\((g,g',\lambda,v)=(1,1,1,1)\) and \((1,2,2,1)\) have the same
internal representations and canonical kinetic form, and both have the
electromagnetic null direction. Their two ratios are respectively
\((8,2)\) and \((16,5)\). They are different declared actions satisfying
the same group/current/kinetic-form premises. They are not an ambiguity
in the output of one fully specified simulator.

With the source's fundamental trace normalization fixed, a color action
\((\kappa/2)\operatorname{tr}_3\mathcal F^2\), with
\(\mathcal A=g_sA\), has \(g_s^2=1/\kappa\). Every positive
\(\kappa\) preserves the same internal symmetry. Choices one and two
give \(\alpha_s=1/(4\pi)\) and \(1/(8\pi)\) in this declared
curvature convention. The group fixes the invariant tensor shape; a
specified source action and its canonical response determine its
coefficient. A physical running coupling additionally needs the scale
and renormalization map.

A common dimensional conversion can be a convention. Rescaling
\(v\mapsto sv\) with all dimensionless coefficients held fixed
multiplies the displayed masses by \(s\); changing the clock unit
inversely describes the same dimensionless predictions. This operation
does not change either ratio above or remove dependence on
\(\lambda/g^2\), \(g'/g\) or \(\kappa\). In a fixed microscopic
program, changing a Higgs scale relative to another fixed energy is an
action change, rather than a mere unit change.

## The remaining physical response map

The exact local quantum-scalar coupling, scalar-charge interaction and
native electric/hopping ratio are computed quantities of their configured
sources. The conditional gauge ray reduces the gauge parameters once its
stated physical response and entropy hypotheses are adopted. A chosen
complete spatial action may also declare those coefficients directly.

To compute its physical outputs, use that action's state and clock to
derive canonical gauge-current and curvature responses, the scalar
effective potential or susceptibility, and the relevant phase. Carry
the same complete detector responses through controlled spatial/cutoff
limits and a specified renormalization map. These calculations decide
whether a finite gap is the radial Higgs mode and whether a finite color
coefficient becomes the coupling measured by a particular scattering
observable. No further architectural uniqueness requirement is imposed.

Parameters chosen from target masses or couplings retain calibration
ancestry. A common energy unit cannot turn fitted dimensionless ratios
into independent predictions. The exact two-choice certificate proves
insufficiency of the group/current premises alone; it does not exclude a
specific source response from selecting and measuring those parameters.

Run `python3 -m model.particle_observer_outputs.gauge_input` for the
dependency certificate and reused readout enclosures. Run
`python3 -m pytest -q model/particle_observer_outputs/test_gauge_input.py`
for the fixed-source, exact-ratio and dependency checks.
