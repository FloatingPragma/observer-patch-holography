**Paper release:** `r1549`
**Released:** July 18, 2026

# Scope

This is a technical companion outside the release bundle. It states the inflation-free branch against the compact, screen-microphysics, dark-sector, and finite-source CMB papers.

# Finite Quotient and Claim-Tier Firewall

#### Finite quotient ensemble theorem surface.

Fix a finite regulator $`r`$. The physical presentation space is $`\Sigma_r`$, the presentation redundancy groupoid is $`\Gamma_r`$, and the finite physical quotient is
``` math
Q_r=\Sigma_r/\Gamma_r,\qquad \pi_r:\Sigma_r\to Q_r.
```
The quotient removes nonphysical presentation data: gauge representatives, port relabelings, mesh labels, shard or worker identifiers, queue order, repair schedule identifiers, retry counters, timestamps unless declared semantic, hidden carrier coordinates, and inert ancillary labels. If only settled configurations carry probability, the probability space is the normal-form subset
``` math
N_r=n_r(Q_r).
```
The map $`n_r`$ is a normal-form map, not a probability law. Any promoted physical branch inherits this firewall: it must declare the quotient-intrinsic source law or action before a normal form can be read as a selection or prediction claim.

#### Observable algebras and reference states.

In the finite classical case the quotient observable algebra is
``` math
\mathcal O_r=\ell^\infty(Q_r).
```
In the finite quantum case the physical algebra is a declared quotient algebra $`\mathcal A^{\rm phys}_r`$ with state
``` math
\omega_r(A)=\operatorname{Tr}(\rho_r A).
```
When the reference object is obtained from a finite lifted carrier, the load-bearing data are not an abstract groupoid cardinality alone. They are a tracially pointed quotient
``` math
\left(\mathcal A^{\rm phys}_{r,b},\tau^0_{r,b}\right),
\qquad
\mathcal A^{\rm phys}_{r,b}
=
z_{r,b}B(\widetilde{\mathcal H}_r)^{G_r}z_{r,b},
```
where $`U_r:G_r\to U(\widetilde{\mathcal H}_r)`$ is the compact gauge action and $`z_{r,b}`$ is the central projection for the declared boundary or superselection sector. The reference trace is
``` math
\tau^0_{r,b}(A)
=
\frac{\operatorname{Tr}_{\widetilde{\mathcal H}_r}(A)}
{\operatorname{Tr}_{\widetilde{\mathcal H}_r}(z_{r,b})}.
```
If
``` math
z_{r,b}\widetilde{\mathcal H}_r
\cong
\bigoplus_\alpha V_\alpha\otimes M_\alpha,
```
with $`d_\alpha=\dim V_\alpha`$ and $`m_\alpha=\dim M_\alpha`$, then
``` math
\mathcal A^{\rm phys}_{r,b}
\cong
\bigoplus_\alpha I_{V_\alpha}\otimes B(M_\alpha),
\qquad
p_{r,\alpha}
=
\frac{d_\alpha m_\alpha}{\sum_\beta d_\beta m_\beta}.
```
These are the induced central-sector weights only after the carrier representation and boundary sector have been fixed.

#### OPH quotient ensemble.

An OPH quotient ensemble is specified by a quotient-intrinsic base weight and action
``` math
m_r:Q_r\to \mathbb R_{>0},
\qquad
S_r:Q_r\to \mathbb R\cup\{+\infty\},
```
and
``` math
w_r(q)=m_r(q)e^{-S_r(q)},\qquad
Z_r=\sum_{q\in Q_r}w_r(q),\qquad
\mu_r(q)=Z_r^{-1}w_r(q).
```
Equivalently, one may state an intrinsic projective prior $`\nu_r`$ on $`Q_r`$ and set $`\mu_r=(n_r)_\#\nu_r`$. Uniform quotient counting, uniform representative counting pushed to the quotient, groupoid weights, and tracial central-sector weights are different physical claims. The paper must declare which one is being used.

#### Normal-form projector non-selection.

For any retraction $`N:Q\to Q_{\rm nf}`$ onto a subset $`Q_{\rm nf}\subseteq Q`$, that is, any map whose restriction to $`Q_{\rm nf}`$ is the identity, the induced map on laws
``` math
\mathcal C_Q(\mu)=N_\#\mu
```
is idempotent:
``` math
\mathcal C_Q^2=\mathcal C_Q.
```
Every law supported on $`Q_{\rm nf}`$ is fixed; both statements use the retraction property. Therefore settlement or canonicalization never selects a unique physical probability law by itself.

#### Selection-gap corollary.

Let $`X\subseteq Q_{\rm nf}`$ be a finite set of quotient-normal candidates distinguished by visible invariants. Normal-form data determine $`X`$ and its quotient-visible invariants, but they do not choose a member of $`X`$. If two laws $`\mu,\nu`$ are supported on $`X`$ and concentrate on different candidates, both are fixed by $`\mathcal C_Q`$. Unique sector selection therefore requires source data: an intrinsic action with a unique minimizer, a declared physical ensemble, or a refinement-stable gap certificate. A defect or holonomy classification can classify possible sectors without choosing the physical sector, and a contraction or repair generator can certify convergence toward a declared target without creating the target law.

#### Finite MaxEnt quotient ensemble.

For finite $`Q`$, positive $`m`$, and quotient observables $`F_1,\ldots,F_k`$, maximizing
``` math
\mathcal H_m(\nu)=-\sum_q\nu(q)\log\frac{\nu(q)}{m(q)}
```
subject to
``` math
\sum_q\nu(q)=1,\qquad \sum_q\nu(q)F_a(q)=c_a
```
has the full-support solution, when the feasible full-support surface is nonempty,
``` math
\mu(q)=
\frac{m(q)\exp[-\sum_a\theta_aF_a(q)]}{Z(\theta)}.
```
Boundary optima obey the same formula after restricting to their support. On a finite noncommutative quotient algebra with faithful reference state $`\sigma_r`$,
``` math
\rho_r=
\frac{\exp(\log\sigma_r-\sum_a\theta_aF_{r,a})}
{\operatorname{Tr}\exp(\log\sigma_r-\sum_a\theta_aF_{r,a})}.
```
The finite constraint ledger must name every $`F_{r,a}`$, its units and support, the target expectation and source, sector or zero-mode treatment, refinement transformation, and proof that no run output or observational output entered the source definition.

#### Refinement compatibility and RG closure.

For $`s\succeq r`$, let $`c_{sr}:Q_s\to Q_r`$ be the physical coarse map. Exact compatibility of weighted ensembles is equivalent to the fiber-sum identity
``` math
\sum_{q':\,c_{sr}(q')=q}m_s(q')e^{-S_s(q')}
=
\alpha_{sr}m_r(q)e^{-S_r(q)}
```
for a constant $`\alpha_{sr}>0`$ independent of $`q`$. Then
``` math
(c_{sr})_\#\mu_s=\mu_r.
```
If the one-step defects are
``` math
\delta_{k+1,k}
=
\left\|(c_{k+1,k})_\#\mu_{k+1}-\mu_k\right\|_{\mathrm{TV}},
```
then
``` math
\left\|(c_{nr})_\#\mu_n-\mu_r\right\|_{\mathrm{TV}}
\le
\sum_{k=r}^{n-1}\delta_{k+1,k}.
```
For exponential-family refinement, exact closure requires the fine conditional free energy
``` math
G_{sr,\theta}(q_r)
=
-\log\mathbb E_{m_s^0}\left[
\exp[-\theta\cdot F_s(Q_s)]\mid c_{sr}(Q_s)=q_r
\right]
```
to equal $`\kappa_{sr}(\theta)+R_{sr}(\theta)\cdot F_r(q_r)`$. If the residual is uniformly bounded by $`\varepsilon`$, the induced total-variation defect is bounded by $`\tanh\varepsilon`$.

#### Implementation invariance and representative lifting.

If implementations $`A,B`$ have quotient bijections $`h_r:Q_r^A\to Q_r^B`$ satisfying
``` math
m_r^B(h_rq)=m_r^A(q),\qquad
S_r^B(h_rq)=S_r^A(q),
\qquad
h_r\circ c_{sr}^A=c_{sr}^B\circ h_s,
```
then
``` math
(h_r)_\#\mu_r^A=\mu_r^B.
```
For tracially pointed quantum quotients the corresponding equivalence is a trace-preserving quotient equivalence. It is invariant under unitary intertwiners preserving the gauge action and sector, and under inert trivial ancillas $`A\mapsto A\otimes I_{\rm anc}`$. It is not invariant under arbitrary changes of gauge-representation multiplicities.

If an implementation stores representatives, a representative-level law must be a conditional lift
``` math
\widetilde\mu_r(x)=\mu_r(\pi_r x)\kappa_r(x\mid \pi_r x),
\qquad
\sum_{x:\pi_r(x)=q}\kappa_r(x\mid q)=1.
```
Then $`(\pi_r)_\#\widetilde\mu_r=\mu_r`$. Uniform representative sampling yields orbit-size weights and is physical only if representative counting is the declared base measure.

#### Quotient-lumpable kernels and sampler correctness.

A representative kernel $`\widetilde P(x,y)`$ descends to $`Q_r`$ only when
``` math
P_Q(q,q')
=
\sum_{y:\pi(y)=q'}\widetilde P(x,y)
```
is independent of the chosen representative $`x\in\pi^{-1}(q)`$. For $`w(q)=m(q)e^{-S(q)}`$ and proposal $`R(q,q')`$ with reciprocal support, the Metropolis–Hastings acceptance rule
``` math
a(q,q')=\min\left\{1,
\frac{w(q')R(q',q)}{w(q)R(q,q')}
\right\}
```
gives detailed balance
``` math
\mu(q)R(q,q')a(q,q')=\mu(q')R(q',q)a(q',q).
```
Repair-informed proposals must include the Hastings asymmetry term; otherwise the stationary law is generically changed.

#### Repair generators are not selectors.

A repair generator of the form
``` math
L_{\rm rep}=\sum_C c_C(I-E_C)
```
is a relaxation or sampling object after a law has been selected. Conditional expectations $`E_C`$ are defined on $`L^2(X_r,\pi_r)`$, so the reference law $`\pi_r`$ is input. On overlapping collars the expectations need not commute. The correct finite gap certificate is the Poincare constant
``` math
\kappa_r
=
\inf_{f\perp 1}
\frac{\sum_C\|(I-E_C)f\|^2}{\|f\|^2}.
```
If local fiber rates have a positive lower bound $`\gamma_*`$, then
``` math
L_{\rm rep}\ge \gamma_*\kappa_r(I-P_0).
```
Finite repair completeness gives $`\kappa_r>0`$ at fixed regulator. A uniform refinement lower bound $`\inf_r\kappa_r>0`$ is a separate theorem or receipt.

#### Finite evidence accuracy.

For bounded coarse observables $`O`$, if
``` math
\|\widehat\mu_s-\mu_s\|_{\mathrm{TV}}\le\epsilon_{\rm samp}
```
and the refinement defects sum to $`\epsilon_{\rm ref}`$, then
``` math
\left|\mathbb E_{\widehat\mu_s}[O\circ c_{sr}]-\mathbb E_{\mu_r}[O]\right|
\le
2\|O\|_\infty(\epsilon_{\rm samp}+\epsilon_{\rm ref}).
```
Continuum-facing observables require a realization map and correlation Cauchy bound in addition to a finite histogram.

#### Vacuum promotion gate.

A stationary sampler is not a physical vacuum. For any faithful target law one can build a positive transfer operator with that law as ground state, so positivity alone is not a selector. Vacuum promotion requires source Euclidean slab data
``` math
\mathfrak S_r^E=(Q_r,m_r^0,J_r,V_r,a_{t,r})
```
whose conductance $`J_r(q,q')=J_r(q',q)\ge0`$, local potential $`V_r`$, and slab thickness $`a_{t,r}`$ are derived without using the target law or sampler output. With connected event graph,
``` math
(H_r^Ef)(q)
=
\frac{1}{m_r^0(q)}
\sum_{q'}J_r(q,q')\bigl(f(q)-f(q')\bigr)
+
V_r(q)f(q)
```
is self-adjoint and bounded below on $`L^2(Q_r,m_r^0)`$; its finite Feynman–Kac semigroup is positivity improving. Perron–Frobenius gives a unique positive normalized ground state $`\Omega_r`$, and the finite vacuum law is
``` math
\mu_r^{\rm vac}(q)=|\Omega_r(q)|^2m_r^0(q).
```
For $`T_r=e^{-a_{t,r}(H_r^E-E_{0,r})}`$, the Doob kernel is stochastic and detailed-balanced with $`\mu_r^{\rm vac}`$. Continuum promotion additionally requires reflection positivity or equivalent reconstruction plus refinement compatibility of the transfer family.

#### Primordial and cosmological prediction firewall.

A screen covariance is not automatically a primordial curvature spectrum. The map
``` math
C_\ell
=
4\pi
\int_0^\infty
\frac{dk}{k}
\Delta_\zeta^2(k)j_\ell^2(k\chi_\star)
```
has a null space unless a radial prior, finite parametric family, or source-stress bridge is declared. OPH primordial promotion requires the source-only stress, single-clock, entropy-repair, curvature-evolution, adiabatic-mode, phase-coherence, screen-to-radial-lift, radial-null-space, and forward-projection receipts. Observable CMB comparison also requires declared source, solver, dataset, covariance, nuisance, data-use, and pooled-reducer provenance.

#### Claim tiers and required receipts.

Every ensemble-facing run records its ensemble id, claim tier, regulator, representative schema, gauge action, canonicalizer, base measure, action coefficients, coarse maps, zero-mode projector, amplitude convention, sampler, smoothing policy, source provenance, and explicit nonclaims. The seed belongs to the run receipt rather than the ensemble definition. The claim tiers are
``` math
\begin{array}{ll}
E0:&\text{seed noise, proposal noise, repair jitter},\\
E1:&\text{conventional reference ensemble},\\
E2:&\text{OPH-native quotient ensemble},\\
E3:&\text{OPH vacuum},\\
E4:&\text{OPH primordial field},\\
E5:&\text{observable cosmological prediction}.
\end{array}
```
The evidence bundle must keep separate receipts for stationary-law schedule invariance, detailed balance of the aggregate kernel, and pathwise partition invariance. Deterministic replay of semantic random streams or a canonical serial chain is useful, but it is not pathwise partition invariance. Smoothing must preserve raw coefficients, raw spectra, smoothing kernels, smoothed coefficients, smoothed spectra, and hashes of each stage; it is not part of $`S_r`$ unless explicitly declared.

# Finite Screen Spectrum Theorem Package

This fragment owns the screen-level spectrum theorem used by the staged cosmology branch. It does not own TT, TE, EE, lensing, likelihoods, or physical CMB promotion. Those remain downstream Boltzmann and data-contract gates.

<div id="def:oph-screen-regulator" class="definition">

**Definition 1** (finite screen regulator). *For regulator $`r`$, a finite screen regulator is
``` math
\mathfrak S_r=(\mathcal T_r,M_r,\Gamma_r,\mathcal Q_r,C_{sr},J_{rs}),
```
where $`\mathcal T_r`$ is a finite cellulation of the screen, $`M_r\succ0`$ is the area/quadrature mass matrix, $`\Gamma_r`$ is the hidden-presentation groupoid, $`\mathcal Q_r=\Sigma_r/\Gamma_r`$ is the physical quotient, and $`C_{sr}`$, $`J_{rs}`$ are the coarse-graining and interpolation maps. Shape regularity and quadrature convergence require
``` math
f^{\mathsf T}M_rg \to \int_{S^2}fg\,d\Omega
```
with an explicit band-limited error bound
``` math
\left|f^{\mathsf T}M_rg-\int f g\,d\Omega\right|
  \le \varepsilon_{M,r}(L)\|f\|_{H^s}\|g\|_{H^s}.
```
Patch count alone is therefore not an angular-resolution certificate.*

</div>

<div id="def:oph-screen-scalar" class="definition">

**Definition 2** (geometric screen scalar). *Let the quotient-visible collar-volume readout be
``` math
J_{X,r}(x)=\lambda_r(x)\sqrt{\det\sigma_{AB,r}(x)}>0,
```
and let $`\bar J_{X,r}`$ be emitted by an independently defined homogeneous or frozen-background branch, not by a CMB fit. Define
``` math
q_{0,r}=\frac13\log\frac{J_{X,r}}{\bar J_{X,r}}.
```
On a certified spherical chart, with $`B_r=[1,n_x,n_y,n_z]`$,
``` math
\Pi_{<2,r}=B_r(B_r^{\mathsf T}M_rB_r)^{-1}B_r^{\mathsf T}M_r,\qquad
  q_r=(I-\Pi_{<2,r})q_{0,r}.
```
On an irregular screen, $`B_r`$ is replaced by the certified generalized eigenprojector onto the finite $`\ell=0,1`$ subspace. Hard-coded feature z-scores or observer-summary weights are diagnostic proxies, not $`q_r`$.*

</div>

<div id="prop:oph-screen-scalar-invariance" class="proposition">

**Proposition 3** (quotient invariance). *If a recharting $`U`$ preserves the screen inner product and geometric readout,
``` math
U^{\mathsf T}M'_rU=M_r,\qquad q'_{0,r}=Uq_{0,r},\qquad B'_r=UB_r,
```
then $`\Pi'_{\ge2,r}U=U\Pi_{\ge2,r}`$, and $`q'_r=Uq_r`$.*

</div>

<div class="proof">

*Proof.* Substitution gives
``` math
\begin{aligned}
\Pi'_{<2,r}U
&=UB_r(B_r^{\mathsf T}U^{\mathsf T}M'_rUB_r)^{-1}B_r^{\mathsf T}U^{\mathsf T}M'_rU\\
&=UB_r(B_r^{\mathsf T}M_rB_r)^{-1}B_r^{\mathsf T}M_r=U\Pi_{<2,r}.
\end{aligned}
```
Subtracting from $`U`$ gives the high-mode identity. ◻

</div>

<div id="def:oph-screen-precision" class="definition">

**Definition 4** (physical scalar precision and action). *The scalar precision $`K_r`$ must have a physical origin on $`V_r=\operatorname{im}\Pi_{\ge2,r}`$. Allowed branches are:
``` math
K_r^{\rm Hess}=\left.\nabla^2\Phi_r\right|_{q=0},
```
for a quotient-visible mismatch or release free energy; a repair Dirichlet form
``` math
K_r^{\rm diss}=H_r-T_r^{\mathsf T}H_rT_r,\qquad T_r^{\mathsf T}H_rT_r\preceq H_r;
```
or the reversible generator form
``` math
K_r^{\rm Dir}=-\frac12(L_r+L_r^{\dagger_H}).
```
A raw nonsymmetric repair matrix or caller-supplied $`\kappa`$ is not a precision operator. The absolute normalization must be fixed independently, for example by
``` math
K_0Y_{\ell m}=\ell(\ell+1)Y_{\ell m}
```
on the $`\theta=0`$ branch. Once normalized,
``` math
S_{{\rm scr},r}[q]=\frac{1}{2A_{q,r}}\langle q,K_rq\rangle_r
  =\frac{1}{2A_{q,r}}q^{\mathsf T}M_rK_rq .
```*

</div>

<div id="thm:oph-screen-maxent" class="theorem">

**Theorem 5** (finite MaxEnt screen covariance). *Let $`K_re_i=\kappa_i e_i`$ with an $`M_r`$-orthonormal basis of $`V_r`$, $`\kappa_i>0`$, and $`d_r=\dim V_r`$. Among continuous densities in coefficients $`q=\sum_i a_ie_i`$ with $`\mathbb E[a_i]=0`$ and
``` math
\frac12\mathbb E\langle q,K_rq\rangle_r=E^{\rm src}_{q,r},
```
the unique entropy maximizer is
``` math
d\mu_r(q)=Z_r^{-1}
  \exp\!\left[-\frac{1}{2A_{q,r}}\langle q,K_rq\rangle_r\right]dq,
  \qquad
  A_{q,r}=\frac{2E^{\rm src}_{q,r}}{d_r},
```
and
``` math
\mathbb E[a_ia_j]=\delta_{ij}\frac{A_{q,r}}{\kappa_i}.
```*

</div>

<div class="proof">

*Proof.* The Euler–Lagrange equation for entropy with normalization and quadratic-energy constraints gives a Gaussian density with inverse temperature $`\beta_r`$. The expected energy is
``` math
E^{\rm src}_{q,r}=\frac12\sum_i\kappa_i\frac1{\beta_r\kappa_i}=\frac{d_r}{2\beta_r}.
```
Thus $`A_{q,r}=\beta_r^{-1}=2E^{\rm src}_{q,r}/d_r`$. Strict concavity of entropy gives uniqueness. ◻

</div>

<div class="remark">

*Remark 6* (discrete quotient branch). For a discrete finite quotient with intrinsic base weight $`m_r(q)`$, the exact Gibbs law is
``` math
\mu_r(q)=Z_r^{-1}m_r(q)e^{-\beta_r E_r(q)}.
```
A Gaussian statement on that branch requires a separate Laplace or central-limit theorem with controlled higher-order terms.

</div>

<div id="def:oph-screen-release-energy" class="definition">

**Definition 7** (source release energy). *Let $`m_r^0`$ be the tracially pointed base weight on the finite physical collar quotient. Let $`F_{r,a}`$ be a finite ledger of primitive collar observables: scalar occupancy, protected-center presence, released volume, conserved source charges, and the release clock. The ledger excludes screen energy, sky spectra, CMB residuals, likelihoods, fitted amplitudes, and measurement-calibrated proxies. For source-emitted constraint values $`c_{r,a}^{\rm col}`$, the selected release law is
``` math
\nu^{\rm rel}_r(x)=\frac{m_r^0(x)
  \exp[-\sum_a\lambda_{r,a}F_{r,a}(x)]}{Z_r},
  \qquad
  \mathbb E_{\nu_r^{\rm rel}}F_{r,a}=c_{r,a}^{\rm col}.
```
Strict concavity of relative entropy gives uniqueness on full support. A boundary optimum uses the same form on its certified support. The constraint ledger, base weight, source DAG, multipliers, support, and refinement maps belong to the receipt.*

*For this independently selected quotient ensemble on released collar normal forms, define
``` math
E^{\rm src}_{q,r}
  =\frac12\int \langle q_r(x),K_rq_r(x)\rangle_r\,d\nu^{\rm rel}_r(x).
```
The receipt must expose $`\nu^{\rm rel}_r`$, the base measure, $`K_r`$, $`d_r`$, $`E^{\rm src}_{q,r}`$, $`A_{q,r}`$, no-observation ancestry, and the same operator normalization used in the action. MaxEnt alone does not determine amplitude.*

</div>

<div id="prop:oph-screen-amplitude-nogo" class="proposition">

**Proposition 8** (microscopic-law necessity). *For fixed positive $`K_r`$, the Gaussian family $`Z(A)^{-1}\exp[-\langle q,K_rq\rangle/(2A)]`$ exists for every $`A>0`$. The quadratic MaxEnt form therefore leaves one positive scale free. The source release law in Definition <a href="#def:oph-screen-release-energy" data-reference-type="ref" data-reference="def:oph-screen-release-energy">7</a>, or an equivalent source-only vacuum law, supplies the energy that fixes $`A_{q,r}=2E^{\rm src}_{q,r}/d_r`$.*

</div>

<div id="thm:oph-screen-angular-spectrum" class="theorem">

**Theorem 9** (rotational screen spectrum). *If the continuum precision $`K_\theta`$ commutes with the $`SO(3)`$ action on scalar functions, then by Schur’s lemma each $`\mathcal H_\ell`$ is an eigenspace. For the exact conformal-shell family,
``` math
\Lambda_\ell(\theta)
  =\frac{\Gamma(\ell+2+\theta/2)}{\Gamma(\ell-\theta/2)},\qquad \ell\ge2.
```
Then $`\Lambda_\ell(0)=\ell(\ell+1)`$, and the MaxEnt covariance gives
``` math
C_\ell^q
  =A_q\,\frac{\Gamma(\ell-\theta/2)}{\Gamma(\ell+2+\theta/2)}.
```
Consequently $`\mathcal D_\ell^q=\ell(\ell+1)C_\ell^q/(2\pi)\sim(A_q/2\pi)\ell^{-\theta}`$.*

</div>

<div id="def:oph-screen-conformal-precision" class="definition">

**Definition 10** (source-derived conformal precision). *Let $`L_r\succeq0`$ be the normalized detailed-balanced scalar collar Dirichlet operator on $`V_r`$, with continuum limit $`-\Delta_{S^2}`$. Set
``` math
B_r=(L_r+\tfrac14I)^{1/2},\qquad
  K_{\theta,r}=
  \frac{\Gamma(B_r+\frac32+\frac\theta2)}
       {\Gamma(B_r-\frac12-\frac\theta2)}.
```
For $`-2<\theta<4`$ on the retained $`\ell\ge2`$ branch, this operator is positive. The gamma recurrence gives $`K_{0,r}=L_r`$ exactly. The finite receipt carries detailed balance, positivity, the $`\ell(\ell+1)`$ normalization, anisotropy splitting, and refinement residuals.*

</div>

<div id="thm:oph-screen-refinement" class="theorem">

**Theorem 11** (finite refinement error). *Suppose that for $`2\le\ell\le L`$
``` math
\left|\frac{\lambda_{\ell m,r}}{\Lambda_\ell(\theta)}-1\right|\le\varepsilon_{K,r}(L),
  \qquad
  \left|\frac{A_{q,r}}{A_q}-1\right|\le\varepsilon_{A,r},
  \qquad \varepsilon_{K,r}<1 .
```
Then $`C_{\ell m,r}^q=A_{q,r}/\lambda_{\ell m,r}`$ satisfies
``` math
\left|
  \frac{C_{\ell m,r}^q}{A_q/\Lambda_\ell(\theta)}-1
  \right|
  \le
  \frac{\varepsilon_{A,r}+\varepsilon_{K,r}}{1-\varepsilon_{K,r}} .
```*

</div>

<div class="proof">

*Proof.* Write $`A_{q,r}=A_q(1+a_r)`$ and $`\lambda_{\ell m,r}=\Lambda_\ell(1+b_{\ell m,r})`$. Then
``` math
\frac{C_{\ell m,r}^q}{A_q/\Lambda_\ell}=\frac{1+a_r}{1+b_{\ell m,r}},
```
and the displayed bound follows from $`|a_r|\le\varepsilon_A`$, $`|b_{\ell m,r}|\le\varepsilon_K`$. ◻

</div>

<div id="thm:oph-screen-tilt" class="theorem">

**Theorem 12** (refinement-semigroup tilt). *Let $`R_b`$ be the scalar refinement map for scale ratio $`b>1`$. If one isolated scalar covariance mode has positive eigenvalue $`\lambda(b)`$, $`\lambda(b_1b_2)=\lambda(b_1)\lambda(b_2)`$, and $`\lambda`$ is continuous, then there is a unique real $`\theta`$ with
``` math
\lambda(b)=b^{-\theta},\qquad
  \theta=-\frac{\log\lambda(b)}{\log b}.
```*

</div>

<div class="proof">

*Proof.* Set $`g(t)=-\log\lambda(e^t)`$. The semigroup law gives $`g(t+s)=g(t)+g(s)`$. Continuity makes $`g(t)=\theta t`$, hence the result. ◻

</div>

<div id="def:oph-screen-reserve-generator" class="definition">

**Definition 13** (edge-center reserve generator receipt). *Write $`s=\log b`$ for logarithmic refinement thickness. Let $`u_{\rm full}(s)`$ be the scalar-conditioned covariance-survival cocycle across a full oriented collar, and let $`u_q(s)`$ be its source-facing half-collar restriction. The receipt requires
``` math
u(s+t)=u(s)u(t),\qquad u(0)=1,
```
strong continuity at zero, the source-derived full-collar density
``` math
-u'_{\rm full}(0)=\frac{P_\star}{24},
```
and the orientation-reversal coarea identity
``` math
-u'_q(0)=\frac12[-u'_{\rm full}(0)]=\frac{P_\star}{48}.
```
These are infinitesimal generator statements. A one-step survival probability has exponent $`-\log u_q(\log b)/\log b`$.*

</div>

<div id="thm:oph-p-over-48-tilt" class="theorem">

**Theorem 14** (edge-center tilt and repair-clock reconciliation). *Under Definition <a href="#def:oph-screen-reserve-generator" data-reference-type="ref" data-reference="def:oph-screen-reserve-generator">13</a>,
``` math
u_q(s)=e^{-\theta s},\qquad
\theta=\frac{P_\star}{48},\qquad
n_s=1-\frac{P_\star}{48}.
```
In the coordinate $`\theta=\kappa_{\rm rep}(P_\star-\varphi)`$, the same branch has
``` math
\kappa_{\rm rep}^{\rm edge}
=\frac{P_\star}{48(P_\star-\varphi)}.
```
The value $`e`$ is a separate diagnostic hypothesis unless an additional identity equates it with this source-derived coordinate.*

</div>

<div class="proof">

*Proof.* For $`g(s)=-\log u_q(s)`$, the cocycle law gives the continuous Cauchy equation. Hence $`g(s)=\theta s`$. Differentiation at zero and the half-collar identity give $`\theta=P_\star/48`$. The formula for $`\kappa_{\rm rep}^{\rm edge}`$ is algebraic. ◻

</div>

<div id="thm:oph-screen-radial-source-family" class="theorem">

**Theorem 15** (homogeneous radial source family). *Assume source-stress closure, a single clock, freezeout, multicenter consistency, and translation/rotation covariance. Let $`C_\zeta=M_{\Delta_\zeta^2}`$ on the common physical $`d\ln k`$ mode basis. If a scale-natural source embedding transports refinement to physical dilation and its finite operator residual converges to
``` math
D_s^{-1}C_\zeta D_s=e^{-\theta s}C_\zeta,
\qquad (D_sf)(k,\hat k)=f(e^{-s}k,\hat k),
```
then every positive measurable representative satisfies
``` math
\Delta_\zeta^2(k)=A_\zeta(k/k_\star)^{-\theta}
```
almost everywhere. Continuity gives the pointwise identity. A single shell has an infinite-dimensional radial kernel, including positive ambiguities; complete radial cross-covariances give the independent spherical-Hankel tomography route.*

</div>

<div class="proof">

*Proof.* Conjugation of the multiplication operator gives $`\Delta_\zeta^2(e^sk)=e^{-\theta s}\Delta_\zeta^2(k)`$ almost everywhere. In logarithmic coordinates, rational-translation invariance on a common full-measure set forces $`e^{\theta t}\Delta_\zeta^2(e^t)`$ to be constant almost everywhere. The one-shell and tomography statements follow from the correlation-restriction and spherical-Hankel theorems in the primordial bridge packet. ◻

</div>

<div id="thm:oph-thin-shell-powerlaw-lift" class="theorem">

**Theorem 16** (thin-shell power-law lift). *Assume $`q(\hat n)=Z_\star\Pi_{\ell\ge2}\zeta_\star(R_\star\hat n)`$ and
``` math
\Delta_\zeta^2(k)=A_\zeta(k/k_\star)^{-\theta}.
```
Then
``` math
C_\ell^q=4\pi Z_\star^2\int d\ln k\,
  \Delta_\zeta^2(k)j_\ell^2(kR_\star),
```
and, for $`-2<\theta<4`$,
``` math
A_q
  =\pi^{3/2}Z_\star^2A_\zeta(k_\star R_\star)^\theta
  \frac{\Gamma(1+\theta/2)}{\Gamma(3/2+\theta/2)} .
```
Thus
``` math
A_\zeta
  =\frac{A_q}{\pi^{3/2}Z_\star^2(k_\star R_\star)^\theta}
  \frac{\Gamma(3/2+\theta/2)}{\Gamma(1+\theta/2)} .
```*

</div>

<div id="prop:oph-finite-window-bound" class="proposition">

**Proposition 17** (finite-window bound). *Let $`\Psi_\ell(k)=\int dr\,W(r)j_\ell(kr)`$, $`W\ge0`$, $`\int Wdr=1`$, with mean $`R_\star`$ and variance $`\sigma_R^2`$. If
``` math
\delta_\ell(k)=\frac{k^2\sigma_R^2}{2}
  \sup_{r\in{\rm supp}\,W}|j_\ell''(kr)|,
```
then
``` math
|C_{\ell,W}^q-C_{\ell,{\rm shell}}^q|
  \le
  4\pi Z_\star^2\int d\ln k\,\Delta_\zeta^2(k)
  \left[2\delta_\ell(k)+\delta_\ell(k)^2\right].
```*

</div>

<div id="prop:oph-radial-null" class="proposition">

**Proposition 18** (radial non-identifiability). *For a finite radial basis, $`C=Ap`$. If $`A\in\mathbb R^{N_\ell\times N_k}`$ has rank $`r`$, then $`\dim\ker A=N_k-r`$, and every $`p+v`$ with $`v\in\ker A`$ gives the same screen spectrum. Radial promotion therefore requires either a source theorem restricting $`p`$, or a declared prior with a published null basis, resolution kernels, positivity checks, and prior-sensitivity report.*

</div>

<div id="thm:oph-screen-spectrum-promotion" class="theorem">

**Theorem 19** (conditional source-spectrum promotion). *One hash-locked source DAG may emit
``` math
(q_r,S_{{\rm scr},r},\theta,A_{q,r},C_{\ell,r}^q,A_{\zeta,r},
\Delta_{\zeta,r}^2)
```
as a conditional primordial packet only when the geometric-scalar, collar-law, release-energy, reserve-generator, conformal-precision, refinement, source-stress, single-clock, freezeout, adiabaticity, isocurvature, phase-coherence, thin-shell, radial-null, finite-window, and forward- residual receipts pass on that DAG. Any measurement, likelihood, fitted parameter, or measurement-calibrated ancestor fails the source packet. Physical temperature and polarization spectra require the independent Boltzmann, recombination, nuisance, covariance, and likelihood gates.*

</div>

<div id="target:oph-screen-spectrum-receipts" class="target">

**Target 20** (screen-spectrum receipt set). *A concrete source-derived spectrum requires the following receipt families before primordial promotion:*

1.  *geometry, scalar quotient, low-mode projector, scalar precision, and operator normalization;*

2.  *quotient ensemble selection, scalar release energy, and MaxEnt screen covariance;*

3.  *the infinitesimal reserve generator, scalar refinement tilt, and angular screen spectrum with a finite error budget;*

4.  *source-stress, clock, freezeout, radial-window, radial-null, forward-residual, and transfer firewall receipts.*

*The theorem and algorithm packet is complete. Construction of one finite source DAG that passes this receipt set is work in progress.*

</div>

# Branch Thesis

The branch claim is:
``` math
\boxed{
\begin{gathered}
\text{On a finite-source branch, OPH may replace the logical jobs of}\\
\text{inflation by direct/CMH spatial-holonomy selection, same-boundary or low-}k\text{ repair coherence,}\\
\text{a geometric screen spectrum, MaxEnt release, and ordinary Boltzmann transfer.}
\end{gathered}}
```
This claim is not part of the recovered SM/GR core. It becomes predictive only when the source objects are fixed before observational comparison. Screen synchronization instantiates C3 of the project consistency stack: modular self-consistency, with internal time generated by the state–algebra pair rather than by an external clock.

Coherent-matter same-channel forcing supplies a local scalar-channel identity for nondegenerate OPH-coherent material sources under Scalar Edge-Center Exhaustion. It does not supply the finite covariant dark-stress parent needed for cosmological $`w_A`$, $`c_{s,A}^2`$, $`\sigma_A`$, $`Q_A^\mu`$, recipient-stress, gauge-variable, CDM-limit, and refinement-convergence receipts.

The anomaly abundance selector is owned by `cosmology/oph_dark_matter_paper.tex`. The finite-source branch may use a source-only $`\bar\rho_A(a)`$ only when
``` math
\texttt{ANOMALY\_ABUNDANCE\_SOURCE\_RECEIPT}
```
passes. Otherwise the dark continuation is either $`\texttt{OFF}`$ or $`\texttt{CONDITIONAL\_SOURCE\_STATE}`$.

# Paper Targets

| Target | Required theorem or receipt |
|:---|:---|
| Finite early state | Quotient normal-form branch with declared base measure and implementation invariance. |
| Flatness | Zero clock-slice spatial Levi–Civita holonomy identifies $`\kappa=0`$. Exact selection requires <span class="smallcaps">DirectTheorem</span>, <span class="smallcaps">ConditionalCMH</span>, or <span class="smallcaps">ExplicitAssumption</span>; damping is only an approximate $`\vert{}K\vert{}`$ or $`\vert{}\Omega_K\vert{}`$ bound. |
| Horizon coherence | Same-boundary unique-extension theorem or a uniform low-$`k`$ repair gap. |
| Near scale invariance | Geometric screen scalar, normalized scalar precision, source release energy, infinitesimal reserve-generator tilt, and angular screen-spectrum receipts. |
| Amplitude | The conditional theorem fixes $`A_q=2E_q^{\rm src}/d`$ and its source-family conversion to $`A_\zeta`$. A finite collar release law and energy receipt must supply their values without CMB data. |
| Hot start | MaxEnt release into the realized radiation branch with a common release clock. |
| Acoustic transfer | Declared source functions and Boltzmann transfer owned by the finite-source CMB and data-contract papers. |

# Simulator Contract

The simulator reference is <https://github.com/muellerberndt/oph-physics-sim>. For this paper the simulator may emit diagnostic screen spectra, observer-screen synchronization receipts, and release-amplitude candidates. It must not label them inflation replacement, OPH primordial field, or physical CMB prediction unless the E4 and E5 receipts in the shared theorem fragment pass in the concrete run bundle.

# Flatness Boundary

The flatness branch is split into identification and selection. On a clocked FLRW continuation with $`u_a=-N\nabla_a\chi`$, $`h_{ab}=g_{ab}+u_au_b`$, and $`K(\tau)=\kappa/a(\tau)^2`$, the area-normalized small-loop holonomy of the spatial Levi–Civita connection identifies $`K=0`$ and hence $`\kappa=0`$. It does not select that sector from the raw FLRW fiber. The simulator must therefore report one of $`\textsc{DirectTheorem}`$, $`\textsc{ConditionalCMH}`$, $`\textsc{ExplicitAssumption}`$, or $`\textsc{OpenTheorem}`$. S3 screen/collar defect decay is a repair diagnostic, not a spatial curvature holonomy. Curvature damping can certify only an $`\textsc{ApproximateBound}`$ on $`|K|`$ or $`|\Omega_K|`$.

# Screen-Spectrum Branch

The load-bearing screen result is the angular screen covariance, not a CMB spectrum:
``` math
C_\ell^q
  =
  A_q\frac{\Gamma(\ell-\theta/2)}{\Gamma(\ell+2+\theta/2)}.
```
This statement requires the geometric $`q_r`$ of Definition <a href="#def:oph-screen-scalar" data-reference-type="ref" data-reference="def:oph-screen-scalar">2</a>, a normalized conformal precision operator $`K_{\theta,r}`$, the primitive collar release law, $`E^{\rm src}_{q,r}`$, the MaxEnt covariance theorem, and a finite-to-continuum refinement error budget. The edge-center theorem gives
``` math
\theta=\frac{P_\star}{48},\qquad
n_s=1-\frac{P_\star}{48},\qquad
\kappa_{\rm rep}^{\rm edge}=\frac{P_\star}{48(P_\star-\varphi)}.
```
Its premise is a source-emitted full-collar generator density $`P_\star/24`$ with the orientation-half identity. A finite one-step survival probability is converted through $`-\log u/\log b`$. The $`e(P_\star-\varphi)`$ clock branch is a diagnostic alternative. The finite generator and clock receipt is work in progress.

Under the declared source dilation cocycle, the primordial power law
``` math
\Delta_\zeta^2(k)=A_\zeta\left(\frac{k}{k_\star}\right)^{-\theta}
```
is the unique member of the one-dimensional source family. The exact thin-shell identity fixes the amplitude conversion
``` math
A_\zeta=\frac{A_q}{\pi^{3/2}}
\frac{\Gamma(3/2+\theta/2)}{\Gamma(1+\theta/2)}
```
at the source pivot. A finite run requires $`R_\star`$, the radial window, the unrestricted null-space report, the forward residual, and the physical scale. The screen theorem proves no TT, TE, EE, lensing, acoustic peak, or likelihood statement.

# Unclosed Claim Boundaries

1.  Prove the same-boundary or low-$`k`$ coherence branch in quotient language.

2.  Emit the geometric collar-volume field $`q_r=\Pi_{\ge2,r}(1/3)\log(J_{X,r}/\bar J_{X,r})`$ from one finite source.

3.  Emit the primitive collar release law and pooled release energy that fix $`A_q`$.

4.  Emit the full-collar infinitesimal reserve generator and orientation-half receipt that instantiate the conditional $`P_\star/48`$ theorem.

5.  Pass the source-stress, clock, freezeout, physical-scale, radial-window, null-space, and forward-residual receipts on the same DAG.

6.  Hand the resulting source inputs to the CMB and likelihood-contract papers.
