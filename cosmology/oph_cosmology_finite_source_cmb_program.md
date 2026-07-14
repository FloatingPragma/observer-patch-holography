# What This Paper Contributes

The OPH core papers recover local observer-facing geometry, records, modular flow, gauge/matter structure, and finite screen microphysics. The cosmology work asks a different question: which large-scale data can be predicted once the finite observer-screen state is released into an Einstein-Boltzmann cosmology?

The program paper has four jobs.

1.  It states the inflation-free OPH cosmology branch as a technical companion outside the release bundle.

2.  It separates diagnostic simulator outputs from physical predictions. A good visual or measurement-facing CMB overlay is not a physical CMB prediction unless its inputs are source-only finite OPH artifacts fixed before likelihood comparison.

3.  It names the finite objects required for physical promotion of simulation outputs in CMB, vacuum, or cosmological uses of the dark sector: source provenance, physical mode calibration, freezeout, covariant stress, anomaly kernels, transfer, likelihood, and quotient ensemble receipts. The dark-sector derivations themselves are owned by the released cosmology paper `cosmology/oph_dark_matter_paper.tex`.

4.  It gives the conditional cosmology branch one theorem surface for the finite-source and likelihood claim gates.

# Claim Boundary

This technical companion sits outside the release bundle and sets the shared cosmology contract for paper and simulator work.

The branch claim is:
``` math
\boxed{
\begin{gathered}
\text{OPH can replace the logical jobs of inflation on a finite-source branch:}\\
\text{flatness by direct/CMH spatial-holonomy selection, horizon coherence by same-boundary repair,}\\
\text{near scale invariance by the geometric screen spectrum, hot initial data by MaxEnt release,}\\
\text{and acoustic peaks by Boltzmann transfer from source-only OPH artifacts.}
\end{gathered}}
```

Here “can” means “can on the declared finite-source branch with the listed source, stress, lift, transfer, and likelihood gates closed.” The recovered Phase I SM/GR core excludes inflation replacement, CMB spectra, and dark-sector cosmology.

The branch supplies claim boundaries and measurement-facing diagnostic evidence. Physical CMB prediction, completed Planck/ACT likelihood runs, OPH-native quantum vacuum, full cosmological dark-matter likelihood, and production particle-formation simulation require the gates below.

<div class="remark">

*Remark 1* (No simulator shortcut). The paper-side theorem gates mark the difference between a simulation that resembles CMB data and a prediction whose input artifacts were fixed without using CMB data. Physical CMB promotion requires those gates as mathematical statements and deterministic checks.

</div>

# Finite Quotient Ensembles and Claim Tiers

#### Finite quotient ensemble theorem surface.

Fix a finite regulator $`r`$. The physical presentation space is $`\Sigma_r`$, the presentation redundancy groupoid is $`\Gamma_r`$, and the finite physical quotient is
``` math
Q_r=\Sigma_r/\Gamma_r,\qquad \pi_r:\Sigma_r\to Q_r.
```
The quotient removes nonphysical presentation data: gauge representatives, port relabelings, mesh labels, shard or worker identifiers, queue order, repair schedule identifiers, retry counters, timestamps unless declared semantic, hidden carrier coordinates, and inert ancillary labels. If only settled configurations carry probability, the probability space is the normal-form subset
``` math
N_r=n_r(Q_r).
```
The map $`n_r`$ is a normal-form map, not a probability law.

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

For any map $`N:Q\to Q_{\rm nf}`$, the induced map on laws
``` math
\mathcal C_Q(\mu)=N_\#\mu
```
is idempotent:
``` math
\mathcal C_Q^2=\mathcal C_Q.
```
Every law supported on $`Q_{\rm nf}`$ is fixed. Therefore settlement or canonicalization never selects a unique physical probability law by itself.

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
has a null space unless a radial prior, finite parametric family, or source-stress bridge is declared. OPH primordial promotion requires the source-only stress, single-clock, entropy-repair, curvature-evolution, adiabatic-mode, phase-coherence, screen-to-radial-lift, radial-null-space, and forward-projection receipts. Observable CMB promotion further requires frozen source, solver, likelihood, no-data-use, pooled-reducer, and falsification-rule receipts.

#### Claim tiers and required receipts.

Every ensemble-facing run records immutable ensemble id, claim tier, regulator id, representative schema hash, gauge action hash, canonicalizer hash, base measure, action coefficients, coarse maps, zero-mode projector, amplitude convention, sampler hash, smoothing policy, source provenance, and explicit nonclaims. The seed belongs to the run receipt, not to the ensemble definition. The claim tiers are
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

This fragment gives the screen-level spectrum theorem used by the conditional cosmology branch. TT, TE, EE, lensing, likelihoods, and physical CMB promotion belong to the downstream Boltzmann and data-contract gates.

<div id="def:oph-screen-regulator" class="definition">

**Definition 2** (finite screen regulator). *For regulator $`r`$, a finite screen regulator is
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

**Definition 3** (geometric screen scalar). *Let the quotient-visible collar-volume readout be
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

**Proposition 4** (quotient invariance). *If a recharting $`U`$ preserves the screen inner product and geometric readout,
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

**Definition 5** (physical scalar precision and action). *The scalar precision $`K_r`$ must have a physical origin on $`V_r=\operatorname{im}\Pi_{\ge2,r}`$. Allowed branches are:
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

**Theorem 6** (finite MaxEnt screen covariance). *Let $`K_re_i=\kappa_i e_i`$ with an $`M_r`$-orthonormal basis of $`V_r`$, $`\kappa_i>0`$, and $`d_r=\dim V_r`$. Among continuous densities in coefficients $`q=\sum_i a_ie_i`$ with $`\mathbb E[a_i]=0`$ and
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

*Remark 7* (discrete quotient branch). For a discrete finite quotient with intrinsic base weight $`m_r(q)`$, the exact Gibbs law is
``` math
\mu_r(q)=Z_r^{-1}m_r(q)e^{-\beta_r E_r(q)}.
```
A Gaussian statement on that branch requires a separate Laplace or central-limit theorem with controlled higher-order terms.

</div>

<div id="def:oph-screen-release-energy" class="definition">

**Definition 8** (source release energy). *The amplitude source is an independently selected quotient ensemble $`\nu^{\rm src}_r`$ on frozen or released collar normal forms:
``` math
E^{\rm src}_{q,r}
  =\frac12\int \langle q_r(x),K_rq_r(x)\rangle_r\,d\nu^{\rm src}_r(x).
```
The receipt must expose $`\nu^{\rm src}_r`$, the base measure, $`K_r`$, $`d_r`$, $`E^{\rm src}_{q,r}`$, $`A_{q,r}`$, no-observation ancestry, and the same operator normalization used in the action. MaxEnt alone does not determine amplitude.*

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

<div id="thm:oph-screen-refinement" class="theorem">

**Theorem 10** (finite refinement error). *Suppose that for $`2\le\ell\le L`$
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

**Theorem 11** (refinement-semigroup tilt). *Let $`R_b`$ be the scalar refinement map for scale ratio $`b>1`$. If one isolated scalar covariance mode has positive eigenvalue $`\lambda(b)`$, $`\lambda(b_1b_2)=\lambda(b_1)\lambda(b_2)`$, and $`\lambda`$ is continuous, then there is a unique real $`\theta`$ with
``` math
\lambda(b)=b^{-\theta},\qquad
  \theta=-\frac{\log\lambda(b)}{\log b}.
```*

</div>

<div class="proof">

*Proof.* Set $`g(t)=-\log\lambda(e^t)`$. The semigroup law gives $`g(t+s)=g(t)+g(s)`$. Continuity makes $`g(t)=\theta t`$, hence the result. ◻

</div>

<div id="target:oph-p-over-48-tilt" class="target">

**Target 12** (edge-center $`P_\star/48`$ tilt). *The identity $`\theta=P_\star/48`$ is a theorem-grade statement only after three lemmas pass: an oriented edge-center reserve $`\rho_{\rm full}=P_\star/24`$, a scalar half-collar projector $`\rho_q=P_\star/48`$, and a reserve-to-RG Ward identity equating this reserve with $`-d\log\mathcal D_\ell^q/d\log\ell`$. Until then $`P_\star/48`$ and $`e(P_\star-\varphi)`$ are distinct diagnostic hypotheses.*

</div>

<div id="thm:oph-thin-shell-powerlaw-lift" class="theorem">

**Theorem 13** (thin-shell power-law lift). *Assume $`q(\hat n)=Z_\star\Pi_{\ell\ge2}\zeta_\star(R_\star\hat n)`$ and
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

**Proposition 14** (finite-window bound). *Let $`\Psi_\ell(k)=\int dr\,W(r)j_\ell(kr)`$, $`W\ge0`$, $`\int Wdr=1`$, with mean $`R_\star`$ and variance $`\sigma_R^2`$. If
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

**Proposition 15** (radial non-identifiability). *For a finite radial basis, $`C=Ap`$. If $`A\in\mathbb R^{N_\ell\times N_k}`$ has rank $`r`$, then $`\dim\ker A=N_k-r`$, and every $`p+v`$ with $`v\in\ker A`$ gives the same screen spectrum. Radial promotion therefore requires either a source theorem restricting $`p`$, or a declared prior with a published null basis, resolution kernels, positivity checks, and prior-sensitivity report.*

</div>

<div id="target:oph-screen-spectrum-receipts" class="target">

**Target 16** (screen-spectrum receipt set). *The theorem-grade screen spectrum requires the following receipt families before any primordial or CMB promotion:*

1.  *geometry, scalar quotient, low-mode projector, scalar precision, and operator normalization;*

2.  *quotient ensemble selection, scalar release energy, and MaxEnt screen covariance;*

3.  *scalar refinement tilt and angular screen spectrum with a finite error budget;*

4.  *radial window, radial null-space, forward residual, and transfer firewall.*

</div>

# Physical Scale Bridge and Mode Calibration

This fragment gives the OPH physical scale bridge used by the conditional cosmology branch. It separates physical comoving wavenumber, source-screen angular sectors, observed-sky transfer multipoles, finite graph labels, cap labels, and repair-cycle labels.

<div id="def:cosmology-claim-tiers" class="definition">

**Definition 17** (cosmology claim tiers). *Every cosmology artifact has one of three claim tiers:
``` math
\texttt{DIAGNOSTIC\_PROXY},\qquad
\texttt{CONDITIONAL\_PHYSICAL},\qquad
\texttt{OPH\_NATIVE\_PHYSICAL}.
```
`DIAGNOSTIC_PROXY` covers cap angles, graph mode numbers, repair cycles, sorted covariance indices, unit strings, and other finite labels without a physical scale bridge. `CONDITIONAL_PHYSICAL` imports a frozen FLRW geometry packet and proves the downstream spectral, clock, freezeout, source, transfer, and likelihood contracts relative to that packet. `OPH_NATIVE_PHYSICAL` requires the geometry packet and source embedding to be read from the OPH quotient carrier.*

</div>

<div id="def:cosmology-notation-firewall" class="definition">

**Definition 18** (notation firewall). *The bridge uses separate symbols for separate objects:
``` math
T_r,\qquad \tau,\qquad \eta,\qquad \chi,\qquad K_{\rm FLRW},\qquad
k,\qquad \ell_{\rm src},\qquad L_{\rm CMB},\qquad \ell_\star,\qquad L_{\rm IR}.
```
$`T_r`$ is a finite source-clock label. $`\tau`$ is proper time, $`\eta`$ is conformal time, $`\chi`$ is comoving radial distance, and $`K_{\rm FLRW}\in\{-1,0,+1\}`$ is the spatial-curvature branch. $`k`$ is physical comoving wavenumber on a reconstructed spatial slice. $`\ell_{\rm src}`$ is an angular representation sector on the finite primordial source screen. $`L_{\rm CMB}`$ is the observed-sky multipole produced by Boltzmann line-of-sight transfer. Cap openings, graph mode numbers, and repair-cycle indices are diagnostic labels. The canonical internal unit for $`k`$ is $`{\rm Mpc}^{-1}`$; $`h\,{\rm Mpc}^{-1}`$ is a display unit only after $`h`$ is frozen and provenance-typed.*

</div>

<div id="def:finite-cosmo-geometry-packet" class="definition">

**Definition 19** (finite cosmological geometry packet). *At regulator $`r`$, a finite cosmological geometry packet is
``` math
\mathfrak G_r=
\left(
Q_r,\{\Sigma_{r,n}\},\{\mathcal K_{r,n}\},\{h_{r,n}\},\{\mu_{r,n}\},
\{N_{r,n}\},\{\beta^i_{r,n}\},\{u^a_{r,n}\},\mathcal L_r,\mathcal B_r,
\{\iota_{sr}\},\mathcal R_{\ell_\star}
\right).
```
It consists of the quotient carrier $`Q_r`$, common-clock spatial slices, finite cell complexes, positive spatial metric or Dirichlet-form data, volume measures, lapse, shift, observer congruence, lineage, topology and boundary data, refinement maps, and a dimensional scale-certificate receipt. The native readout has type
``` math
\mathsf{CosmoGeomRead}_r:
\operatorname{nf}(\Sigma_r/\Gamma_r)
\longrightarrow
[\mathfrak G_r]_{\rm diff/iso}.
```
The codomain is modulo the declared discrete diffeomorphism, frame, and isometry redundancies. It is a quotient-natural OPH readout only when
``` math
x\sim_{\Gamma_r}y
\quad\Longrightarrow\quad
\mathsf{CosmoGeomRead}_r(x)\simeq\mathsf{CosmoGeomRead}_r(y),
```
and when port relabeling, worker partition, queue order, repair schedule, retry history, and hidden carrier coordinates leave the packet invariant. Positivity and nondegeneracy require
``` math
0<c^-_r |v|^2\le h_{r,n}(v,v)\le c^+_r |v|^2,\qquad
N_{r,n}\ge N_{\min}>0.
```
Metric, volume, area, Jacobian, frame, connection, and transport readouts must agree within declared residuals; for instance $`d\mu_{h_r}=\sqrt{\det h_r}\,d^3x`$ must match $`\mathsf{AreaVolumeRead}`$ and $`\mathsf{VolumeJacobianRead}`$. Physical length is tied to the scale certificate by $`L_{\rm phys}=\ell_\star L_{\rm dimensionless}`$. Refinement naturality requires
``` math
\mathcal R_{sr}\circ\mathsf{CosmoGeomRead}_s
\simeq
\mathsf{CosmoGeomRead}_r\circ n_{sr},
```
with discrepancy tending to zero on the cofinal refinement family.*

</div>

<div id="def:physical-mode-basis" class="definition">

**Definition 20** (physical mode basis). *For each sector $`s=0,1,2`$, let $`\{\phi_i^{(s)}\}`$ be a conforming finite basis on the reconstructed physical geometry. Define
``` math
M^{(s)}_{r,ij}=\int_{\Sigma_r}
\langle\phi_i^{(s)},\phi_j^{(s)}\rangle\,d\mu_{\gamma_r},
```
``` math
K^{(s)}_{r,ij}=
\mathfrak a_r^{(s)}(\phi_i^{(s)},\phi_j^{(s)}).
```
For scalars,
``` math
\mathfrak a_r^{(0)}(f,g)=
\int_{\Sigma_r}\gamma_r^{ab}\nabla_a f\nabla_b g\,d\mu_{\gamma_r}.
```
Vector and tensor sectors use declared Hodge and Lichnerowicz-type operators after gauge and constrained zero modes are removed. The physical basis solves
``` math
K_r^{(s)}V_r^{(s)}=M_r^{(s)}V_r^{(s)}\Lambda_r^{(s)},\qquad
V_r^{(s)\dagger}M_r^{(s)}V_r^{(s)}=I.
```
The evidence bundle must report
``` math
\epsilon_{\rm orth}=\left|V^\dagger M V-I\right|
```
and, for every eigenpair,
``` math
\epsilon_{{\rm eig},j}
=
\frac{
\left|Kv_j-\lambda_jMv_j\right|_{M^{-1}}
}{
\left(|K|+|\lambda_j||M|\right)|v_j|
}.
```
Degenerate sectors are projector-valued. For a spectral interval $`I`$,
``` math
P_{r,I}=\frac{1}{2\pi i}
\oint_{\partial I}\left(z-M_r^{-1}K_r\right)^{-1}dz,
```
and every physical result must be invariant under
``` math
V_{r,I}\longmapsto V_{r,I}U_I,\qquad
U_I\in U(\operatorname{rank}P_{r,I}).
```*

</div>

<div id="def:finite-spectral-measure" class="definition">

**Definition 21** (finite spectral measure and quadrature). *The finite density of states is intrinsic to the physical operator:
``` math
\nu_r(f)=\frac1{V_r}\operatorname{Tr}
f\!\left(\sqrt{M_r^{-1}K_r}\right).
```
Equivalently,
``` math
d\nu_r(k)=\frac1{V_r}\sum_I
\operatorname{rank}(P_{r,I})\,\delta(k-k_{r,I})\,dk.
```
On a large flat domain the resolved-band target is
``` math
d\nu(k)=\frac{k^2}{2\pi^2}\,dk.
```
The $`d\ln k`$ quadrature weights used by the primordial bridge must be derived from this normalization or from another declared quadrature theorem.*

</div>

<div id="def:source-angular-sector" class="definition">

**Definition 22** (source angular sector). *For a closed source screen $`S_r`$, set
``` math
R_{A,r}=\sqrt{\frac{\operatorname{Area}(S_r)}{4\pi}},
\qquad
\overline{\mathcal L}_{S,r}=-R_{A,r}^{2}\Delta_{S_r}.
```
A spectral projector $`P_{\ell,r}`$ may be labeled by $`\ell_{\rm src}`$ only if
``` math
\operatorname{rank}P_{\ell,r}=2\ell+1
```
and
``` math
\epsilon_{\ell,r}=
\frac{
\left|(\overline{\mathcal L}_{S,r}-\ell(\ell+1))P_{\ell,r}\right|
}{
1+\ell(\ell+1)
}
\le \bar\epsilon_\ell.
```
It must also refine by
``` math
\left|J_{sr}P_{\ell,r}-P_{\ell,s}J_{sr}\right|\to0.
```
The scalar primordial contract normally starts at $`\ell_{\rm src}\ge2`$. A cap or masked screen requires either a certified full-screen extension or a certified mode-coupling operator
``` math
\widetilde a_\alpha=\sum_{\ell m}\mathcal K_{\alpha,\ell m}a_{\ell m}.
```
Without that operator the label is $`\texttt{PSEUDO\_ELL}`$ or $`\texttt{CAP\_MODE\_INDEX}`$. The observed multipole is produced only after transfer:
``` math
C_{L_{\rm CMB}}^{XY}
=4\pi\int d\ln k\,
\mathcal P_\zeta(k)\Delta_{L_{\rm CMB}}^X(k)\Delta_{L_{\rm CMB}}^Y(k).
```*

</div>

<div id="def:source-screen-embedding-packet" class="definition">

**Definition 23** (source-screen embedding packet). *A screen mode has radial meaning only after an embedding packet
``` math
\mathfrak E_r=
\left(
S_r,\jmath_r:S_r\hookrightarrow\Sigma_{r,\star},O_r,\chi,W_r(\chi),J_r
\right)
```
has been supplied. $`S_r`$ is the source screen, $`\jmath_r`$ embeds it in the initial slice, $`O_r`$ is the observer or radial origin, $`\chi`$ is comoving geodesic distance, $`W_r`$ is the radial release window, and $`J_r`$ is the area/volume Jacobian. The window obeys
``` math
W_r(\chi)\ge0,\qquad \int W_r(\chi)\,d\chi=1.
```*

</div>

<div id="prop:scale-bridge-rescaling-nogo" class="proposition">

**Proposition 24** (dimensional rescaling no-go). *If the finite data used by a proposed calibration are unchanged under $`h\mapsto s^2h`$, then $`\Delta_{s^2h}=s^{-2}\Delta_h`$ and physical wavenumber transforms as $`k\mapsto k/s`$. Cap angle, graph incidence, patch count, sorted covariance index, and normalized covariance cannot determine nonzero physical $`k`$ without a dimensional metric receipt.*

</div>

<div id="prop:angular-radial-nogo" class="proposition">

**Proposition 25** (angular projection cannot determine radial $`k`$). *At finite angular cutoff the map
``` math
T:f(k)\mapsto\left\{\int d\ln k\,\Psi_\ell(k)f(k)\right\}_{\ell\le L}
```
has finite-dimensional codomain and infinite-dimensional domain. Hence it has a nontrivial kernel. A single angular screen supplies a forward radial projection plus a null-space ledger, not a unique radial spectrum.*

</div>

<div id="prop:clock-scale-nogo" class="proposition">

**Proposition 26** (clock and unit-string no-goes). *A rank-one source family $`\mathcal F=\widehat{\mathcal F}\circ T_r`$ fixes a common ordering, not proper elapsed time, because $`f\circ T_r`$ gives the same fibers for any strictly increasing $`f`$. Likewise $`h=a^2\gamma`$ is invariant under $`a\mapsto ca,\ \gamma\mapsto c^{-2}\gamma`$. A lapse or operational clock receipt and the normalization $`a(\tau_0)=1`$ are separate physical inputs. Relabeling $`\texttt{inverse\_cap\_opening\_angle\_proxy}`$ as $`\texttt{Mpc}^{-1}`$ changes no numerical object and proves no calibration.*

</div>

<div id="thm:physical-comoving-k" class="theorem">

**Theorem 27** (physical comoving $`k`$). *Assume a valid geometry packet, FLRW reduction, dimensional ruler, topology and boundary sector, solver convention, physical mode operator, projector-valued basis, density of states, and mode normalization. On the flat comoving branch,
``` math
-\Delta_{\gamma_r}v_{rj}=\lambda_{\gamma,rj}v_{rj},\qquad
k_{rj}=\sqrt{\lambda_{\gamma,rj}}.
```
If the operator is assembled from $`h_{ij}^{\rm phys}(\tau)=a(\tau)^2\gamma_{ij}`$, then
``` math
\lambda_{h,rj}=\frac{k_{rj}^2}{a_r^2},\qquad
k_{rj}=a_r\sqrt{\lambda_{h,rj}}.
```
For dimensionless eigenvalue $`\widehat\lambda_{h,rj}=\ell_\star^2\lambda_{h,rj}`$,
``` math
k_{rj}=\frac{a_r}{\ell_\star}\sqrt{\widehat\lambda_{h,rj}}.
```
The verifier recomputes this equation and reports
``` math
\epsilon_{k,rj}
=
\frac{|k_{rj}^2-a_r^2\lambda_{h,rj}|}{k_{rj}^2+k_{\rm floor}^2}.
```
For $`K_{\rm FLRW}\ne0`$, the solver convention must specify
``` math
k_{\rm solver}=\kappa_{K_{\rm FLRW}}(\lambda)
```
and use the matching hyperspherical kernels.*

</div>

<div id="thm:screen-to-spatial-mode-association" class="theorem">

**Theorem 28** (screen-to-spatial mode association). *There is no physical theorem of the form $`\ell_{\rm src}\mapsto k`$. The spatial geometry defines the physical $`k`$-basis, and the source embedding plus radial window define a forward projection from that basis to the screen:
``` math
\Psi_{\ell j}^{(K)}=
\int d\chi\,W_r(\chi)\,\Phi_{\ell,k_j}^{(K)}(\chi),
```
where $`\Phi_{\ell,k}^{(0)}=j_\ell(k\chi)`$ on the flat branch and the curved branches use the declared hyperspherical functions. Then
``` math
C_{\ell_{\rm src},r}^{q}
=4\pi\sum_j w_{rj}^{\ln k}
\Delta_{\zeta,r}^{2}(k_{rj})
\left|\Psi_{\ell_{\rm src},j}^{(K)}\right|^2.
```
If an inverse reconstruction $`C^q=Ap`$ is attempted, the evidence reports prior and support, positivity constraints, singular values, effective rank, null basis, resolution kernels, prior sensitivity, and forward residuals. The Limber estimate $`L+1/2\simeq k\chi_\star`$ is a projection check on a completed bridge; it does not define $`k`$, $`L_{\rm CMB}`$, or $`\chi_\star`$.*

</div>

<div id="thm:physical-time-scale" class="theorem">

**Theorem 29** (proper time, scale factor, conformal time, and redshift). *Let the selected congruence be
``` math
u^a=N^{-1}\left[(\partial_t)^a-\beta^a\right].
```
Proper time is operationally $`d\tau=N\,dt`$ along the congruence, with shift and threading handled by the lineage map. For a lineage-tracked comoving domain $`D_r(\tau)`$,
``` math
V_{D,r}(\tau)=\int_{D_r(\tau)}d\mu_{h_r},\qquad
a_{V,r}(\tau)=
\left[\frac{V_{D,r}(\tau)}{V_{D,r}(\tau_0)}\right]^{1/3}.
```
Independently,
``` math
\log\frac{a_{\Theta,r}(\tau)}{a_{\Theta,r}(\tau_0)}
=
\int_{\tau_0}^{\tau}
\frac{\langle\Theta\rangle_{D,r}}{3}\,d\tau',
\qquad
\Theta=\nabla_a u^a.
```
The scale-factor consistency residual is
``` math
\epsilon_{a,r}=
\sup_\tau
\left|\log a_{V,r}(\tau)-\log a_{\Theta,r}(\tau)\right|.
```
The Hubble variables are $`H=a^{-1}da/d\tau`$ and $`\mathcal H=aH`$, while
``` math
\eta(\tau)-\eta(\tau_0)=\int_{\tau_0}^{\tau}\frac{d\tau'}{a(\tau')}.
```
On a certified comoving FLRW branch $`1+z=1/a`$. Outside that branch,
``` math
1+z=\frac{(u_\mu p^\mu)_{\rm emit}}{(u_\mu p^\mu)_{\rm obs}}.
```
The FLRW residual vector contains normalized shear, vorticity, acceleration, spatial-curvature gradient, domain backreaction, Hamiltonian constraint, momentum constraint, and energy-conservation components. Friedmann and conservation equations are consistency checks, not definitions chosen to force a desired $`a(\tau)`$.*

</div>

<div id="thm:mode-freezeout-common-surface" class="theorem">

**Theorem 30** (mode freezeout and common initial surface). *For every resolved mode $`k_j`$, define
``` math
\mathcal E_j(\tau)=
\left(
\epsilon_{\zeta,j},
\epsilon_{{\rm nad},j},
\epsilon_{{\rm iso},j},
\epsilon_{{\rm dec},j},
\epsilon_{{\rm rep},j},
\epsilon_{{\rm phase},j},
\epsilon_{{\rm constr},j}
\right).
```
The physical mode-freezeout time is
``` math
\tau_{\star,j}=
\inf\left\{
\tau:\mathcal E_j(\tau')\le\overline{\mathcal E}
\quad\forall \tau'\in[\tau,\tau+\Delta\tau_{\rm pers}]
\right\}.
```
This proves a $`\texttt{PHYSICAL\_MODE\_FREEZEOUT\_MAP\_RECEIPT}`$. A Boltzmann initial-value problem additionally requires a common spacelike Cauchy slice. For $`\Sigma_{\star,r}=\{x:T_r(x)=T_{\star,r}\}`$ and signature $`(-+++)`$, spacelikeness requires
``` math
g^{ab}\nabla_aT_r\nabla_bT_r<0
```
on the resolved surface. A common surface is certified either by coincident freezeout,
``` math
\Delta_\star=
\sup_{j\in\mathcal B_{\rm safe}}
\left|H_\star(\tau_{\star,j}-\bar\tau_\star)\right|
\le\bar\delta_\star,
```
or by evolution to a later common slice $`\tau_c`$ with no independent stochastic source and a bounded conserved-state residual. The finite initial packet contains
``` math
\mathfrak I_{\star,r}=
\left(
\Sigma_{\star,r},h_{ij}^\star,K_{ij}^\star,
\{X_{\star,rj}\},\{n^a\nabla_aX_{\star,rj}\},
\mathcal C_{\rm H},\mathcal C_{\rm M},
\hbox{phase convention},\hbox{basis hashes}
\right).
```
A mode-dependent map alone cannot satisfy $`\texttt{PHYSICAL\_FREEZEOUT\_SURFACE\_RECEIPT}`$.*

</div>

<div id="thm:common-primordial-anomaly-basis" class="theorem">

**Theorem 31** (common primordial and anomaly basis). *Primordial covariance and dark/anomaly response use the same physical mode projectors. For a finite anomaly response operator $`\mathcal K_{A,r}^{(\rho)}(a)`$,
``` math
K_{A,r;IJ}^{(\rho)}(a)=
P_{r,I}\mathcal K_{A,r}^{(\rho)}(a)P_{r,J}.
```
The physical contrast response is initially a matrix,
``` math
B_{A,r;IJ}(a)=
\frac{\bar\rho_b(a)}{\bar\rho_{A,\rm eq}(a)}
K_{A,r;IJ}^{(\rho)}(a).
```
It may be reduced to a scalar $`B_A(k,a)`$ only when
``` math
\epsilon_{\rm off}(a)=
\frac{\left(\sum_{I\ne J}|P_I\mathcal K_AP_J|^2\right)^{1/2}}
{|\mathcal K_A|+\epsilon_0}
\le\bar\epsilon_{\rm off}
```
and, inside each degenerate sector,
``` math
P_I\mathcal K_AP_I=b_I(a)P_I+E_I,\qquad
\frac{|E_I|}{|b_I|+\epsilon_0}\le\bar\epsilon_{\rm iso}.
```
The same rule applies to repair rates. A transition eigenvalue per repair cycle is a diagnostic number. With certified physical interval $`\Delta\tau_r`$ and a contractive active response block,
``` math
\Gamma_{{\rm rec},I}
=-\frac1{\Delta\tau_r}\log\rho(P_IL_{\perp,r}P_I),
```
or a declared Lyapunov/semigroup norm bound for a nonnormal block.*

</div>

<div id="thm:scale-bridge-regulator-convergence" class="theorem">

**Theorem 32** (regulator convergence). *Let $`r`$ be a directed regulator containing UV mesh size, IR volume, time step, sampling size, and boundary approximation. Physical scale-bridge receipts require Cauchy or convergence statements:
``` math
\left|J_{sr}P_{r,I}-P_{s,I}J_{sr}\right|\le\epsilon_{sr}^{\rm spec},\qquad
\epsilon_{sr}^{\rm spec}\to0,
```
``` math
\left|\int f(k)\,d\nu_s(k)-\int f(k)\,d\nu_r(k)\right|\to0
```
for every bounded test function $`f`$,
``` math
|a_s-a_r|_\infty+|H_s-H_r|_\infty+|\eta_s-\eta_r|_\infty\to0,
```
and
``` math
d_H(\Sigma_{\star,s},J_{sr}\Sigma_{\star,r})\to0,\qquad
\left|J_{sr}^{*}X_{\star,s}-X_{\star,r}\right|\to0.
```
The directed limit must be independent of the order in which UV, IR, sampling, and time refinements are taken. Every residual $`e_\alpha(r)`$ has a frozen tolerance $`\bar e_\alpha`$, a proof that $`e_\alpha(r)\le\bar e_\alpha`$, and a refinement statement $`e_\alpha(r)\to0`$.*

</div>

<div id="thm:finite-physical-scale-bridge" class="theorem">

**Theorem 33** (finite physical scale bridge). *Assume valid receipts for geometry or imported conditional geometry, spatial physical $`k`$, screen-to-$`k`$ association, source angular sector, proper time and scale-factor history, mode-freezeout map, common primordial initial surface, common primordial/anomaly mode basis, regulator convergence, cross-receipt identity, source-only no-data-use provenance, and no post-hoc calibration. Then
``` math
\mathcal S_{{\rm OPH},r}:
\left(
P^{\rm screen}_{r,\ell},P^{\rm spatial}_{r,I},T_{r,n}
\right)
\longmapsto
\left(
k_{r,I},\ell_{\rm src},a_r(\tau_n),z_r(\tau_n),\eta_r(\tau_n),
\Sigma_{\star,r},\mathfrak I_{\star,r}
\right)
```
is unique modulo the declared spatial diffeomorphism/isometry quotient, unitary rotations inside degenerate spectral projectors, frozen topology and boundary sector, curvature convention, and solver convention. It converges on the declared safe physical band. The tier is $`\texttt{CONDITIONAL\_PHYSICAL}`$ when $`\mathfrak G_r`$ is imported and frozen. The tier is $`\texttt{OPH\_NATIVE\_PHYSICAL}`$ when $`\mathsf{CosmoGeomRead}_r`$ and the source embedding are proved from the OPH quotient carrier.*

</div>

<div id="target:physical-scale-receipts" class="target">

**Target 34** (exact physical-scale receipt conjunctions). *The code and papers use recomputable receipt conjunctions. Physical spatial $`k`$ is
``` math
\begin{aligned}
R_{k,\rm spatial}={}&
R_{\rm geometry}\wedge R_{\rm scale}\wedge R_{\rm FLRW}\wedge
R_{\rm operator}\wedge R_{\rm eigenresidual}\wedge R_{\rm normalization}\\
&\wedge R_{\rm k\ equation}\wedge R_{\rm lineage}\wedge
R_{\rm refinement}\wedge R_{\rm nofit}.
\end{aligned}
```
Then
``` math
\texttt{physical\_k\_units\_calibrated}
\iff
\texttt{PHYSICAL\_SPATIAL\_K\_RECEIPT}
\iff
R_{k,\rm spatial}.
```
The screen association is
``` math
R_{\rm screen\to k}
=
R_{k,\rm spatial}\wedge R_{\rm embedding}\wedge R_{\rm radial\ kernel}
\wedge R_{\rm geometry\ compatibility}\wedge R_{\rm nullspace}
\wedge R_{\rm forward\ residual}.
```
The angular sector is
``` math
R_{\ell_{\rm src}}=
R_{\rm closed/extension}\wedge R_{\rm angular\ operator}\wedge
R_{\rm cluster}\wedge R_{\rm multiplicity}\wedge R_{\rm angular\ refinement}.
```
The scale-factor history is
``` math
\begin{aligned}
R_a={}&
R_{\rm proper\ time}\wedge R_{\rm lineage}\wedge R_{\rm volume}
\wedge R_{\Theta}\wedge R_{a_V=a_\Theta}\\
&\wedge R_{\rm FLRW}\wedge R_{\rm constraints}\wedge
R_{\rm refinement}\wedge R_{\rm nofit}.
\end{aligned}
```
The mode-freezeout map is
``` math
R_{\rm modefreeze}=
R_{\zeta}\wedge R_{\rm adiabatic}\wedge R_{\rm isocurvature}
\wedge R_{\rm growing}\wedge R_{\rm repair}\wedge R_{\rm phase}
\wedge R_{\rm constraints}\wedge R_{\rm persistence}.
```
The common surface is
``` math
R_{\Sigma_\star}=
R_{\rm modefreeze}\wedge R_{\rm spacelike}
\wedge (R_{\rm coincident}\vee R_{\rm transport\ to\ common})
\wedge R_{\rm initial\ data}\wedge R_{\rm refinement}.
```
Only $`R_{\Sigma_\star}`$ instantiates $`\texttt{PHYSICAL\_FREEZEOUT\_SURFACE\_RECEIPT}`$. The aggregate scale bridge is
``` math
R_{\rm scale}=R_{k,\rm spatial}\wedge R_{\rm screen\to k}\wedge
R_{\ell_{\rm src}}\wedge R_a\wedge R_{\Sigma_\star}\wedge
R_{\rm common\ basis}\wedge R_{\rm cross}\wedge R_{\rm nofit}.
```
Every child receipt repeats the same regulator family, generation, geometry, background, clock, scale certificate, source embedding, mode basis, mode lineage, boundary condition, solver convention, freezeout-surface, and source-DAG hashes. An empty cross-receipt manifest fails.*

</div>

<div id="target:cosmo-geometry-extraction" class="target">

**Target 35** (native geometry extraction). *$`\texttt{OPH\_NATIVE\_PHYSICAL}`$ requires
``` math
\mathsf{CosmoGeomRead}_r:\operatorname{nf}(\Sigma_r/\Gamma_r)
\to(\mathfrak G_r,\mathfrak E_r)
```
as a quotient-natural, positive, nondegenerate, causal-order-compatible, observer-clock-compatible, executor-independent, refinement-convergent, physically scaled source embedding. Without that construction, the physical route has tier $`\texttt{CONDITIONAL\_PHYSICAL}`$ with an imported frozen FLRW geometry packet.*

</div>

<div id="target:physical-scale-falsifiers" class="target">

**Target 36** (falsifiers). *The physical scale gate fails under any of the following conditions: dimensional rescaling failure; unit-string relabeling; direct cap-angle calibration; direct $`\ell_{\rm src}\to k`$ inversion from one angular shell without null-space and radial-prior receipts; CMB peak alignment; conflating $`\ell_{\rm src}`$ with $`L_{\rm CMB}`$; dependence on rotations inside degenerate eigenspaces; hidden mode mixing in scalar $`B_A(k,a)`$ tables; repair cycles used as proper time; excessive scale-factor, FLRW, or constraint residuals; a mode-dependent freezeout map reported as a common surface; non-spacelike release; lack of persistence; regulator drift; limit-order dependence; vacuous cross receipts; or imported geometry labeled as native.*

</div>

# Compatibility With the Core Stack

The compatibility table fixes how the inflation-free draft is read against the core papers.

| Imported draft idea | Core-compatible reading |
|:---|:---|
| Inflation-free replacement | Phase III conditional branch. The compact paper’s recovered core covers the SM/GR package; inflation replacement, physical CMB likelihood, the $`H_0/S_8`$ branch, and dark-sector cosmology sit outside that core. |
| Finite normal form | Use quotient-level fixed-cutoff normal forms on the declared branch, unique only modulo boundary redundancy, implementation hiding, inert ancillary stabilization, and any required same-boundary unique-extension condition. Do not read the old draft as selecting a unique microscopic representative. |
| Screen spectrum and $`n_s=1-P/48`$ | Screen-level theorem only after geometric $`q_r`$, normalized $`K_r`$, source release energy, refinement tilt, and angular-spectrum receipts pass. The special $`P/48`$ tilt is an analytic candidate, not lattice-derived, with the $`\kappa_{\mathrm{rep}}`$ certificate pending; two candidate tilt formulas are on record ($`1-P_\star/48`$ and $`1-\mathrm e\,\alpha\sqrt{\pi}`$) with the selection uncounted, and the $`0.27\sigma`$ Planck comparison is a basin location, not a hit. It is theorem-grade only after the edge-center reserve, half-collar projection, and reserve-to-RG receipts pass. Physical $`A_s`$, $`n_s`$, running, isocurvature, phase coherence, and TT/TE/EE spectra require the separate source-stress, single-clock, radial-prior, null-space, forward-residual, transfer, and likelihood receipts. |
| Dark/anomaly slot | Imported from `cosmology/oph_dark_matter_paper.tex`. For a source-only primordial certification run, the dark continuation is `OFF` unless the dark paper's source-only anomaly abundance selector passes. A supplied abundance is explicitly typed `CONDITIONAL_SOURCE_STATE`. Transported $`Q_A`$ without the selector has label `PHYSICAL_PARENT_WITH_CONDITIONAL_ABUNDANCE`. |
| Scalar anomaly rows | Rows for $`\bar\rho_A`$, $`\bar\rho_{A,\mathrm{eq}}`$, and $`B_A`$ are diagnostics unless the finite covariant parent also emits $`w_A`$, $`c_{s,A}^2`$, $`\sigma_A`$, $`Q_A^\mu`$, exchange-current closure, recipient stress for nonzero exchange, gauge-independent variables, CDM-limit recovery, and refinement convergence. A transition-matrix number is $`\gamma_{\mathrm{repair\ step}}`$ until promoted to $`\Gamma_{\rm rec}`$ by the physical clock, active-fiber, and common-parent response receipts. |
| Boltzmann and likelihood comparison | Physical promotion requires immutable source, solver, and likelihood hashes, official likelihood execution, and global CDM-limit reductions. Shard-local `any()` rollups or nonlinear averages before global pooling do not satisfy the contract. |
| Vacuum language | Finite reference states and free-field/lattice baselines remain reference ensembles. OPH-native vacuum promotion requires the quotient-ensemble selector plus source Euclidean slab data and a transfer/reflection-positive reconstruction gate from the compact/screen-microphysics theorem surface. |
| Simulator receipts | Any physical-promotion receipt cited here must carry exact source/config/solver/likelihood hashes. The 256k run used for diagnostics records `git_commit=unknown`, so it is not used as exact replay provenance. |

# Source Corpus

This paper consolidates these active surfaces.

| Surface | Imported role |
|:---|:---|
| `cosmology/oph_inflation_without_inflaton_observer_screen_synchronization.tex` | Main inflation-free branch: finite screen federation, flatness, horizon coherence, geometric screen spectrum, screen amplitude, radial lift, MaxEnt release, anomaly slot, transfer, and likelihood gates. |
| `cosmology/physical_cmb_theorem_program.md` | Theorem-contract file identifying the physical CMB source, scale, stress, kernel, transfer, and likelihood contracts. |
| `paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex` | Compact recovered-core claim boundary and the first surface for later propagation. |
| `paper/tex_fragments/PAPER.tex` | Main technical paper surface for global OPH theorem language. |
| `paper/tex_fragments/PRIMORDIAL_BRIDGE_THEOREMS.tex` | Primordial bridge theorem boundary and screen-to-primordial finite-source conditions. |
| `paper/screen_microphysics_and_observer_synchronization.tex` | Finite observer-screen carriers, records, patch ports, evidence bundles, and implementation invariance. |
| `cosmology/oph_dark_matter_paper.tex` | Collar-remainder dark/anomaly stress, galaxy equilibrium, dynamic stress, and cosmological kernel contracts. This paper imports those contracts by reference; it does not duplicate the dark-sector derivations. |
| `oph-physics-sim/runs/gcp_large_vis_256k_modkernel_20260623` | 256k-patch simulator evidence: strong CMB diagnostics, physical-prediction receipt with failed gates, and full visualizer export data. The run manifest records `git_commit=unknown`; the canonical simulator repository is <https://github.com/muellerberndt/oph-physics-sim>, and the public source reference available for this consolidation is `aeea9e502491a277b5e21a80bb89df6d089a074d`. |

<div class="remark">

*Remark 37* (Simulator receipt provenance). Whenever this paper mentions simulator receipts, the receipt artifacts should be traced to the responsible commit in the canonical simulator repository, <https://github.com/muellerberndt/oph-physics-sim>. For the 256k run used here, the local artifact `runs/gcp_large_vis_256k_modkernel_20260623/manifest.json` reports `git_commit=unknown`. Therefore the numerical run is cited as diagnostic evidence with incomplete commit provenance. The available public source context for this consolidation is <https://github.com/muellerberndt/oph-physics-sim/commit/aeea9e502491a277b5e21a80bb89df6d089a074d>, but that link is not asserted to be an exact run-replay hash for the 256k artifacts. Physical-promotion receipts require exact source, config, solver, and likelihood hashes.

</div>

<div class="remark">

*Remark 38* (Ownership boundary for the dark-matter paper). The released dark-matter paper `cosmology/oph_dark_matter_paper.tex` keeps ownership of the dark-sector theory: collar-remainder stress, galaxy equilibrium, dynamic repair-stress transport, cluster behavior, and the detailed perturbation-kernel derivations. This cosmology program paper only records the interface that paper must expose to the physical CMB/large-scale-structure pipeline. When it names $`\bar\rho_A(a)`$, $`B_A(k,a)`$, stress closure, repair-stress transport, or a certified $`\Gamma_{\rm rec}(k,a)`$, it is naming imported source functions and receipts, not rederiving the dark sector here.

</div>

# Imported Theorem Map

The inflation-free draft contains many theorem statements. The table gives them a public-core landing surface and keeps their dependency roles visible.

| Imported theorem family | Role in this paper |
|:---|:---|
| Finite synchronization normal form | Defines the early observer-facing quotient state whose repair normal form replaces an assumed smooth primordial substrate, with uniqueness read modulo the compact-paper quotient equivalences. |
| Clocked FLRW curvature holonomy, flat-sector selection, and curvature bounds | Owns the flatness branch: zero clock-slice spatial Levi-Civita holonomy identifies $`\kappa=0`$; exact selection requires a direct theorem, conditional CMH selector, or explicit flat-branch assumption. Curvature damping supplies only $`|K|`$ or $`|\Omega_K|`$ bounds. |
| Diffusive horizon no-go, same-boundary coherence, and low-$`k`$ repair gap | Separates a failed purely local diffusion story from the allowed same-boundary or uniform low-$`k`$ synchronization mechanisms. |
| Geometric screen scalar, precision operator, MaxEnt screen covariance, Mellin lift, edge-center scalar opportunities, $`\mathbb Z_6`$ reserve, and half-collar sampling | Owns the screen-level near-scale-invariant covariance. The red tilt $`P/48`$ is an analytic candidate ($`\kappa_{\mathrm{rep}}`$ certificate pending, not lattice-derived) and is theorem-grade only after the reserve-to-RG receipts pass; physical primordial promotion waits for the source-only radial bridge receipts. |
| Finite release collar state, scalar release code, and scalar amplitude theorem | Owns $`A_\zeta`$ as a finite source artifact rather than a fitted amplitude. |
| MaxEnt release, entropy transport, BBN-safe release, recombination inheritance, and adiabaticity | Connects the synchronized finite screen/collar state to a hot radiation branch with ordinary low-energy thermal history. |
| Freezeout, scalar quotient record, gauge-invariant $`\zeta`$, isocurvature decay, and superhorizon conservation | Imported only through the stricter primordial bridge: total stress closure, single clock, repair gap, freezeout, growing mode, isocurvature, phase coherence, radial prior, null-space, and forward residual. |
| Finite homogeneous anomaly charge, no double counting, packet moments, cold anomaly, Bianchi exchange, no-slip lensing, and linear anomaly transfer | Imported from `cosmology/oph_dark_matter_paper.tex`; this paper only records the CMB/LSS handoff interface and receipt gates. |
| Universe Simulation certificate bundle, no-data-use theorem, finite-certificate branch, Boltzmann stability, likelihood promotion, and CMB admissibility | Owns the promotion rule: source-only finite artifacts first, frozen transfer and likelihood second, physical CMB claim last. |
| Tensor repair gap and finite-collar non-Gaussianity | Records secondary predictions to be tested only after the primary scalar/source/stress gates close. |

# Extended Cosmology Agenda

The local cosmology notes identify several OPH-native prediction surfaces beyond the immediate CMB gate. They are program tracks, not independent claims.

1.  **CMB as observer-consensus fossil.** The highest-value target is the CMB: largest observable scales should carry relics of finite observer synchronization rather than a fitted inflationary initial spectrum.

2.  **$`\Lambda`$ and observer capacity.** The screen-capacity branch relates the de Sitter scale to finite horizon capacity. Under principle SL-4 the capacity $`N`$ is estimated from measured $`\Lambda`$ with the scale bridge (`../STRANGE_LOOP_PRINCIPLES.md`); $`\Lambda`$ is an input, not an OPH output. A full cosmology paper should derive the equation of state, possible time evolution, and any deviation from $`w=-1`$. The missing object is the readback map $`F`$ with a certified contraction interval — generator G2 of `../CONSISTENCY_STACK.md`, closing CL-7 and moving CL-3/CL-4; its formal specification lives at `code/capacity_readback/F_READBACK_SPEC.md`.

3.  **Full cosmological dark sector.** The dark paper’s galaxy layer is not enough. CMB, BAO, weak lensing, growth, and $`S_8`$ need the transported repair-stress perturbation branch.

4.  **$`H_0`$ and $`S_8`$ predictions.** With the scale bridge and anomaly kernels fixed, the branch should emit blind values or intervals for $`H_0`$ and $`S_8`$, without tuning toward observed tensions.

5.  **Large-scale topology.** Finite observer-accessible screens make global topology and matching-circle style signatures natural tests.

6.  **Consensus hierarchy of structure.** Galaxy, cluster, and supercluster formation should be investigated as a hierarchy of stable consensus attractors, with density perturbation growth in a fixed background as one layer of the story.

7.  **Cosmological arrow of time.** OPH has native records, checkpoints, and continuation laws; these should be used to formulate the low-entropy past condition as a finite-record boundary condition.

8.  **Cosmological neutrino background.** The particle/mass branch should eventually feed neutrino background and free-streaming predictions into the same frozen likelihood contract.

# Finite Cosmological State

<div class="definition">

**Definition 39** (Early OPH screen federation). *An early OPH cosmological state at regulator $`r`$ is a finite patch federation
``` math
\mathcal F_r=
\left(
V_r,E_r,\{\mathcal A_i,\rho_i,\mathcal R_i\}_{i\in V_r},
\{\mathcal I_e,\pi_{i,e}\}_{e\in E_r},
\mathcal U_r,\mathrm{Chk}_r
\right),
```
with finite accessible algebras, record algebras, visible overlap interfaces, repair instruments, and checkpoint data. Its physical content is the quotient under hidden carrier coordinates that do not affect visible interface records, repair maps, or checkpoint continuation.*

</div>

<div class="definition">

**Definition 40** (Source-only finite artifact). *A finite artifact $`X`$ is source-only for cosmology if it is a deterministic functional of declared OPH source data, release-branch constants, and standard non-CMB physical constants, and if its dependency manifest excludes CMB, BAO, supernova, weak-lensing, RSD, SPARC, cluster, or other observational likelihood values used to evaluate the prediction.*

</div>

<div class="definition">

**Definition 41** (Physical CMB promotion). *An OPH CMB run is physically promoted only when the primordial source, scale bridge, dark/anomaly stress kernels, Boltzmann solver input files, solver version/tolerances, and likelihood code/data hashes are all frozen before likelihood evaluation, and all source functions are source-only finite artifacts.*

</div>

<div class="definition">

**Definition 42** (Cosmology artifact promotion ladder). *Cosmology artifacts first carry the three claim tiers of Definition <a href="#def:cosmology-claim-tiers" data-reference-type="ref" data-reference="def:cosmology-claim-tiers">17</a>. Within a tier, source artifacts also carry an explicit semantic type
``` math
\textsc{DiagnosticProxy}
\prec
\textsc{FiniteScreenSource}
\prec
\textsc{FinitePrimordialSource}
\prec
\textsc{PhysicalBoltzmannSource}
\prec
\textsc{PhysicalCMBPrediction}.
```
No receipt may skip a boundary in this order. A close TT curve, a scalar CMI row, a cap-angle mode number, a repair-clock eigenvalue, or a visual freezeout surface is therefore `DIAGNOSTIC_PROXY` until the required source, physical-scale, stress, transfer, and likelihood receipts promote it step by step. Imported FLRW geometry can promote downstream labels only to `CONDITIONAL_PHYSICAL`; `OPH_NATIVE_PHYSICAL` additionally requires Target <a href="#target:cosmo-geometry-extraction" data-reference-type="ref" data-reference="target:cosmo-geometry-extraction">35</a>.*

</div>

<div class="definition">

**Definition 43** (Proof-carrying source artifact). *A promotable artifact is a record
``` math
\mathfrak A_v=
(id_v,\tau_v,x_v,\mathcal X_v,\mathcal U_v,\mathrm{Par}(v),
\mathcal D_v,\mathcal E_v,h_v^{\rm theorem},h_v^{\rm witness},
h_v^{\rm build},h_v^{\rm artifact}),
```
where $`\tau_v`$ is the semantic type, $`x_v`$ is the value or tensor/function, $`\mathcal X_v`$ records axes, domain, grid, and support, $`\mathcal U_v`$ records units, gauge, frame, pivot, and conventions, $`\mathrm{Par}(v)`$ is the complete parent set, $`\mathcal D_v`$ is the data-class ledger, $`\mathcal E_v`$ is the numerical and theorem error envelope, and the hashes commit to theorem version, witness, build environment, and artifact payload.*

</div>

<div class="theorem">

**Theorem 44** (Sound promotion by source DAG induction). *Let $`G=(V,E)`$ be a finite acyclic artifact graph. For each node $`v`$, let $`P_v(\mathfrak A_v)`$ be the mathematical predicate claimed by the artifact and let $`V_v(\mathfrak A_v,w_v)`$ be a sound verifier, so
``` math
V_v(\mathfrak A_v,w_v)=1\Longrightarrow P_v(\mathfrak A_v).
```
Define
``` math
\operatorname{Pass}(v)
=V_v(\mathfrak A_v,w_v)
\wedge\bigwedge_{p\in\operatorname{Par}(v)}\operatorname{Pass}(p)
\wedge\operatorname{SourceOnly}(v).
```
Then $`\operatorname{Pass}(v)=1`$ implies $`P_v(\mathfrak A_v)`$, and every required parent predicate holds.*

</div>

<div class="proof">

*Proof.* Order the finite DAG topologically. Root nodes satisfy their predicates by local verifier soundness. If every parent of $`v`$ has passed, then the local verifier proves $`P_v`$ and the parent hypotheses are available. The claim follows by induction. ◻

</div>

<div class="target">

**Target 45** (Source ancestry and freeze receipts). *Physical promotion requires
``` math
\texttt{TRANSITIVE\_SOURCE\_ANCESTRY\_RECEIPT},\quad
\texttt{HERMETIC\_READ\_SET\_RECEIPT},\quad
\texttt{SOURCE\_MODEL\_FREEZE\_RECEIPT}.
```
The first receipt walks every ancestor of every promoted quantity and rejects any CMB, BAO, supernova, weak-lensing, growth, cluster, or similar evaluation-data dependency. The second requires an intercepted, content-addressed read set for files, network, environment, command-line configuration, caches, randomness, code, dependencies, and human-selected branch identifiers. The third freezes branch, model family, basis, pivot, cutoff, packetization, smoothing, thresholds, freezeout tolerances, and hyperparameter ranges before comparison.*

</div>

# Inflation-Free Branch

The inflation-free branch replaces inflation’s jobs rather than reproducing an inflaton. Its load-bearing theorem stack is:

1.  quotient-level observer-facing normal form for the early screen/collar federation, with the compact-paper quotient equivalences;

2.  zero clock-slice spatial Levi-Civita holonomy to identify the flat branch, plus a direct flatness theorem, conditional CMH selector, or explicit flat-branch assumption to select it;

3.  same-boundary scalar normal form or a low-$`k`$ repair gap for horizon coherence;

4.  geometric screen spectrum for near scale invariance;

5.  protected $`\mathbb Z_6`$ half-collar reserve for the red spectral exponent;

6.  finite scalar release code for $`A_\zeta`$;

7.  MaxEnt release and common release clock for a hot adiabatic start;

8.  imported dark/anomaly stress variables for the pre-recombination CDM-like slot, when the `cosmology/oph_dark_matter_paper.tex` source contract permits them;

9.  Boltzmann transfer and frozen likelihood for observational comparison.

<div class="theorem">

**Theorem 46** (Conditional OPH replacement of inflation). *Assume a cofinal family of finite early OPH screen federations has a quotient-level observer-facing normal-form projection on a declared branch, with uniqueness understood modulo boundary redundancy, implementation hiding, inert ancillary stabilization, and any same-boundary unique-extension condition required by the branch. Assume the clocked FLRW continuation boundary $`B_{\cos,r}`$ supplies $`u_a=-N\nabla_a\chi`$, $`h_{ab}=g_{ab}+u_au_b`$, $`K(\tau)=\kappa/a^2`$, and topology policy, and that one of the exact flat-sector gates is present: <span class="smallcaps">DirectTheorem</span>, <span class="smallcaps">ConditionalCMH</span>, or <span class="smallcaps">ExplicitAssumption</span>. Assume every CMB-scale scalar mode shares the same boundary normal form or satisfies a uniform low-$`k`$ repair gap, the finite screen-spectrum receipts derive a geometric $`q_r`$, normalized $`K_r`$, source release energy $`A_q`$, and refinement tilt, the source-only radial bridge receipts pass, and the imported dark/anomaly sector supplies a permitted CDM-like source state through recombination. Then the branch supplies flatness, horizon coherence, near scale invariance, hot adiabatic initial data, and acoustic-transfer input without an inflaton degree of freedom.*

</div>

<div class="proof">

*Proof.* The proof is a dependency composition. The quotient normal-form projection gives a finite observer-facing early state on the declared branch and leaves hidden representatives unselected. The spatial Levi-Civita holonomy theorem identifies the $`\kappa=0`$ branch; the direct, CMH, or explicit-assumption gate selects it. Curvature damping is not used as a selector. Same-boundary repair or a low-$`k`$ gap gives coherence on the observed CMB band. The screen-spectrum theorem package first gives
``` math
C_\ell^q
=A_q\frac{\Gamma(\ell-\theta/2)}{\Gamma(\ell+2+\theta/2)}.
```
The value $`\theta=P_\star/48`$ follows only if the edge-center reserve, half-collar scalar projection, and reserve-to-RG receipts pass; otherwise $`\theta`$ is measured from the actual scalar refinement operator or kept as a diagnostic hypothesis. The source-only bridge receipts then type that screen scalar as a primordial curvature source and convert $`A_q`$ to $`A_\zeta`$ only after the exact thin-shell or finite-window radial lift, null-space, and forward-residual receipts pass. Without them it remains a screen theorem. MaxEnt release gives a hot radiation state and the common clock gives adiabatic leading perturbations. If the imported anomaly stress branch is an admissible source state through recombination, the acoustic peaks are computed by ordinary Einstein-Boltzmann transfer from those finite inputs. ◻

</div>

<div class="proposition">

**Proposition 47** (Curvature suppression is not sector selection). *On a fixed clocked FLRW sector,
``` math
K(\tau)=\kappa/a(\tau)^2,\qquad \dot K=-2HK .
```
Any repair or expansion estimate that proves a finite bound on $`|K|`$ or $`|\Omega_K|`$ is an <span class="smallcaps">ApproximateBound</span>. It does not change $`\kappa`$, does not distinguish globally inequivalent flat topologies, and may not be reported as <span class="smallcaps">DirectTheorem</span> or <span class="smallcaps">ConditionalCMH</span>.*

</div>

<div class="remark">

*Remark 48* (Theorem boundary). The theorem is conditional at the CMB-prediction tier unless the finite source package, physical scale bridge, stress parent, Boltzmann transfer, and frozen likelihood receipts all pass for the same generation. Imported frozen geometry can support the `CONDITIONAL_PHYSICAL` route. `OPH_NATIVE_PHYSICAL` also requires the quotient-derived $`\mathsf{CosmoGeomRead}_r`$ and source embedding theorem.

</div>

# Physical CMB Finite-Source Contract

The physical CMB source contract requires theorem-grade finite artifacts:
``` math
A_\zeta,\qquad q_{\rm IR},\qquad \ell_{\rm IR},\qquad
\bar\rho_A(a),\qquad \bar\rho_{A,\mathrm{eq}}(a),\qquad B_A(k,a),\qquad
\gamma_{\mathrm{repair\ step}}(k,a)\ \hbox{or certified}\ \Gamma_{\rm rec}(k,a),
```
plus full stress variables, freezeout, physical mode calibration, and $`N_{\rm CRC}`$ provenance.

<div class="definition">

**Definition 49** (Finite primordial source package). *At regulator $`r`$, a physical primordial-source claim is carried by
``` math
\mathfrak S_r^{\rm prim}=
(\mathcal Q_r,n_r,\mu_r,q_r,J_r,\Pi_r,K_r,M_r,\Xi_r,\Sigma_r,E_r^{\rm rel},
A_{q,r},\theta_r,q_{{\rm IR},r},t_{{\rm IR},r},\Sigma_{\star,r},W_r,
\Delta_{\zeta,r}^2,\mathcal P_{IJ,r},\mathcal E_r).
```
Here $`\mathcal Q_r`$ is the physical quotient state space, $`n_r`$ the normal-form map, $`\mu_r`$ the source probability law, $`q_r`$ the quotient-visible scalar source, $`J_r`$ the finite carrier-to-harmonic readout, $`\Pi_r`$ the gauge/monopole/dipole projector, $`K_r`$ the positive scalar precision operator, $`M_r`$ the mobility or repair-response operator, $`\Xi_r`$ the stochastic forcing covariance, $`\Sigma_r=\operatorname{Cov}_{\mu_r}(q_r)`$, $`E_r^{\rm rel}`$ the source release-energy artifact, $`A_{q,r}`$ the finite screen amplitude, $`\theta_r`$ the spatial/RG spectral exponent, $`q_{{\rm IR},r}`$ the IR suppression depth, $`t_{{\rm IR},r}`$ the IR semigroup time, $`\Sigma_{\star,r}`$ the freezeout cut, $`W_r`$ the radial window, $`\Delta_{\zeta,r}^2`$ the lifted primordial spectrum, $`\mathcal P_{IJ,r}`$ the species-level initial covariance, and $`\mathcal E_r`$ the complete error and refinement envelope.*

</div>

<div class="target">

**Target 50** (Finite-source input contract). *The paper stack must define pass/fail receipts for:*

1.  *a source-provenance dependency DAG for all CMB source inputs;*

2.  *globally pooled sufficient-statistic reducers before nonlinear source estimates;*

3.  *the role of $`N_{\rm CRC}`$ as a consensus invariant, additive capacity, or theorem-side constant;*

4.  *finite derivations of $`A_\zeta`$, $`q_{\rm IR}`$, and $`\ell_{\rm IR}`$;*

5.  *total stress closure, single-clock normal form, entropy repair gap, curvature evolution, adiabatic growing mode, isocurvature bound, phase coherence, screen-to-radial lift, radial null-space, forward residual, and finite freezeout receipts;*

6.  *$`\texttt{physical\_k\_units\_calibrated}`$, derived exactly from $`\texttt{PHYSICAL\_SPATIAL\_K\_RECEIPT}`$;*

7.  *$`\texttt{screen\_to\_physical\_k\_association\_calibrated}`$, derived from $`\texttt{SCREEN\_TO\_PHYSICAL\_K\_ASSOCIATION\_RECEIPT}`$;*

8.  *$`\texttt{source\_angular\_sector\_calibrated}`$, derived from $`\texttt{SOURCE\_ANGULAR\_SECTOR\_RECEIPT}`$, with $`\ell_{\rm src}`$ kept distinct from $`L_{\rm CMB}`$;*

9.  *$`\texttt{calibrated\_a\_evolution}`$, derived from $`\texttt{CALIBRATED\_A\_EVOLUTION\_RECEIPT}`$;*

10. *$`\texttt{physical\_mode\_freezeout\_map\_calibrated}`$, derived from $`\texttt{PHYSICAL\_MODE\_FREEZEOUT\_MAP\_RECEIPT}`$;*

11. *$`\texttt{common\_primordial\_initial\_surface\_calibrated}`$, derived from $`\texttt{PHYSICAL\_FREEZEOUT\_SURFACE\_RECEIPT}`$;*

12. *a nonempty cross-receipt manifest with matching regulator-family, generation, geometry, background, clock, scale-certificate, source-embedding, mode-basis, mode-lineage, boundary, solver-convention, freezeout-surface, and source-DAG hashes;*

13. *source-only imported dark-sector functions $`\bar\rho_A(a)`$, $`\bar\rho_{A,\mathrm{eq}}(a)`$, $`w_A(a)`$, $`c_{s,A}^2(k,a)`$, $`\sigma_A(k,a)`$, $`Q_A^\mu`$, $`B_A(k,a)`$, and the $`\gamma_{\mathrm{repair\ step}}`$ diagnostic or certified $`\Gamma_{\rm rec}(k,a)`$, where the chosen dark-continuation mode permits them.*

*A mode-dependent freezeout map does not satisfy $`\texttt{PHYSICAL\_FREEZEOUT\_SURFACE\_RECEIPT}`$; the receipt denotes a common spacelike initial surface with initial data and normal derivatives in one frozen generation.*

</div>

<div class="target">

**Target 51** (Scalar, harmonic, and IR source receipts). *The scalar-source branch must supply
``` math
\texttt{SCALAR\_SOURCE\_MAP\_RECEIPT},\quad
\texttt{A5\_TO\_SO3\_INTERTWINER\_RECEIPT},\quad
\texttt{SPHERICAL\_QUADRATURE\_RECEIPT},\quad
\texttt{LOW\_MODE\_PROJECTOR\_RECEIPT}.
```
The map $`\mathsf{ScalarSourceMap}_r:\operatorname{nf}(\mathcal Q_r)\to
V_r^{\rm scalar}`$ must be quotient-invariant, refinement-natural, and normalized without per-run standardization. Harmonic reconstruction uses a weighted Gram matrix $`G_L=Y_L^\dagger W Y_L`$ and the low-mode projector
``` math
\Pi_{\ge2}=I-Y_{01}(Y_{01}^\dagger W Y_{01})^{-1}Y_{01}^\dagger W.
```
The IR kernel is parameterized by the semigroup time
``` math
S_{\rm IR}=I-q_{\rm IR}e^{-t_{\rm IR}L},\qquad
F_{\rm IR}(\ell)=1-q_{\rm IR}e^{-t_{\rm IR}\ell(\ell+1)}.
```
The validator must enforce $`0\le q_{\rm IR}\le1`$, $`t_{\rm IR}\ge0`$, and derive the display scale
``` math
\ell_{\rm IR}=\frac{-1+\sqrt{1+4/t_{\rm IR}}}{2}
```
when $`t_{\rm IR}>0`$. The old dimension count $`12+20+1=33\Rightarrow
\ell_{\rm IR}=32`$ is not a source theorem.*

</div>

<div class="target">

**Target 52** (Fluctuation covariance and tilt receipts). *Deterministic settling does not choose a fluctuation ensemble. The finite quotient ensemble must supply $`\mu_r(q)=m_r(q)e^{-S_r(q)}/Z_r`$ with projective compatibility, sampler correctness, and implementation invariance. Repair drift does not determine covariance: $`D\Sigma+\Sigma D^{\mathsf T}=\Xi`$ requires the noise covariance $`\Xi`$. A fluctuation-dissipation branch may use $`D=M_rK_r`$, $`\Xi=2A_rM_r`$, and $`\Sigma_r=A_rK_r^{-1}`$ when it also supplies
``` math
\texttt{QUOTIENT\_ENSEMBLE\_RECEIPT},\quad
\texttt{LYAPUNOV\_COVARIANCE\_RECEIPT},\quad
\texttt{FLUCTUATION\_DISSIPATION\_RECEIPT},\quad
\texttt{MODE\_RESOLVED\_PRECISION\_OPERATOR\_RECEIPT}.
```
The temporal repair gap, spatial RG exponent, primordial tilt, and physical exchange-rate kernel are distinct types:
``` math
\texttt{TEMPORAL\_REPAIR\_GAP},\quad
\texttt{SPATIAL\_RG\_EIGENEXPONENT},\quad
\texttt{PRIMORDIAL\_SPECTRAL\_TILT},\quad
\texttt{PHYSICAL\_EXCHANGE\_RATE\_KERNEL}.
```
Only after the radial theorem identifies $`\Delta_\zeta^2(k)\propto k^{-\theta}`$ may the branch write $`n_s=1-\theta`$.*

</div>

<div class="proposition">

**Proposition 53** (Why the 256k CMB curve is diagnostic). *The 256k OPH-FPE run has measurement-facing diagnostic status only. Its promotion report marks the finite-source input contract false and lists missing source provenance, missing pooled reducers, missing $`N_{\rm CRC}`$ invariant role, non-finite $`A_\zeta`$, missing screen-to-primordial lift, non-finite $`q_{\rm IR}`$ and $`\ell_{\rm IR}`$, non-finite $`\bar\rho_A`$ and $`B_A`$, missing finite covariant parent receipt, missing stress closure, missing full fluid/exchange variables, missing gauge/causal/refinement certificates, missing freezeout, missing official likelihood and CDM-limit reductions, and missing frozen solver and likelihood hashes.*

</div>

# Physical Scale Bridge

Section <a href="#sec:physical-scale-bridge" data-reference-type="ref" data-reference="sec:physical-scale-bridge">4</a> owns the bridge definitions, no-go propositions, and conditional physical-scale theorem. The practical rule is:
``` math
\texttt{inverse\_cap\_opening\_angle\_proxy},\quad
\texttt{repair\_cycle},\quad
\texttt{graph\_mode\_index}
```
are nonpromotable diagnostic labels unless the structured receipt conjunctions in Target <a href="#target:physical-scale-receipts" data-reference-type="ref" data-reference="target:physical-scale-receipts">34</a> pass. A frozen imported FLRW geometry packet can exercise the complete spectral/clock/freezeout/transfer pipeline at tier `CONDITIONAL_PHYSICAL`. The `OPH_NATIVE_PHYSICAL` tier additionally requires a quotient-derived $`\mathsf{CosmoGeomRead}_r`$ and source embedding theorem with the same downstream interfaces.

# Dark/Anomaly Interface Imported from the Dark-Matter Paper

The dark/anomaly slot is the largest cosmology-facing bridge between the OPH dark-matter paper and the CMB program. The released dark-matter paper `cosmology/oph_dark_matter_paper.tex` owns the derivation: collar remainder stress, settled galaxy response, dynamic transport, cluster behavior, and perturbation kernels. This paper only states the interface the CMB program requires from that paper. Cosmology needs the imported dark sector to expose a transported stress component with background density, perturbations, exchange terms, gauge-invariant variables, and the CDM-limit switch-off check.

<div class="definition">

**Definition 54** (Imported finite covariant collar-packet parent interface). *For the purpose of this cosmology program, the finite covariant collar-packet parent is an imported `cosmology/oph_dark_matter_paper.tex`-owned finite artifact interface
``` math
\mathcal P_A=
(\mathcal C_r,Z_r,A_r,R_r,G_r,\pi_r,L_r,Q_r,D_r)
```
where $`Z_r`$ splits anomaly and recipient packet states, $`\pi_r`$ is the finite equilibrium packet functional, $`L_r`$ is the repair generator, $`Q_r`$ records energy-momentum reaction channels, and $`D_r`$ records any causal auxiliary response. The interface emits packet stress moments
``` math
\rho_A,\quad P_A,\quad q_A,\quad \pi_A,
```
and source variables
``` math
\bar\rho_A(a),\quad \bar\rho_{A,\mathrm{eq}}(a),\quad w_A(a),\quad c_{s,A}^2(k,a),\quad
\sigma_A(k,a),\quad Q_A^\mu,\quad B_A(k,a),\quad
\gamma_{\mathrm{repair\ step}}(k,a)\ \hbox{or certified}\ \Gamma_{\rm rec}(k,a),
```
with explicit conservation or exchange equations, recipient stress for nonzero exchange, gauge projection, causal response, regulator-refinement convergence, and $`B_A`$ built from the anomaly-frame baryon density $`n_b^{(A)}=-u_{A\mu}J_b^\mu`$.*

</div>

<div class="definition">

**Definition 55** (Issue \#319 conditional source receipt). *The imported dark-sector parent first has to pass
``` math
\texttt{ISSUE\_319\_CONDITIONAL\_SOURCE\_RECEIPT}.
```
This receipt is a conjunction of source-side receipts, not a scalar CMI row. It requires exactly one declared route,
``` math
\texttt{SOURCE\_ROUTE\_RECEIPT},
```
unit metadata
``` math
\texttt{ENTROPY\_UNIT\_RECEIPT},
```
the exact carried modular identity
``` math
\langle K_{AB}+K_{BD}-K_B-K_{ABD}\rangle=I(A:D|B),
```
and the first-variation classification that prevents raw CMI from being silently identified with a linear fixed-reference modular-energy variation. The source route is either `FIXED_REFERENCE_MODULAR_ENERGY` or `NONLINEAR_CMI_STRESS`. The second route also requires
``` math
\texttt{CMI\_TO\_MODULAR\_SOURCE\_MATCHING\_RECEIPT}.
```
Both routes require collar-localization, BW normalization, physical diamond scale, $`\ell^4`$ scaling, cover-independence, and stress-tomography receipts before the proper-radius normalization
``` math
T_A^{ab}u_a u_b
  =
  \frac{15\hbar c}{8\pi^2\ell^4}R^{\rm src}_{r,\ell}
```
may be imported. Classical diagonal CMI rows, cap-angle proxies, duplicated collar counts, or absolute-value rectified density rows do not satisfy this receipt by themselves.*

</div>

<div class="target">

**Target 56** (Imported stress-parent receipt). *The simulator receipt
``` math
\texttt{FINITE\_COVARIANT\_PARENT\_RECEIPT}
```
must fail unless the owning dark-sector paper and its evidence bundle certify all of the following:*

1.  *$`\texttt{ISSUE\_319\_CONDITIONAL\_SOURCE\_RECEIPT}`$;*

2.  *$`\texttt{CMI\_TO\_MODULAR\_SOURCE\_MATCHING\_RECEIPT}`$ on the nonlinear-CMI route, or $`\texttt{FIXED\_REFERENCE\_MODULAR\_ENERGY\_RECEIPT}`$ on the fixed-reference route;*

3.  *$`\texttt{STRESS\_TOMOGRAPHY\_RECEIPT}`$ and $`\texttt{SM\_CURRENT\_NULL\_RECEIPT}`$;*

4.  *exact finite stress-energy closure for anomaly plus recipient stresses;*

5.  *explicit recipient stress and equal-and-opposite exchange current whenever the exchange branch is nonzero;*

6.  *detailed balance or a declared nonequilibrium exchange law;*

7.  *gauge independence of the emitted perturbation variables, including the anomaly-frame baryon-density definition of $`B_A`$;*

8.  *causal response with subluminal finite characteristics;*

9.  *convergence under regulator refinement;*

10. *recovery of the CDM branch when exchange, pressure, sound speed, and anisotropic stress are switched off.*

*The finite parent receipt is source-side. It does not include official CMB likelihood execution. The later $`\texttt{PHYSICAL\_CMB\_PROMOTION\_RECEIPT}`$ additionally requires the physical scale bridge, physical $`\rho_A`$, $`B_A`$, and $`\Gamma_{\rm rec}`$ kernels, Boltzmann transfer, frozen source/solver/likelihood hashes, and the official likelihood run.*

</div>

<div class="remark">

*Remark 57* (Dark continuation modes). For a source-only primordial certification run, the dark-sector paper recommends `dark_continuation = OFF`. A supplied dark abundance may enter only as `dark_continuation = CONDITIONAL_SOURCE_STATE`. A source-selected abundance may enter as `SOURCE_ONLY_ANOMALY_ABUNDANCE` only when `ANOMALY_ABUNDANCE_SOURCE_RECEIPT` passes. A parent with physical transport and no selector is labelled `PHYSICAL_PARENT_WITH_CONDITIONAL_ABUNDANCE`. The physical cosmology branch may set `dark_continuation = PHYSICAL_PARENT` only after the finite covariant parent receipt, physical scale bridge, response/kernel receipts, and CDM-limit recovery all pass for the same frozen generation.

</div>

# Physical Kernel Interface

The physical kernels are not scalar rows in a table. They are imported finite functionals with units, gauge-invariant definitions, paired controls, and refinement margins. Their derivation stays in `cosmology/oph_dark_matter_paper.tex`; this paper only states what the CMB/LSS pipeline must receive.

<div class="target">

**Target 58** (Physical anomaly kernel interface). *The released dark-matter paper must define deterministic receipts*

<div class="center">

|                              |
|:-----------------------------|
| *`RHO_A_SOURCE_RECEIPT`*     |
| *`RHO_A_EQ_SOURCE_RECEIPT`*  |
| *`B_A_KERNEL_RECEIPT`*       |
| *`GAMMA_REC_SOURCE_RECEIPT`* |

</div>

*These receipts pass only when $`\bar\rho_A(a)`$, $`B_A(k,a)`$, and $`\Gamma_{\rm rec}(k,a)`$ are emitted as certified physical functions with $`\bar\rho_{A,\mathrm{eq}}(a)`$, $`w_A(a)`$, $`c_{s,A}^2(k,a)`$, $`\sigma_A(k,a)`$, and $`Q_A^\mu`$ from the finite parent, physical scale calibration, gauge consistency, paired perturbation controls, active-fiber/physical-clock response receipts, CDM-limit recovery, and refinement convergence. The cosmology program consumes the emitted functions but does not redefine the dark-sector source law.*

The imported `RHO_A_SOURCE_RECEIPT` is the conjunction of `RHO_A_TRANSPORT_RECEIPT` and `ANOMALY_ABUNDANCE_SOURCE_RECEIPT`. The second receipt is supplied by the source-only anomaly abundance selector in `cosmology/oph_dark_matter_paper.tex`.

*The physical kernel receipt imports these inner receipts from the dark/anomaly paper:*

<div class="center">

|                                                 |
|:------------------------------------------------|
| *`COMMON_SOURCE_FUNCTIONAL_RECEIPT`*            |
| *`ADMISSIBLE_SOURCE_TANGENT_RECEIPT`*           |
| *`CONSTRAINT_PRESERVING_RETRACTION_RECEIPT`*    |
| *`SOURCE_VECTOR_SUFFICIENCY_RECEIPT`*           |
| *`B_A_SOURCE_LIFT_INDEPENDENCE_RECEIPT`*        |
| *`SOURCE_DESIGN_IDENTIFIABILITY_RECEIPT`*       |
| *`FINITE_DIFFERENCE_ORDER_RECEIPT`*             |
| *`C1_REFINEMENT_RECEIPT`*                       |
| *`ORDER_OF_LIMITS_RECEIPT`*                     |
| *`NATIVE_REPAIR_GENERATOR_RECEIPT`*             |
| *`PHYSICAL_RATE_UNIT_RECEIPT`*                  |
| *`QUOTIENT_LUMPABILITY_RECEIPT`*                |
| *`STATIC_DYNAMIC_RESPONSE_CONSISTENCY_RECEIPT`* |

</div>

*A diagnostic candidate may contain paired finite differences, controls, and refinement rows. The physical receipt is false unless the source functional, delivered source tangent, source-vector sufficiency, lift-independence, native generator, and static-dynamic consistency receipts pass in the same frozen generation.*

</div>

Diagnostic rows remain useful for debugging. They should be plotted, stress-tested, and compared to controls. They should not be accepted by the physical CMB contract unless they are promoted by the kernel receipts above. Until the physical clock, active fiber, conserved-sector decomposition, and common-parent response pole are certified, the simulator should emit the transition number as $`\gamma_{\mathrm{repair\ step}}`$, not as $`\Gamma_{\rm rec}`$. Likewise finite rows for $`\rho_A`$, $`\rho_{A,\mathrm{eq}}`$, and $`B_A`$ have label `DIAGNOSTIC_SCALAR_RESPONSE` until promoted by the parent and kernel receipts.

<div class="target">

**Target 59** (Finite covariant parent import). *The imported dark-sector parent is not a readiness boolean. The receipt consumed by this paper is the conjunction of primitive parent receipts owned by `cosmology/oph_dark_matter_paper.tex`:
``` math
\begin{gathered}
\texttt{SOURCE\_ROUTE\_RECEIPT},\quad
\texttt{ENTROPY\_UNIT\_RECEIPT},\quad
\texttt{MODULAR\_NONADDITIVITY\_IDENTITY\_RECEIPT},\\
\texttt{CMI\_FIRST\_VARIATION\_CLASSIFICATION\_RECEIPT},\quad
\texttt{SOURCE\_LOCALIZATION\_SATURATION\_RECEIPT},\\
\texttt{BW\_BALL\_NORMALIZATION\_RECEIPT},\\
\texttt{PHYSICAL\_DIAMOND\_SCALE\_RECEIPT},\\
\texttt{ELL4\_SCALING\_PLATEAU\_RECEIPT},\\
\texttt{COVER\_INDEPENDENCE\_RECEIPT},\quad
\texttt{STRESS\_TOMOGRAPHY\_RECEIPT},\\
\texttt{FINITE\_PACKET\_KINEMATICS\_RECEIPT},\\
\texttt{PACKET\_MASS\_SHELL\_RECEIPT},\\
\texttt{TRANSPORT\_COVARIANCE\_RECEIPT},\\
\texttt{CHANNEL\_FOUR\_MOMENTUM\_RECEIPT},\\
\texttt{FINITE\_PACKET\_STRESS\_READOUT\_RECEIPT},\\
\texttt{VARIATIONAL\_MOMENT\_STRESS\_AGREEMENT\_RECEIPT},\\
\texttt{LOCAL\_FRAME\_COVARIANCE\_RECEIPT},\\
\texttt{CARRIER\_QUOTIENT\_INVARIANCE\_RECEIPT},\\
\texttt{TOTAL\_STRESS\_CLOSURE\_RECEIPT},\\
\texttt{EXPLICIT\_RECIPIENT\_STRESS\_RECEIPT},\\
\texttt{EXCHANGE\_CURRENT\_CLOSURE\_RECEIPT},\\
\texttt{SM\_CURRENT\_NULL\_RECEIPT},\\
\texttt{COSMOLOGICAL\_GAUGE\_INVARIANCE\_RECEIPT},\\
\texttt{FINITE\_DOMAIN\_OF\_DEPENDENCE\_RECEIPT},\\
\texttt{SUBLUMINAL\_CHARACTERISTICS\_RECEIPT},\\
\texttt{RETARDED\_RESPONSE\_RECEIPT},\quad
\texttt{RESPONSE\_STABILITY\_RECEIPT},\\
\texttt{REFINEMENT\_CONVERGENCE\_RECEIPT},\quad
\texttt{CDM\_LIMIT\_RECOVERY\_RECEIPT}.
\end{gathered}
```
The source parent does not contain solver or likelihood hashes. Those belong to the frozen generation and likelihood layer below.*

</div>

# Boltzmann Transfer and Frozen Likelihood

For finite physical source functions, the transfer problem is ordinary cosmology. OPH does not need a special CMB plotting rule. It needs a frozen handoff to CAMB, CLASS, or a declared independent Einstein-Boltzmann solver.

<div class="target">

**Target 60** (Frozen transfer and likelihood protocol). *A physical OPH CMB prediction requires:*

1.  *immutable hashes for source artifacts;*

2.  *immutable hashes for solver source, version, tolerances, and input files;*

3.  *immutable hashes for likelihood code, datasets, covariances, masks, and nuisance priors;*

4.  *a no-data-use manifest showing that observational likelihood values did not enter source artifact generation;*

5.  *global pooled reducers for nonlinear source estimates and global CDM-limit checks;*

6.  *official likelihood execution readiness beyond a diagnostic binned-table comparison;*

7.  *a falsification rule fixed before the likelihood run.*

</div>

Passing `ANOMALY_ABUNDANCE_SOURCE_RECEIPT` promotes the dark abundance inside the source species. It does not promote CMB spectra. The later `PHYSICAL_CMB_PREDICTION_RECEIPT` requires frozen transfer and official likelihood execution.

<div id="def:model-generation" class="definition">

**Definition 61** (Model generation). *A cosmology prediction generation is
``` math
G_n=(H_{\rm source},H_{\rm scale},H_{\rm parent},H_{\rm kernel},H_{\rm init},
H_{\rm solver},H_{\rm build},H_{\rm config},H_{\rm likelihood},H_{\rm data},
H_{\rm priors},\mathcal A_n),
```
where $`\mathcal A_n`$ is the registry of observational datasets, residuals, posterior summaries, and diagnostic plots exposed during source development. Hashing the final artifact is not enough: if a Planck residual influenced source or model selection, the generation is calibrated on Planck even when the final bytes are frozen.*

</div>

<div id="thm:physical-oph-cosmology-promotion" class="theorem">

**Theorem 62** (Physical OPH cosmology promotion). *Let $`R_{\rm src}`$ be source-only provenance and pooled-reducer receipt, $`R_{\rm scale}`$ the physical mode/clock/lift receipt, $`R_{\rm parent}`$ the finite covariant parent receipt, $`R_{\rm kernel}`$ the physical response/kernel receipt, $`R_{\rm init}`$ the regular initial-mode and Einstein-constraint receipt, $`R_{\rm transfer}`$ the Boltzmann well-posedness and numerical-convergence receipt, $`R_{\rm freeze}`$ the model-generation freeze receipt, and $`R_{\rm like}`$ the official likelihood execution receipt. Define
``` math
R_{\rm spectrum}=
R_{\rm src}\wedge R_{\rm scale}\wedge R_{\rm parent}\wedge R_{\rm kernel}
\wedge R_{\rm init}\wedge R_{\rm transfer}\wedge R_{\rm freeze},
```
and for an observational dataset family $`D`$,
``` math
R_{\rm tested}^{(D)}=R_{\rm spectrum}\wedge R_{\rm like}^{(D)}.
```
If every receipt in $`R_{\rm spectrum}`$ is labelled $`\texttt{PASS}`$, then the resulting TT/TE/EE, lensing, matter power, growth, BAO, or RSD spectra are $`\texttt{FROZEN\_PHYSICAL\_SPECTRUM}`$ outputs of generation $`G_n`$. If $`R_{\rm tested}^{(D)}`$ also passes for a frozen dataset and likelihood family $`D`$, the output is a $`\texttt{LIKELIHOOD\_EVALUATED\_PHYSICAL\_PREDICTION}`$. If any receipt is $`\texttt{FAIL}`$, $`\texttt{OPEN\_GATE}`$, $`\texttt{INCOMPLETE}`$, or $`\texttt{INVALIDATED}`$, the spectra remain diagnostic or conditional according to the parameter ledger.*

</div>

<div class="proof">

*Proof.* Each receipt states a typed premise needed by the next map: source fields, physical coordinates, parent stress, response kernels, initial state, transfer evolution, freeze integrity, and likelihood execution. Composition through $`R_{\rm spectrum}`$ gives a frozen physical spectrum. Composing that frozen spectrum with $`R_{\rm like}^{(D)}`$ evaluates it against $`D`$. A failed or open receipt removes at least one premise, so downstream software cannot promote the output. ◻

</div>

<div class="target">

**Target 63** (Outcome and prediction-class ledger). *Likelihood outcomes are
``` math
\begin{gathered}
\texttt{INVALIDATED},\quad
\texttt{NUMERICALLY\_INCONCLUSIVE},\quad
\texttt{FALSIFIED\_ON\_DECLARED\_TEST},\\
\texttt{PASSED\_DECLARED\_TEST}.
\end{gathered}
```
`INVALIDATED` is used when provenance or freeze integrity fails; it is not a physical falsification. Inputs are typed as
``` math
\begin{gathered}
\texttt{OPH\_DERIVED},\quad
\texttt{EXTERNAL\_FIXED},\quad
\texttt{CONDITIONAL\_COSMOLOGY},\quad
\texttt{NUISANCE},\\
\texttt{MEASUREMENT\_CALIBRATED},
\end{gathered}
```
and outputs as
``` math
\begin{gathered}
\texttt{SOURCE\_ONLY\_OPH\_PREDICTION},\quad
\texttt{CONDITIONAL\_OPH\_PREDICTION},\\
\texttt{FITTED\_OPH\_MODEL},\quad
\texttt{DIAGNOSTIC\_BASELINE}.
\end{gathered}
```
The 256k curve has label $`\texttt{CMB0\_LCDM\_DIAGNOSTIC\_BASELINE}`$.*

</div>

<div class="remark">

*Remark 64* (Global reducers only). Shard-local nonlinear averages and shard-local `any()` rollups are not promotion rules. The source artifact must pool additive sufficient statistics globally, with validated units, coordinate grids, coverage, duplicate policy, interpolation policy, and covariance, before any nonlinear quantity such as amplitude, rank, condition number, isocurvature leakage, $`B_A`$, $`\bar\rho_A`$, or $`\Gamma_{\rm rec}`$ is evaluated; a transition diagnostic is recorded separately as $`\gamma_{\mathrm{repair\ step}}`$.

</div>

<div class="remark">

*Remark 65* (Comparison versus prediction). The best OPH diagnostic model in the 256k run is scientifically useful because it is close: shape correlation $`0.9951542364`$, normalized RMSE $`0.0984455548`$, amplitude-fit $`\chi^2/{\rm bin}=41.6684770597`$, zero mean peak-location error, and mean peak-height fractional error about $`0.0704905`$. Those numbers justify fixing the theorem gates. They do not replace the gates.

</div>

# Vacuum and Quantum-Fluctuation Boundary

The same finite-source discipline applies to vacuum claims. Seed noise, repair jitter, and reference free-field baselines are not an OPH-native quantum vacuum.

<div class="target">

**Target 66** (OPH-native vacuum gate). *An OPH-native vacuum simulation needs a quotient ensemble measure or density operator on the observable quotient algebra, refinement compatibility, implementation invariance, a sampler correctness proof, and source Euclidean slab data $`(Q_r,m_r^0,J_r,V_r,a_{t,r})`$ whose transfer operator is derived independently of the target law and sampler. Without this gate, vacuum-like fields in the simulator are reference baselines or diagnostics.*

</div>

# Claim-Gate Summary

The paper-side gates are quotient ensemble measure, OPH-native vacuum promotion, Boltzmann transfer, frozen likelihood closure, finite-source CMB contract hardening, physical scale bridge, finite covariant collar-packet stress parent, physical anomaly kernels, flat-sector selection, screen-spectrum derivation, conditional collar source routing, repair-stress transport, and linear repair-stress perturbations. GitHub issues carry the mutable ownership, dependency, and completion state for those gates.

# Physical Promotion Rule

Physical promotion requires exact receipt semantics for the objects the simulator is being asked to produce:

1.  implement the finite-source dependency DAG and no-data-use firewall;

2.  implement global pooled reducers for all nonlinear source estimates;

3.  implement the `CONDITIONAL_PHYSICAL` imported-FLRW reference backend for the spectral, clock, freezeout, and transfer pipeline;

4.  implement the structured physical scale bridge and freezeout certificates from Target <a href="#target:physical-scale-receipts" data-reference-type="ref" data-reference="target:physical-scale-receipts">34</a>;

5.  implement the imported finite covariant parent and stress-closure receipts from the `cosmology/oph_dark_matter_paper.tex` dark-sector theory;

6.  implement physical $`\bar\rho_A`$, $`B_A`$, and $`\Gamma_{\rm rec}`$ kernel receipts as imported dark-sector source functions;

7.  freeze source, solver, and likelihood hashes;

8.  run official likelihood comparisons and global CDM-limit reductions;

9.  label the result `CONDITIONAL_PHYSICAL` when all imported-geometry physical receipts pass, or `OPH_NATIVE_PHYSICAL` when Target <a href="#target:cosmo-geometry-extraction" data-reference-type="ref" data-reference="target:cosmo-geometry-extraction">35</a> also passes.

The simulator may run diagnostics before physical promotion: visualization payloads, CMB overlays, vacuum baselines, defect worldlines, observer cameras, and H3 readouts. Those outputs carry diagnostic labels until the relevant receipt passes.

# Conclusion

The 256k simulator almost-fits are not random noise. They expose stable structure close enough to justify the theorem gates, and they show why those gates are necessary. A physical cosmology prediction needs finite source-only objects, physical scale calibration, covariant stress closure, physical anomaly kernels, and frozen likelihoods. This paper is the staging surface for that work and the simulator contract for physical promotion.

<div class="thebibliography">

99

B. Müller, A. Osika, M. Poneder, K. Xue, P. Nguyen, M. A. Visser, and D. Matscheko, *Recovering Relativity and the Standard Model from Observer Overlap Consistency*, OPH compact paper, GitHub main PDF. <https://github.com/FloatingPragma/observer-patch-holography/blob/main/paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf>

B. Müller, A. Osika, M. Poneder, K. Xue, B. Cassie, P. Nguyen, M. A. Visser, K. A. Anirudha, D. Matscheko, and J. Hill, *Observers Are All You Need*, OPH synthesis paper, GitHub main PDF. <https://github.com/FloatingPragma/observer-patch-holography/blob/main/paper/observers_are_all_you_need.pdf>

B. Müller, K. Xue, K. A. Anirudha, D. Matscheko, and J. Hill, *Reality as a Consensus Protocol*, OPH consensus paper, GitHub main PDF. <https://github.com/FloatingPragma/observer-patch-holography/blob/main/paper/reality_as_consensus_protocol.pdf>

B. Müller, A. Osika, K. Xue, B. Cassie, M. A. Visser, and D. Matscheko, *Federated Echosahedral Screen Microphysics: Patch Hardware, Records, and Observer Synchronization in OPH*, OPH microphysics paper, GitHub main PDF. <https://github.com/FloatingPragma/observer-patch-holography/blob/main/paper/screen_microphysics_and_observer_synchronization.pdf>

B. Müller and D. Matscheko, *Observer-Patch Holography and the Dark Matter Phenomenon*, OPH cosmology paper, 2026. <https://github.com/FloatingPragma/observer-patch-holography/blob/main/cosmology/oph_dark_matter_paper.pdf>

B. Müller, *Inflation Without an Inflaton: Observer-Screen Synchronization as an OPH Cosmology Branch*, workspace cosmology draft, 2026.

B. Müller, *OPH-FPE finite screen-consensus and cosmology diagnostics*. <https://github.com/muellerberndt/oph-physics-sim>

Planck Collaboration, *Planck 2018 results. VI. Cosmological parameters*, Astron. Astrophys. 641, A6, 2020, arXiv:1807.06209. <https://arxiv.org/abs/1807.06209>

A. Lewis, A. Challinor, and A. Lasenby, *Efficient computation of CMB anisotropies in closed FRW models*, Astrophys. J. 538, 473–476, 2000, arXiv:astro-ph/9911177. <https://arxiv.org/abs/astro-ph/9911177>

D. Blas, J. Lesgourgues, and T. Tram, *The Cosmic Linear Anisotropy Solving System (CLASS). Part II: Approximation schemes*, JCAP 07, 034, 2011, arXiv:1104.2933. <https://arxiv.org/abs/1104.2933>

</div>
