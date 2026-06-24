# Status

This is a staged template, not a release-pipeline paper. It coordinates the vacuum/fluctuation part of `cosmology/oph_cosmology_finite_source_cmb_program.tex`, the screen-microphysics finite carrier, the compact theorem stack, and OPH-FPE simulator outputs in <https://github.com/muellerberndt/oph-physics-sim>.

# Shared Finite Quotient Theorem Surface

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
The finite constraint ledger must name every $`F_{r,a}`$, its units and support, the target expectation and source, sector or zero-mode treatment, refinement transformation, and proof that no simulator or observational output entered the source definition.

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

If a simulator stores representatives, a representative-level law must be a conditional lift
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

#### Simulator accuracy.

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
Continuum-facing observables require a realization map and correlation Cauchy bound, not just a finite histogram.

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
The simulator must keep separate receipts for stationary-law schedule invariance, detailed balance of the aggregate kernel, and pathwise partition invariance. Deterministic replay of semantic random streams or a canonical serial chain is useful, but it is not pathwise partition invariance. Smoothing must preserve raw coefficients, raw spectra, smoothing kernels, smoothed coefficients, smoothed spectra, and hashes of each stage; it is not part of $`S_r`$ unless explicitly declared.

# Vacuum Claim Boundary

The strongest allowed simulator labels are:
``` math
E0:\hbox{ seed noise or repair jitter},\qquad
E1:\hbox{ conventional reference ensemble}.
```
An OPH-native vacuum label requires an E3 receipt. The source data must include
``` math
\mathfrak S_r^E=(Q_r,m_r^0,J_r,V_r,a_{t,r})
```
with $`J_r`$, $`V_r`$, and $`a_{t,r}`$ derived from OPH source structure rather than from a target law or a desired visual spectrum. Accepted repair moves are generally directed and cannot simply be reused as reversible Euclidean conductances.

# Structure-Formation Targets

| Object | Required source-side label |
|:---|:---|
| Fluctuating screen field | E0 or E1 unless an E2 quotient ensemble receipt passes. |
| OPH vacuum field | E3 only after source Euclidean slab, ground-state, Doob, reflection-positivity, and transfer-refinement receipts. |
| Proto-object worldlines | Diagnostic until production matter/particle and observer-facing object population receipts pass. |
| Growth seeds | E4 only after the primordial bridge, radial lift, null-space, positivity, and forward-projection receipts pass. |
| LSS comparison | E5 only after frozen likelihood/data contracts pass in the data-contract paper. |

# Simulator Contract

The simulator should emit the following receipt names even when they fail:

- `TRACIALLY_POINTED_QUOTIENT_RECEIPT`

- `SOURCE_EUCLIDEAN_SLAB_RECEIPT`

- `TRANSFER_REFINEMENT_RECEIPT`

- `STATIONARY_LAW_SCHEDULE_INVARIANCE_RECEIPT`

- `PATHWISE_PARTITION_INVARIANCE_RECEIPT`

False receipts are part of the scientific result. They prevent a reference baseline or rendering field from being mislabeled as an OPH vacuum.

# Open Work

1.  Derive the reversible source conductances $`J_r`$, potential $`V_r`$, and slab thickness $`a_{t,r}`$ from finite OPH source data.

2.  Prove reflection positivity or an equivalent reconstruction theorem on the selected branch.

3.  Relate passed E3 vacuum receipts to E4 primordial growth seeds without using CMB or LSS data.

4.  Specify proto-object/worldline receipts that remain distinct from production particle claims.
