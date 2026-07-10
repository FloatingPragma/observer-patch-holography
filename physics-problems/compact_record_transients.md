# Compact Record Transients

## Motivating Result

This note entered the queue because compact transients exposed two separate
surprises at once. First, several fast radio bursts point to old compact
environments rather than only young magnetar channels: FRB 20200120E was
localized to a globular cluster in M81
([Nature 602, 585-589, 2022](https://www.nature.com/articles/s41586-021-04354-w)),
and FRB 20240209A was traced to the outskirts of an old elliptical host
([Keck Observatory, 2025](https://keckobservatory.org/fbr/)). Second,
gravitational-wave catalogs contain high-mass merger candidates such as
GW190521 where hierarchical black-hole recycling is a serious possibility
([LIGO/Virgo, 2020](https://www.ligo.org/detections/GW190521.php)).

The OPH question is not "can OPH explain FRBs or black-hole mergers" in
general. The question is whether compact record surfaces can define a
quotient-correct transient law for repeaters, old-host compact sources,
hierarchical black-hole genealogies, packet emissions, detector thinning,
censoring, controls, and frozen likelihoods without silently promoting a
visual or simulator artifact into a physical prediction.

**Status:** solved as a conditional compact-transient theorem package and
simulator contract. The default simulator output is conditional
phenomenology. A frozen physical prediction requires the full receipt
conjunction below. A source-only OPH prediction remains work in progress until
the compact source action, emission microphysics, physical compact-repair
clock, old-host FRB source law, and black-hole genealogy prior are derived from
OPH source data.

Date: 2026-07-08

## Claim Boundary

Compact record surfaces are not enough by themselves. The lane must keep these
non-identifications explicit:

```text
compact record surface is not an FRB rate law
normal form is not a burst timing law
repair eigenvalue is not a physical duration without a clock map
packet label is not emitted fluence, gamma energy, or GW strain
old host is not an OPH source theorem
hierarchical-generation tag is not allowed to read the ringdown residual
detected catalog rows are not the full likelihood without exposure and censoring
```

A compact-transient simulator may generate plausible bursts, ringdown tails, or
catalog rows only as diagnostics until it declares the quotient, source law,
repair-emission kernel, packet parent, physical clock, propagation map,
detector thinning, censoring, point-process likelihood, controls, refinement
audit, frozen hashes, and promotion receipt.

## Compact-Transient Objects

The compact-record surface uses

```math
Q_r^{\rm CR},\quad
B_r,\quad
\Phi_R,\quad
\mu_r^{\rm CR},\quad
K_{\Gamma,r},\quad
E_c,\quad
\mathcal P_{\rm prop},\quad
\mathcal R_{\rm det},\quad
\Lambda_c.
```

The complete compact-transient lane also requires

```math
\mathsf{Hist}_r,\quad
\mathsf{ObsWin},\quad
\mathsf{Thin}_c,\quad
\mathsf{Cens}_c,\quad
\mathsf{Ctrl}_c,\quad
\mathsf{Ref}_{sr},\quad
\mathsf{Err}_r,\quad
\mathsf{Promote}_r.
```

`Hist_r` stores compact-object histories across multiple events. `ObsWin`
records survey time, cadence, dead time, sky coverage, and non-detection
intervals. `Thin_c` is the detector selection kernel. `Cens_c` records upper
limits, vetoes, and exposure-only non-detections. `Ctrl_c` is the frozen control
family. `Ref_sr` compares regulators. `Err_r` is the simulator error ledger.
`Promote_r` maps receipts to CR0-CR4 claim tiers.

## History Space

For regulator $r$, define the compact-transient history space

```math
\mathsf{Hist}_r^{\rm CR}
=
\left\{
H=(q_0,\gamma_1,\Pi_1,\tau_1,\ldots,\gamma_N,\Pi_N,\tau_N,q_N)
\right\}.
```

Here $q_k\in Q_r^{\rm CR}$ is a compact record-surface quotient state,
$\gamma_k:q_{k-1}\to q_k$ is an accepted compact repair path, $\Pi_k$ is the
emitted packet ledger, and $\tau_k$ is the physical or dimensionless commit
time. A history is source-valid when $q_0$ is drawn from $\mu_r^{\rm CR}$, each
transition is drawn from $K_{\Gamma,r}$, and each emitted packet ledger obeys

```math
C(B(q_{k-1}))-C(B(q_k))
=
C_{\rm pkt}(\Pi_k)+C_{\rm env,k}+\rho_k.
```

This object is load-bearing for FRB repeaters and hierarchical black holes:
they are histories, not independent catalog rows.

## Observing Windows, Thinning, and Censoring

An observing window is

```math
\mathsf{ObsWin}_c
=
(T_{\rm start},T_{\rm stop},\Omega_{\rm sky},\mathcal E_c,\mathcal D_c,\mathcal S_c).
```

For FRBs this includes beam pattern, time-on-source, fluence threshold,
bandwidth, dispersion-measure search range, and downtime. For gravitational
waves this includes detector network, PSD, antenna pattern, duty cycle, trigger
threshold, and parameter-estimation prior.

The detector thinning kernel is

```math
\mathsf{Thin}_c(d\mathcal O\mid y_c,\mathsf{ObsWin}_c)
=
p_{\rm det}(y_c;\mathsf{ObsWin}_c)
\mathcal R_{{\rm det},c}(d\mathcal O\mid y_c).
```

The censoring kernel is

```math
\mathsf{Cens}_c(dU\mid y_c,\mathsf{ObsWin}_c)
=
\left[1-p_{\rm det}(y_c;\mathsf{ObsWin}_c)\right]
\mathcal U_c(dU\mid y_c).
```

The observed record space is

```math
\mathcal Z_c=\mathcal O_c\sqcup\mathcal U_c.
```

Detections and non-detections must be scored together.

## Compact-Transient Law

The full compact-transient law is the pushforward

```math
\mathbb P_{c,r}^{\rm CT}
=
\mu_r^{\rm CR}
\to
K_{\Gamma,r}^{\rm hist}
\to
E_{c,r}
\to
\mathcal P_{\rm prop,c}
\to
(\mathsf{Thin}_c\sqcup\mathsf{Cens}_c).
```

Expanded:

```math
\begin{aligned}
\mathbb P_{c,r}^{\rm CT}(dZ)
&=
\int
\mu_r^{\rm CR}(dq_0)
\int
K_{\Gamma,r}^{\rm hist}(dH\mid q_0)
\int
E_{c,r}(dy\mid H) \\
&\quad\times
\mathcal P_{\rm prop,c}(dy'\mid y,z,h)
\left[
\mathsf{Thin}_c(dZ\mid y',\mathsf{ObsWin}_c)
+
\mathsf{Cens}_c(dZ\mid y',\mathsf{ObsWin}_c)
\right].
\end{aligned}
```

This is the object the simulator must sample or approximate.

## Theorem Package

### Theorem 1: Compact Record Realization

**Statement.** A compact object is an OPH compact record surface when its finite
presentation realizes observer-like material:

```math
\mathsf O_i=
\left(
\mathcal A_i,\rho_i,\mathcal R_i,
\{(\mathcal I_e,\pi_{i,e})\}_{e\ni i},
\mathcal U_i,\mathrm{Chk}_i
\right)
```

with compact support and boundary data.

**Proof.** A compact record surface adds compact support, boundary state,
mismatch ledger, and packet ledger to the abstract observer-patch tuple. OPH
identifies patches by accessible algebra, visible interface maps, record
readout, repair instruments, and checkpoint continuation. Hidden coordinates
and presentation labels are silent under the quotient. $\square$

### Theorem 2: Compact Quotient Invariance

**Statement.** If every simulator readout factors through
$Q_r^{\rm CR}=\Sigma_r^{\rm CR}/\Gamma_r^{\rm CR}$, compact-transient
observables are invariant under gauge representatives, coordinate charts, mesh
labels, port labels, worker IDs, repair schedules, hidden carrier coordinates,
and inert ancillas.

**Proof.** Equivalent presentations map to the same quotient state. If source,
transition, emission, propagation, detector, and likelihood maps depend only on
that quotient state, they assign the same observed-record law. $\square$

### Theorem 3: Normal Form Is Not a Transient Source Law

**Statement.** The compact normal-form map
$n_r:Q_r^{\rm CR}\to N_r^{\rm CR}$ does not determine FRB rates, burst timings,
black-hole generation rates, host distributions, or waveform residual
amplitudes.

**Proof.** The pushforward $\mathcal C_Q(\mu)=n_{r\#}\mu$ is idempotent, and
every law supported on $N_r^{\rm CR}$ is fixed. Therefore many source laws are
compatible with the same normal-form map. $\square$

### Theorem 4: Compact Repair Normal Form

**Statement.** If compact repair transitions are quotient-local, preserve
protected boundary and sector data except for declared packet emissions,
strictly lower a well-founded mismatch ledger, satisfy local diamond, and are
repair-complete, then accepted repair induces a total, idempotent,
boundary-preserving, schedule-independent normal-form map
$\operatorname{Rep}_{R,r}:Q_r^{\rm CR}\to Q_r^{\rm CR}$.

**Proof.** Exact descent gives termination. Local diamond plus repair
completeness gives confluence. Termination plus confluence gives a unique
normal form independent of schedule. $\square$

### Theorem 5: Packet-Conserving Repair Discharge

**Statement.** A compact repair path is physical only if it emits a packet
ledger $\Pi(\gamma)$ satisfying

```math
C(B(q_0))-C(B(q_N))
=
C_{\rm pkt}(\Pi(\gamma))+C_{\rm env}(\gamma)+\rho_{\rm cons}(\gamma),
\qquad
|\rho_{\rm cons}|\le\varepsilon_{\rm cons}.
```

**Proof.** Local conservation holds at each accepted step. Summing over steps
telescopes the intermediate boundary terms. $\square$

### Theorem 6: Finite Packet Parent

**Statement.** A compact emission channel is physical only if its packet parent
supplies finite covariant packet data sufficient for stress, exchange,
response, propagation, and detector readout.

**Proof.** A scalar event row does not determine stress, propagation, energy
exchange, polarization, chirality, waveform content, or detector response. The
finite packet parent supplies those missing quotient-visible records. $\square$

### Theorem 7: Packetized Kernel Disintegration

**Statement.** A compact repair-emission kernel
$K_{\Gamma,r}(dq',d\Pi,d\ell,d\tau\mid q)$ has state, packet, receipt, and
clock marginals. A simulator that samples those pieces separately is faithful
only when its joint sampler reconstructs the same coupled law.

**Proof.** Marginals exist for finite or standard Borel kernels. Equal marginals
do not determine the joint law, so independence shortcuts require a declared
factorization theorem. $\square$

### Theorem 8: Physical Clock

**Statement.** A compact repair eigenvalue becomes a physical duration only
after a physical clock map
$\Delta\tau_{\rm phys}=\mathcal T_{\rm clock}(q,\ell;b,\theta,h)$ is supplied.
For a reversible finite kernel,

```math
\Gamma_{\rm phys}
=
-\frac{\log r_\star}{\Delta\tau_{\rm phys}},
\qquad
r_\star=\max_{\lambda_j\ne1}|\lambda_j(K_{\Gamma,r}^Q)|.
```

**Proof.** Eigenvalues measure contraction per repair step. Physical rates are
logarithmic contraction divided by physical step time. $\square$

### Theorem 9: Compact Event Pushforward

**Statement.** The observer-facing compact-transient law is the composed
Markov law

```math
(\mathsf{Thin}_c\sqcup\mathsf{Cens}_c)
\circ\mathcal P_{\rm prop,c}
\circ E_{c,r}
\circ K_{\Gamma,r}^{\rm hist}
\circ\mu_r^{\rm CR}.
```

**Proof.** Composition of Markov kernels gives a probability law on detection
and censoring records. Quotient invariance follows when each map is defined on
the compact quotient. $\square$

### Theorem 10: Detector Thinning and Non-Detection

**Statement.** If $\lambda_{\rm src}(y)$ is the source-event intensity and
$p_{\rm det}(y)$ is the selection probability, detected events have intensity
$p_{\rm det}(y)\lambda_{\rm src}(y)$ and censored exposure has complementary
intensity $(1-p_{\rm det}(y))\lambda_{\rm src}(y)$.

**Proof.** Independent thinning retains each event with probability
$p_{\rm det}$ and rejects it with probability $1-p_{\rm det}$. $\square$

### Theorem 11: Marked Catalog Likelihood

**Statement.** A marked compact-transient catalog in observing window
$\mathsf{ObsWin}_c$ has inhomogeneous point-process likelihood

```math
\log\mathcal L(D\mid M)
=
\sum_{i=1}^N\log \Lambda_c(\mathcal O_i)
-
\int_{\mathsf{ObsWin}_c}\Lambda_c(d\mathcal O).
```

**Proof.** This is the likelihood of an inhomogeneous marked Poisson process:
observed events contribute the log intensity and the compensator scores the
expected event count over exposure. $\square$

### Theorem 12: Repeater History Likelihood

**Statement.** For a repeating compact source $s$ with marked conditional
intensity $\lambda_s(t,m\mid\mathcal H_t)$, detected bursts have likelihood

```math
\log\mathcal L_s
=
\sum_k\log\lambda_s(t_k,m_k\mid\mathcal H_{t_k})
-
\int_{T_{\rm obs}}\int\lambda_s(t,m\mid\mathcal H_t)p_{\rm det}(t,m)\,dm\,dt.
```

**Proof.** This is the conditional-intensity likelihood with detector thinning.
The compensator subtracts expected detectable bursts across the exposure
window. $\square$

### Theorem 13: FRB Repair-Reload

**Statement.** If an old compact record source has a mismatch reservoir whose
hazard increases with reservoir level, and if fluence is monotone in discharge,
then the waiting time to the next high-fluence burst is stochastically
increasing in the previous high fluence after conditioning on source identity,
cadence, exposure, and host class.

**Proof.** A larger burst discharges more reservoir, leaving a lower
post-burst reservoir. Under monotone reload and monotone hazard, a lower
starting reservoir reaches the same hazard threshold later. $\square$

The first FRB prediction is therefore:

```text
old/GC repeaters show fluence-conditioned recovery after cadence and exposure correction.
```

### Theorem 14: Host-Mixture Identifiability

**Statement.** In a host model
$S_{\rm FRB}=A_y{\rm SFR}+A_oM_{\star,\rm old}+A_{\rm GC}M_{\rm GC}$, the
coefficients are identifiable only when the exposure-corrected features are
not collinear under the survey selection function.

**Proof.** If a nonzero vector annihilates the exposure-weighted feature vector
almost everywhere, two coefficient vectors give the same source density. Full
rank removes that degeneracy locally. $\square$

### Theorem 15: Record-Surface Recycling

**Statement.** A hierarchical black-hole merger is a compact record-surface
recycling morphism that consumes two compact boundary ledgers and an orbital
ledger, emits a GW packet, and updates a genealogy DAG:

```math
\mathcal G_3=\operatorname{Join}(\mathcal G_1,\mathcal G_2,y_{\rm GW}),
\qquad
g_3=1+\max(g_1,g_2).
```

**Proof.** Conservation assigns the remnant ledger to the two progenitors plus
orbital ledger minus the emitted GW packet. The genealogy DAG records the
public descent relation. $\square$

### Theorem 16: No Generation Leakage

**Statement.** A ringdown repair-tail test is valid only if the generation prior
$p(g_{\rm gen}\ge2\mid M_1,M_2,\chi_1,\chi_2,q,h_{\rm env})$ is computed
without using the ringdown residual being tested.

**Proof.** If the residual defines the generation prior, the predictor contains
the target variable. If mass, spin, mass ratio, and environment define the
prior, the residual remains held out. $\square$

### Theorem 17: Linear Repair Tail

**Statement.** If a remnant repair perturbation satisfies $\dot x=A_{\rm rep}x$
and an active eigenpair is $-\Gamma_{\rm rep}\pm i\omega_{\rm rep}$, then the
held-out repair template is a damped sinusoid plus faster modes. If its
amplitude is proportional to an independently computed generation probability,
hierarchical candidates should have larger expected repair-tail amplitude than
first-generation controls.

**Proof.** Linear systems decompose into eigenmodes. A complex eigenpair gives
a damped sinusoid, and an independently supplied initial projection controls
its expected amplitude. $\square$

### Theorem 18: Control-Model Dominance

**Statement.** An OPH compact-transient claim is physical only if the frozen OPH
model improves held-out likelihood relative to declared controls by a
predeclared threshold.

**Proof.** A model that does not outperform relevant controls adds no
predictive content. Post-hoc controls or thresholds do not test the frozen
law. $\square$

### Theorem 19: Refinement Stability

**Statement.** If source laws, repair kernels, and emission maps cohere under
coarse maps $c_{sr}:Q_s^{\rm CR}\to Q_r^{\rm CR}$ with defects
$\delta^\mu_{sr}$, $\delta^K_{sr}$, and $\delta^E_{sr}$, the output event laws
differ by at most the accumulated defects plus detector and propagation error.

**Proof.** Couple source draws, transition paths, and emitted packets under the
declared defects. Propagation and detector kernels are Markov contractions in
total variation. $\square$

### Theorem 20: Simulator Accuracy

**Statement.** If the simulator approximates source, kernel, emission,
propagation, detector, canonicalization, clock, and Monte Carlo steps with
declared errors, then

```math
|\mathbb P_{\rm sim}-\mathbb P_{\rm paper}|_{\rm TV}
\le
\varepsilon_\mu+\mathbb E[N]\varepsilon_K+\varepsilon_E+\varepsilon_P
+\varepsilon_D+\varepsilon_{\rm canon}+\varepsilon_{\rm clock}
+\varepsilon_{\rm MC}.
```

**Proof.** Couple each stage. The final observer-record law differs only if at
least one coupling fails. $\square$

### Theorem 21: Promotion

**Statement.** A compact-transient simulator output is a physical OPH
prediction only if all promotion gates pass:

```text
QUOTIENT SOURCE KERNEL CLOCK PACKET PROP DETECTOR POINTPROCESS CENSORING
CONTROLS REFINEMENT FREEZE LIKE
```

**Proof.** Each gate is a required map in the source-to-record chain. Missing
one gate leaves the law undefined, unphysical, unscored, or post hoc. $\square$

### Theorem 22: Falsification

**Statement.** A compact-record transient channel is falsified at claim tier CR3
if frozen source, exposure, detector, controls, and likelihood receipts exclude
the distinguishing held-out predictions below predeclared amplitudes.

**Proof.** A CR3 model is a frozen probability law. Held-out data that exclude
its predeclared distinguishing amplitudes falsify that law on the declared
domain. A law changed after inspection was not tested. $\square$

## Claim Tiers

Use exactly these labels:

```text
CR0_VOCABULARY_ONLY
CR1_QUOTIENT_DIAGNOSTIC
CR2_CONDITIONAL_PHENOMENOLOGY
CR3_FROZEN_PHYSICAL_PREDICTION
CR4_SOURCE_ONLY_OPH_PREDICTION
```

`CR0` means record-surface vocabulary exists but no event law exists. `CR1`
means quotient and normal forms exist but no physical event law exists. `CR2`
means source law, repair-emission kernel, packet parent, detector model,
censoring, and likelihood exist but at least one map is phenomenological. `CR3`
means all maps are frozen, controls are declared, censoring is included,
refinement is checked, and likelihood receipts pass. `CR4` means source
action, emission microphysics, physical clock, old-host source law, and
genealogy prior are OPH-derived without target leakage.

## Required Receipts

The compact-transient receipt bundle must include:

```text
COMPACT_HISTORY_RECEIPT
COMPACT_QUOTIENT_RECEIPT
COMPACT_SOURCE_LAW_RECEIPT
PACKETIZED_KERNEL_RECEIPT
PHYSICAL_CLOCK_RECEIPT
FINITE_PACKET_PARENT_RECEIPT
PACKET_CONSERVATION_RECEIPT
PROPAGATION_RECEIPT
DETECTION_THINNING_RECEIPT
CENSORING_AND_UPPER_LIMIT_RECEIPT
POINT_PROCESS_LIKELIHOOD_RECEIPT
REPEATER_HISTORY_LIKELIHOOD_RECEIPT
FRB_SOURCE_IDENTITY_RECEIPT
FRB_CADENCE_EXPOSURE_RECEIPT
BH_GENEALOGY_DAG_RECEIPT
NO_GENERATION_LEAKAGE_RECEIPT
CONTROL_MODEL_RECEIPT
REFINEMENT_STABILITY_RECEIPT
SIMULATOR_ACCURACY_RECEIPT
FROZEN_HASHES_RECEIPT
HELDOUT_LIKELIHOOD_RECEIPT
PROMOTION_AUDIT_RECEIPT
```

## Simulator Fail-Closed Rules

The simulator must refuse physical promotion if any are true:

1. quotient schema missing;
2. canonicalizer non-deterministic;
3. likelihood reads representative labels;
4. source law missing;
5. source law fitted but labeled source-only;
6. packetized kernel missing;
7. physical clock missing;
8. packet parent missing;
9. conservation residual exceeds bound;
10. detector thinning missing;
11. censoring or upper-limit model missing;
12. point-process compensator missing;
13. controls missing;
14. refinement stability missing;
15. frozen hashes missing;
16. FRB source identity missing;
17. FRB cadence or exposure missing;
18. black-hole genealogy DAG missing;
19. generation prior uses ringdown residuals;
20. waveform template is tuned after residual inspection.

## First Implementation Targets

Implement FRBs first with three frozen controls:

```text
M0 = young-only
M1 = young+old/GC with Poisson or Weibull timing
M2 = young+old/GC with repair-reload timing
```

The promotion condition is

```math
\log\mathcal L_{\rm heldout}(M_2)
-
\max[
\log\mathcal L_{\rm heldout}(M_0),
\log\mathcal L_{\rm heldout}(M_1)
]
>
\Delta_{\rm min}.
```

Implement black-hole recycling second:

```text
genealogy DAG
-> generation prior
-> frozen repair-tail template
-> injection tests
-> stacked ringdown likelihood
```

Never use:

```text
ringdown residual -> generation label -> claim success
```

That is target leakage.

## Scope

This note closes the compact-transient lane at conditional simulator level.
The simulator samples or approximates

```text
mu_CR -> K_hist -> E -> propagation -> thinning/censoring -> catalog likelihood
```

and emit the receipt ladder. It should not claim that OPH explains compact
transients without rate, host, timing, waveform, exposure, censoring, control,
refinement, and frozen-likelihood receipts.
