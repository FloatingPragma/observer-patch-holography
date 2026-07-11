# Compact Record Transients

## Motivating Result

This note entered the queue because compact transients exposed two separate
surprises at once. First, several fast radio bursts point to old compact
environments rather than only young magnetar channels: FRB 20200120E was
localized to a globular cluster in M81
([Nature 602, 585-589, 2022](https://www.nature.com/articles/s41586-021-04354-w)),
and FRB 20240209A was traced to the outskirts of an old elliptical host
([Eftekhari et al., 2025](https://doi.org/10.3847/2041-8213/ad9de2)). Second,
gravitational-wave catalogs contain high-mass merger candidates such as
GW190521 where hierarchical black-hole recycling is a serious possibility
([Abbott et al., 2020](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.101102)).

These observations establish population diversity and motivate delayed FRB
channels and hierarchical-merger models; they do not select OPH. FRBs and
black-hole mergers are separate physical case studies. They share an
observer-facing source-to-record audit architecture here, not a demonstrated
common repair mechanism.

The OPH question is not "can OPH explain FRBs or black-hole mergers" in
general. The question is whether an observer-like compact-source patch can be
instantiated with local state, boundaries or ports, readback, persistent
records, feedback or repair moves, emitted packets, and public evidence, and
whether that object can be connected to a normalized transient law without
silently promoting a visual or schema artifact into a physical prediction.

**Status:** completed visual/schema contract at
`CR1_QUOTIENT_DIAGNOSTIC`, not validated compact-transient phenomenology. The
objects and receipts below specify what an end-to-end model must contain, but
this note supplies no normalized channel-specific source law, fitted detector
model, catalog analysis, held-out likelihood result, or OPH-derived emission
microphysics. `CR2_CONDITIONAL_PHENOMENOLOGY` begins only after concrete toy or
empirical kernels, units, parameters, exposure, data, and reproducible outputs
are attached. A source-derived OPH prediction remains work in progress.

Date: 2026-07-11

## Standard Baseline, Remaining Problem, and OPH Contribution

**Established physics and inference.** Young magnetars are an observationally
supported FRB channel, while globular-cluster and quiescent hosts motivate
additional delayed formation routes. Standard repeater analyses use renewal,
Poisson, Weibull, clustered, periodic, and state-space timing models with survey
selection. Standard hierarchical black-hole models use cluster or AGN
population dynamics, numerical-relativity remnant fits, detector selection,
and genealogical mixture priors. General-relativistic ringdown is already a
sum of damped quasinormal modes.

**Unresolved part.** The open work is to identify channel-specific source and
emission microphysics, infer rates and host mixtures without selection bias,
and make distinguishing held-out predictions. For black holes, an additional
tail must be distinguished from ordinary Kerr modes, overtones, nonlinear
modes, precession, eccentricity, waveform error, and detector noise.

**OPH contribution in this note.** OPH supplies the explicit observer-patch and
receipt boundary: source state and history, ports, packet readback, record-
conditioned updates, propagation, exposure, detector response, censoring,
controls, and promotion evidence are kept in one auditable chain. The Markov
and point-process machinery is standard. OPH becomes physically distinctive
only if its source/repair structure yields a frozen prediction that standard
controls do not already contain.

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

A compact-transient visualization may generate plausible bursts, ringdown
tails, or catalog rows only as diagnostics until a model declares the quotient,
source law, repair-emission kernel, packet parent, physical clock, propagation
map, detector thinning, censoring, point-process likelihood, controls,
refinement audit, frozen hashes, and promotion receipt. A receipt name or
schema field is not evidence that the corresponding physical map has been
implemented or validated.

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
time. Physical and dimensionless clocks must never be mixed without an explicit
clock map. A history is source-valid when $q_0$ is drawn from
$\mu_r^{\rm CR}$, each
transition is drawn from $K_{\Gamma,r}$, and each emitted packet ledger obeys

```math
C(B(q_{k-1}))-C(B(q_k))
=
C_{\rm pkt}(\Pi_k)+C_{\rm env,k}+\rho_k.
```

Here $C$ is a declared vector of conserved quantities, such as four-momentum,
angular momentum, and charge, with a sign convention, common units, and a
componentwise or normed residual bound $\|\rho_k\|\le\varepsilon_{\rm cons}$.
Four-momenta at separated events must be compared in a declared common tetrad
or by a specified parallel-transport/boundary-flux convention. Angular
momentum additionally needs an origin or center-of-mass convention and a stated
spin/orbital split. In a generic curved spacetime, a global conserved vector
exists only with the required symmetry or asymptotic/boundary construction;
otherwise the ledger must use the appropriate local covariant flux law. The
equation is a required model assumption until a channel-specific emission
calculation supplies those quantities.

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

The complement

```math
p_{\rm miss}(y_c;\mathsf{ObsWin}_c)=1-p_{\rm det}(y_c;\mathsf{ObsWin}_c)
```

describes a latent missed event. It is not automatically an observed censored
row. A censoring or upper-limit kernel is appropriate only when a trial,
follow-up observation, or externally triggered event is itself known:

```math
\mathsf{Cens}_c(dU\mid y_c,w)
=p_{\rm cens}(y_c;w)\mathcal U_c(dU\mid y_c,w).
```

For a blind survey, a window with no detected events is scored by the
point-process compensator or zero-count probability, not by fabricating one row
for every latent missed event. The observed record space, when both kinds of
records genuinely exist, is

```math
\mathcal Z_c=\mathcal O_c\sqcup\mathcal U_c.
```

Detected events, known upper limits, and zero-detection exposure must be scored
together without double counting.

## Compact-Transient Law

**Model contract.** Let $\zeta$ collect source context such as host,
environment, redshift, distance, and sky position. A channel-specific model
must define a normalized source law $\mu_r^{\rm CR}(dq_0,d\zeta)$ and measurable
kernels for history, emission, propagation, and observation. Its trial-level
pushforward is

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
\mathsf{Obs}_c.
```

Expanded, with every latent variable bound:

```math
\begin{aligned}
\mathbb P_{c,r}^{\rm CT}(dZ)
&=
\int
\mu_r^{\rm CR}(dq_0,d\zeta)
\int
K_{\Gamma,r}^{\rm hist}(dH\mid q_0,\zeta)
\int
E_{c,r}(dy\mid H,\zeta) \\
&\quad\times
\int \mathcal P_{\rm prop,c}(dy'\mid y,\zeta)
\mathsf{Obs}_c(dZ\mid y',\mathsf{ObsWin}_c).
\end{aligned}
```

$\mathsf{Obs}_c$ may include detection and genuine censoring outcomes for a known
trial. A blind catalog additionally requires a source-event intensity and
exposure compensator. The arrows and integral above are a schema until all
spaces, normalized kernels, units, parameters, assumptions, and data products
are attached; no such attachment is claimed in this note.

## Conditional Definitions and Lemmas

### Operational Definition 1: Compact Record Realization

An astrophysical compact object is modeled as an OPH compact record surface
only when its finite presentation instantiates observer-like material:

```math
\mathsf O_i=
\left(
\mathcal A_i,\rho_i,\mathcal R_i,
\{(\mathcal I_e,\pi_{i,e})\}_{e\ni i},
\mathcal U_i,\mathrm{Chk}_i
\right)
```

with compact support and boundary data.

This is an operational modeling definition, not a proof that every compact
object is self-reading. A channel must identify the physical boundary, ports,
readback, persistent record, feedback/update dynamics, and emitted public
packet. A compact record surface adds compact support, boundary state,
mismatch ledger, and packet ledger to the abstract observer-patch tuple. OPH
identifies patches by accessible algebra, visible interface maps, record
readout, repair instruments, and checkpoint continuation. Only coordinates and
presentation labels proven inert under all retained dynamics may be quotiented.

### Conditional Lemma 2: Compact Quotient Invariance

**Statement.** If every simulator readout factors through
$Q_r^{\rm CR}=\Sigma_r^{\rm CR}/\Gamma_r^{\rm CR}$, compact-transient
observables are invariant under gauge representatives, coordinate charts, mesh
labels, inert port relabelings, worker IDs, and inert ancillas.

**Proof.** Equivalent presentations map to the same quotient state. If source,
transition, emission, propagation, detector, and likelihood maps depend only on
that quotient state, they assign the same observed-record law. A repair
schedule or hidden state that changes timing, emission, or likelihood is
physical history and must not be quotiented. $\square$

### Lemma 3: Normal Form Is Not a Transient Source Law

**Statement.** The compact normal-form map is typed as an idempotent
$n_r:Q_r^{\rm CR}\to Q_r^{\rm CR}$ with image
$N_r^{\rm CR}\subseteq Q_r^{\rm CR}$. It does not determine FRB rates, burst timings,
black-hole generation rates, host distributions, or waveform residual
amplitudes.

**Proof.** The pushforward $`\mathcal C_Q(\mu)=n_{r\#}\mu`$ is idempotent, and
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

### Evidence Requirement 5: Packet-Conserving Repair Discharge

A compact repair path is promoted as physical only if each accepted local step
obeys a declared vector conservation law and the path emits a packet ledger
$\Pi(\gamma)$ satisfying

```math
C(B(q_0))-C(B(q_N))
=
C_{\rm pkt}(\Pi(\gamma))+C_{\rm env}(\gamma)+\rho_{\rm cons}(\gamma),
\qquad
\|\rho_{\rm cons}\|\le\varepsilon_{\rm cons}.
```

**Proof.** Local conservation holds at each accepted step. Summing over steps
telescopes the intermediate boundary terms. $\square$

### Evidence Requirement 6: Finite Packet Parent

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

### Conditional Lemma 8: Kernel Relaxation Clock

**Model assumptions.** A finite, time-homogeneous, irreducible, aperiodic,
reversible kernel has a stationary law, and each transition uses the same
physical step duration supplied by
$\Delta\tau_{\rm phys}=\mathcal T_{\rm clock}(q,\ell;b,\theta,h)$.
When $0<r_\star<1$, its slowest spectral relaxation rate is

```math
\Gamma_{\rm phys}
=
-\frac{\log r_\star}{\Delta\tau_{\rm phys}},
\qquad
r_\star=\max_{\lambda_j\ne1}|\lambda_j(K_{\Gamma,r}^Q)|.
```

**Proof.** Under the assumptions above, $r_\star$ controls asymptotic
$L^2$ contraction per step, and physical relaxation rates are logarithmic
contraction divided by the physical step time. $\square$

This is a relaxation or mixing rate, not automatically an event duration,
emission clock, or compact-object repair time. Those identifications require a
channel-specific physical derivation or calibration.

### Conditional Lemma 9: Compact Event Pushforward

**Statement.** The observer-facing compact-transient law is the composed
Markov law

```math
\mathsf{Obs}_c
\circ\mathcal P_{\rm prop,c}
\circ E_{c,r}
\circ K_{\Gamma,r}^{\rm hist}
\circ\mu_r^{\rm CR}.
```

**Proof.** Composition of normalized measurable Markov kernels gives a
probability law on the declared observation records. Quotient invariance
follows when each map is defined on the compact quotient. This lemma applies
only after the schematic maps have been instantiated. $\square$

### Standard Lemma 10: Detector Thinning and Latent Misses

**Statement.** If $\lambda_{\rm src}(y)$ is the source-event intensity and
$p_{\rm det}(y)$ is the selection probability, detected events have intensity
$p_{\rm det}(y)\lambda_{\rm src}(y)$ and the latent rejected process has
complementary intensity $(1-p_{\rm det}(y))\lambda_{\rm src}(y)$.

**Proof.** Independent thinning retains each event with probability
$p_{\rm det}$ and rejects it with probability $1-p_{\rm det}$. $\square$

The rejected process is not observed censoring. Zero-detection exposure is
scored through the compensator; a censoring row requires a known trial or
external trigger.

### Standard Lemma 11: Conditional Poisson Catalog Likelihood

**Model assumption.** Conditional on model $M$, suppose the detected catalog is
an inhomogeneous marked Poisson process with intensity density $\lambda_c$
relative to base measure $\nu_c$. Then

```math
\log\mathcal L(D\mid M)
=
\sum_{i=1}^N\log \lambda_c(\mathcal O_i)
-
\int_{\mathsf{ObsWin}_c}\lambda_c(\mathcal O)\,\nu_c(d\mathcal O).
```

**Proof.** This is the likelihood of an inhomogeneous marked Poisson process:
observed events contribute the log intensity and the compensator scores the
expected event count over exposure. Event-level measurement uncertainty
requires integrating $\lambda_c$ against the measurement likelihood. Repeating
sources with history-dependent intensity are not Poisson catalogs. $\square$

### Standard Lemma 12: Repeater History Likelihood with Selection

If $\lambda_s(t,m\mid\mathcal H_t)$ is the source intensity and the complete
latent history $\mathcal H_t$ is known, independent detector thinning gives

```math
\log\mathcal L_s
=
\sum_k\log\!\left[
p_{\rm det}(t_k,m_k)\lambda_s(t_k,m_k\mid\mathcal H_{t_k})
\right]
-
\int_{T_{\rm obs}}\int\lambda_s(t,m\mid\mathcal H_t)p_{\rm det}(t,m)\,dm\,dt.
```

**Proof.** This is the conditional-intensity likelihood with detector thinning.
The compensator subtracts expected detectable bursts across the exposure
window. $\square$

The $\log p_{\rm det}$ event terms may be omitted only as a stated
parameter-independent constant. If missed bursts change a reservoir or future
hazard, the latent history is not known. The observed conditional intensity
must then be filtered,

```math
\lambda_{\rm obs}(t,m\mid\mathcal H_t^{\rm obs})
=p_{\rm det}(t,m)\,
\mathbb E\!\left[
\lambda_s(t,m\mid\mathcal H_t^{\rm latent})
\mid\mathcal H_t^{\rm obs}
\right],
```

and $\lambda_{\rm obs}$ is used in both event and compensator terms.

### Design Hypothesis 13: FRB Reservoir Reload

Let $R^+=g(R^-,E)$ be the post-burst reservoir after intrinsic emitted energy
$E$. Assume, rather than infer, that larger $E$ leaves an ordered lower
post-burst state after conditioning on $R^-$, that the reload flow preserves
that ordering, and that the hazard of a high-energy event is nondecreasing in
$R$. Under those assumptions, the waiting time to the next high-energy burst is
stochastically increasing with the previous emitted energy.

**Derivation.** Ordered post-burst states remain ordered under the assumed
reload flow, so the high-energy hazard along the trajectory following the
larger discharge is pointwise no greater. The survival function
$\exp[-\int_0^t\lambda_H(R_u)du]$ is therefore larger. Without the conditioned
pre-burst state and order-preserving reset assumption, a larger observed burst
need not leave a lower reservoir. $\square$

The first conditional FRB hypothesis is therefore:

```text
old/GC repeaters show emitted-energy-conditioned recovery after source,
bandpass, cadence, exposure, propagation, and missed-burst correction.
```

This is not unique to OPH. Frozen controls must include renewal/Weibull,
periodic, Hawkes or self-inhibiting, self-organized-criticality, and ordinary
magnetar reservoir models.

### Necessary Rank Condition 14: Host Mixture

In a host model
$S_{\rm FRB}=A_y{\rm SFR}+A_oM_{\star,\rm old}+A_{\rm GC}M_{\rm GC}$, the
coefficients have units that convert the corresponding star-formation or mass
tracer to an event intensity. A necessary local identifiability condition is
that the exposure-corrected features are not collinear under the survey
selection function.

**Proof.** If a nonzero vector annihilates the exposure-weighted feature vector
almost everywhere, two coefficient vectors give the same source density. Full
rank removes that linear degeneracy locally, but host-association uncertainty,
missing hosts, tracer error, and a non-identifiable observation model can still
prevent global identification. $\square$

### Operational Genealogy 15: Record-Surface Recycling

A hierarchical black-hole population model may represent a merger as a
record-surface recycling morphism that consumes two compact boundary ledgers
and an orbital ledger, emits a GW packet, and updates a genealogy DAG:

```math
\mathcal G_3=\operatorname{Join}(\mathcal G_1,\mathcal G_2,y_{\rm GW}),
\qquad
g_3=1+\max(g_1,g_2).
```

This is genealogy bookkeeping. A physical implementation must use
numerical-relativity remnant mass, spin, radiated energy, angular momentum, and
recoil fits, together with environmental escape velocity and retention
probability. Conservation assigns the remnant ledger to the two progenitors plus
orbital ledger minus the emitted GW packet. The genealogy DAG records the
public descent relation.

### Evidence Rule 16: No Generation Leakage

**Statement.** A ringdown repair-tail test is valid only if the generation prior
$p(g_{\rm gen}\ge2\mid M_1,M_2,\chi_1,\chi_2,q,h_{\rm env})$ is computed
without using the ringdown residual being tested. Because full-waveform mass
and spin posteriors can contain ringdown information, the generation predictor
must use inspiral-only information, independent environmental data, or an
explicit cross-fitting scheme.

**Proof.** If the residual defines the generation prior, the predictor contains
the target variable. If mass, spin, mass ratio, and environment define the
prior, the residual remains held out. $\square$

### Undemonstrated Ansatz 17: Additional Linear Repair Tail

General relativity already predicts a superposition of damped quasinormal modes
for a perturbed remnant. OPH has not derived an additional mode in this note.
As a provisional ansatz only, suppose an extra repair perturbation obeys
$\dot x=A_{\rm rep}x$, has an active eigenpair
$-\Gamma_{\rm rep}\pm i\omega_{\rm rep}$, and every other retained extra mode
has real part strictly below $-\Gamma_{\rm rep}$. Then that extra component has
a damped-sinusoid tail plus faster extra modes.

**Derivation.** Linear systems decompose into eigenmodes, and a complex
eigenpair gives a damped sinusoid. This algebra does not establish that the mode
exists in nature. A physical OPH claim must derive and freeze its frequency,
damping, dimensionless strain amplitude, start time, and dependence on an
independently computed generation probability, then show that it is not
absorbed by Kerr fundamentals, overtones, nonlinear modes, precession,
eccentricity, waveform error, calibration error, or noise. $\square$

### Promotion Rule 18: Control-Model Dominance

A frozen OPH compact-transient model gains distinguishing predictive support
only if it improves a preregistered held-out proper score relative to declared
controls by a calibrated threshold and also passes absolute calibration,
injection recovery, and multiple-testing checks.

**Proof.** A model that does not outperform relevant controls adds no
predictive content. Post-hoc controls or thresholds do not test the frozen
law. $\square$

### Conditional Lemma 19: Refinement Stability

If source laws, repair kernels, and emission maps cohere in total variation
under
coarse maps $c_{sr}:Q_s^{\rm CR}\to Q_r^{\rm CR}$ with defects
$\delta^\mu_{sr}$, $\delta^K_{sr}$, and $\delta^E_{sr}$, the output event laws
differ by at most the accumulated dimensionless defects plus detector and
propagation kernel error, provided the expected history length is finite.

**Proof.** Couple source draws, transition paths, and emitted packets under the
declared defects. Propagation and detector kernels are Markov contractions in
total variation. $\square$

### Conditional Error Bound 20: Model Approximation

If an implementation approximates source, kernel, emission,
propagation, detector, canonicalization, clock, and Monte Carlo steps with
declared errors, then

```math
\|\mathbb P_{\rm approx}-\mathbb P_{\rm model}\|_{\rm TV}
\le
\varepsilon_\mu+\mathbb E[N]\varepsilon_K+\varepsilon_E+\varepsilon_P
+\varepsilon_D+\varepsilon_{\rm canon}+\varepsilon_{\rm clock}
+\varepsilon_{\rm MC}.
```

**Proof.** Couple each stage. The final observer-record law differs only if at
least one coupling fails. This bound requires uniform per-kernel total-
variation errors, compatible state spaces, finite $\mathbb E[N]$, and all clock
or canonicalization errors already converted to total variation. Until those
conditions and errors are supplied, the equation is a target bound. $\square$

### Evidence Rule 21: Promotion

A compact-transient visualization or model output is a physical OPH prediction
only if all promotion gates are implemented and evidence-backed:

```text
QUOTIENT SOURCE KERNEL CLOCK PACKET PROP DETECTOR POINTPROCESS CENSORING
CONTROLS REFINEMENT FREEZE LIKE
```

**Proof.** Each gate is a required map in the source-to-record chain. Missing
one gate leaves the law undefined, unphysical, unscored, or post hoc. $\square$

### Evidence Rule 22: Falsification

A compact-record transient model at claim tier CR3 is falsified on its declared
domain when it predicts a distinguishing effect $A\ge A_{\min}$ and a frozen
source, exposure, detector, control, and likelihood analysis places a
predeclared upper bound $A<A_{\min}$, or when its held-out proper score misses a
predeclared acceptance region.

**Proof.** A CR3 model is a frozen probability law. Held-out data that exclude
its predeclared distinguishing region falsify that law on the declared domain.
A law changed after inspection was not tested. $\square$

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

The present Markdown/schema package is `CR1`. Derivation provenance and
empirical validation should also be reported separately: an OPH-derived but
untested model is not empirically stronger than a phenomenological model that
passes held-out and replication tests.

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

Each receipt has status `MISSING`, `PLACEHOLDER`, `DECLARED`, or `VALIDATED`.
A schema key, prose equation, or generated Boolean is not a validated receipt;
`VALIDATED` requires an evidence path, content hash, assumptions, units, and a
passing check.

## Model Fail-Closed Rules

Any visualization or implementation must refuse physical promotion if any are
true:

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

An FRB implementation should begin with at least these frozen controls:

```text
M0 = young-only
M1 = young+old/GC with Poisson, renewal, or Weibull timing
M2 = young+old/GC with periodic/Hawkes/self-inhibiting timing
M3 = young+old/GC with the frozen reservoir-reload timing law
```

The promotion condition is

```math
S_{\rm heldout}(M_3)
-
\max_{j<3}S_{\rm heldout}(M_j)
>
\Delta_{\rm min}.
```

$S$ is a preregistered proper score with its per-event normalization,
uncertainty, and calibrated $\Delta_{\min}$ stated before evaluation. Passing
this comparison supports the timing law; it does not by itself identify OPH as
the unique mechanism.

Implement black-hole recycling second:

```text
genealogy DAG
-> generation prior
-> frozen GR/numerical-relativity ringdown baseline
-> optional predeclared extra-tail ansatz
-> injection tests
-> stacked ringdown likelihood
```

Never use:

```text
ringdown residual -> generation label -> claim success
```

That is target leakage.

## Scope

This note closes a visual/schema-level audit at `CR1`; it does not close either
compact-transient physics case. A future channel-specific model must instantiate
and sample or approximate

```text
mu_CR -> K_hist -> E -> propagation -> observation/exposure -> catalog likelihood
```

and emit evidence-backed receipts. It must not claim that OPH explains compact
transients without rate, host, timing, waveform, exposure, selection,
censoring, control, refinement, held-out likelihood, and replication evidence.

## References and Baselines

- F. Kirsten et al., "A repeating fast radio burst source in a globular cluster
  in the nearby galaxy M81", Nature, 2022.
  https://www.nature.com/articles/s41586-021-04354-w
- V. Shah et al., "A Repeating Fast Radio Burst Source in the Outskirts of a
  Quiescent Galaxy", Astrophysical Journal Letters, 2025.
  https://doi.org/10.3847/2041-8213/ad9ddc
- T. Eftekhari et al., "The Massive and Quiescent Elliptical Host Galaxy of the
  Repeating Fast Radio Burst FRB 20240209A", Astrophysical Journal Letters,
  2025. https://doi.org/10.3847/2041-8213/ad9de2
- C. D. Bochenek et al., "A fast radio burst associated with a Galactic
  magnetar", Nature, 2020.
  https://www.nature.com/articles/s41586-020-2872-x
- CHIME/FRB Collaboration, "The First CHIME/FRB Fast Radio Burst Catalog",
  Astrophysical Journal Supplement Series, 2021.
  https://doi.org/10.3847/1538-4365/ac33ab
- R. Abbott et al., "GW190521: A Binary Black Hole Merger with a Total Mass of
  150 Solar Masses", Physical Review Letters, 2020.
  https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.101102
- M. Capano et al., "Multimode Quasinormal Spectrum from a Perturbed Black
  Hole", Physical Review Letters, 2023.
  https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.131.221402
- I. Mandel, W. M. Farr, and J. R. Gair, "Extracting distribution parameters
  from multiple uncertain observations with selection biases", Monthly Notices
  of the Royal Astronomical Society, 2019.
  https://doi.org/10.1093/mnras/stz896
