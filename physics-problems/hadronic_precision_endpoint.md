# Hadronic Precision Endpoint

## Motivating Result

This note entered the queue because several precision results kept pointing at
the same hadronic boundary. The immediate trigger was the final Fermilab muon
`g-2` result: the combined Run 1--6 value reached 127 ppb precision, while the
new experimental world average reached 124 ppb
([Muon g−2 Collaboration, 2025](https://muon-g-2.fnal.gov/result2025.pdf)).
The same hadronic boundary also appears in rare-decay audits: LHCb reported the
rare `Sigma+ -> p mu+ mu-` observation in 2024
([CERN](https://home.cern/lhcb-investigates-rare-s-pmm-decay/)), while updated
rare `B`-decay lepton-universality measurements showed how earlier anomaly
claims can move when the full dataset and calibration change
([CERN](https://home.cern/lhcb-brings-leptons-line/)). The OPH question is
whether these apparently separate tensions share one missing source object: a
same-scheme hadronic precision backend.

**Status:** conditional backend specification. It identifies the QCD objects
and receipts that a source-derived calculation would need, but OPH does not yet
derive the QCD parameter map, physical measure, continuum/confined theory, or
precision observables. Existing finite outputs are prototypes or conventional
QCD comparisons, not source-only hadronic predictions.

Date: 2026-07-08

## Origin

This note records the OPH hadronic-precision audit prompted by the motivating
precision results. It asks whether the fine-structure endpoint, muon $g-2$,
rare $B$ decays, and rare $\Sigma$ decays share one OPH bottleneck.

They share nonperturbative-QCD inputs, but they are not one scalar problem.
HVP, HLbL, and long-distance rare-decay amplitudes require different
correlators, operators, and analytic continuations. A useful common interface
is therefore a typed hadronic precision backend rather than one informal
`rho_had(s)`:

```math
\mathscr H_{\rm had}^{\rm OPH}(P)
=
\left(
\mathfrak Q_r^{\rm QCD}(P),
\mathfrak S_r^{E,{\rm QCD}}(P),
\mathcal H_{\rm had}(P),
J_Q^{R,\mu}(P),
d\rho_Q^{(2)}(s;P),
\Gamma_{QQQQ}^{(4)}(P),
\mathcal A_{\rm tr}^{B,\Sigma}(P),
\Xi_Q(P),
\mathcal E_{\rm sys}(P)
\right).
```

## Shared Backend, Distinct Observables

The two-current spectral measure $d\rho_Q^{(2)}(s;P)$ is enough for the
fine-structure hadronic endpoint and hadronic vacuum polarization in muon
$g-2$. It is not enough for hadronic light-by-light or rare-decay
long-distance amplitudes. Those require tensor-valued four-current correlators
and generally complex transition amplitudes:

```math
\Gamma_{QQQQ}^{(4)}(P),
\qquad
\mathcal A_{\rm tr}^{B,\Sigma}(P).
```

The scalar two-point measure is therefore one marginal of a larger hadronic
precision functor, not the whole backend.

## What Is Standard, What Is Open, And What OPH Adds

Lattice QCD, dispersive methods, effective field theory, Ward identities,
spectral representations, and experiment already provide the standard
framework. The hard unresolved OPH task is to derive the QCD parameters and
physical source measure from upstream OPH data with controlled regulator,
continuum, renormalization, and analytic-continuation errors.

OPH's present contribution is organizational: represent the calculation as an
observer-like self-reading patch with local gauge state, physical boundaries,
current readback, refinement and Ward-repair moves, immutable correlator
records, and a public evidence DAG. That can make cross-observable scheme and
provenance errors easier to expose. It does not uniquely solve nonperturbative
QCD or replace established methods.

Claims below are either definitions, explicit assumptions, or conditional
derivations. A receipt records a recomputed property; its name alone is not
evidence that the property holds.

## OPH-QCD Quotient Ensemble

For regulator $r$, the source QCD ensemble is

```math
\mathfrak Q_r^{\rm QCD}(P)
=
\left(
\Sigma_r^{\rm QCD},
\Gamma_r^{\rm QCD},
Q_r^{\rm QCD},
m_r^{\rm QCD},
S_r^{\rm QCD}(P),
c_{sr},
\mathcal A_r^{\rm phys},
\mathcal O_r^{\rm QCD}
\right),
\qquad
Q_r^{\rm QCD}=\Sigma_r^{\rm QCD}/\Gamma_r^{\rm QCD}.
```

$\Sigma_r^{\rm QCD}$ contains finite gauge links, quark fields, boundary
conditions, source masses, couplings, current insertions, operator basis,
regulator metadata, and current-normalization metadata.
$\Gamma_r^{\rm QCD}$ removes gauge representatives and demonstrably inert
presentation metadata such as mesh names, worker IDs, queue order, and random
number bookkeeping. Physical boundary conditions, current insertions, port
locations, quark masses, topology, and any variable that changes an observable
must remain visible. Quotient invariance has to be proved for each proposed
redundancy; it cannot be obtained by listing a physical label as hidden.

The source law is

```math
\mu_r^{\rm QCD}(q;P)
=
Z_r^{-1}m_r^{\rm QCD}(q)\exp[-S_r^{\rm QCD}(q;P)].
```

This positive Gibbs notation is appropriate only when the Euclidean weight is
real, nonnegative, normalizable, and dimensionless--for example the standard
zero-density, $\theta_{\rm QCD}=0$ branch after fermions are handled in a
positivity-preserving formulation. At generic $\theta$, nonzero baryon
chemical potential, or with a discretization/flavor content whose determinant
is not nonnegative, the measure can have a sign or phase problem; then a
positive $\mu_r^{\rm QCD}$ cannot simply be assumed. The base measure,
boundary conditions, determinant treatment, and normalization must be
declared before precision comparison. A normal form is not a probability law.

## Source-Law Non-Substitution Theorem

**Statement.** A simulator output cannot be promoted to a source-only OPH
hadronic prediction unless it is sampled from, or deterministically approximates,
a declared quotient-intrinsic law $\mu_r^{\rm QCD}$ with refinement maps
$c_{sr}$, current definitions, scheme metadata, and no-target-leakage
receipts.

**Proof sketch.** Let $c_r:Q_r\to Q_r$ be a canonicalizer with
$c_r^2=c_r$ and image $N_r\subseteq Q_r$. It does not determine a
probability law: every law supported on $N_r$ is fixed by its pushforward.
Therefore settled or canonical status can classify
candidate hadronic sectors, but cannot choose the physical QCD law. If the law
is inferred from $\alpha(0)$, $g-2$, $e^+e^-\to{\rm hadrons}$, hadron
masses, rare-decay data, or PDG QCD fits, the result is empirical closure or
calibration, not source-only OPH.

## Source QCD Parameter Map

The missing source map is

```math
\mathcal P_{\rm QCD}^{\rm src}:
P
\mapsto
\left(
g_3(P),
\theta_{\rm QCD}(P),
m_u(P),m_d(P),m_s(P),m_c(P),m_b(P),m_t(P),
Z_{\rm scheme}(P)
\right).
```

Without this map, a lattice or simulator backend can be a high-quality
Standard Model/QCD comparison engine, but not a fully OPH-native source theorem.
If quark masses or $\Lambda_{\overline{\rm MS}}$ are imported from PDG or
hadron spectroscopy, the row is not source-only.

## Illustrative Finite Transfer Slab

For receipt testing, one may introduce the following finite graph transfer
surrogate:

```math
\mathfrak S_r^{E,{\rm QCD}}(P)
=
\left(
Q_r^{\rm QCD},
m_r^0,
J_r,
V_r,
a_{t,r},
\Theta_r^{\rm RP}
\right).
```

With $J_r(q,q')=J_r(q',q)\ge0$, the finite Euclidean Hamiltonian has form

```math
(H_r^E f)(q)
=
\frac{1}{m_r^0(q)}
\sum_{q'}J_r(q,q')\bigl(f(q)-f(q')\bigr)
+
V_r(q)f(q).
```

If $J_r$ is symmetric nonnegative, the event graph is connected,
$m_r^0(q)>0$, and $V_r:Q_r\to\mathbb R$ is real, this finite operator
is self-adjoint and bounded below on $L^2(Q_r,m_r^0)$. Normalize a
nondegenerate strictly positive ground state by
$\sum_q|\Omega_r(q)|^2m_r^0(q)=1$. It defines the finite ground-state
Born measure

```math
\mu_r^{\rm vac}(q;P)=|\Omega_r(q;P)|^2m_r^0(q),
```

The Euclidean transfer operator is

```math
T_r=e^{-a_{t,r}(H_r^E-E_{0,r})}.
```

It is not generally a Markov kernel because $T_r1\ne1$. If the semigroup
is positivity preserving, the ground-state/Doob transform
$(P_tf)(q)=\Omega_r(q)^{-1}\bigl(e^{-t(H_r^E-E_{0,r})}(\Omega_rf)\bigr)(q)$
is Markov and reversible with stationary law $\mu_r^{\rm vac}$. The Born
measure, Euclidean transfer, and Doob dynamics are three related but distinct
objects.

This is a generic graph-Laplacian construction, not a derivation of the lattice
QCD transfer matrix. Calling the law a physical QCD vacuum additionally
requires a local gauge action with fermion content, Osterwalder--Schrader
reflection positivity, the correct symmetries and continuum limit, and
refinement compatibility. An HMC stationary distribution or this surrogate
alone is not automatically an OPH vacuum.

## Confined Hadronic Hilbert Quotient

Assuming a valid OPH-QCD continuum and confinement certificate, define

```math
\mathcal H_{\rm had}(P)
=
\overline{\mathcal H_{\rm phys}^{\rm QCD}(P)}.
```

Here $\mathcal H_{\rm phys}^{\rm QCD}$ denotes gauge-invariant physical
states satisfying the declared constraints and carrying the relevant
color-singlet hadronic spectrum. Confinement is not a known elementary
projector, and Ward normalization acts on currents and correlators rather than
as a generally commuting state-space projector. This definition is conditional
on unquenched QCD with quarks, confinement, physical currents, and spectral
matrix elements. A Yang--Mills continuum route is not enough by itself.

Required continuum certificates include:

- `QCD_CONTINUUM_CERTIFICATE`
- `REFLECTION_POSITIVITY_CERTIFICATE`
- `TRANSFER_CONVERGENCE_CERTIFICATE`
- `CONFINEMENT_SINGLET_CERTIFICATE`

## Ward-Projected Renormalized Current

The electromagnetic current must be a declared conserved or renormalized
current on the same $U(1)_Q$ branch as the Maxwell endpoint. For a local
continuum-like lattice current, the schematic definition is

```math
J_Q^{R,\mu}(x;P)
=
Z_V(P)
\left[
\sum_f Q_f\,\bar q_f(x)\gamma^\mu q_f(x)
\right]
+
J_{\rm improvement}^{\mu}(x;P).
```

An exactly conserved discretized current may instead have $Z_V=1$. Contact
terms belong to time-ordered correlator Ward identities, not generically to
the one-current operator. After their declared subtraction, the correlator
must satisfy

```math
\partial_\mu J_Q^{W,R,\mu}=0,
\qquad
q_\mu\Pi_Q^{\mu\nu}(q;P)=0.
```

The current-normalization ledger is

```math
\mathcal N_Q(P)
=
\left(
Z_V,
J_{\rm conserved},
J_{\rm local},
J_{\rm contact},
\Pi_{\rm Ward},
\Omega_Q,
\text{scheme}
\right).
```

Without this ledger, the spectral measure has an arbitrary normalization
ambiguity.

## Two-Point Spectral Measure

Let $M_{\rm had}^2=P_\mu P^\mu$ be the declared nonnegative invariant-mass
operator on the physical hadronic Hilbert space, after fixing metric signature
and subtracting the vacuum energy. In a finite-volume rest-frame construction
this is the square of the excitation energy, with the momentum projection and
continuum conversion stated. Write $E_{\rm had}(ds)$ for the projection-valued
spectral measure of $M_{\rm had}^2$; it is not the spectral measure of an
otherwise unspecified graph Hamiltonian.

Let $J_Q^{R,\lambda}(f)$ be a spatial current smeared with a declared test
function or finite-volume operator and set
$\psi_\lambda=J_Q^{R,\lambda}(f)\Omega_P$. The positive two-current
spectral marginal is the quadratic spectral measure

```math
d\rho_Q^{(2)}(s;P)
=
\frac13\sum_{\lambda=1}^{3}
\langle\psi_\lambda,
E_{\rm had}(ds;P)\psi_\lambda\rangle.
```

If $\mathcal H_{\rm had}(P)$ exists and the smeared current states lie in
it, then
$d\rho_Q^{(2)}$ is a positive locally finite measure supported on the
confined hadronic spectrum. The smearing/subtraction is essential because a
point current acting on the vacuum is an operator-valued distribution and can
carry ultraviolet divergences.

## Stieltjes/Jacobi Export

The simulator should not fit a free continuous spectrum as the primary object.
It should export a positive Stieltjes representation. Raw positive moments
need not exist in an ultraviolet QFT, so use either a declared cutoff:

```math
\mu_n^{(\Lambda)}(P)=\int_{s_{\rm th}}^{\Lambda^2}
s^n\,d\rho_Q^{(2)}(s;P),
```

or finite inverse moments such as
$\nu_n(Q_0^2)=\int(s+Q_0^2)^{-n-1}d\rho(s)$, with the subtraction and
kernel stated. For any resulting finite positive moment sequence, the
corresponding Hankel matrices must be positive semidefinite.
The endpoint export is

```math
\left(
J_{24,Q}(P),
\omega_Q(P),
\Xi_Q(P)
\right).
```

$J_{24,Q}$ must come from the moment sequence or a Lanczos/Stieltjes
procedure. If its nodes and weights are chosen after seeing $\alpha(0)$, the
object is a back-solve.

## Same-Scheme Endpoint Remainder

The remainder is a mutually exclusive same-scheme ledger:

```math
\Xi_Q(P)
=
\Xi_{\rm sub}
+
\Xi_{\rm EW}
+
\Xi_{\rm QED}
+
\Xi_{\rm isospin}
+
\Xi_{\rm heavy}
+
\Xi_{\rm FV}
+
\Xi_{\rm cont}
+
\Xi_{\rm contact}.
```

All terms must be in the same electromagnetic convention as $A_Z(P)$ and
$A_{\rm Th}(P)$. `QED` and `EW` here may contain only corrections not
already included in $\Delta_{\rm lep}$, $\Delta_{\rm EW}$, the
spectral kernel, or current renormalization; otherwise the endpoint is double
counted. Each row needs an inclusion/exclusion definition and units.

## Fine-Structure Endpoint Closure

If OPH emits $d\rho_Q^{(2)}(s;P)$ and $\Xi_Q(P)$, then the following
declared spacelike kernel is a conditional endpoint map:

```math
\Delta_{\rm had}(P)
=
\frac{m_Z(P)^2}{3\pi}
\int_{s_{\rm th}}^\infty
\frac{d\rho_Q^{(2)}(s;P)}
{s[s+m_Z(P)^2]}
+
\Xi_Q(P)
```

is source-side under the normalization used to define $d\rho$. The plus sign
in $s+m_Z^2$ makes this a spacelike Euclidean kernel. A timelike
$Z$-pole quantity requires the corresponding analytic continuation,
principal-value/resonance treatment, and scheme conversion; it cannot be
silently identified with this integral. If

```math
A_{\rm Th}(P)
=
A_Z(P)
+
\Delta_{\rm lep}(P)
+
\Delta_{\rm had}(P)
+
\Delta_{\rm EW}(P)
```

and $G(P)=\varphi+\sqrt\pi/A_{\rm Th}(P)$ is a contraction on the declared
interval, then the pixel fixed point $P_\star$ is unique and
$\alpha(0)=1/A_{\rm Th}(P_\star)$.

## Full Hadronic Precision Functor

The common interface must preserve object type. Define schematically

```math
\mathscr R_{\rm had}^{\rm OPH}
:
(J_1,\ldots,J_n;\mathcal O_i;H_i\to H_f)
\mapsto
\left(
\Pi^{(2)},
\Gamma^{(4)},
\mathcal A^{\rm tr}_{H_i\to H_f;i},
\Xi
\right).
```

Required special cases:

```math
\mathscr R_{\rm had}^{\rm OPH}(J_Q,J_Q)=\Pi_Q^{\mu\nu}
\quad\leftrightarrow\quad d\rho_Q^{(2)},
```

```math
\mathscr R_{\rm had}^{\rm OPH}(J_Q,J_Q,J_Q,J_Q)=\Gamma_{QQQQ}^{\mu\nu\rho\sigma},
```

```math
\mathscr R_{\rm had}^{\rm OPH}(\mathcal O_i^{b\to s},J_Q;B\to K^{(*)})
\quad=\mathcal A_{B\to K^{(*)};i}^{\rm tr},
```

```math
\mathscr R_{\rm had}^{\rm OPH}(\mathcal O_i^{\Delta S=1},J_Q;\Sigma\to p)
\quad=\mathcal A_{\Sigma\to p;i}^{\rm tr}.
```

The four-current and transition objects are tensor-valued correlators or
complex amplitudes, not generally one positive scalar spectral measure. The
scalar insufficiency statement is immediate: the two-current measure is
sufficient for running $`\alpha`$ hadronic transport and HVP $g-2$, but it is
insufficient for HLbL $g-2$ and rare-decay long-distance amplitudes.

## Simulator Claim Tiers

| Tier | Meaning |
| --- | --- |
| H0 | raw diagnostic run, no physical claim |
| H1 | conventional QCD reference ensemble |
| H2 | OPH-native quotient ensemble declared |
| H3 | OPH finite vacuum promoted |
| H4 | continuum/confined QCD certificate candidate |
| H5 | Ward two-point spectral source interval |
| H6 | full hadronic precision functor interval |
| H7 | source-only precision prediction beating direct comparison interval |

These names are promotion gates, not evidence that a tier has been reached.
This note itself establishes no tier above an H2-style declared source
prototype. The first implementation milestone is `HVP_ALPHA_SOURCE_PROTOTYPE`, not rare
decays. It should output $C_{QQ}(t)$, $d\rho_Q^{(2)}$ or Stieltjes bounds,
$J_{24,Q}(P)$, $\omega_Q(P)$, $\Xi_Q(P)$, and $\Delta_{\rm had}(P)$,
then compare only after freezing.

## Minimum Receipt Bundle

The simulator backend must emit a directory with:

- `manifest.json`
- `source_dag.json`
- `qcd_ensemble/quotient_schema.json`
- `qcd_ensemble/gamma_groupoid.json`
- `qcd_ensemble/base_measure.json`
- `qcd_ensemble/source_action.json`
- `qcd_ensemble/source_parameter_map.json`
- `qcd_ensemble/coarse_maps.json`
- `vacuum/euclidean_slab.json`
- `vacuum/transfer_operator.json`
- `vacuum/reflection_positivity.json`
- `vacuum/vacuum_promotion.json`
- `currents/ward_current_definition.json`
- `currents/current_normalization_ZV.json`
- `currents/contact_terms.json`
- `currents/ward_residuals.csv`
- `correlators/vector_current_2pt_raw.json`
- `correlators/vector_current_2pt_covariance.json`
- `correlators/disconnected_diagrams.json`
- `correlators/autocorrelation_report.json`
- `spectral/moments.json`
- `spectral/hankel_positivity.json`
- `spectral/stieltjes_bounds.json`
- `spectral/J24Q.json`
- `spectral/omegaQ.json`
- `spectral/spectral_interval.json`
- `endpoint/kernel_definition.json`
- `endpoint/Xi_same_scheme.json`
- `endpoint/Delta_had_interval.json`
- `endpoint/ATh_interval.json`
- `endpoint/pixel_contraction_interval.json`
- `higher_point/Q4_HLbL_receipt.json`
- `higher_point/transition_B_to_K_receipt.json`
- `higher_point/transition_Sigma_to_p_receipt.json`
- `controls/no_target_leak_dag.json`
- `controls/empirical_data_exclusion_manifest.json`
- `controls/frozen_code_hashes.json`
- `controls/replay_receipts.json`
- `controls/comparison_data_manifest.json`
- `claim.md`

`claim.md` must say exactly one of:

- `CONVENTIONAL_QCD_REFERENCE`
- `SOURCE_PROTOTYPE_NOT_PROMOTED`
- `SOURCE_INTERVAL_PROMOTED`
- `EMPIRICAL_CLOSURE_ONLY`
- `COMPARISON_ONLY`

## Forbidden Source Inputs

The source lane must exclude:

- `CODATA_ALPHA`
- `MUON_G_MINUS_2`
- `EE_TO_HADRONS`
- `RARE_DECAY_DATA`
- `HADRON_MASS_TARGETS`
- `PDG_QCD_FITS`

If any forbidden target enters the source DAG, the claim must fail closed as
`TARGET_LEAK_DETECTED`.

## Open Theorem

Construct, from OPH source data alone, a refinement-compatible family

```math
\left(
\mathfrak Q_r^{\rm QCD}(P),
\mathfrak S_r^{E,{\rm QCD}}(P),
J_{Q,r}^{W,R,\mu}(P),
c_{sr}
\right)
```

such that the quotient ensemble law is declared, the QCD parameter map is
source-emitted, the Euclidean slab promotes to a physical vacuum, the continuum
limit gives an unquenched confined QCD Hilbert quotient, the electromagnetic
current is Ward-projected and endpoint-normalized, $d\rho_Q^{(2)}$ and
$\Xi_Q$ are emitted with computable errors, $J_{24,Q}$ is derived from
source correlators rather than fitted to $\alpha(0)$, the higher-point and
transition extensions are emitted, and the no-target-leak DAG has no path from
the forbidden comparison data.

Until this theorem is proved or constructively certified, the hadronic
precision endpoint remains source-open.
