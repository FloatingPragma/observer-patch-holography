# Fractional Excitons As OPH Quotient-Sector Readouts

## Motivating Result

This note entered the queue after the 2026 anyon-trion photoluminescence report
in twisted MoTe2. The surprising claim gave photoluminescence access to a trion
apparently binding to a fractional charge inside a fractional Chern-insulator
state
([Nature, 2026](https://www.osti.gov/pages/biblio/3021156)). A second optical
report identified signatures of a fragile `nu=-1/3` fractional quantum
anomalous Hall state in the same material family
([Phys. Rev. Lett. 136, 056601, 2026](https://arxiv.org/abs/2602.04561)).
The OPH question is how a photon-visible bound state can be attached to a
fractional topological sector without turning every optical line into a phase
identification.

**Status:** solved as a conditional fractional-material sandbox theorem package
and simulator receipt contract. Material-specific promotion remains blocked at
the frozen Hamiltonian/source-law and phase-certificate receipts.

Date: 2026-07-08

## What This Note Contributes

The OPH fractional Hall writeup proves the Hall-collar quotient, the Abelian $(K,t)$ normal form, the non-Abelian repair-sector criterion, and the normal-form non-selection theorem. That is enough for fractional Hall transport as an OPH normal-form statement. It is not enough for fractional excitons in a real material, because the optical experiment asks a sharper question: which fractional topological sector is coupled to which photon-visible bound state?

This note adds the missing objects. A material presentation $\mathcal F_{x,r}$ defines the finite quotient. A source law or Hamiltonian supplies the selector that the normal form cannot supply. A topological ledger records charge, statistics, fusion, and edge data. An optical module records how neutral and charged optical complexes sit over that ledger. The line fan is the public optical readout. The simulator certificate says which finite numerical artifact may be used as evidence.

## Claim Boundary

The result is a theorem package for a sandbox branch. It says how OPH reads a fractionalized material once the material presentation, source law, and certificates are supplied. It does not claim that an arbitrary moire sample, Chern band, fractional Chern insulator, fractional quantum anomalous Hall device, or optical line belongs to a specified topological phase.

A material-specific proof requires public evidence for the source Hamiltonian, active band projector, Chern number, band geometry, many-body gap, ground-sector degeneracy, flux-insertion pump, Hall conductance, edge spectrum, sector ledger, refinement stability, and no-target-leak audit. Without those receipts, the branch remains a diagnostic analogy and not a material theorem.

## Fractional Material Presentation

**Definition 1** (Fractional material presentation). *At regulator $r$, a fractional material presentation for a sample or simulation target $x$ is
``` math
\mathcal F_{x,r}
=
\left(
\Sigma^{\mathrm{mat}}_{x,r},
\Gamma_{x,r},
Q_{x,r},
\mathcal R_{x,r},
\mathcal U_{x,r},
\mathrm{Chk}_{x,r}
\right).
```
$\Sigma^{\mathrm{mat}}_{x,r}$ is the finite presentation space: lattice, geometry, active bands, filling, interactions, disorder or strain fields, electromagnetic ports, optical ports, edge/collar charts, and finite Hilbert-space cutoffs. $\Gamma_{x,r}$ quotients gauge representatives, basis choices, mesh labels, orbital labels, edge charts, repair schedules, hidden carrier coordinates, and inert ancillas. The physical quotient is
``` math
Q_{x,r}=\Sigma^{\mathrm{mat}}_{x,r}/\Gamma_{x,r}.
```
$\mathcal R_{x,r}$ is the public record algebra, $\mathcal U_{x,r}$ is the allowed repair/update family, and $\mathrm{Chk}_{x,r}$ is the refinement and receipt data.*

**Theorem 2** (Material quotient normal form). *Suppose the accepted repair relation on $Q_{x,r}$ is terminating, quotient-descended, locally confluent, and repair-complete on the declared branch. Then the material presentation has a schedule-independent normal form
``` math
n_{x,r}:Q_{x,r}\to N_{x,r}.
```
The components of $N_{x,r}$ are exactly the quotient-visible edge/collar algebra, holonomy or higher-holonomy data, sector ledger, optical-module data, and public record readout.*

*Proof.* The finite consensus theorem applies to the physical quotient, not to hidden representatives. The four stated repair hypotheses give a unique normal form in each repair-connected quotient class. Hidden labels collapse under $\Gamma_{x,r}$. Edge/collar records, holonomy defects, fusion data, optical couplings, and public records survive exactly when they change the quotient readout. $\square$

## Source Law And Non-Selection

The canonicalizer does not choose a material phase. A material branch must declare a source law
``` math
\mu_{x,r,T}(q)
=
Z^{-1}_{x,r,T}m_{x,r}(q)\exp[-S_{x,r,T}(q)]
```
or a quantum Gibbs or ground-projector source built from a frozen Hamiltonian. The required tags are
``` math
\mathrm{SOURCE\_LAW\_REQUIRED},
\qquad
\mathrm{NORMAL\_FORM\_IS\_NOT\_SELECTOR}.
```

**Theorem 3** (Normal-form non-selection). *Let $n_{x,r}:Q_{x,r}\to N_{x,r}$ be a material normal-form map, and let $\mathcal X$ be a finite set of candidate topological sectors inside $N_{x,r}$. The map $n_{x,r}$ does not select a unique material sector unless the branch supplies a quotient-intrinsic source law, Hamiltonian, transfer operator, or vacuum certificate.*

*Proof.* The pushforward $n_{\#}\mu$ is idempotent on laws supported on normal forms. Every law concentrated on a candidate sector is fixed by the canonicalizer. Therefore the canonicalizer can classify candidate sectors, but it cannot choose their material weights. Selection requires the source object. $\square$

## Minimal Hamiltonian And Promotion Certificate

A practical fractional-Chern sandbox may use a projected interacting Chern-band Hamiltonian
``` math
H^{\mathrm{FCI}}_{x,r}
=
P_C\left[
\sum_{\bm k} \varepsilon_{\bm k} c^\dagger_{\bm k}c_{\bm k}
+
\frac12\sum_{\bm q}V(\bm q)\rho_{-\bm q}\rho_{\bm q}
+
H_{\mathrm{dis/strain/gate}}
\right]P_C ,
```
with optical extension
``` math
H^{\mathrm{opt}}_{x,r}
=H^{\mathrm{FCI}}_{x,r}+H_{\mathrm{hole}}+H_{\mathrm{light}}+H_{\mathrm{bind}}.
```
This is a schema, not a material proof. The branch must freeze the source before comparison and emit a promotion certificate
``` math
\mathrm{PhaseCert}_{x,r}
=
\left(
C,\nu,\Delta,G,\sigma_{xy},\mathrm{Pump},
\mathrm{Edge},\mathrm{Ent},\mathrm{Anyon},\mathrm{Modular},
\mathrm{Refine}
\right).
```

The receipt names are:
``` math
\begin{gathered}
\mathrm{SOURCE\_HAMILTONIAN\_FROZEN},\quad
\mathrm{ACTIVE\_BAND\_PROJECTOR},\quad
\mathrm{CHERN\_NUMBER},\\
\mathrm{BAND\_GEOMETRY},\quad
\mathrm{MANYBODY\_GAP},\quad
\mathrm{GROUND\_SECTOR\_DEGENERACY},\\
\mathrm{FLUX\_INSERTION\_PUMP},\quad
\mathrm{HALL\_CONDUCTANCE},\quad
\mathrm{EDGE\_SPECTRUM},\\
\mathrm{TOPOLOGICAL\_SECTOR\_LEDGER},\quad
\mathrm{REFINEMENT\_STABILITY},\quad
\mathrm{NO\_TARGET\_LEAK}.
\end{gathered}
```

**Theorem 4** (Hamiltonian-to-ledger promotion). *Let $H_{x,r}$ be frozen before comparison, and suppose $\mathrm{PhaseCert}_{x,r}$ supplies the receipts above. Let $\mathcal X_{x,r}$ be the candidate sector set and let
``` math
\pi_{x,r}:\mathcal X_{x,r}\to \mathcal L_{x,r}
```
map each candidate to its public topological ledger. If $\pi_{x,r}$ is injective on the candidate set, then the Hamiltonian promotes the material branch to the unique ledger $\mathcal L_{x,r}$ selected by the certificate. If $\pi_{x,r}$ is not injective, the correct simulator output is $\mathrm{SECTOR\_AMBIGUOUS}$.*

## Topological Ledger And Fractional Readout

On the Abelian branch, the ledger contains an integral pair $(K,t)$. For a quasiparticle vector $\ell$,
``` math
\nu=t^T K^{-1}t,\qquad
Q_\ell/e=t^TK^{-1}\ell,
```
``` math
\theta_\ell=\pi\ell^TK^{-1}\ell,\qquad
\theta_{\ell,\ell'}=2\pi\ell^TK^{-1}\ell' .
```
On the non-Abelian branch, the ledger carries
``` math
\mathcal C_r=(\mathrm{Irr},N,F,R,S,T,d,\theta),
```
together with the local electron object and electromagnetic charge grading.

**Theorem 5** (Fractional charge and statistics are quotient-sector readouts). *For a certified fractional Hall or fractional Chern branch, fractional charge, exchange phase, mutual braiding, fusion, chiral central charge, and edge response are functions of the quotient-sector ledger. They are not properties of a bare local representative in $\Sigma^{\mathrm{mat}}_{x,r}$.*

*Proof.* Gauge, basis, mesh, orbital, and hidden carrier relabelings can change local representatives without changing the public sector. The fractional observables above are invariant under those relabelings and are read at the edge, through holonomy, interferometry, transport, flux pump, or modular data. Therefore they factor through $Q_{x,r}$ and its normal form. $\square$

## Optical Module And Line Fan

**Definition 6** (Optical module). *An optical fractional-exciton ledger is a module category
``` math
\mathcal M_{\mathrm{opt}}^{x,r}
```
over the topological category $\mathcal C_r$. Each optically active sector $m$ has a topological shadow $\tau(m)\in \mathrm{Irr}(\mathcal C_r)$, total electromagnetic charge $Q_{\mathrm{tot}}(m)$, oscillator data, polarization data, and a binding energy term.*

A neutral fractional exciton may have
``` math
Q_{\mathrm{tot}}(m)=0,
\qquad
\tau(m)\ne 1 .
```
It is neutral in the electromagnetic slope channel but still fractional in the topological ledger.

**Definition 7** (Optical line fan). *The optical line fan is
``` math
\mathcal L_{\mathrm{opt}}^{x,r}
=
\{(E,I,\partial_g E,\mathrm{pol},\tau,Q_{\mathrm{tot}},\eta)\}.
```
It records line energy, intensity, gate or field slope, polarization, the topological shadow, total charge, and residual identification metadata.*

**Theorem 8** (Line-fan decomposition). *Let the optical operators be certified as quotient-descended maps from $\mathcal M_{\mathrm{opt}}^{x,r}$ to public optical records. Then the sector spectral function decomposes as a sum over optical-module sectors:
``` math
A_{\mathrm{opt}}(\omega,g)
=
\sum_{m\in \mathrm{Irr}(\mathcal M_{\mathrm{opt}}^{x,r})}
I_m(g)\,\delta_\eta(\omega-E_m(g)).
```
The public line fan is the quotient readout of this decomposition.*

**Theorem 9** (Fractional optical slope boundary). *A gate or field slope identifies total charge only when the derivative of the binding energy is independently bounded. In symbols, if
``` math
\left|\partial_g E_{\mathrm{bind}}(m)\right|\le b_m ,
```
then a measured slope determines $Q_{\mathrm{tot}}(m)$ only up to that bound. For a neutral fractional exciton, the leading charge slope can vanish while $\tau(m)\ne 1$.*

*Proof.* The measured slope is the derivative of the whole optical energy, including terms beyond the electromagnetic charge contribution. Binding drift can mimic or mask a charge slope. When its derivative is bounded, the residual interval gives the certified charge range. When the total charge is zero, the charge slope can vanish even though the sector shadow remains topologically nontrivial. $\square$

**Theorem 10** (Optical identifiability). *Let
``` math
\mathrm{ID}_{x,r}:\mathrm{Irr}(\mathcal M_{\mathrm{opt}}^{x,r})\to
(E,I,\partial_g E,\mathrm{pol},\tau,Q_{\mathrm{tot}},\eta)
```
be the public optical identifier after all stated uncertainties. The optical line selects a unique sector iff $\mathrm{ID}_{x,r}$ is injective on the candidate sectors compatible with the frozen source and topological ledger. If it is not injective, the correct output is $\mathrm{OPTICAL\_SECTOR\_AMBIGUOUS}$.*

## Simulator Quotient Correctness

The simulator is evidence only if it implements the declared quotient. It must not define missing objects after the target comparison has been inspected.

**Theorem 11** (Simulator quotient correctness). *Let a simulator store representatives $s\in\Sigma^{\mathrm{mat}}_{x,r}$, a quotient map $q:\Sigma^{\mathrm{mat}}_{x,r}\to Q_{x,r}$, a canonicalizer $c$, observables $O_a$, transition kernels $K$, source law $\mu$, and refinement maps $c_{sr}$. The simulation is quotient-correct for the fractional branch only if it emits:
``` math
\begin{gathered}
\mathrm{CANONICALIZER\_IDEMPOTENCE},\quad
\mathrm{REPRESENTATIVE\_INVARIANCE},\quad
\mathrm{QUOTIENT\_LUMPABILITY},\\
\mathrm{DETAILED\_BALANCE\_OR\_DECLARED\_NONEQUILIBRIUM},\quad
\mathrm{REFINEMENT\_COMPATIBILITY},\\
\mathrm{NO\_ORBIT\_SIZE\_BIAS},\quad
\mathrm{NO\_TARGET\_LEAK}.
\end{gathered}
```
If any required receipt fails, the material conclusion is not promoted.*

*Proof.* Idempotence makes the canonicalizer a normal-form map. Representative invariance makes observables functions on $Q_{x,r}$. Lumpability makes the stochastic or deterministic transition law descend to the quotient. Source freezing and the no-target-leak audit prevent posterior fitting. Refinement compatibility prevents a finite-regulator artifact from being reported as a stable material sector. Orbit-size control prevents the simulator from selecting a sector because it has more hidden representatives. $\square$

Implementation receipts exist on both simulator surfaces. The active simulator helpers live in `oph-physics-sim/oph_fractional/`, with the generated sandbox bundle at `oph-physics-sim/runs/fractional/quotient_sector_sandbox/`. The paper-stack mirror lives in `reverse-engineering-reality/code/particles/fractional/`, with the generated bundle at `reverse-engineering-reality/code/particles/runs/fractional/quotient_sector_sandbox/`. Both bundles stop at `FRACTIONAL_QUOTIENT_SANDBOX_DIAGNOSTIC` and block promotion at `MATERIAL_SPECIFIC_HAMILTONIAN_PROOF_RECEIPT`.

## Failure States

| State | Meaning |
| --- | --- |
| `SOURCE_NOT_FROZEN` | comparison occurred before the material source was fixed |
| `NOT_QUOTIENT_INVARIANT` | an observable changed under hidden representative relabeling |
| `CANONICALIZER_NOT_IDEMPOTENT` | repeated canonicalization changed the normal form |
| `KERNEL_NOT_LUMPABLE` | transition probabilities failed to descend to quotient sectors |
| `ORBIT_SIZE_BIAS_DETECTED` | hidden representative multiplicity affected sector weights |
| `NO_GAP_CERTIFICATE` | the many-body gap was not certified |
| `CHERN_NUMBER_UNSTABLE` | the active-band Chern number failed stability checks |
| `PHASE_CERTIFICATE_NONINJECTIVE` | the material certificate did not select a unique ledger |
| `SECTOR_AMBIGUOUS` | transport or topological data admit more than one sector |
| `OPTICAL_OPERATOR_UNCERTIFIED` | optical operators were not shown to descend to the quotient |
| `BINDING_DRIFT_UNBOUNDED` | optical slopes cannot be interpreted as charge readouts |
| `OPTICAL_SECTOR_AMBIGUOUS` | optical line-fan identifiers are not injective |
| `TARGET_LEAK_DETECTED` | measured target data entered the source or selector |
| `REFINEMENT_DEFECT_TOO_LARGE` | finite-regulator results failed refinement stability |
| `DIAGNOSTIC_ONLY` | the result is useful for visualization or debugging but not a material claim |

## First Experimental Reading

Twisted transition-metal dichalcogenide bilayers are a natural first sandbox because the same finite material platform can expose moire bands, Chern-band transport, edge or pump signatures, and optical excitonic spectra. The OPH prediction is not one fitted peak. It is a sector-indexed line fan: the optical peaks should organize by the same quotient-sector ledger that explains charge, edge, and transport. Charged anyon-trion branches can have nonzero slopes after binding drift is bounded. Neutral fractional-exciton branches can have $Q_{\mathrm{tot}}=0$ and $\tau\ne1$.

The direct falsifier is equally simple. If no single frozen source and no single ledger $(\mathcal C,\chi,\Theta)$ can index the transport, edge, and optical records without leakage or ambiguity, this OPH material branch fails.

## Outcome

This note closes at note level as a conditional theorem package. The OPH answer is that fractional excitons and related optical peaks are quotient-sector readouts when, and only when, the material source, topological ledger, optical module, line fan, identifiability map, and simulator certificate all descend to the same finite quotient. The simulation's job is receipt production: gaps, Chern numbers, sector data, optical spectra, refinement defects, no-target-leak audits, quotient lumpability, and identifiability. It is not a substitute for the material proof.
