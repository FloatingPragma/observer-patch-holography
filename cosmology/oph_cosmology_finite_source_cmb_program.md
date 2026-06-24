---
abstract: |
  This paper stages the OPH cosmology branch inside the public core repository without yet promoting it to the published paper set. It consolidates the inflation-free observer-screen synchronization draft, the dark/anomaly stress paper, the screen-microphysics paper, the compact OPH theorem stack, and the latest OPH-FPE simulator diagnostics into one finite-source prediction program. The central claim is deliberately conditional. OPH may replace the logical jobs of inflation by finite observer-screen synchronization: flatness by visible holonomy selection, horizon coherence by a same-boundary or low-$k$ repair gate, near scale invariance by a screen Green spectrum, hot initial data by MaxEnt release, and acoustic peaks by ordinary Boltzmann transfer from source-only finite OPH artifacts. The latest 256k-patch OPH-FPE run strengthens the program because its best diagnostic CMB curve has shape correlation $0.9951542364$, normalized RMSE $0.0984455548$, and zero mean peak-location error against the local binned TT table. It also hardens the claim boundary: the physical CMB prediction receipt remains false until finite source provenance, physical scale calibration, a covariant collar-packet stress parent, physical $\rho_A(a)$, $B_A(k,a)$, and $\Gamma_{\rm rec}(k,a)$ kernels, Boltzmann transfer, and frozen likelihood receipts pass.
author:
- B. Müller
title: |
  OPH Cosmology as a Finite-Source Prediction Program:\
  Observer-Screen Synchronization, Dark Stress, and the Physical CMB Gate
---

# What This Paper Contributes {#what-this-paper-contributes .unnumbered}

The OPH core papers recover local observer-facing geometry, records, modular flow, gauge/matter structure, and finite screen microphysics. The cosmology work asks a different question: which large-scale data can be predicted once the finite observer-screen state is released into an Einstein-Boltzmann cosmology?

This staging paper contributes four things.

1.  It consolidates the inflation-free OPH cosmology branch into the core repository under `cosmology/`, while keeping it outside the release build until the theorem gates close.

2.  It separates diagnostic simulator outputs from physical predictions. A good visual or measurement-facing CMB overlay is not a physical CMB prediction unless its inputs are source-only finite OPH artifacts fixed before likelihood comparison.

3.  It names the finite objects that must exist before simulation work should resume on physical CMB, vacuum, or cosmological uses of the dark sector: source provenance, physical mode calibration, freezeout, covariant stress, anomaly kernels, transfer, likelihood, and quotient ensemble receipts. The dark-sector derivations themselves remain in `extra/`.

4.  It gives the open GitHub theorem issues a single paper-side target.

# Claim Boundary

This is a continuation paper, not a release checkpoint and not a published OPH paper bundle member. It is part of the public core repository so that later paper and simulator work can target one coherent cosmology contract.

The branch claim is: $$\boxed{
\begin{gathered}
\text{OPH can replace the logical jobs of inflation on a finite-source branch:}\\
\text{flatness by visible-holonomy selection, horizon coherence by same-boundary repair,}\\
\text{near scale invariance by the screen Green spectrum, hot initial data by MaxEnt release,}\\
\text{and acoustic peaks by Boltzmann transfer from source-only OPH artifacts.}
\end{gathered}}$$

Here "can" means "can on the declared Phase III finite-source continuation branch if the listed source, stress, lift, transfer, and likelihood gates close." It does not mean that inflation replacement, CMB spectra, or dark-sector cosmology are part of the recovered Phase I SM/GR core.

The branch does not currently claim a physical CMB prediction, a completed Planck/ACT likelihood run, an OPH-native quantum vacuum, a full cosmological dark-matter likelihood, or a production particle-formation simulation. The latest simulator run is measurement-facing and useful, but its own promotion receipts keep the physical claim closed.

::: remark
*Remark 1* (No simulator shortcut). The paper-side theorem gates are not bureaucratic checks. They are the distinction between a simulation that happens to resemble CMB data and a prediction whose input artifacts were fixed without using CMB data. The simulator should not resume physical CMB promotion until those gates are mathematically specified here and then implemented as deterministic receipts.
:::

# Compatibility With the Current Core Stack

This section records the compatibility pass against the current core papers after the recent issue closures. It is the guardrail for importing the older inflation-free draft.

  Imported draft idea                      Current core-compatible reading
  ---------------------------------------- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Inflation-free replacement               Phase III continuation only. The compact paper's recovered core does not include an inflation replacement, physical CMB likelihood, $H_0/S_8$ branch, or dark-sector cosmology claim.
  Finite normal form                       Use quotient-level fixed-cutoff normal forms on the declared branch, unique only modulo boundary redundancy, implementation hiding, inert ancillary stabilization, and any required same-boundary unique-extension condition. Do not read the old draft as selecting a unique microscopic representative.
  Screen Green spectrum and $n_s=1-P/48$   Screen-level theorem or target until the source-only primordial bridge passes. Physical $A_s$, $n_s$, running, isocurvature, phase coherence, and TT/TE/EE spectra require the source-stress, single-clock, repair-gap, freezeout, radial-prior, null-space, and forward-residual receipts.
  Dark/anomaly slot                        Imported from `extra/`. For a source-only primordial certification run, the dark continuation is `OFF` unless a supplied abundance is explicitly typed `CONDITIONAL_SOURCE_STATE`; transported $Q_A$ is not by itself an OPH derivation of the homogeneous anomaly abundance.
  Scalar anomaly rows                      Rows for $\bar\rho_A$, $\bar\rho_{A,\mathrm{eq}}$, and $B_A$ are diagnostics unless the finite covariant parent also emits $w_A$, $c_{s,A}^2$, $\sigma_A$, $Q_A^\mu$, $\Gamma_{\rm rec}$, recipient stress when $\Gamma_{\rm rec}>0$, gauge-independent variables, CDM-limit recovery, and refinement convergence.
  Boltzmann and likelihood comparison      Physical promotion requires immutable source, solver, and likelihood hashes, official likelihood execution, and global CDM-limit reductions. Shard-local `any()` rollups or nonlinear averages before global pooling do not satisfy the contract.
  Vacuum language                          Finite reference states and free-field/lattice baselines remain reference ensembles. OPH-native vacuum promotion requires the quotient-ensemble selector and positive-transfer or reflection-positive gate from the compact/screen-microphysics theorem surface.
  Simulator receipts                       Any future receipt cited here must carry exact source/config/solver/likelihood hashes. The 256k run used for diagnostics records `git_commit=unknown`, so it is not used as exact replay provenance.

# Source Corpus

This paper consolidates these active surfaces.

  Surface                                                                                                    Imported role
  ---------------------------------------------------------------------------------------------------------- ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  `cosmology/inflation_without_inflaton_observer_screen_synchronization.tex`                                 Main inflation-free branch: finite screen federation, flatness, horizon coherence, screen Green spectrum, amplitude, MaxEnt release, anomaly slot, transfer, and likelihood gates.
  `cosmology/physical_cmb_theorem_program.md`                                                                Recent theorem-roadmap file identifying the missing physical CMB source, scale, stress, kernel, transfer, and likelihood contracts.
  `paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex`   Compact recovered-core claim boundary and the first surface for later propagation.
  `paper/tex_fragments/PAPER.tex`                                                                            Main technical paper surface for global OPH theorem language.
  `paper/tex_fragments/PRIMORDIAL_BRIDGE_THEOREMS.tex`                                                       Primordial bridge theorem status and screen-to-primordial finite-source conditions.
  `paper/screen_microphysics_and_observer_synchronization.tex`                                               Finite observer-screen carriers, records, patch ports, evidence bundles, and implementation invariance.
  `extra/oph_dark_matter_paper.tex`                                                                          Collar-remainder dark/anomaly stress, galaxy equilibrium, dynamic stress, and cosmological kernel contracts. This paper imports those contracts by reference; it does not duplicate the dark-sector derivations.
  `oph-physics-sim/runs/gcp_large_vis_256k_modkernel_20260623`                                               Latest 256k-patch simulator evidence: strong CMB diagnostics, explicit false physical-prediction receipt, and full visualizer export data. The run manifest records `git_commit=unknown`; the canonical simulator repository is <https://github.com/muellerberndt/oph-physics-sim>, and the current public source reference available during this consolidation is `aeea9e502491a277b5e21a80bb89df6d089a074d`.

::: remark
*Remark 2* (Simulator receipt provenance). Whenever this paper mentions simulator receipts, the receipt artifacts should be traced to the responsible commit in the canonical simulator repository, <https://github.com/muellerberndt/oph-physics-sim>. For the 256k run used here, the local artifact `runs/gcp_large_vis_256k_modkernel_20260623/manifest.json` reports `git_commit=unknown`. Therefore the numerical run is cited as diagnostic evidence with incomplete commit provenance. The available public source context at the time of this consolidation is <https://github.com/muellerberndt/oph-physics-sim/commit/aeea9e502491a277b5e21a80bb89df6d089a074d>, but that link is not asserted to be an exact run-replay hash for the 256k artifacts. Future physical-promotion receipts must carry exact source, config, solver, and likelihood hashes.
:::

::: remark
*Remark 3* (Ownership boundary for extra papers). Published and staging papers under `extra/` keep ownership of the dark-sector theory: collar-remainder stress, galaxy equilibrium, dynamic repair-stress transport, cluster behavior, and the detailed perturbation-kernel derivations. This cosmology program paper only records the interface those papers must expose to the physical CMB/large-scale-structure pipeline. When it names $\bar\rho_A(a)$, $B_A(k,a)$, $\Gamma_{\rm rec}(k,a)$, stress closure, or repair-stress transport, it is naming imported source functions and receipts, not rederiving the dark sector here.
:::

# Imported Theorem Map

The long inflation-free draft contains many theorem statements. This paper does not copy every proof verbatim. Instead it gives them one public-core landing surface and keeps their dependency roles visible.

  Imported theorem family                                                                                                                                    Role in this paper
  ---------------------------------------------------------------------------------------------------------------------------------------------------------- ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Finite synchronization normal form                                                                                                                         Defines the early observer-facing quotient state whose repair normal form replaces an assumed smooth primordial substrate, with uniqueness read modulo the current compact-paper quotient equivalences.
  FLRW curvature holonomy, flat-sector selection, and curvature damping                                                                                      Owns the flatness branch: exact $k=0$ when zero visible holonomy is selected, or a bounded residual when the curvature-repair margin stays positive.
  Diffusive horizon no-go, same-boundary coherence, and low-$k$ repair gap                                                                                   Separates a failed purely local diffusion story from the allowed same-boundary or uniform low-$k$ synchronization mechanisms.
  Screen Green spectrum, Mellin lift, edge-center scalar opportunities, $\mathbb Z_6$ reserve, and half-collar sampling                                      Owns the screen-level near-scale-invariant covariance and the red tilt target $n_s=1-P/48$; physical primordial status waits for the source-only bridge receipts.
  Finite release collar state, scalar release code, and scalar amplitude theorem                                                                             Owns $A_\zeta$ as a finite source artifact rather than a fitted amplitude.
  MaxEnt release, entropy transport, BBN-safe release, recombination inheritance, and adiabaticity                                                           Connects the synchronized finite screen/collar state to a hot radiation branch with ordinary low-energy thermal history.
  Freezeout, scalar quotient record, gauge-invariant $\zeta$, isocurvature decay, and superhorizon conservation                                              Imported only through the stricter primordial bridge: total stress closure, single clock, repair gap, freezeout, growing mode, isocurvature, phase coherence, radial prior, null-space, and forward residual.
  Finite homogeneous anomaly charge, no double counting, packet moments, cold anomaly, Bianchi exchange, no-slip lensing, and linear anomaly transfer        Imported from the dark-sector papers in `extra/`; this paper only records the CMB/LSS handoff interface and receipt gates.
  Universe Simulation certificate bundle, no-data-use theorem, finite-certificate branch, Boltzmann stability, likelihood promotion, and CMB admissibility   Owns the promotion rule: source-only finite artifacts first, frozen transfer and likelihood second, physical CMB claim last.
  Tensor repair gap and finite-collar non-Gaussianity                                                                                                        Records secondary predictions to be tested only after the primary scalar/source/stress gates close.

# Extended Cosmology Agenda

The local cosmology notes identify several OPH-native prediction surfaces beyond the immediate CMB gate. They are not independent claims yet; they are the next paper tracks once the finite-source contract is stable.

1.  **CMB as observer-consensus fossil.** The highest-value target is still the CMB: largest observable scales should carry relics of finite observer synchronization, not merely a fitted inflationary initial spectrum.

2.  **$\Lambda$ from observer capacity.** The screen-capacity branch relates the de Sitter scale to finite horizon capacity. A full cosmology paper should derive the equation of state, possible time evolution, and any deviation from $w=-1$.

3.  **Full cosmological dark sector.** The dark paper's galaxy layer is not enough. CMB, BAO, weak lensing, growth, and $S_8$ need the transported repair-stress perturbation branch.

4.  **$H_0$ and $S_8$ predictions.** Once the scale bridge and anomaly kernels are fixed, the branch should emit blind values or intervals for $H_0$ and $S_8$, rather than tuning toward the current tensions.

5.  **Large-scale topology.** Finite observer-accessible screens make global topology and matching-circle style signatures natural tests.

6.  **Consensus hierarchy of structure.** Galaxy, cluster, and supercluster formation should be investigated as a hierarchy of stable consensus attractors, not only as density perturbation growth in a fixed background.

7.  **Cosmological arrow of time.** OPH has native records, checkpoints, and continuation laws; these should be used to formulate the low-entropy past condition as a finite-record boundary condition.

8.  **Cosmological neutrino background.** The particle/mass branch should eventually feed neutrino background and free-streaming predictions into the same frozen likelihood contract.

# Finite Cosmological State

::: definition
**Definition 4** (Early OPH screen federation). *An early OPH cosmological state at regulator $r$ is a finite patch federation $$\mathcal F_r=
\left(
V_r,E_r,\{\mathcal A_i,\rho_i,\mathcal R_i\}_{i\in V_r},
\{\mathcal I_e,\pi_{i,e}\}_{e\in E_r},
\mathcal U_r,\mathrm{Chk}_r
\right),$$ with finite accessible algebras, record algebras, visible overlap interfaces, repair instruments, and checkpoint data. Its physical content is the quotient under hidden carrier coordinates that do not affect visible interface records, repair maps, or checkpoint continuation.*
:::

::: definition
**Definition 5** (Source-only finite artifact). *A finite artifact $X$ is source-only for cosmology if it is a deterministic functional of declared OPH source data, release-branch constants, and standard non-CMB physical constants, and if its dependency manifest excludes CMB, BAO, supernova, weak-lensing, RSD, SPARC, cluster, or other observational likelihood values used to evaluate the prediction.*
:::

::: definition
**Definition 6** (Physical CMB promotion). *An OPH CMB run is physically promoted only when the primordial source, scale bridge, dark/anomaly stress kernels, Boltzmann solver input files, solver version/tolerances, and likelihood code/data hashes are all frozen before likelihood evaluation, and all source functions are source-only finite artifacts.*
:::

# Inflation-Free Branch

The inflation-free branch replaces inflation's jobs rather than reproducing an inflaton. Its load-bearing theorem stack is:

1.  quotient-level observer-facing normal form for the early screen/collar federation, with the current compact-paper quotient equivalences;

2.  zero visible spatial holonomy or positive curvature-holonomy damping for flatness;

3.  same-boundary scalar normal form or a low-$k$ repair gap for horizon coherence;

4.  screen Green spectrum for near scale invariance;

5.  protected $\mathbb Z_6$ half-collar reserve for the red spectral exponent;

6.  finite scalar release code for $A_\zeta$;

7.  MaxEnt release and common release clock for a hot adiabatic start;

8.  imported dark/anomaly stress variables for the pre-recombination CDM-like slot, when the `extra/` dark-sector source contract permits them;

9.  Boltzmann transfer and frozen likelihood for observational comparison.

::: theorem
**Theorem 7** (Conditional OPH replacement of inflation). *Assume a cofinal family of finite early OPH screen federations has a quotient-level observer-facing normal-form projection on a declared branch, with uniqueness understood modulo boundary redundancy, implementation hiding, inert ancillary stabilization, and any same-boundary unique-extension condition required by the branch. Assume visible FLRW curvature holonomy is selected to zero or dynamically damped, every CMB-scale scalar mode shares the same boundary normal form or satisfies a uniform low-$k$ repair gap, the scale-free scalar screen repair cascade realizes the half-collar spectral reserve, the source-only primordial bridge receipts pass, the finite release code emits $A_\zeta$, and the imported dark/anomaly sector supplies a permitted CDM-like source state through recombination. Then the branch supplies flatness, horizon coherence, near scale invariance, hot adiabatic initial data, and acoustic-transfer input without an inflaton degree of freedom.*
:::

::: proof
*Proof.* The proof is a dependency composition. The quotient normal-form projection gives a finite observer-facing early state on the declared branch, not a unique hidden representative. Visible holonomy selection or damping gives the flat sector. Same-boundary repair or a low-$k$ gap gives coherence on the observed CMB band. The screen Green spectrum with half-collar reserve gives the screen-level power law $$\Delta_\zeta^2(k)=A_\zeta\left(\frac{k}{k_\star}\right)^{-P/48},
\qquad
n_s=1-\frac{P}{48}.$$ The source-only bridge receipts type that screen scalar as a primordial curvature source; without them it remains a screen theorem. The finite release code fixes the amplitude before comparison. MaxEnt release gives a hot radiation state and the common clock gives adiabatic leading perturbations. If the imported anomaly stress branch is an admissible source state through recombination, the acoustic peaks are computed by ordinary Einstein-Boltzmann transfer from those finite inputs. ◻
:::

::: remark
*Remark 8* (Current status of the theorem). The theorem is conditional because several premises are still open theorem issues. The branch is sharp enough to simulate, but not sharp enough to promote diagnostic CMB agreement to a physical prediction.
:::

# Physical CMB Finite-Source Contract

The simulator's false physical CMB receipt identified a precise theorem contract. The following objects must be theorem-grade finite artifacts: $$A_\zeta,\qquad q_{\rm IR},\qquad \ell_{\rm IR},\qquad
\bar\rho_A(a),\qquad \bar\rho_{A,\mathrm{eq}}(a),\qquad B_A(k,a),\qquad \Gamma_{\rm rec}(k,a),$$ plus full stress variables, freezeout, physical mode calibration, and $N_{\rm CRC}$ provenance.

::: target
**Target 9** (Finite-source input contract). *The paper stack must define pass/fail receipts for:*

1.  *a source-provenance dependency DAG for all CMB source inputs;*

2.  *globally pooled sufficient-statistic reducers before nonlinear source estimates;*

3.  *the status of $N_{\rm CRC}$ as a consensus invariant, additive capacity, or theorem-side constant;*

4.  *finite derivations of $A_\zeta$, $q_{\rm IR}$, and $\ell_{\rm IR}$;*

5.  *total stress closure, single-clock normal form, entropy repair gap, curvature evolution, adiabatic growing mode, isocurvature bound, phase coherence, screen-to-radial lift, radial null-space, forward residual, and finite freezeout receipts;*

6.  *physical $k$, $\ell$, scale-factor, and redshift calibration;*

7.  *source-only imported dark-sector functions $\bar\rho_A(a)$, $\bar\rho_{A,\mathrm{eq}}(a)$, $w_A(a)$, $c_{s,A}^2(k,a)$, $\sigma_A(k,a)$, $Q_A^\mu$, $B_A(k,a)$, and $\Gamma_{\rm rec}(k,a)$, where the chosen dark-continuation mode permits them.*
:::

::: proposition
**Proposition 10** (Why the 256k CMB curve is not yet a physical prediction). *The 256k OPH-FPE run may be used as a measurement-facing diagnostic, but not as a physical CMB prediction, because its promotion report marks the finite-source input contract false and lists missing source provenance, missing pooled reducers, missing $N_{\rm CRC}$ invariant status, non-finite $A_\zeta$, missing screen-to-primordial lift, non-finite $q_{\rm IR}$ and $\ell_{\rm IR}$, non-finite $\bar\rho_A$ and $B_A$, missing finite covariant parent receipt, missing stress closure, missing full fluid/exchange variables, missing gauge/causal/refinement certificates, missing freezeout, missing official likelihood and CDM-limit reductions, and missing frozen solver and likelihood hashes.*
:::

# Physical Scale Bridge

The screen simulator has natural finite labels: patches, caps, overlaps, spherical harmonics, repair depths, and collar coordinates. A physical CMB prediction needs more. It needs a theorem that maps those labels to physical comoving $k$, multipole $\ell$, scale factor $a$, redshift $z$, and freezeout surface.

::: target
**Target 11** (OPH scale bridge). *Define a finite bridge $$\mathcal S_{\rm OPH}:
(r,\hbox{screen mode},\hbox{collar mode},\hbox{repair clock})
\longmapsto
(k,\ell,a,z,\eta)$$ that is invariant under hidden carrier representatives, compatible with regulator refinement, and independent of post-hoc peak alignment.*
:::

Proxy variables such as inverse cap opening angle can remain useful diagnostics. They cannot clear the physical $k$-gate unless the scale bridge proves that the proxy converges to the physical comoving mode label used by the transfer solver.

# Dark/Anomaly Interface Imported from Extra

The dark/anomaly slot is the largest cosmology-facing bridge between the OPH dark-matter paper and the CMB program. The dark-sector papers in `extra/` own the derivation: collar remainder stress, settled galaxy response, dynamic transport, cluster behavior, and perturbation kernels. This paper only states the interface the CMB program requires from those papers. Cosmology needs the imported dark sector to expose a transported stress component with background density, perturbations, exchange terms, gauge-invariant variables, and the CDM-limit switch-off check.

::: definition
**Definition 12** (Imported finite covariant collar-packet parent interface). *For the purpose of this cosmology program, the finite covariant collar-packet parent is an imported `extra/`-owned finite artifact interface $$\mathcal P_A=
(\mathcal C_r,Z_r,A_r,R_r,G_r,\pi_r,L_r,Q_r,D_r)$$ where $Z_r$ splits anomaly and recipient packet states, $\pi_r$ is the finite equilibrium packet functional, $L_r$ is the repair generator, $Q_r$ records energy-momentum reaction channels, and $D_r$ records any causal auxiliary response. The interface emits packet stress moments $$\rho_A,\quad P_A,\quad q_A,\quad \pi_A,$$ and source variables $$\bar\rho_A(a),\quad \bar\rho_{A,\mathrm{eq}}(a),\quad w_A(a),\quad c_{s,A}^2(k,a),\quad
\sigma_A(k,a),\quad Q_A^\mu,\quad B_A(k,a),\quad \Gamma_{\rm rec}(k,a),$$ with explicit conservation or exchange equations, recipient stress for nonzero exchange, gauge projection, causal response, regulator-refinement convergence, and $B_A$ built from the anomaly-frame baryon density $n_b^{(A)}=-u_{A\mu}J_b^\mu$.*
:::

::: target
**Target 13** (Imported stress-parent receipt). *The simulator receipt for the imported parent must fail unless the owning dark-sector paper and its evidence bundle certify all of the following:*

1.  *exact finite stress-energy closure for anomaly plus recipient stresses;*

2.  *explicit recipient stress whenever $\Gamma_{\rm rec}\ne0$;*

3.  *detailed balance or a declared nonequilibrium exchange law;*

4.  *gauge independence of the emitted perturbation variables, including the anomaly-frame baryon-density definition of $B_A$;*

5.  *causal response with subluminal finite characteristics;*

6.  *convergence under regulator refinement;*

7.  *recovery of the CDM branch when exchange, pressure, sound speed, and anisotropic stress are switched off.*
:::

::: remark
*Remark 14* (Dark continuation modes). For a source-only primordial certification run, the current dark-sector paper recommends `dark_continuation = OFF`. A supplied dark abundance may enter only as `dark_continuation = CONDITIONAL_SOURCE_STATE` until the `extra/` theorem surface derives the homogeneous anomaly abundance and finite kinetic stress construction.
:::

# Physical Kernel Interface

The physical kernels are not scalar rows in a table. They are imported finite functionals with units, gauge-invariant definitions, paired controls, and refinement margins. Their derivation stays in `extra/`; this paper only states what the CMB/LSS pipeline must receive.

::: target
**Target 15** (Physical anomaly kernel interface). *The `extra/` dark-sector papers must define deterministic receipts $$\texttt{RHO\_A\_SOURCE\_RECEIPT},\qquad
\texttt{B\_A\_KERNEL\_RECEIPT},\qquad
\texttt{GAMMA\_REC\_SOURCE\_RECEIPT}.$$ These receipts pass only when $\bar\rho_A(a)$, $B_A(k,a)$, and $\Gamma_{\rm rec}(k,a)$ are emitted with $\bar\rho_{A,\mathrm{eq}}(a)$, $w_A(a)$, $c_{s,A}^2(k,a)$, $\sigma_A(k,a)$, and $Q_A^\mu$ from the finite parent, physical scale calibration, gauge consistency, paired perturbation controls, CDM-limit recovery, and refinement convergence. The cosmology program consumes the emitted functions but does not redefine the dark-sector source law.*
:::

Diagnostic rows remain useful for debugging. They should be plotted, stress-tested, and compared to controls. They should not be accepted by the physical CMB contract unless they are promoted by the kernel receipts above.

# Boltzmann Transfer and Frozen Likelihood

Once source functions are finite and physical, the transfer problem becomes ordinary cosmology. OPH does not need a special CMB plotting rule. It needs a frozen handoff to CAMB, CLASS, or a declared independent Einstein-Boltzmann solver.

::: target
**Target 16** (Frozen transfer and likelihood protocol). *A physical OPH CMB prediction requires:*

1.  *immutable hashes for source artifacts;*

2.  *immutable hashes for solver source, version, tolerances, and input files;*

3.  *immutable hashes for likelihood code, datasets, covariances, masks, and nuisance priors;*

4.  *a no-data-use manifest showing that observational likelihood values did not enter source artifact generation;*

5.  *global pooled reducers for nonlinear source estimates and global CDM-limit checks;*

6.  *official likelihood execution readiness, not only a diagnostic binned-table comparison;*

7.  *a falsification rule fixed before the likelihood run.*
:::

::: remark
*Remark 17* (Global reducers only). Shard-local nonlinear averages and shard-local `any()` rollups are not promotion rules. The source artifact must pool additive sufficient statistics globally, with validated units, coordinate grids, coverage, duplicate policy, interpolation policy, and covariance, before any nonlinear quantity such as amplitude, rank, condition number, isocurvature leakage, $B_A$, $\bar\rho_A$, or $\Gamma_{\rm rec}$ is evaluated.
:::

::: remark
*Remark 18* (Comparison versus prediction). The current best OPH diagnostic model is scientifically useful because it is close: shape correlation $0.9951542364$, normalized RMSE $0.0984455548$, amplitude-fit $\chi^2/{\rm bin}=41.6684770597$, zero mean peak-location error, and mean peak-height fractional error about $0.0704905$. Those numbers justify fixing the theorem gates. They do not replace the gates.
:::

# Vacuum and Quantum-Fluctuation Boundary

The same finite-source discipline applies to vacuum claims. Seed noise, repair jitter, and reference free-field baselines are not an OPH-native quantum vacuum.

::: target
**Target 19** (OPH-native vacuum gate). *An OPH-native vacuum simulation needs a quotient ensemble measure or density operator on the observable quotient algebra, refinement compatibility, implementation invariance, and a sampler correctness proof such as exact sampling or detailed balance. Without this gate, vacuum-like fields in the simulator are reference baselines or diagnostics.*
:::

# Issue Ledger

The paper-side open issues should target this manuscript before simulation promotion resumes.

  Issue   Paper-side theorem gate                                                 Primary section here
  ------- ----------------------------------------------------------------------- -----------------------------------------------
  #360    Quotient ensemble measure and OPH-native vacuum gate                    Vacuum and quantum-fluctuation boundary
  #363    Boltzmann transfer and frozen likelihood closure                        Boltzmann transfer and frozen likelihood
  #371    Physical CMB finite-source contract hardening                           Physical CMB finite-source contract
  #372    OPH physical scale bridge and mode calibration                          Physical scale bridge
  #373    Finite covariant collar-packet stress parent                            Dark/anomaly interface imported from `extra/`
  #374    Physical $B_A(k,a)$, $\rho_A(a)$, and $\Gamma_{\rm rec}(k,a)$ kernels   Physical kernel interface
  #329    Flat-sector selection for the inflation-free branch                     Inflation-free branch
  #330    Screen-spectrum derivation for OPH inflation alternative                Inflation-free branch and source contract
  #319    Collar remainder stress source                                          Dark/anomaly interface imported from `extra/`
  #322    Repair-stress transport and relaxation                                  Dark/anomaly interface imported from `extra/`
  #323    Linear OPH repair-stress perturbations                                  Physical kernel interface
  #375    Consolidation under `/cosmology`                                        Closed by this staging paper

# Simulator Resumption Rule

The simulator should resume physical-promotion work only after the paper defines the exact receipt semantics for the objects it is being asked to produce. A good order is:

1.  implement the finite-source dependency DAG and no-data-use firewall;

2.  implement global pooled reducers for all nonlinear source estimates;

3.  implement the physical scale bridge and freezeout certificate;

4.  implement the imported finite covariant parent and stress-closure receipts from the `extra/` dark-sector theory;

5.  implement physical $\bar\rho_A$, $B_A$, and $\Gamma_{\rm rec}$ kernel receipts as imported dark-sector source functions;

6.  freeze source, solver, and likelihood hashes;

7.  run official likelihood comparisons and global CDM-limit reductions;

8.  only then label the result a physical CMB prediction.

The simulator can still run diagnostics earlier, including visualization payloads, CMB overlays, vacuum baselines, defect worldlines, observer cameras, and H3 readouts. Those outputs should carry diagnostic labels until the relevant receipt passes.

# Conclusion

The current OPH cosmology program is promising because the simulator almost-fits are not random noise: they expose stable structure close enough to make the theorem gates worth closing. The same data also explain why the gates are necessary. A physical cosmology prediction needs finite source-only objects, physical scale calibration, covariant stress closure, physical anomaly kernels, and frozen likelihoods. This paper is the staging surface for that work. Once its issue gates close, the simulator can resume with a precise target instead of trying to infer missing theorem semantics from diagnostic data.

::: thebibliography
99

B. Müller, A. Osika, K. Xue, and P. Nguyen, *Recovering Relativity and the Standard Model from Observer Overlap Consistency*, OPH compact paper, release , 2026. <https://github.com/FloatingPragma/observer-patch-holography/blob/main/paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf>

B. Müller, A. Osika, K. Xue, B. Cassie, P. Nguyen, M. Poneder, and K. A. Anirudha, *Observers Are All You Need*, OPH synthesis paper, release , 2026. <https://github.com/FloatingPragma/observer-patch-holography/blob/main/paper/observers_are_all_you_need.pdf>

B. Müller, A. Osika, K. Xue, and P. Nguyen, *Reality as Consensus Protocol*, OPH consensus paper, release , 2026. <https://github.com/FloatingPragma/observer-patch-holography/blob/main/paper/reality_as_consensus_protocol.pdf>

B. Müller and K. Xue, *Federated Echosahedral Screen Microphysics: Patch Hardware, Records, and Observer Synchronization in OPH*, OPH microphysics paper, release , 2026. <https://github.com/FloatingPragma/observer-patch-holography/blob/main/paper/screen_microphysics_and_observer_synchronization.pdf>

B. Müller, *Observer-Patch Holography and the Dark Matter Phenomenon*, OPH extra paper, 2026. <https://github.com/FloatingPragma/observer-patch-holography/blob/main/extra/oph_dark_matter_paper.pdf>

B. Müller, *Inflation Without an Inflaton: Observer-Screen Synchronization as an OPH Cosmology Branch*, workspace cosmology draft, 2026.

B. Müller, *OPH-FPE finite screen-consensus and cosmology diagnostics*. <https://github.com/muellerberndt/oph-physics-sim>

Planck Collaboration, *Planck 2018 results. VI. Cosmological parameters*, Astron. Astrophys. 641, A6, 2020, arXiv:1807.06209. <https://arxiv.org/abs/1807.06209>

A. Lewis, A. Challinor, and A. Lasenby, *Efficient computation of CMB anisotropies in closed FRW models*, Astrophys. J. 538, 473--476, 2000, arXiv:astro-ph/9911177. <https://arxiv.org/abs/astro-ph/9911177>

D. Blas, J. Lesgourgues, and T. Tram, *The Cosmic Linear Anisotropy Solving System (CLASS). Part II: Approximation schemes*, JCAP 07, 034, 2011, arXiv:1104.2933. <https://arxiv.org/abs/1104.2933>
:::
