# Row 31. Orbital levels: theorems, premises, gap (worker TD, 2026-11-14)

`RER/` = `reverse-engineering-reality/`; `EOC` = `RER/paper/executable_observer_construction.tex`; `FL` =
`RER/flagship/from_observer_consensus_to_standard_physics.tex`; packet = `openwave-push/mf/O_rows_28_31/ROW_31/STATUS.md`.
Pin package (RER_MAP.json): `RER/model/openwave_atomic_refinement`.

## 1. Chain

| # | Statement | Where | Kind |
| --- | --- | --- | --- |
| 1 | One-unit native kernel: under law C, uniform seams, fair coin, a single unit's step law is P = I - L/60 on the carrier, for every background | `openwave-tuning/measurement_first/gravity_clocks/METHOD_NOTE.md` "Exact one-unit theorem" | proof + exec |
| 2 | Unread completed transfer F = D^{-1}B is self-adjoint in the metric G with 1/20 <= F <= 1; H = -log F is an effective Euclidean generator per completed tick | EOC 2744-2746 | proof |
| 3 | A positive state-dependent weight D gives the Feynman-Kac operator DP, similar to sqrt(D) P sqrt(D): eigenvalues real and nonnegative, no oscillatory eigenphase | `openwave-tuning/measurement_first/record_foundations_review/METHOD.md` lines 91-99 | proof (negative for "standing waves") |
| 4 | Agreement of b independent coin prefixes gives the Gibbs family of the same cost at integer multiplier b; a Lean instance builds the repair and calibration data | FL 980-985 | proof (+ lean per FL) |
| 5 | Records: weighted statistic W_tau = E[exp(-(gamma/2M) sum u)] has rates in shells 1, 4, 9; 2p 1.0506, 3p 1.0548, 3d 1.0174, 1.0628 of Ry = gamma^2/(18 pi^2); weight table from 320,000 releases, exponent 1.019 | packet section 4 | exec (N+, weight supplied) |
| 6 | s levels 1.554, 1.394, 1.276; excess falls with a_B (0.089 at 2.36 cells, 0.048 at 3.14); first-order coefficient off by 17 to 21 % (G10) | packet section 4 | exec |
| 7 | Gauss atom on the Dirichlet Green function, mass 1, coupling^2 = 4 pi: nested-cube maximum principle gives G_large >= G_small; min-max makes every ordered eigenvalue nonincreasing under domain extension; Collatz bounds on the ground level | `RER/model/openwave_atomic_refinement/PROOF.md` | proof + exec (declared field) |
| 8 | Walk reduction on an installed walk W = S C_mu P_v: ker F(omega) <-> eigenspace of W at exp(-i omega), F'(omega) >= 6 I; lambda0 = abs(k)^2/8 + O(k^4) | `openwave-push/reviews/ROW_31_WALK_REDUCTION_FABLE_B/REVIEW.md` statements 1a-1c, 3 | proof (reviewed, installed source) |
| 9 | Record Coulomb table: mass-tape mean field 1.00125 x 1/(6 pi), exponent 1.0085 | `openwave-push/blockers/X1g_mass_field_lapse/VERDICT.md` | exec (N+) |

## 2. Premises

| Id | Text | State |
| --- | --- | --- |
| UNREG-TD-WEIGHT | the reader's weight exp(-(gamma/2M) sum u) and its strength gamma | open; supplied. Three record normalizations give gamma = 1, 3.14, 6 and no record selects (packet section 3). Theory route: step 4, "selection owed" (SS 941-958) |
| PR-04 | declared phase-sensitive effect | declared input |
| PR-21 | oscillator lift of the seam-current generator | constructed on declared inputs; "passive native repair has no wave sector (OP7)" |
| PR-43 | star-automorphism flow | obstructed: "the executed quantum channel strictly attenuates a coherent mode" |
| PR-54 | source gauge field, current and action | obstructed on unrecorded native channels |
| PR-65 | source-produced phase preparation | obstructed (native coherence floor) |
| UNREG-TD-PHASE | phase access (DAG F4): an unread choice register (candidate D') or a pinned reflection as memory law | open; owner decision. Step 3 proves the committed schedule has no phase |
| PR-18, PR-19 / PR-29 to PR-31 | the 1/r interaction (row 20, Gauss shell route) | open; step 9 is a record route with supplied decoders |
| UNREG-TD-ROUTING | cubic routing of the 3-dim block (`model/native_spatial_refinement`) | declared (decision); PR-16 for a physical continuum |
| PR-15, PR-22 | energy and clock: naming -log lambda an energy, rates per sweep as frequencies (F6) | declared input; packet keeps rates only |
| PR-05, PR-20 | uniform schedule, fair coin | declared |

## 3. The gap

The decisive open component is UNREG-TD-WEIGHT, the source origin of the
attractive interaction. UNREG-TD-PHASE must be interpreted with the
observer-first reconstruction in section7: positive Euclidean transfer does
not force an unread primitive oscillator as the only route to standing modes.

1. **One optional event-time phase route (UNREG-TD-PHASE).** Steps 1 and 3 prove that with every register read the spectrum is real.
   An event-clock wave phase needs a protocol outside the passive fixed-IID record law. An adaptive classical
   controller can also oscillate, so oscillation alone establishes no quantum phase. The "choice unread"
   schedule of D' (`openwave-push/onecandidate/DP1_retention_schedule/CANDIDATE_DPRIME.md`) plus a controlled
   phase exp(-iV(x)) is one declared design, touching PR-21, PR-43, PR-54. Its retained/dephased operational
   contrast and the source origin of that phase must be demonstrated on the same candidate.
2. **Theorem (after 1).** On the D' block of side L with V = -(gamma'/12) u, u the record table, the unitary step
   W has eigenvalues exp(-i omega_n) with omega_n = omega_inf - R/n^2 + O(a_B^-2) for l >= 1, shells 1, 4, 9,
   R fixed by the walk's kinetic constant (step 8 gives the reduction; the ordinal count is not certified).
3. **Execution (UNREG-TD-WEIGHT).** Replace the supplied weight by recorded acceptance: at each own tick the
   observer draws b prefixes and keeps the history when all agree, b read from row 20's tape at its port, every
   rejection kept (design on file). Gates: >= 4 levels, (E2 - E1)/(E3 - E1) within 5 % of 27/32. With step 4
   this turns gamma into a record quantity; a theorem linking b to the tape mean is owed with it.
4. **Execution.** Blocks of side >= 96 at a_B >= 2.2 cells for the s levels (G5), and the second-order
   coefficient behind G10.

## 4. Shared premises

With row 20: PR-18, PR-19, PR-29 to PR-31, UNREG-TD-WEIGHT (the interaction). With rows 3, 29, 30 and F4:
UNREG-TD-PHASE, PR-04, PR-21, PR-43, PR-54, PR-65. With rows 2, 27 (F6): PR-15, PR-22. With all rows on A:
PR-05, PR-20.

## 5. Finite-window extension and operational-time scope (Codex MFR)

[New proof](work/2026-11-14_codex_record_poles/PROOF.md): for every finite graph and
nonnegative fixed total, uniform or fixed-IID seam repair has only eigenvalues
in [0,1]. Ordering configurations by Q=sum n² makes the matrix block triangular;
each same-Q seam block is an orthogonal projection. Any passive finite-window
function of committed letters has gap correlation μT_f K^nT_g1 and no new poles.
Positive one-state Feynman-Kac weights retain nonnegative-real poles; arbitrary
history-dependent weights and adaptive schedules are outside that corollary.

[Executed instance](work/2026-11-14_codex_record_poles/RECEIPT.json): 71,160 pinned
native calls, totals1,2,3, all27,240 branches and43,200 nonlinear windows retained.
Native poles are1,9/10,(55±sqrt(5))/60; the changed-law decoder detects four nonreal
poles and the independent planted control detects two. This extends E-03 and
E-08; it does not discharge PR-21, PR-43, PR-22 or the common-source bridge.

Event count is not assumed to be physical time. A relational time or an Euclidean
reconstruction must establish its time/algebra/state/readout transfer. The existing
`model/openwave_native_rest_energy/PROOF.md` already proves adaptive current-only
phase blindness from partial-injection Kraus histories, while allowing such an
alternate reconstruction. No theorem here excludes all observer-first physics.
The written new proof awaits another group's review; no new Lean build is claimed.

## 6. Positive reconstructed-observable route (finite sources inspected)

A committed diagonal record algebra and the algebra reconstructed from record
statistics are different objects. `RER/model/native_quantum_reconstruction/PROOF.md`
derives all4,356 matrix units of M66(R) from nonzero native path words (length<=4),
then its canonical complexification and the positive functional
rho_L=P^L/Tr(P^L) from closed-return probabilities. The current-only census and
signed tape attach these ideal counts with an explicit finite apparatus error;
the stated effect-statistic systematic bound is3.76596e-6. This is a source-derived
positive matrix functional, not a supplied density matrix.

`RER/model/native_modular_dynamics/PROOF.md` reconstructs ordered modular
correlations as finite signed/complex combinations of native event probabilities,
with exact binomial truncation bounds, and fixes the comparison t=Ls between
loop-length references. `RER/model/native_effective_fields/README.md` derives
G=sum_e(I-S_e), P=I-G/60 and H=-logP, including a nonquadratic interaction witness.
The spatial metric of its finite12-point support comes from current responses.
Thus linear mean-load superposition does not imply absence of every reconstructed
interaction, and a commutative-algebra no-flow theorem does not apply to M66(C).

The packages keep their remaining boundary explicit: signed reconstructed effect
weights need not be a single native measurement, the native returned working
state is diag(rho_L), and its evolution is not the modular unitary. Native attempts,
Euclidean parameter, modular parameter and physical duration are distinct. The
same-history operational algebra/state/clock/readout intertwiner remains owed.
This positive reconstruction chain should be used before any new field equation
is chosen; the fixed-IID window-pole theorem in section5 does not contradict it.

## 7. Positive Euclidean transfer and standing modes

The absence of primitive event-lag oscillation is not by itself an obstruction
to observer-first standing waves. For a finite connected one-unit source with
symmetric positive-definite stochastic transfer P and a positive diagonal weight D,
T=sqrt(D)P sqrt(D) is positive. Let its Perron eigenpair be(t0,psi0),
psi0>0. The normalized transfer T/t0 is a positive contraction, and
Q_ij=T_ij psi0_j/(t0 psi0_i) is a stochastic reversible kernel with reference
psi0_i² after normalization. Its Euclidean histories are reflection positive.
The reconstructed H=-log(T/t0) gives stationary spatial eigenvectors and
their unitary standing phases. No primitive unread register is required for
this mathematical reconstruction. A possible zero eigenvalue requires
restriction to the positive spectral support; the stated strict-P case avoids it.

This corrects the earlier inference from “real nonnegative eigenvalues” to
absence of every standing-wave interpretation. It does not derive D. The
existing atomic packet still supplied its attractive weight and gamma, so
its interaction claim remains conditional. The new free-field and first-hit
clock results make the source of that interaction the next concrete target.
Current work tests a native first-hit feedback controller before assigning
any potential or computing its orbital spectrum. No row is promoted here.


## 8. Native feedback interaction and controlled orbital shells

The sealed `work/2026-11-14_codex_probe_binding` packet executes305157
native calls, including all relocation, inner first-hit, auxiliary-hold and
outer response events. Its observer never reads a numerical potential or
gamma table. Target-hit feedback derives
L_eff=L_external+(1-h_B)L_internal, hence a band operator
S=a^-2 L_space+lambda(1-h_B). The constant internal band is unaffected;
positive native bands have a stable attractive well. Disabling feedback
removes it. The fully capped implementation has a retained preparation error;
the exact Markov matrix refers to exact/almost-sure preparation.

The sealed written proof `work/2026-11-14_codex_coulomb_limit` supplies the
previously missing source-derived shell limit. With a=1/m, B=m³ and y=a²z,
(S-lambda)/a² converges below threshold to
-Delta-lambda c_*/|y|, c_*=1/(4pi G_cubic(0)). Discrete Sobolev control
removes a singular-cell artifact; compact Birman–Schwinger convergence
preserves every fixed shell's eigenvalue count. The spatial multiplicity is
n², and source squared-level shifts are−a²(lambda c_*)²/(4n²)+o(a²).
Internal degeneracies multiply that spatial count and remain declared.

The same positive scalar record reconstruction yields standing energies
E_a(s)=2asinh(a sqrt(s)/2)/a. Binding is measured from E_a(lambda), not
sqrt(lambda), retaining the same-order vacuum lattice correction. A time
period2m³+1 gives exponentially small positive-band clock error, o(a²).
Finite native caps/costs and a computable increasing spectral menu are
specified. No large-shell production calculation or physical charge/species
identification is asserted. This written proof has same-group review; owner
pin, another-group review and the registered platform task remain open.


## 9. Actual global production source and corrected band decoder

`work/2026-11-14_codex_production_field` executes347533 native calls:82008
census,768 full-vector checks,262144 outer/auxiliary and2613 inner calls.
Every occupation, band/curl readout and periodic first-hit clock is retained.
The exact source-mode residual is below2.347e-14, compressed covariance error
below7.958e-5 and clock error below4.180e-14. Feedback-disabled histories
remove the well; the constant internal sector changes by at most6.555e-16.

Independent review found the initial NumPy view-normalization decoder bug.
The original producer-pass/peer-fail attempt remains intact. A preregistered
copy/projector-check correction changes only the eight band readouts; all
raw native histories, occupations, clocks, constant readouts, curls and
population predictions remain identical. The corrected replay and same-
group review pass. This closes the actual production join, not a large-
shell numerical accuracy gate; the full shell conclusion remains section8's
controlled written limit. No external acceptance is inferred.
