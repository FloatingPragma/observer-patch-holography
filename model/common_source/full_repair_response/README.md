# Full repair response and mass-observable tests

This package follows the retained
[hierarchy diagnosis](../hierarchy_diagnosis/README.md), keeping the native
complement, quantum coherence sectors and observer clock distinctions.
The [research report](../../../docs/research/OPH_PARTICLE_MASS_RESPONSE.md#full-repair-response)
contains the derivation, tested attachments and open physical identifications.
No measured Standard Model spectrum is claimed.

| Module | Calculation |
| --- | --- |
| `memory.py` | Complete stopped population histories; exact complementary memory; minimal observable realization; incident-clock stopping map |
| `quantum_response.py` | Actual C/E Kraus response in all signed-difference sectors; terminal failures; charged-coherence and modular-flow diagnostics |
| `coupling.py` | Harmonic elimination and frequency derivative; actual gauge-consistent Higgs insertion of a declared rank-one family kernel; uniform-repair control |
| `coherent_arc.py` | Existing coherent continuation with mixed-cap weights, full arc complement, source-executed coin reflection and terminal Ramsey outputs |
| `mass_observable.py` | Calibration of charged rest poles, dispersion and Euclidean decay using the actual declared chiral candidate and clock |
| `current_response.py` | Exact native current fluctuations, charge continuity, complete record-preparation response and its finite preparation budget |
| `run.py` | Reproducible report, source custody and retained coherent-work arrays |

## Passive stationary schedule selection

The [stationary selection receipt](stationary_selection_results.json) checks
all 7,920 native input/coin addresses. For the same faithful core-biased
reference, any nonnegative state-independent passive seam schedule is
stationary only if all eight core-to-outer weights vanish. The other
22 seams, including the core-internal seam, are admitted. Every resulting
recorded history conserves core occupation. The maximal schedule's
static Schur family seed has rank two; inserted into the same declared
matter/Higgs route, it leaves seven charged Dirac fields massless, counting
color. A crossing edge of weight 1/64 restores seed rank three but causes
the exact reference defect 157/5760000. All-positive repair instead selects
the uniform diagonal reference and a scalar one-body family feature.
This limits passive schedule tuning at this reference. Active correction,
other reference choices and other charged mass couplings are separate.

## A measured native current response

The current is the integer actually recorded by the native rule:
$x'=x-B_e j$, with $B_e=e_a-e_b$ for the addressed seam. The calculation uses
the same faithful protected-core reference. All 2,640 source state/seam/coin
branches determine the mean, fluctuation and conservation identities; no
coupling or observed particle value enters.

For the mixed complement schedule, let $D=\mathrm{diag}(p_e)$ be the
admitted seam probabilities, $B$ its incidence matrix and
$M=I-BDB^T/2$. The recorded currents have zero mean, and their exact
stationary covariances are

```math
\Gamma_0=dD,\qquad
\Gamma_k=-\frac d4 DB^TM^{k-1}BD\quad(k\geq1),\qquad d=\frac{19}{150}.
```

Independent contractions through all 66 native states verify the first three
lags. The parameter is admitted native work-proposal count. The retained
one-round selector has survival $q=27/32$ and terminal timeout $5/32$.
For its actual stopped records, assigning zero after termination gives
$\mathbb E[J_kJ_0^T]=q^{k+1}\Gamma_k$; a timeout is not restarted.

The normalized admitted-proposal chain has the exact zero-frequency
current-noise form

```math
S(0)=dC,\qquad
C=D-DB^T(BDB^T)^+BD,\qquad
BC=0,\quad C D^{-1} C=C.
```

Thus $v^TCv=(Cv)^TD^{-1}(Cv)\geq0$. Its rank is eleven, the cycle dimension
of the twenty-edge connected ten-port complement. The finite Green–Kubo
sums retain their exact remainder $dDB^TM^N(BDB^T)^+BD$.
This supplies a concrete transition-current counterpart to the
[finite reversible-kernel theorem](../../../Lean/Thermodynamics/GreenKubo.lean).
The infinite sum belongs to the normalized diagnostic chain; it is not an
infinite surviving history of the terminal-timeout instrument or a physical
conductivity measurement.

The signed records also recover a finite occupation-field readout without an
initial occupation measurement. For an actual window of $N$ events, let
$z^N=1$ and $z\ne1$. Summing the pathwise continuity equation gives

```math
\sum_{n=0}^{N-1}z^n x_n
=\frac{B\left(\sum_{n=0}^{N-1}J_n-z\sum_{n=0}^{N-1}z^nJ_n\right)}{1-z}.
```

This holds on open paths: the first sum retains the endpoint difference
$B\sum J_n=x_0-x_N$. No periodic endpoint is assumed, and the final current
cancels from the result. A complete four-event native branch audit checks
all three nonzero Fourier coordinates on 1,056 histories, including open
paths and zero-current outcomes. Identical zero-current words with different
occupation sums demonstrate the remaining DC ambiguity. A terminated window
uses its actual event count and retains its terminal flag. These frequencies
refer to the native event index; no physical-time or coherent phase service
is inferred from performing Fourier arithmetic on records.

Summing the exact covariance series gives the whole stationary current-noise
spectrum. With $c=\cos\omega<1$,

```math
S_J(\omega)=d\left[D-\frac12 DB^T(cI-M)
(I-2cM+M^2)^{-1}BD\right],\qquad
BS_JB^T=2(1-c)S_x,
```

where $S_x=d(I-M^2)(I-2cM+M^2)^{-1}$ on the dynamic density modes.
The frequency Ward identity attaches the current and density measurements
to the same source. Conserved density components have separate DC weight.
On a density-gradient direction in the coordinates $D^{-1/2}J$, the
spectrum is $d(1+m)(1-c)/[(1-m)^2+2m(1-c)]$, with the corresponding
transfer eigenvalue $19/27\leq m<1$. Thus every nonconstant pole of an
instantaneous linear current readout is a density relaxation pole.

For any cycle $Bc_{\mathrm{loop}}=0$, the scalar record with weights
$h=D^{-1}c_{\mathrm{loop}}$ instead has zero conditional mean for every
native state. Under the fixed independent selector it is a martingale
difference and has a white two-point spectrum. An explicit bounded readout
is $J_{23}-J_{24}+J_{34}/2$, taking unselected seam currents as zero; its
variance is $19/1620$. Its mean response to every earlier record is zero.
This does not make its outcomes independent or Gaussian. The eleven cycle
directions therefore supply measured transverse record noise, without an
oscillatory photon pole or a Maxwell kinetic normalization. This conclusion
concerns these linear native records; nonlinear fields and additional
coherent operations have different response contracts. In particular, the
existing [coherent arc walk](coherent_arc.py) uses an uncommitted event register,
an explicit reflection, and ideal cap weights different from this finite
selector's admitted weights. Its eigenphases cannot be substituted into the
current-noise measurement without deriving a common instrument and clock.

The [collective field's local charges](../hierarchy_diagnosis/collective_response.py)
also have a current in this same source. For $Q_j=n_0n_j$ when $j\ne0$ and
$Q_0=0$, the binary total-two sector gives $\sum_jQ_j=n_0$.
Complement repair preserves $n_0$, hence its current is $J_e=n_0j_e$ and
$Q'=Q-BJ$ pathwise. The same covariance formula holds with
$d_Q=1/36$, verified through the full native state histories. Reading this
weighted current requires access to the protected $n_0$ occupation as well
as the raw seam current. That occupation read is a separate intervention
from the collective field's coherence-sensitive Ramsey measurement.

An addressed two-event experiment makes the raw-current response directly
operational. After recording $j_0$ on a complement seam, two native coin
reads create a flag with probability $(1+j_0/2)/2$. Both flag values remain
in the record, and the same work carrier undergoes the second event. All
1,056 histories give

| Readout | Exact mean |
| --- | ---: |
| Second current, conditioned on the accepted flag | $-19/600$ |
| Second current, conditioned on the rejected flag | $19/600$ |
| Signed flag times second current, all outcomes | $-19/600$ |
| Same signed read with an independent fair flag | $0$ |
| Reversed current-dependent flag | $19/600$ |

Zero-current records remain in the denominators. There is no work reset;
auxiliary coin reads continue on their own carrier. The exact relation is
$\mathbb E[j_1\mid\mathrm{accept}]=\epsilon\mathbb E[j_1j_0]$ with
$\epsilon=1/2$. It is a response to an implemented record-based preparation,
not to an assumed electromagnetic force.

The existing finite reference preparation bounds the signed-current mean
error by less than $1.50\times10^{-6}$, compared with signal $19/600$.
Counting each preparation failure and timeout as a trial with zero readout
gives ideal attempted-trial mean $-19/3960$. A fixed two-million-attempt
comparison with the unbiased control has each error below $1/100$ under
the exact source and its uniform preparation bound. A bounded-increment
exponential argument needs no exact reset or actual trial independence.
The receipt gives the large worst-case native-attempt budget; this campaign
is not executed, and it supplies no physical duration or apparatus-error
certificate.

The calculation also fixes a useful distinction for electromagnetic readout.
The native current instrument has mean-effect operator
$Q_e=(n_a-n_b)/2$ but second-moment effect $2Q_e^2$.
Its record variance is $19/150$ while the operator variance is $19/300$.
Every $Q_e$ commutes with the diagonal reference's modular generator, so
its modular retarded commutator is exactly zero. The nonzero transition
current cannot therefore be replaced by a modular charge derivative.
Charge continuity fixes relative source charge/current units; it does not
fix a Maxwell kinetic coefficient. A photon identification needs a derived
propagating transverse field, a conserved Lorentzian current and their
physical-clock pole/residue normalization. Wissner-Gross's
[entangled-pair scalar](../../observer_particles/results/entangled_alpha.json)
is not identified with either of these current variances.

Use the repository root's Python environment and `requirements.txt`.
The retained calculation uses Python 3.12 and the pinned NumPy/SciPy/SymPy
dependencies. In a clean checkout, create an environment if necessary:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q model/common_source/full_repair_response
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.common_source.full_repair_response.run --output /tmp/fresh-full-repair-response.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.common_source.full_repair_response.run --verify model/common_source/full_repair_response/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.common_source.full_repair_response.stationary_selection --check model/common_source/full_repair_response/stationary_selection_results.json
```

The producer refuses to overwrite JSON or its companion NPZ. Verification
recomputes this package's report and coherent-work arrays, and checks source
bindings. It verifies the earlier candidate's retained weights and custody;
it does not rerun that earlier full experiment. Floating source residuals
and spectral diagnostics are distinguished from exact rational identities
and analytic bounds. Replay compares exact claims and rational strings
exactly, and floating diagnostics within the explicit tolerances in the
receipt; the retained NPZ itself is checked by its byte digest.
No physical service-error or finite-shot budget is
inferred from numerical agreement.

To resume, first reproduce this result. Then address the operator map from
native repair to the charged left/right field response. The most specific
obligations are charge transport, the field algebra, the same-history clock,
kinetic normalization, interaction vertices and the reconstructed spacetime
domain. The full matrix algebra generated by compressed native family
features permits many couplings; by itself it selects none. The
[full-mode calculation](../../yukawa_source_dynamics/PROOF.md) proves that
the twelve-port mean-occupation and completion matrices commute: their
compressed commutator comes entirely from intermediate complement ports.
Its canonically normalized zero-frequency Schur response has rank one;
its nine positive relaxation poles are not three fermion masses. A proposed
continuation must produce its ratios before the empirical comparison and
retain every complementary, timeout and detector output relevant to its
claimed scope.

The [native Higgs dual](../../yukawa_source_operator/PROOF.md) extracts all
family coefficients once their actual Higgs derivative and kinetic
normalization are supplied. The
[complete Higgs join](../../yukawa_source_running/PROOF.md) derives its
no-hopping Pauli-blocking costs and reconstructs `abs(g)*sigma_i(Y)` from
both spectral lines and all blocked outcomes. The distinct
[ordered feedback controller](../../yukawa_source_feedback/PROOF.md)
uses the [complete-memory theorem](../../particle_memory_response/PROOF.md)
for the coherent correction and noise, and derives the exact next
family coefficient on its full 432-dimensional channel. These results keep their
specified preparations, pulses and clocks; none identifies its diagnostic
coefficient with a physical Yukawa coupling.

All earlier retained evidence is included in the sibling packages. Their
replay and verification commands are documented separately:
[mass response](../mass_response/README.md),
[particle identification](../particle_identification/README.md),
[chiral candidate](../chiral_candidate/README.md), and
[hierarchy diagnosis](../hierarchy_diagnosis/README.md).
