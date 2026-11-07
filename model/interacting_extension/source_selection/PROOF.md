# Native records and conditional references

## Scope and result

This extension uses the pinned r2062 endpoint law, commit
`3cf3c14810671a16e61b7915f69079a1773a12a2`. It does not use any measured
mass, mixing matrix, energy calibration, or empirical fitted parameter.

There are two constructive results.

1. A native recorded event selects a finite orbit of faithful three-dimensional
   reference states through an explicit complete chart transducer. Its
   orientation is selected by the observed seam and record, rather than a
   separately chosen reset port. The no-current companions and every
   transducer outcome remain in the calculation. Successive native records
   select further references with a fully normalized joint law; some of
   those finer conditional references have smaller rank.
2. A positive tilt of the **actual recorded move count** gives a one-parameter
   family of completely positive, trace-preserving quantum Doob instruments.
   Every coefficient and stationary state is then determined by the native
   law and that one parameter. A native-plus-filter implementation includes
   every failed outcome. A three-dimensional quotient determines the
   family's Perron data and its clock statistics exactly.

Neither result supplies an autonomous relativistic vacuum or a physical
time/energy calibration. The source does not identify the tilt with public
`P`, and public record capacity `N` does not occur in the native endpoint
transition. Their selection requires additional source dynamics.

The source distinction is deliberate. Paper 11 describes the native repair
and supplied logical services and says that the closure value of `N` is not
consumed by an execution component. Paper 17 distinguishes a stationary
reference and declared variational problem from selection of realized
Hamiltonian flow. Paper 27 explicitly requires a source Euclidean slab,
including conductance, potential, and slab duration, before assigning a
sampler the status of a physical vacuum. The present calculation does not
remove those qualifications.

## 1. The complete native instrument

Let `X` be the 66 binary twelve-port load configurations with total load
two. This is an exact invariant sector of the literal endpoint rule, not
a truncation based on a small amplitude. Let `E` be the thirty source seams.
For seam `e=(a,b)` and fair placement coin `c`, the source returns

\[
j={x_a-x_b-(2c-1)|x_a-x_b|\over2},\qquad
f_{e,c}(x)=x-j\mathbf e_a+j\mathbf e_b.
\]

The source C record resolves `j`; E resolves the initial contrast
`d=x_a-x_b`. Source coding uses `CURRENT_CODES={0:0,1:1,-1:2}`.
C has a four-dimensional raw pointer. E adds the parity bit and has an
eight-dimensional raw pointer. Thus the complete chosen policy has

\[
30\times2\times4=240\quad\hbox{C labels},\qquad
30\times2\times8=480\quad\hbox{E labels}.
\]

Only 120 C labels or 180 E labels are attained on the binary sector. All
other raw labels remain as exactly zero Kraus operators. The implementation
reads the raw-code convention from the pinned source; it does not silently
replace a four- or eight-dimensional source pointer with a three-bin
Fourier pointer. It is an exact finite instrument reconstruction using the
source's phase-exact convention; it does not reexecute every compiler pulse.

For the full label `r=(e,c,raw)`, define

\[
K_r={1\over\sqrt{60}}
\sum_{x:\,r(x)=r}|f_{e,c}(x)\rangle\langle x|.
\]

Within a fixed record the endpoint map is injective: the output and recorded
current reconstruct the input. On the binary sector E's contrast and coin
also reconstruct that current. Therefore each `K_r^† K_r` is diagonal. Each
input has exactly sixty nonzero event columns, proving exactly

\[
\sum_rK_r^\dagger K_r=I_{66}.
\]

Every complete outcome is the subnormalized matrix
`K_r rho K_r^†`. Its trace is its Born probability and its normalized form,
when nonzero, is its state update. The `complete_outcomes` API returns every
label. The maps apply to arbitrary native/reference entangled inputs by
tensoring with the untouched reference identity; no classical-input
assumption is needed for this completeness statement. Native measurement
backaction can change those correlations.

The same API with `binary=False` retains the whole 78-state native carrier,
including all twelve double-load states. It has 4680 nonzero event columns,
with 180 attained C labels or 300 attained E labels within the same padded
240/480 alphabets. Those full-carrier maps are saved as well. If `P_D` is the
double-load projector, the literal endpoint columns prove

\[
\Phi^\dagger(P_D)=\frac56 P_D,\qquad
\operatorname{Tr}(P_D\Phi^n(\rho))
=\left(\frac56\right)^n\operatorname{Tr}(P_D\rho).
\]

Thus the binary shell is both invariant and approached from the full
carrier; double-load records are not assumed nonexistent. The exact
self-adjoint transfer and Doob formulas below are stated on the binary
invariant sector and are not incorrectly extended to the nonunital
78-state transient map.

Uniform independent seam choice is the unique memoryless **seam-only** policy
invariant under the transitive native seam action. This is a conditional
selection statement: covariance does not select memorylessness, whether C
or E is measured, a feedback controller, or a waiting-time process. The
placement coin, seam choice and raw repair record remain distinct records.

### The native stationary state

The population channel is the lazy two-particle exclusion walk on the
connected icosahedral graph. It is symmetric and irreducible, and has the
uniform stationary distribution. The full quantum channel has the exact
signed-load-difference blocks used by the source's `quantum_mixing.py`.
Their sizes are

\[
1\times66,\qquad132\times10,\qquad2970\times1,
\]

which account for all `4356=66^2` operator directions. For a nonzero signed
difference `delta`, the C row sum is

\[
1-\frac1{60}\sum_{(a,b)\in E}|\delta_a-\delta_b|<1,
\]

and the E row sum is

\[
1-\frac1{30}\#\{(a,b)\in E:\delta_a\ne\delta_b\}<1.
\]

Thus no off-diagonal direction is stationary. The unique native stationary
state on this invariant sector is `I_66/66`. It is a conclusion from the
complete quantum instrument, including backaction. Its modular Hamiltonian
is a scalar. A nontrivial coherent vacuum is not concealed in this particular
unconditioned finite instrument.

The full 78-state stationary state must have zero double-load population by
the identity above; positivity then kills its cross-support matrix elements.
It is therefore the same `I_66/66` embedded in the full carrier. This does not
make it faithful on all 78 states.

An exact factorization of the complete 66-state population characteristic
polynomial, saved in the receipt, gives the second eigenvalue
`r=(55+sqrt(5))/60`. The off-diagonal blocks have smaller spectral radii, so
for every binary input state

\[
\frac12\|\Phi^n(\rho)-I_{66}/66\|_1
\le\frac{\sqrt{65}}2 r^n.
\]

The reference-stable half-diamond distance to the replacement channel is at
most `min(1,33 r^n)`. Hence this native reference has an actual finite
preparation approximation through the same complete recorded dynamics.
The approximation cost is an attempt count and its retained tape capacity,
not a physical temperature or an energy unit. The exact stationary state is
a mathematical limit, not an exactly attained finite-tape preparation.

## 2. A reference selected by an observed native record

Let `v_i` be the normalized source frame rays and
`P_i=|v_i><v_i|`. The exact frame gives

\[
\sum_{i=0}^{11}P_i=4I_3,\qquad
\operatorname{tr}(P_aP_b)=\frac15\quad\hbox{on every seam}.
\]

Define a complete chart measurement/preparation transducer by

\[
L_{x,i}=\sqrt{\frac{x_i}{2}}\,|v_i\rangle\langle x|,
\qquad x\in X,\quad i\in\{0,\ldots,11\}.
\]

Since the total load is two,

\[
\sum_{x,i}L_{x,i}^\dagger L_{x,i}=I_{66}.
\]

The complete composed instrument has Kraus operators `L_(x,i) K_r` and
records `(r,x,i)`, including all zero possibilities: 190,080 C labels or
380,160 E labels. A version retaining the measured native load replaces
`|v_i>` by `|x> tensor |v_i>` and has the same family marginal. This is an
explicit declared chart readout/preparation service. It is not claimed to
be passive native time evolution, and it destroys native coherences. Its
completeness is an identity for arbitrary inputs.

Conditioning only on the native raw record `r`, and reducing the transducer
outputs, yields

\[
F_r={1\over p_r}\sum_x
\langle x|K_r\rho K_r^\dagger|x\rangle
\left({1\over2}\sum_i x_iP_i\right).
\]

The full `(r,x,i)` branch family state is pure `P_i`. The faithful mixed
`F_r` is the reduced state conditional on the native record alone; the
additional transducer outcomes belong to the complete law. This prevents
confusing a reduced observer reference with the state conditional on every
fine-grained apparatus outcome.

### First record from the native stationary state

Take `rho=I_66/66`. A nonzero current from `a` to `b` leaves `b` occupied and
the other occupied port uniformly among the ten ports other than `a,b`.
Consequently

\[
F_{a\to b}
={1\over2}P_b+{1\over20}\sum_{i\notin\{a,b\}}P_i
={I_3\over5}+{9P_b-P_a\over20}.
\]

Its exact spectrum is

\[
\boxed{\left\{\frac{20-\sqrt{145}}{50},\frac15,
                  \frac{20+\sqrt{145}}{50}\right\}.}
\]

The probability of a particular moving C record is `10/(66*60)=1/396`;
the sixty moving labels have total probability `5/33`. No moving outcome is
silently renormalized into a certain event.

For a C zero-current record, let `a` be the endpoint whose lone occupation
would have moved and `b` the coin's preferred endpoint. The 56 admitted
unchanged loads give

\[
F^{C,0}_{a,b}={40I_3-9P_a+P_b\over112},
\qquad
\operatorname{spec}F^{C,0}
=\left\{\frac{90-\sqrt{145}}{280},\frac5{14},
                  \frac{90+\sqrt{145}}{280}\right\}.
\]

Their combined probability is `28/33`. For an E even-contrast record, 46
loads are admitted: both endpoints empty or both occupied. It yields

\[
F^{E,0}_{a,b}={9I_3-2(P_a+P_b)\over23},
\qquad
\operatorname{spec}F^{E,0}
=\left\{\frac{7-2/\sqrt5}{23},\frac{7+2/\sqrt5}{23},\frac9{23}\right\}.
\]

Their total probability is `23/33`. The E odd-contrast branches have exactly
one occupied endpoint; their post-repair reference has the same spectrum
as `F_(a->b)`, whether that particular coin moved it or left it in place.
Their combined probability is `10/33`. Thus an E no-move event need not
belong to the even-contrast reference class.

Every matrix above is faithful. Their exact characteristic polynomials,
not merely rounded eigenvalues, are checked in `results.json`. The native
rotations act transitively on oriented seams. If `R` takes `(a,b)` to
`(a',b')`, it gives

\[
F_{a'\to b'}=R F_{a\to b}R^T.
\]

The code explicitly checks all sixty oriented-edge rotations, including
their action on all twelve source rays. Hence orientation is selected by
the record while the first-record spectral orbit is fixed. There is no
continuous choice of eigenvalues hidden in a chosen seam direction.

Summing **all** first-record subnormalized references gives `I_3/3`.
Source symmetry and completeness are consistent with nonuniform conditional
references. They do not select one conditional branch as a global vacuum.

### Successive records and loss of faithfulness

For a two-record history `(r,s)`, every surviving initial load has amplitude
`1/60`. Its probability is therefore exactly

\[
p_{r,s}={m_{r,s}\over66\,60^2}
={m_{r,s}\over237600},
\]

where `m_(r,s)` is the integer number of surviving initial configurations.
The code saves every integer numerator for all 57,600 C record pairs and
230,400 E record pairs, including zeros. Their numerator sum is exactly
237600. Both one-record marginals are checked exactly. Each history's
reference is an explicit positive sum of the same source `P_i`, with the
integer port numerators saved in the archive.

The two-record packet applies the chart transducer at the final time. For
the stated tracial native input, every intermediate conditional native state
is diagonal. Retained-copy transducers inserted between events therefore
leave this particular raw-record law unchanged. Their reference outputs can
be correlated; the packet does not replace a common history with
independent freshly drawn references.

These references are generally not all unitarily equivalent. The first
record selects one spectral orbit; subsequent observations refine the
posterior and can select a different spectrum. Some histories determine
only one or two ray directions, making the reference rank one or two.
The six antipodal ray lines have positive three-line Gram determinants
`2/5 +/- 2/(5 sqrt(5))`, so the exact rank is the smaller of three and the
number of distinct antipodal lines with positive port weight. Exact rank
counts and probability numerators are reported. Small negative numerical
eigenvalues of order `1e-17` are roundoff of these proved positive sums;
they are not clipped into a false statement of full rank.

A full-space modular logarithm cannot be applied to those rank-deficient
conditional histories. One must explicitly restrict the observation
algebra to their support or supply a source-selected faithful refresh
process. Taking finer observations does not preserve faithfulness for free.

The conditional references and their histories are source-derived within
the specified transducer. They do not supply a Yukawa operator, a vacuum
preparation rule for the interacting gauge–matter theory, or an energy
normalization merely because their matrices have three eigenvalues.

## 3. Conditioning the actual repair activity

The recorded activity is `a_r=|j_r|`, equal to zero or one. It is determined
by the actual raw record and coin. Reserved impossible raw codes have zero
Kraus operators, so assigning them cost zero cannot affect any probability.
Choosing this observable as a variational cost is an explicit additional
premise. It has not been identified with the compiler's gate count, a
source spacetime action, or the public coupling/pressure parameter. In the
unmodified hidden-state event law every one of the sixty seam/coin choices
has probability `1/60`; its bare per-event log probability alone does not
penalize moves differently from no-move events.

For explicit `0<z<=1`, define

\[
\mathcal E_z(X)=\sum_r z^{a_r}K_rXK_r^\dagger,
\qquad z=e^{-s}.
\]

This is a tilted CP map, not a trace-preserving time evolution. Its positive
Perron effect is diagonal, `h=diag(h_x)`, and obeys
`E_z(h)=E_z^†(h)=lambda h`. The native population operator is symmetric;
the signed-difference blocks are symmetric too. Positive tilting preserves
the primitive Kraus-product support. The strictly positive Perron effect
therefore has a simple dominant eigenvalue.

The source carrier's pair distances are 1, 2 and 3, with multiplicities
30, 30 and 6. The exact population quotient, acting on orbit probabilities,
is

\[
M_z={1\over60}
\begin{pmatrix}
52+4z&4z&0\\
4z&50+4z&10z\\
0&2z&50
\end{pmatrix}.
\]

In the orthonormal orbit basis it is symmetric:

\[
B_z={1\over60}
\begin{pmatrix}
52+4z&4z&0\\
4z&50+4z&2\sqrt5z\\
0&2\sqrt5z&50
\end{pmatrix}.
\]

Writing `u=60 lambda-50`, the largest root satisfies the exact equation

\[
u^3-(2+8z)u^2+(8z-20z^2)u+40z^2+80z^3=0.
\]

The transformed complete Kraus family is

\[
\boxed{\widetilde K_r
=\lambda^{-1/2}z^{a_r/2}h^{1/2}K_rh^{-1/2}.}
\]

It is trace preserving because

\[
\sum_r\widetilde K_r^\dagger\widetilde K_r
={h^{-1/2}\mathcal E_z^\dagger(h)h^{-1/2}\over\lambda}=I.
\]

Its unique stationary state is

\[
\rho_z={h^2\over\operatorname{Tr}h^2}.
\]

Substitution proves stationarity without relying on an eigenvalue tolerance.
The generic discrete-time quantum Doob construction is an established
trajectory-conditioning method; this calculation supplies its complete
native instrument, quotient, reference, clock and failure-record realization.
See Cilluffo et al., *Microscopic biasing of discrete-time quantum trajectories*,
[arXiv:2007.15659](https://arxiv.org/abs/2007.15659).

At `z=1`, `h=I`, `lambda=1` and the instrument is exactly native. The code
uses this exact identity instead of allowing a numerical Perron vector to
manufacture spurious tiny filtering failures. At `z=1/2`, the orbit masses
are approximately `(0.5961232608,0.3470516788,0.0568250604)`; at `z=1/4`
they are `(0.8018503763,0.1795448370,0.01860478675)`. These are forward
results for declared dimensionless controls, not fitted physical states.

The nonuniformity concerns pair-distance correlations. Every `rho_z` remains
invariant under the native rotations, so every port has mean load
`2/12`. Applying the unconditioned one-body chart transducer therefore gives
`I_3/3` for every displayed tilt, not a flavor hierarchy. The anisotropic
references in section 2 require the actual oriented record to be retained
in the observer's conditioning. A scalar activity tilt alone does not
replace that operational conditioning or select an interacting flavor law.

### Complete native implementation with failures

The changed Kraus coefficients are not passed off as the unchanged repair
instrument. For `h_min<=h_x<=h_max`, use the input filter

\[
B=\sqrt{h_{\min}}\,h^{-1/2},
\]

then the unchanged native repair, the activity acceptance factor
`z^(a_r/2)`, and the output filter

\[
C=h^{1/2}/\sqrt{h_{\max}}.
\]

The failure operators are the positive complementary square roots at each
stage. The complete outcomes are:

- input-filter failure;
- for every native raw label, activity-filter failure;
- for every native raw label, output-filter failure;
- for every native raw label, full success.

There are 721 C or 1441 E outcomes, all retained. The total success map is

\[
\rho\longmapsto\alpha\,\widetilde{\mathcal E}_z(\rho),
\qquad
\alpha=\lambda h_{\min}/h_{\max},
\]

with state-independent success probability `alpha`. The whole instrument
is trace preserving. Its conditional successful channel is the Doob channel.
The unconditional native-plus-filter channel includes the failures and is
not silently identified with that successful channel. A deterministic
implementation of the Doob coefficients would be a further controlled
instrument with those explicit coefficients.

### What a variational principle would select

On classical stationary paths retaining every native event label, let
`q_r(x)` be an alternative conditional event law and `pi` its stationary
load distribution. The native probability for each allowed event is `1/60`.
For the Doob law,

\[
\widetilde q_r(x)={z^{a_r}\over60\lambda}
                  {h_{f_r(x)}\over h_x}.
\]

Stationarity makes the mean increment of `log h` vanish. Hence the exact
relative-entropy-rate identity is

\[
\mathcal D(q\Vert q_0)+s\langle a\rangle_q
=\mathcal D(q\Vert\widetilde q)-\log\lambda.
\]

This selects the Doob law uniquely once `s` and the allowed event grammar are
specified. It is a genuine path-law optimization; it is stronger than
choosing a state by one-time relative-entropy projection. It does not
select `s`, identify a cost multiplier with `P`, distinguish C/E coherence
backaction, or fix a physical event duration.

For finite histories the stationary Doob law includes Perron boundary
factors. The exact operator product is

\[
\widetilde K_{r_n}\cdots\widetilde K_{r_1}
=\lambda^{-n/2}z^{A_n/2}h^{1/2}
 K_{r_n}\cdots K_{r_1}h^{-1/2}.
\]

Thus it must not be confused with the unmodified finite-horizon native
distribution multiplied only by `z^A` and renormalized; endpoint factors
matter. Every attempt count and every boundary state is specified in this
extension.

## 4. Native clocks and the Euclidean transfer

Let `A_n` be the count of actual nonzero currents in `n` attempts. At the
native stationary state,

\[
\mathbb E A_n={5n\over33}.
\]

The count is **not iid Bernoulli**. The tilted Perron eigenvalue gives the
scaled cumulant function `theta(s)=log lambda(e^-s)`. Implicit differentiation
of the exact cubic gives

\[
-\theta'(0)={5\over33},\qquad
\boxed{\theta''(0)={1597\over11979}.}
\]

The iid variance would be `140/1089`, which is different. The same variance
is checked independently by the exact three-orbit Poisson equation. The
implementation also computes all six move-count probabilities after five
attempts, with their conditional quantum outputs, for each displayed tilt
and both C/E instruments. C and E have the same population/count statistics
and different coherence decay.

The three displayed Perron values are enclosed by exact rational real-root
isolation with width below `10^-60`. Each reported floating value is compared
as its exact binary rational to that enclosure. The receipt therefore
reports its actual absolute-error upper bound rather than treating a rounded
eigenvalue as an exact root. Conditional reference characteristic polynomials,
two-record probabilities, and preparation contraction formulas are checked
in exact arithmetic separately from their displayed numerical values.

All un-tilted C superoperator blocks have exact Gershgorin lower bound
`40/60`; E has bound `22/60`. For `0<z<=1`,

\[
\mathcal E_z=(1-z)\mathcal E_0+z\mathcal E_1.
\]

Here `E_0` is a positive diagonal operator on the matrix-unit Hilbert space,
and `E_1` is strictly positive there. Therefore `E_z` is strictly positive
and self-adjoint on Hilbert–Schmidt operator space. Its normalized transfer
has a well-defined finite self-adjoint logarithm,

\[
H_z^E=-\log(\mathcal E_z/\lambda_z).
\]

Its ground vector is `h/sqrt(Tr h^2)`. Equivalently, the Heisenberg Doob map
is self-adjoint in the inner product

\[
\langle A,B\rangle_{\rho_z}
=\operatorname{Tr}(\rho_z^{1/2}A^\dagger\rho_z^{1/2}B),
\]

and is intertwined with `E_z/lambda` by
`A -> h^(1/2) A h^(1/2)/sqrt(Tr h^2)`.

This is a concrete positive transfer and a canonical logarithm once the
native instrument and tilt are fixed. It is an operator-space Euclidean
construction. The arbitrary choice to regard one repair attempt as one
Euclidean time slab, or to realize `exp(-it H_z^E)` as coherent physical time,
is not an identity with the original recorded repair process. The latter
has irreversible state updates and a growing tape. Paper 27's additional
physical slab and continuum conditions are not satisfied simply by taking
this logarithm. The finite transfer gap is therefore explicitly not labeled
a particle mass in the receipt.

The reference spectra in section 2 similarly determine dimensionless
modular gaps `log(lambda_i/lambda_j)`. They do not determine a coefficient
converting modular time to observed proper time or energy. Supplying that
coefficient only when convenient would reintroduce the energy-unit problem.

## 5. Retained quantum records do not make a phase invisible

Rephasing a raw Kraus operator, `K_r -> exp(i phi_r)K_r`, leaves the committed
CP instrument and every function of its committed records unchanged. This
includes the complete activity clock. It is a different coordinate choice
on an unobserved environment until a particular coherent read of that
environment is specified.

The supplied witness keeps seam and coin choices committed. It prepares a
superposition of the two native inputs with one shared spectator load and
opposite occupancy of the first seam. For that seam and coin, the unchanged
and moving repair outputs coincide while the raw pointer differs. A declared
Fourier read on the actual four- or eight-dimensional repair pointer before
its commitment detects a sign on the moving raw branch. All pointer outputs
remain. The total-variation difference is `1/120` for C and approximately
`0.0100592231766` for E. If the repair pointer is committed first, both laws
coincide exactly. There is no interference between committed seam
or coin choices.

Thus full coherent access can test a specified dilation convention. Diagonal
clock and cost constraints cannot pick that convention. The source's actual
phase-exact pulse word fixes a convention for that configured implementation;
it does not prove it is the unique physical continuation of the committed
native law. This extension does not relabel an added Fourier read or an
environment phase as a source-selected interaction.

## 6. Attempt count, record capacity, and public P

Use three separate objects:

- `n`: number of attempted native repairs;
- `N_rec=log D_rec`: logarithmic capacity of the retained record register;
- `P`: the separate public coupling/pressure parameter in the source papers.

For any fixed binary initial configuration, the sixty possible seam/coin
choices give exactly `60^n` distinct complete record sequences: the choices
are part of the record. Across all 66 inputs, the number `R_n` of attained
raw histories therefore satisfies

\[
60^n\le R_n\le\min(m^n,66\,60^n),
\]

where `m=120` for C and `m=180` for E counts attained one-step raw labels.
Lossless retention requires `N_rec>=log R_n>=n log 60`. A fixed-alphabet
compressed tape of size `m^n` is sufficient; retaining the literal padded
pointer register instead uses `240^n` or `480^n` states. The latter is an
implementation choice, not the information-theoretic minimum. The system,
purifier, additional transducer records and any other apparatus require
separately accounted capacity. None of these bounds identifies the public
closure value of `N` with a chosen attempt count.

**Prefix no-deformation theorem.** Suppose a finite tape can hold at least
`k` events, and no capacity-dependent control acts before it fills. Its
`k`-record instrument is exactly the same composition
`K_(r_k)...K_(r_1)` for every larger capacity. Appending any trace-preserving
future instrument and summing its outcomes preserves that prefix law,
including any spectator reference. This is the repeated identity
`sum_r K_r^†K_r=I`, not an asymptotic estimate. Capacity alone therefore
does not deform any fixed-prefix probability or native transfer spectrum.

Finite-capacity effects can occur after a specified halt, overwrite,
postselection, resource-conditioned control, or capacity-dependent endpoint
law is introduced. Different such laws need not have the same correction.
In particular, conditioning on a future count or a nonoverflow event changes
the ensemble; it cannot be treated as an unconditional prediction of `N`.

The implemented native instrument contains neither `P` nor `N`. Identifying
`z=exp(-P)` is an available **candidate law**, but no source derivation
in this packet equates the cost multiplier with that parameter. Reporting
the family as conditional exposes exactly what a future stationarity or
synchronization theorem must determine.

## 7. Stationary information projection is a different question

For a declared Hamiltonian `H`, faithful prior `tau`, commutation constraint
`[rho,H]=0`, and specified mean energy, minimizing `D(rho||tau)` gives

\[
\rho_*=Z^{-1}\exp\bigl(\mathcal P_H(\log\tau)-\beta H\bigr).
\]

The spectral pinching is applied to `log tau`, not replaced by `log P_H(tau)`.
This is a valid constrained information projection and its multiplier is
fixed by the chosen moment when the constraint is feasible. It does not
select `H`, the prior, the moment, the physical time scale, or an autonomous
preparation mechanism. Many complete instruments have the same stationary
state. Here the actual native path law is used first, and the activity
relative-entropy-rate problem selects a changed instrument only after its
cost multiplier has been supplied. The two variational constructions answer
different operational questions and should not be exchanged.

## Receiving files and reproducibility

- `native_record_dynamics.py`: full finite instruments, count outputs,
  conditioned references, history laws and complete filter implementation.
- `results.json`: source hashes, exact formulas, zero-label counts,
  normalization checks and source-status declarations.
- `complete_native_instruments.npz`: full raw alphabets, all 3960 native
  event columns per candidate, Doob amplitudes, all 4356 transfer eigenvalues
  for each displayed control, full transducer probabilities, exact two-record
  probability/port numerators and exact conditional ranks.

Run from the project workspace:

```bash
python oph_interacting_extension/source_selection/native_record_dynamics.py --repo executable_audit_repo
```

All finite matrix dimensions, padding outcomes and conditional services are
part of the recorded receipt. No numerical small-mode cut is used. This
packet supplies an operationally selected reference **after observation**
and a fully specified native-conditioned dynamical family. A common physical
state-selection, coupling, chiral continuum and energy-clock construction
is required before any of its eigenvalues can be interpreted as the
measured particle spectrum.

## Executable contract

The live source bytes are checked against the immutable source manifest before
reconstruction. Equal-event counting functions require the unchanged native
amplitudes and reject a weighted or Doob-transformed instrument. The latter
remains available to the full channel and record-count APIs.

The committed certificate contains exact integer counts, symbolic identities,
rational Perron enclosures and verified numerical bounds. Full transfer spectra
and conditioned arrays are retained separately. Their floating entries are
checked with absolute tolerance `3e-12`; integer entries are checked exactly.
Small eigensolver differences do not change the exact certificate. The explicit
Perron approximation bound is `1e-14`, verified against its rational enclosure.
These numerical controls do not constitute interval arithmetic for every
entry of the retained array packet.
