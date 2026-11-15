# A native reference state supplies a modular frequency

This packet connects the existing native event statistics to a programmed,
observer-relative modular readout, including a source-compiled flow on the
full native `M66` algebra. The faithful Palm reference is
`sigma = diag(degree/600)`, derived in
[`event_circulation`](../event_circulation/README.md). Its two eigenvalues have
ratio `5/4`. The ratio is a degree census of the carrier: an occupied pair
leaves 8 moving branches when its ports are adjacent, since their shared seam
cannot move, and 10 otherwise, so `5/4 = 10/8` and the frequency is a
combinatorial identity of the icosahedral incidence, not a dynamical result. The reference therefore supplies a nonzero modular frequency
`log(5/4)` without choosing a Gibbs state from an externally supplied
Hamiltonian. This is a finite reference-state attachment, not a recovery of
interacting spacetime or a new comparison with experimental data.

For the active-corner Ramsey experiment below, the reference is **not the
probe preparation**. That probe begins at
`I66/66`, uses the existing fixed-seam signed-current herald, and prepares a
coherent mixture with the compiled native balanced mixer. Each herald sign
has probability `5/66`; failed preparation has probability `28/33` and remains
in the raw record. The ten orientation pairs have degree census
`(8,10):2`, `(10,8):2`, `(8,8):2`, `(10,10):4`. Thus four pairs see a
nontrivial relative modular phase and six do not.

## Exact relation and executed computation

With Heisenberg convention `sigma^(is) A sigma^(-is)`, the dual state action
uses `sigma^(-is)`. On the full native algebra,

```
-log(sigma) = log(60) I + log(5/4) P_degree8.
```

The twenty-dimensional active corner commutes with the reference. There the
implemented diagonal target is `D(theta)=exp(i theta P_degree8)`, up to the
common scalar phase, with `theta=s log(5/4)`. Normalizing the restricted
reference changes only that scalar. Its full-reference weight is `23/75`.
That experiment's extension assigns zero phase to the remaining native addresses;
it omits 22 other degree-eight states. **It does not implement the modular
flow on the full 66-dimensional native algebra.**

Each observer program executes:

1. Native signed-current herald, retaining failure and both signs.
2. On success, the 34-pulse compiled balanced mixer `B`.
3. A fresh compiled diagonal source word `D` on seven work bits, lifted by
   two spectator bits to all 512 encoded and padding coordinates.
4. The same 34-pulse mixer `B`, then the native pointer with sign-aware
   decoding. Failure receives its explicit read record and no phase cost.

The diagonal compiler acts on eight declared degree-eight addresses among
the twenty active addresses. All source matrices are raw binary pulse
products, with no polar projection. Padding remains in the complete runtime;
native herald padding is failure and pointer padding has an unbiased bit.
There are six joint herald/read records, including the zero-probability
failed-preparation/read-one record. The C and E candidate laws are both run.
The reference preparation itself reuses the separately executed Palm-herald
construction; this packet does not rerun that experiment as a probe.

All settings consume the actual `runtime.Program` events, controls and
Kraus arrays with a factor executor: the 66 initial basis columns retain
their original binary weight `1/66`, and every hidden Kraus alternative
adds a column. Final densities are weighted sums of outer products on all
512 coordinates. This avoids subtractive cancellation at exactly zero
ideal probabilities without discarding small amplitudes, projecting states
or normalizing outcomes. At phase `pi/2`, the generic dense observer runtime
independently reproduces the full output and serialized checkpoint
continuation. Its relative positivity guard can reject cancellation noise
at nominally zero outcomes; it is not relaxed for this packet.

The exact ideal conditional sign-match law is

```
P(match | successful herald, s) = 4/5 + cos(s log(5/4))/5.
```

| Programmed phase | Modular parameter | Conditional match | Joint match, each sign |
|---|---|---|---|
| 0 | 0 | 1 | 5/66 |
| pi/2 | pi/(2 log(5/4)) | 4/5 | 2/33 |
| pi | pi/log(5/4) | 3/5 | 1/22 |
| 2 pi | 2 pi/log(5/4) | 1 | 5/66 |

These settings use exact rational-complex target entries `1`, `i` and `-1`.
No binary logarithm or exponential is treated as an exact source target.
The displayed decimal frequency is diagnostic only. Uniform reference
`I66/66` gives identity modular action and conditional match one at every
setting. Conversely, the native unread dissipative attempt has a different,
law-dependent contraction. Programmed modular control is not asserted to be
that natural waiting evolution.

## Error boundary and operational meaning

The certificate converts every original binary component of the local
128-by-128 phase map to an exact rational, compares every coordinate against
the exact diagonal target, and encloses its Frobenius error `epsilon`.
Consequently its operator error is at most `epsilon`, its half-diamond error
at most `delta=epsilon(2+epsilon)/2`, and its channel norm at most
`g=(1+epsilon)^2`. The spectator lift and arbitrary untouched quantum
references preserve these operator bounds. The existing eight-dimensional
mixer certificate is reused. For the whole raw `B D B` sequence the bound is

```
delta_B (1 + g_D g_B) + delta_D g_B.
```

Both mixer uses and numerical amplification are charged. Exact ideal native
instruments and adaptive failed branches contract this comparison. This
bound covers the compiled gates against their stated complete extensions;
it does **not** certify binary native Kraus production, subsequent floating
state propagation, or hardware service errors. Full terminal density and
probability comparisons and checkpoint replay are separately labeled
numerical diagnostics. Physical service discrepancy remains unbounded.
For portable replay the fresh exact phase comparison must pass the fixed
operator-error envelope `1e-10`; the mixer must pass `1e-11`. The exported
whole-sequence half-diamond allowance is outward rounded from those fixed
envelopes, rather than depending on platform-specific last bits.

The law is an ensemble clock candidate on the monotone window
`0 <= s log(5/4) <= pi`; one bit is not an exact clock reading, and sign and
period aliases remain. A noncommuting prepared probe is necessary: the
reference itself is stationary under its modular flow. Coherent services,
the probe design and selected phase settings are supplied. Calling `s`
an observer's candidate time needs no universal clock or conversion to SI
seconds. What remains is to establish consistency with continuing clock
records, other accessible correlations and interventions, and effective
field dynamics. This packet supplies a concrete reference-derived frequency
and testable finite readout law, rather than that larger identification.

The finite modular-inner attachment is the existing
[`PrivateInner.lean`](../../Lean/Dynamics/PrivateInner.lean); the current
construction uses explicit finite spectral calculus. The prior native
interface supplies the exact herald and pointer laws, while
[`compiled_native_interface`](../compiled_native_interface/README.md) supplies
the complete source mixer and its numerical comparison theorem. No additional
axiom or further abstract theorem is needed.

## Finite preparation of the faithful reference

[`preparation.py`](preparation.py) removes exact equilibrium preparation as
a premise for this reference approximation. For any initial binary-sector
state, including one entangled with an untouched reference, the existing
[native warm-up theorem](../record_laws/README.md) gives
`D(rho_nR, I66/66 tensor rho_R) <= delta = min(1, 33(19079/20000)^n)`.
The preparation tape and controller are unavailable to subsequent work;
the bound retains the theorem's exogenous seam-coverage hypotheses.

One subsequent **global random-seam movement herald**, the same CP success
map for C and E, sends `I66/66` to `q sigma`, with `q=5/33`. The complete
success/failure instrument gives `|p-q|<=delta` and conditional work-reference
distance at most `delta/q`. Its source effect also bounds `2/15<=p<=1/6`.
Failure is part of the experiment and its resource cost. This full
Palm herald is distinct from the fixed-seam probe herald, whose successful
support has dimension 20 and cannot prepare a faithful full-algebra state.

Put `t=66 delta`. Positivity gives the operator sandwich
`(1-t)sigma <= M(rho_n)/q <= (1+t)sigma`. For `t<1`, the actual successful
work state `sigma_n=M(rho_n)/p` is faithful, with minimum eigenvalue at least
`(1-t)q/[75 min(1/6,q+delta)]`. Operator monotonicity of the logarithm and
Duhamel's identity then give

```
||log(sigma_n)-log(sigma)+log(p/q) I|| <= t/(1-t),
half-diamond(modular_n(s), modular(s)) <= min(1, |s| t/(1-t)).
```

The scalar phase cancels; no commutation assumption is needed. At **369
native warm-up attempts plus one herald**, the eigenvalue floor exceeds
`0.01333244` and the channel error is below `0.000971` uniformly for
`|s|<=16`, including arbitrary probe references. Repeating the complete
preparation on failure costs at most `2442.015` expected native attempts.
Faithfulness concerns the successful work algebra, not a possibly singular
external reference or every classical record block. The explicit replay
also retains exact finite diagonal histories and failure masses.

Each such state admits the existing finite inner-flow theorem and
[`gibbsState_kms`](../../Lean/QFT/StructuralInheritance.lean) at
`H=-log(sigma_n), beta=1`; its Gibbs Heisenberg parameter is `-s` in the
modular convention used here. These are finite algebraic conclusions.
The full-support construction below approximates the actual `sigma_n` flow
without computing its unknown logarithm. Passive repair is not identified
with modular evolution, and neither attempt counts nor `s` are calibrated
to physical seconds.

## Full native modular flow from the same source services

[`full_flow.py`](full_flow.py) compiles
`U(s)=exp(i s log(5/4) P_degree8)` on **all 30 degree-eight native addresses**,
including the 22 omitted by the active-corner experiment. This is the dual
Palm modular action on every native state and coherence, up to the common
scalar `exp(i s log(60))`. The declared extension acts as identity on all
446 nonnative/padding addresses; it supplies no faithful modular state there.

The complete nine-bit Walsh expansion is verified at every one of the 512
addresses. Closed Gray parity walks execute it using 512 existing `I`/`Z`
phase requests and 510 existing phase-exact CNOTs: **2042 source pulses**.
Each CNOT uses the existing `H J H` word from the same proper-code source
controls. No new global phase primitive is required. Intermediate words
retain every encoded address: projecting onto the native 66-dimensional
shell during execution loses amplitude and is explicitly rejected by a
control. The certificate covers the complete terminal operation, not an
experiment that measures a partial word.

For any requested exact rational `|s|<=16`, rational series enclose
`log(5/4)` and its complex phase. Every entry of the actual raw 512-by-512
source map is compared to this target; the Frobenius enclosure must pass a
portable operator-error gate `epsilon=1e-9`. Thus the half-diamond bound is
`delta_word=epsilon(2+epsilon)/2`, with arbitrary untouched references.
The retained executions use `s=0,1,16`; other requested words must pass their
own numerical source certificate. The actual prepared-reference comparison is

```
half-diamond(source(s), modular_sigma_n(s))
    <= delta_word + min(1, |s| 66 delta_n / (1-66 delta_n)).
```

At 369 native warm-up attempts and one full movement herald this is below
`0.000971` throughout the declared modular window, for each certified word.
The complete preparation-and-source experiment, retaining failure, instead
has half-trace error at most `delta_n + p_max delta_word` against its ideal
Palm experiment. These bounds concern different comparisons; successful
preparation is never treated as free or certain.

The Ramsey witness uses native addresses 16 and 144, outside the active
corner, and predicts `(1+cos(s log(5/4)))/2`. Its initialization and readout
are complete native protocols:

1. A [current-only controller](../native_clock/configuration.py) resets the
   work to code 16. A separate native auxiliary is reset to the
   [pointer's ready configuration](../native_interface/pointer.py).
   Each reset permits arbitrary total-two inputs, including collisions;
   both public timeout flags remain. Private reset tapes are discarded.
2. Two existing source `Y` pulses surround the 2042-pulse modular word.
3. Six source pulses couple work bit 7 to the auxiliary orientation. One
   native coin/current repair reads it. The ideal reduced target is the
   continuing work-bit Lüders instrument on all 512 work coordinates;
   the literal source word retains its separate numerical error bound.

The maximum cost is **5357 native attempts, 2050 source pulses and 4100
proper-code service events**. The native work and pointer preparation needs
no external blank or random-seam service. With four seed scans and 32 retries per
mixed edge, the two complete reset errors sum to `140/2^32`. The ideal
coupling/readout theorem covers arbitrary work/reference inputs. Its raw
six-pulse map has operator error at most `1e-11`, half-diamond error
`delta_pointer` and norm bound `g_pointer=(1+1e-11)^2`. Thus

```
complete history error <= 140/2^32 + g_pointer*delta_Ramsey + delta_pointer
                       < 3.4e-8.
```

The norm factor accounts for possible amplification by a raw numerical
map. This compares the complete history with the ideal Palm Ramsey
experiment, including empty ideal timeout blocks. It is distinct from the
`0.000971` comparison to an unknown finitely prepared reference above.
The retained execution averages the exact native reset laws over independent
`I78/78` inputs; the uniform bound does not require that input or independence.
On a correlated input, successful reset can steer the reference. The proof
uses the complete failure branch and never assumes otherwise.

All four public outcomes are retained for both C and E: work-reset timeout,
auxiliary-reset timeout, bit zero and bit one. At `s=16`, the successful
match law is about `0.04524740944`. Padding is retained in both bit outcomes;
nonnative weight is a diagnostic, not a measured membership flag. Raw
outputs are not renormalized. The independent source-array, rank-one reset
and complete pointer-instrument checks establish composition on this common
finite apparatus.

The retained `full_flow_results.json` witness contains the four probability
marginals and trace/leakage diagnostics for the declared `I78/78` inputs.
The complete quantum-output and arbitrary-reference bound follows by
composition of the independently checked reset, source-array and pointer
instruments. The pointer also has separate reference-entangled executions.
Failure quantum channels are specified by the finite native controller;
the full-flow receipt does not enumerate their dense output matrices.
Replay checks the probability and diagnostic domains and the complete
record mass as well as the comparison bound. Portable numerical tolerances
do not permit negative probabilities or relax execution error gates.

Coherent proper-code load/export, their service ancillas and routing remain
architecture services; the native resets do not construct those service blanks.
Their physical errors and durations are not authenticated by source-array
arithmetic. Committed native repair from a joint diagonal accessible input
with fresh independent diagonal ancillas cannot generate the needed coherent
controls. No passive-repair/modular-time or
laboratory-clock identification follows from implementing this experiment.

## Regional algebras under one native modular action

[`regional.py`](regional.py) constructs regional algebras directly on the
native total-two binary shell. For a set of ports `R`, the bare algebra is
`A_R = direct_sum_q M(local q-configurations) tensor I(outside)`. The local
occupation records lie in its diagonal part. Every such classical record
commutes with the Palm reference, so its modular action is trivial. A moving
modular clock needs a coherent probe and the admitted control services.

The full Palm action does not generally preserve `A_R`: the phase on a local
matrix unit depends on the outside configuration. For the adjacent ports
`R=(0,1)`, one coherence has degree ratio `4/5` in one outside context and
`1` in another. Its reduced faithful state is tracial in each charge block,
so its own modular flow is identity. Faithful regional restriction alone
does not establish restriction of the ambient clock.

There is an exact, canonical completion. Group outside configurations when
their positive degree vectors over local configurations are proportional.
Each resulting block has `sigma=p mu tensor nu`, with both factors faithful.
The algebra `C_R = direct_sum_blocks M(local) tensor I(outside context)` is
modular invariant, and its state-preserving expectation is

```
E_R(X) = direct_sum_blocks [(id tensor nu)(P_block X P_block) tensor I].
```

Its Kraus formula gives complete positivity; the displayed factorization
gives unitality, the bimodule identity, idempotence and preservation of
`sigma`. It need not preserve the ordinary matrix trace. The Schrödinger
restriction to all local factors is CPTP and retains every context.

This is the **smallest modular-invariant algebra containing `A_R`**. Spectral
components of a local matrix unit separate its distinct degree ratios;
products with the reverse unit give the corresponding outside-context
projectors. Their joint partition is exactly proportionality of the degree
columns. Minimality proves `C_R subset C_S` whenever `R subset S`; the one
ambient modular group agrees on every common observable. For the cap
`R=(0,1,5,7)`, the algebra dimension grows from 53 to 149. This supplies
finite regional modular compatibility, with explicit additional context.
Its implementation by local records and coherent controls is an access
obligation. Disjoint completed regions can fail to commute; the receipt
exhibits a nonzero commutator for `(0,1)` and `(2,8)`.

Nested expectations compose: `E_R E_S=E_S E_R=E_R` for `R subset S`.
The bimodule identity and state preservation make each expectation the
orthogonal projection for the faithful inner product `Tr(sigma X*Y)`;
ordinary nested-projection composition then applies.

These explicit finite maps instantiate the relevant algebraic formulas,
not the entire [`CPRegionalNet`](../../Lean/QFT/CPRestrictionNet.lean)
interface: that interface additionally requires trace-preserving expectations,
locality and its other net laws. There is no claimed spacetime causal net.
The existing full source word implements the ideal action simultaneously on
all `C_R`. Regional partial trace preserves the finite-preparation operator
sandwich, giving the same modular-error bound for the actual regional states;
exact invariance of the actual finite-prepared ambient state is not assumed.

## One native entropy completion and its regional clock

[`completed_clock.py`](completed_clock.py) joins the
[native complement-repair entropy completion](../native_equilibration/cut_completion.py)
to a clock on the **same prepared state and bare cap algebra**. For the cap
above, run the admitted native C/E repairs on complementary seams. Their
derived tracial completion `E_R^trace` preserves every cap observable and
selects `bar_sigma=E_R^trace(sigma)`. The source degree counts give its five
eigenvalues `5/336, 19/1200, 37/2400, 1/75, 1/60`. It is the unique entropy
maximum at the retained Palm cap marginal. Its distance from Palm is `8/175`;
the tracial entropy completion and Palm-preserving expectation are distinct.

Here `bar_sigma` factors with a uniform complement in each charge block, so
its modular flow **does preserve the original 53-dimensional cap algebra**.
The same counted state gives its cap modular gaps:

| Cap charge | Nonzero modular frequency magnitude |
| --- | --- |
| 0 | None |
| 1 | `log(38/37)` |
| 2 | `log(5/4)` |

The dimensionless frequency ratio `log(38/37)/log(5/4)` has no fitted rate.
It is an internal prediction of this preparation and access protocol.
No physical temperature, atomic transition or SI time conversion is attached.

Two complete 2042-pulse words use the existing source controls: the direct
state action `U_bar(s)=bar_sigma^(-is)` up to the common native scalar, and
the relative action `V_s=bar_sigma^(-is) sigma^(is)`. The exact diagonal
identity is `U_bar(s)=V_s U_Palm(s)` in the same gauge. The Heisenberg relative
cocycle is `V_s*`; this relates two reference-dependent clocks without
identifying their frequencies. Rational logarithm/phase enclosures certify
every original binary source-array entry against its target. Target centers
are rounded to dyadic rationals with both component errors charged; the
executed arrays are unchanged. Each requested rational `|s|<=16` word must
pass operator error `1e-9`; the retained words use `s=16`.

Finite preparation also joins. Put `t=66 delta`, where `delta` is the Palm
warm-up error, and `u=75 eta`, where `eta` is the complement-repair error.
Complete positivity and `bar_sigma>=I/75` give the actual unnormalized
successful work state a relative sandwich of radius `r=t+u+tu`, up to its
herald normalization scalar. Logarithm stability gives

```
half-diamond(source_bar(s), actual_completed_modular(s))
    <= delta_word + min(1, |s| r/(1-r)).
```

At **369 warm-up attempts, one global movement herald, then 7440 complement
attempts on success**, this is below `0.001` for each certified word in the
window. The cap-only comparison uses `t` instead of `r`: every complement
attempt preserves the complete cap marginal. The complete flagged
single-trial preparation/source error is `delta + (5/33) eta + p_max delta_word`.
All failed herald outcomes remain records. Repeating failed Palm trials and
then running one successful complement stage costs fewer than 9882.015
expected native attempts. Preparation records/controllers are unavailable to
later work; attempts and `s` retain their different meanings.

The same clock also respects the continuing cut protocol. Every actual C/E
complement Kraus branch preserves the cap configuration, and therefore
commutes with `bar_sigma` and its modular unitary. This proves covariance
of the **complete labelled instrument**, including arbitrary untouched
references and any finite concatenation of its records. Ideal modular words
can be moved between complete complement attempts without changing those
current-record probabilities. Numerical source words carry their certified
control error; measurements during a partial word are outside this claim.

An executed cap Ramsey probe uses native codes 328 and 73, with the same
outside configuration and cap charge one. Seven existing source pulses
prepare their coherence and seven analyze it. The complete work-register
bit-0 PVM reads `(1+cos(s log(38/37)))/2`; all 512 coordinates contribute to
its two outcomes. At `s=16` the match probability lies between `0.95` and
`0.96`. The 2056-pulse experiment retains padding throughout; neither a
membership projection nor outcome normalization is inserted. Configuration
preparation and the projective bit read remain admitted apparatus, with
their physical costs outside this pulse count.

This supplies the finite faithful-state/inner-flow hypotheses used in
[`PrivateInner`](../../Lean/Dynamics/PrivateInner.lean) and the finite Gibbs/KMS
identity in [`StructuralInheritance`](../../Lean/QFT/StructuralInheritance.lean)
for the selected native entropy completion. The cut schedule, discarded
records, coherent services and cap choice are explicit. Passive cut repair
is not the modular unitary. Exact regional invariance belongs to the ideal
completed reference; actual finite-prepared flows have the stated error
bound. No boost, area, smooth spacetime or laboratory-clock identification
is inferred.

## One reference for overlapping caps

[`common_caps.py`](common_caps.py) supplies one native-selected reference
for a specified family of overlapping bare cap algebras. Fairly select
`R=(0,1,5,7)` or `S=(0,1,2,8)`, then a uniform complementary native seam
and a fair placement coin. The 20-seam union is connected on the ten ports
outside `I=R intersection S=(0,1)`. For both C and E, the exact native
channel has fixed algebra `A_I` and quantum gap at least `1/600` per
attempt. Its unread marginal converges to the tracial completion `E_I`.
Records may remain stored if subsequent work does not access them or use
them for feedback.

Applied to the same globally heralded Palm reference, this protocol selects
`sigma_common=E_I(Palm)`. Its per-configuration weights are

| Occupation of the protected core I | Weight | Configurations |
| --- | --- | --- |
| 0 | `68/4500` | 45 |
| 1 | `69/4500` | 20 |
| 2 | `60/4500` | 1 |

This single density lies in `A_I`, hence its modular group preserves every
bare `A_R` containing `I`; all 1024 such cap algebras are checked. It is
central in `A_I`, so **the shared intersection's own modular action is
the identity**. Nontrivial clocks require larger-cap coherences transferring
occupation between the protected core and another cap port. The mixed
preparation preserves the core marginal, not the individual larger-cap
marginals. The state and clock are common across these caps without
independently selecting a reference for each observer.

One 2042-pulse source word implements the common modular action. Two
cap-supported selected-pair Ramsey probes use it unchanged: codes 328/73
for `R`, and codes 0/73 for `S`. With their preparation/analyzer mixers,
they cost 2056 and 2068 pulses and both read
`(1+cos(s log(69/68)))/2` through the complete 512-coordinate bit-0 PVM.
At the retained `s=16`, the ideal probability is approximately
`0.986421907373`. The complete binary source arrays and both mixer errors
are certified; padding is retained without projection or renormalization.
Configuration preparation, coherent source services and the projective bit
read remain admitted apparatus. The selected pairs have the same cap
exterior; this does not assert that each source mixer is a local operation
on arbitrary exterior inputs.

With 369 Palm warm-up attempts, one global movement herald, and 18,600
mixed-cap attempts after success, the same relative-sandwich proof bounds
the source/common-actual-state modular half-diamond discrepancy by
`0.000989302` at `s=16`. The preparation costs fewer than 21,042.015 expected
native attempts. The complete flagged single-trial bound retains every
herald failure. Each other requested rational `|s|<=16` requires its own
word certificate; a continuous family of floating implementations is not
certified by the retained word. Regional invariance is exact for the ideal
common reference; finite-prepared modular actions have the stated error.

The all-region boundary is also exact. On the binary two-occupation shell over
`N>=3` ports, a faithful diagonal reference with weights `w_ij` preserves
**every bare two-port algebra** iff `w_ij=lambda*r_i*r_j` with positive
factors. Equivalently, every four distinct ports satisfy
`w_ij*w_kl=w_ik*w_jl=w_il*w_jk`. Local matrix units require their modular
ratio to be independent of the spectator port; factorization gives the
converse. The generator is a scalar minus `sum_i log(r_i) N_i`.
Native Palm, the one-cap completion and this common-core completion fail
that all-pair condition. A connected symmetric native repair preserves a
factorized reference only when all factors agree. A concrete unprotected
two-cap schedule has a connected 28-seam union and uniquely selects
`I66/66`, whose stationary modular action is trivial.

These are a constructive protected-core family and an exact boundary for
the indicated regional algebras. The protected core and preparation protocol
are specified choices. No arbitrary spacetime region net, physical boost,
laboratory time or continuum geometry is inferred.

## A common observer from two local stopping clocks

[`local_preparation.py`](local_preparation.py) connects the regional stopping
theorem to the common-cap completion. On the same binary total-two carrier,
count incidents at the shared core `I=(0,1)`, including zero currents. After
the declared warm-up count, retain both outcomes of the **next specified**
incident movement flag. This is not the first successful movement. On
success, run the mixed complementary repair of the two caps above, stopping
after the declared number of incidents at the first cap. The core has no
incidents during this complementary repair, so it cannot serve as that
second counter.

The first observer has incident fraction `3/10`; the second has fraction
`5/24` during mixed completion. The complete ideal herald succeeds with
probability `5/33`. Its completion selects one common density:

| Core occupation | Per-configuration weight | Multiplicity |
| --- | --- | --- |
| 0 | `2/225` | 45 |
| 1 | `1/36` | 20 |
| 2 | `2/45` | 1 |

The nonzero elementary modular gaps are `log(25/8)` and `log(8/5)`. These
follow from this local preparation; the global-Palm gaps are not substituted.
The density lies in the protected-core algebra and is faithful. Its one
inner modular action therefore preserves all 1024 core-containing bare cap
algebras. The common core's own action remains trivial. This instantiates
the finite faithful-state, inner-action and Gibbs/KMS contracts with
`H=-log(sigma)` and inverse temperature one in modular units; it does not
turn the unread repair semigroup into that unitary action.

Both stopping clocks have finite controllers. Every tick has a bounded
remote-wait deadline. The existing reusable native auxiliary supplies each
scheduler choice through its coin-only instrument, whose effect is `I/2`
even when the auxiliary is correlated with the work. No external random
service, auxiliary reset or independent auxiliary state is assumed. A
bounded rejection selector retains its failure. The two complementary
menus have different sizes, so their finite acceptance factors differ.
They are included in the actual channel rather than canceled by conditioning.

For either phase let `A` be the accepted incident subchannel, `R` the
accepted remote subchannel, and `t` the state-independent selector failure
probability. For a per-tick deadline `h`, the executed success and timeout
maps are

```
T_h = A (Id + R + ... + R^(h-1)),
F_h = t (Id + R + ... + R^(h-1)) + R^h.
```

They form a complete instrument. Repeating `n` ticks gives success `T_h^n`
and terminal failure `F_h sum_(j=0)^(n-1) T_h^j`, preserving the actual
disturbed failure state. The code evaluates these maps by binary powering,
without enumerating exponentially many fine histories. In the uncapped
comparison, the positive-contraction similarity proof applies on the
orthogonal complement of the common fixed algebra, replacing the tracial
projection by `E_I`. Finite selector and waiting losses enter the complete
error **before** any conditioning or logarithm stability bound.

The retained flags are warm-up timeout, herald failure, completion timeout,
and success. Preparation seam/coin/current tapes may be stored but are not
accessible to later feedback in this comparison. The controller needs
private bounded loop/deadline counters; the public preparation clocks count
local incidents. This removes the global warm-up timer and separately
selected cap references, not the need for controller hardware.

The retained precision uses **496 warm-up ticks, one specified herald tick
and 19,200 completion ticks**, with a 256-slot deadline per tick and 32
selector rounds. Its complete preparation distance is below `5.682e-8`;
conditional success distance is below `7.500e-7`. Expected work plus auxiliary
cost is below **105,445 native calls**. The conservative finite maximum is
659,460,352 calls. The reusable controller-register subtotal is 51 bits;
program tables, interpreter storage, private tapes and quantum registers
are excluded from that subtotal. These are bounds for the stated controller,
not optimal costs or physical durations.

## Continuing source-clock reads of that prepared state

[`local_history.py`](local_history.py) executes two coherent interventions
and native current reads on the **same prepared mixed carrier**. It does
not reset to code 328 or code 0 between probes. Each exact two-address
rotation fixes every other native and padding address; the cheaper
unconditional XOR mixers valid for selected pure inputs would not implement
this mixed-input instrument. Both source words use the same local reference
and modular parameter. The second rotation's sign depends only on the first
retained signed current. Both actual C/E coin/current instruments retain
their zero outcomes and continuing quantum states.

The complete 512-address source arrays are certified entry by entry against
rational enclosures of their exact targets. Every native export also has
an explicit complement abort, on which the work is discarded and an
untouched reference is retained. No escaped probability is renormalized
away. The receipt joins the finite preparation maps to this history and
compares the complete quantum outputs, including all three preparation
failures, both export aborts and every native read label. It reports the
full source-pulse/service budget separately from preparation's work and
auxiliary attempts. The all-input/reference guarantee is compositional;
the retained numerical preparation witness starts from one explicit basis
state and does not replace that proof.

At retained modular parameter one, the full preparation/history bound is
below `6.32e-8`. The two-probe source word uses **31,768 pulses and 63,536
proper-code service events**, plus two native current reads and two complete
native exports. Its single-clock comparison with the actual finitely prepared
state's modular action has a separate bound below `8.437e-5`. The source-array
enclosures and analytic preparation bounds are certificates; intermediate
floating state propagation and its reported residuals are diagnostics.

Coherent proper-code load/export and source controls are specified
architecture. Selecting these two address pairs requires exterior context;
this is not a derivation of cap-local coherent control on arbitrary
spectators. No projective work-bit read or independently selected pure probe
is needed. The comparison supplies neither physical service fidelity nor a
laboratory time, temperature, atom or field dictionary.

## Bare-cap clock operations on arbitrary spectators

[`cap_local_clock.py`](cap_local_clock.py) implements a coherent occupation
transfer uniformly across all ten spectator configurations of each read seam,
`(1,7)` and `(1,8)`. The generator
`J_ab = sum_k (|a,k><b,k| - |b,k><a,k|)` belongs to the bare two-port algebra;
its `pi/4` rotation therefore requires no selected exterior configuration.
The same prepared reference has weights `2/225`, `1/36`, `2/45` in the three
protected-core charge sectors. Its modular action belongs to both cap
algebras, and each native C/E read instrument belongs to its seam algebra.
These exact memberships establish locality of the completed ideal operations
for arbitrary exterior states and untouched references. The unread operation
preserves the outside marginal. A retained local outcome may conditionally
steer an initially correlated outside system.

Each isolated read of this common state has the same signed-current law:

```text
E[j] = (17/200) cos(s log(25/8)) + (1/120) cos(s log(8/5))
P(j != 0) = 181/900.
```

The nine spectators outside the protected core supply the first frequency;
the remaining core spectator supplies the second. No extra frequency or
initial pure probe is fitted. Including the ideal herald probability `5/33`,
the complete-trial expected current is strictly decreasing for
`1/2 <= s <= 1`, with derivative magnitude greater than `1/200`.
The inverse expected-current calibration is thus 200-Lipschitz on that
interval. This is a calibration of the mean: individual random records are
not monotone, opposite modular orientations give the same isolated law, and
empirical sampling and physical service fidelity are separate requirements.
The inner action has the finite [`PrivateInner`](../../Lean/Dynamics/PrivateInner.lean)
contract on the derived cap algebra. The inverse mean-response result does
not supply the strictly monotone individual-record map required by
[`ClockReadout`](../../Lean/Time/ClockReadout.lean).

The continuing protocol uses both caps on the same carrier, retains every
C/E outcome, and uses the first signed current to choose the second analyzer
rotation; its preparation rotation stays fixed. The second marginal is computed from its disturbed state;
it is not asserted to obey the isolated-reference formula. All preparation
failures and both native-export aborts remain in the complete output law.
At retained modular parameter one, the joint preparation/history bound is
below `6.32e-8`. A closed-Gray conditional-Y compiler certifies all original
512-address arrays, including padding, using unmodified full source-factor arrays
and **45,656 source pulses and 91,312
proper-code service events**, plus two native current reads and two complete
exports. Preparation's native calls are counted separately.

Exterior-state selection is removed from the ideal completed instrument.
Its encoded coherent controls and load/export services remain admitted
architecture; the intermediate pulses are not proved spatially local.
The source-array errors bound the departure from the ideal local action.
Neither a signal cone nor a laboratory-second, temperature or particle
identification follows from this operational clock calibration.

## Learn the common reference from one native current setting

The complete local preparation followed by one unread core repair supplies
three charge weights with multiplicities `(45,20,1)`, in the two-moment
maximum-entropy family of
[`common_variation.py`](../native_equilibration/common_variation.py).
The finite preparation theorem bounds coherences and initially correlated
references. Complement equilibration alone would leave additional core
information and would not justify this reduction.

[`reference_tomography.py`](reference_tomography.py) identifies this family
from one signed-current setting on `(1,7)`. For complete-trial frequencies
`u,v` of positive/negative successful current and success frequency `t`,

```text
b0 = (22 v + 2 u - t)/54
b1 = (t - 10 v - 2 u)/6
b2 = 5 u + 15 v - 3 t/2.
```

The scalar success mass cancels from modular ratios. The signed-score inverse
is statistically expensive: its retained `2^24`-trial clock has a two-read
inference allowance about `0.589`. Its confidence-inclusive complete bound
below `0.096` includes preparation failures and is not a successful-clock
error of that size. These remain reproducible alternate readout contracts.

### A finite occupation census from repair records

[`native_census.py`](native_census.py) derives a more informative readout
without adding an occupation projector or consulting hidden loads. It uses
eleven fixed native seams forming a connected tree. On each edge, read the
existing coin/current record until the coin first changes, with deadline
`h=32`; a timeout retains the disturbed carrier and its failure flag.

For binary endpoint loads, repetition of the first coin produces no further
current. At the first flip, its current reveals whether the endpoints were
unequal. If the first coin is `k`, first current `j`, and flip current `j'`,
the pre-edge contrast is `2j+(2k-1)|j'|`. Subtracting the difference of the
accumulated recorded transfers gives the **original** endpoint contrast.
The eleven contrasts and conserved total occupation two uniquely recover all
twelve original occupations. The controller receives only records; the
simulator's hidden configuration is used to check its answer. This training
readout spans all twelve ports under the declared seam-access policy; it is
not a cap-local measurement. The later completed clock reads retain their
separate bare-cap locality certificate.

Each coin has effect `I/2` for both C and E, including on a reference-entangled
input. Thus an edge completes with state-independent probability
`a=1-2^(1-h)`. Every successful fine record has only one possible initial
basis state, so its effect is rank one. Summing records decoded to original
core charge `q` gives exactly

```text
E_q = a^11 Pi_q,       E_timeout = (1-a^11) I.
```

These are **effects**, not an assertion that the postmeasurement operation is
Lueders or nondisturbing. Successful census destroys occupation coherences;
timeout branches retain their actual native disturbance. The executable
checks the C/E branch identities on all 4,356 matrix units, all 135,168
initial-state/compressed-coin paths, complete mass and untouched references.
The maximum read cost is 352 native attempts; expected cost is at most 33.
An explicit logical register layout uses 139 bits for the census, 168 for
eight training counters and 21 for the trial counter. These are subtotals,
excluding the preparation controller, quantum carriers, program/tree tables,
interpreter and private tapes; they are not measured Python memory. A fixed
tree schedule and finite storage are controller architecture inputs.

### Confidence and a successful clock history

Every preparation trial, including failures and census timeouts, contributes
to the counts. For the three successful charge frequencies `f_q`, the unbiased
ideal weights are estimated by `b_q=f_q/g_q`, where `g=(45,20,1)`.
Their common unknown factor is preparation success times `a^11`; neither it
nor the known reference weights center the estimate. This direct nonnegative
estimator avoids cancellation of large current counts to estimate a rare
charge. It identifies the prepared family, not arbitrary coherent states.

The reference-uniform preparation bound supplies a conditional probability
bias `e` after every prior public trial history. For fixed `N` and confidence
parameter `k`, bounded Bernoulli martingale Bernstein gives

```text
|f-p| <= e + 2k/(3N) + sqrt(2k(p+e)/N).
A = e + 2k/(3N), c = 2k/N,
r(f) = A + c/2 + sqrt(c(f+e+A) + c^2/4).
```

Substituting `p<=f+|f-p|` proves `|f-p|<=r(f)`, with a rational upper square
root in the implementation. A union bound covers all three charges with
failure probability at most `6*2^-k`. No independent copies, observed-variance
plug-in assumption or optional stopping is used. Finite-preparation error and the
separately audited dyadic-kernel bias enter `e`; intervals reaching zero halt.

[`census_clock.py`](census_clock.py) fixes `N=2^20`, `k=10` and a successful
history error cap `1/10` before reading data. For a faithful confidence box,
it compiles the inferred ratios with the existing source services. If
`h_i=log(w_hat_i/w_i)`, the generator error modulo identity is bounded by
`max_{i!=j}(upper(h_i)-lower(h_j))/2`. One coordinate cannot take its two
endpoints simultaneously. The resulting channel bound covers arbitrary
work/reference inputs. It composes through both continuing cap reads,
record-selected feedback, raw 512-address source arrays and export aborts.

For preparation error `B`, ideal herald probability `q`, raw-source norm cap
`G` and two-read channel allowance `D`, the distinct budgets are

```text
successful prepared history: G*2B/(q-B) + D
complete flagged history:    G*B + q*D.
```

Admission tests the first quantity, without multiplying by the herald
probability. Source certification failure and failed admission retain explicit
stops and the same carrier. The complete joint protocol uses the fixed cap
`G*B+q*(1/10-G*2B/(q-B))` plus `6*2^-k*(G+1)/2`, below `0.021011`;
both comparison protocols keep the same training
transcript. The confidence failure probability is not conditioned on admission
or successful preparation. Preparation conditioning is charged by
`2B/(q-B)`; it does not grant a postselected statistical guarantee.

The retained census run gives a successful two-read history allowance below
`0.075` using sixteen times fewer preparations than the signed-current run.
It continues C/E histories from the final training carrier, retaining every
failure. Its worst native training budget is below `6.915e14` attempts and
its expected budget below `1.107e11`, including auxiliary selectors. These
remain large physical costs. The demonstration samples an exactly audited
53-bit finite record/state kernel; it does not execute every microscopic
preparation attempt or authenticate the seeded PRNG as physical randomness.
All terminal counts, terminal-state counts and the stream hash are retained.

The positive inferred state instantiates faithful-reference, inner-flow and
finite Gibbs/KMS interfaces on the same core-containing cap algebras. The
source-derived family supplies the state/algebra compatibility premise; its
inferred moments specify the entropy family. These are analytic finite
connections with executable certificates, not new Lean formalizations.
Identification needs native records; clock enactment uses declared
coherent services. Physical service fidelity, intermediate-pulse locality,
laboratory seconds and a field or particle dictionary are separate.

## Infer an oriented clock parameter from complete native records

[`phase_readout.py`](phase_readout.py) turns the mean calibration above into
a finite-record inverse with an orientation. After each cap Ramsey stage,
read seam `(1,b)`, with `b=7` or `8`, recording coin/current `(c,j)`;
then read the core `(0,1)` on that disturbed carrier, recording `(d,k)`.
The public score is

```text
X = 2 j 1[d=c, k=-c],        X in {-2,0,2}.
```

Both reads execute even when the first current is zero. No selected outcome
is divided out: every unmatched coin, zero current, preparation failure,
export abort and numerical abort remains in the trial count, with score
zero. The two native reads have the exact score effect on their input
`(1-N0)(N1-Nb)/2`; the effect of its squared score is
`(1-N0)|N1-Nb|`. The native second read therefore removes the unwanted
core-spectator frequency from this statistic without a separate occupation
measurement or a selected exterior configuration.

The two analyzer settings use the existing real seam rotation and its
quadrature `exp(-i pi X_ab/4)`, compiled from the existing pair-Y and
pivot-Z source controls. For the same prepared three-sector reference,
the two complete-trial means obey

```text
E[X_cos] = A cos(s Delta),    E[X_sin] = A sin(s Delta),
A = (9/2)(b1-b0),            Delta = log(bhat1/bhat0).
```

Here `b_q=q a_q` are the unknown ideal subnormalized reference weights,
while `bhat_q` are inferred from the earlier native training records.
The source word enacts that learned center ratio. The training confidence
box certifies `A>0`; the decoder neither receives nor fits the amplitude
or the programmed parameter `s`. The common amplitude cancels in the
ratio of the two means, and its positive sign fixes the orientation.
The exact native effects also give
`E[X^2]=9q(a0+a1)<=9q/20`, uniformly in the phase and normalized family.

[`phase_inference.py`](phase_inference.py) uses complete, fixed-length
quadrature batches and exact rational arithmetic. A uniform conditional
complete-trial half-trace allowance `T` gives mean bias at most `4T`
and second moment at most `V=9q/20+8T`. The complete source budget is
`GB+qE`, with preparation error `B`, raw source error `E` and source
norm cap `G`; numerical sampling requires its separately audited allowance.
The decoder rejects an allowance smaller than the source budget.
For `N` trials per quadrature, it constructs an outward rational radius
`r` satisfying

```text
N r^2 >= 2 k (V + 4r/3).
```

Martingale Bernstein then gives simultaneous quadrature coverage failure
at most `4*2^-k`, plus the training confidence failure. The radius adds
`4T` once as conditional mean bias; it neither multiplies that bias by
the trial count nor suppresses it by conditioning on successful heralds.
The conditional source theorem and exact sampler-column audit apply after
past public records, so independent input copies are unnecessary.
Two caps share one training event; their simultaneous failure allowance
adds the four quadrature tails without counting training twice.

A declared parameter window must lie inside a single right phase
half-plane, and the observed cosine confidence interval must be strictly
positive. Rational alternating-series arctangent enclosures give the
signed phase interval; dividing by the enacted learned gap gives the
parameter interval. Nonpositive reference orientation, unresolved cosine,
incompatible amplitude, excessive winding windows and intervals missing
the declared window produce an explicit unresolved result. The selected
two means are periodic under a shift `2 pi/Delta`; they cannot determine
unrestricted winding. This score alias does not claim that the entire
fine-record instruments have identical laws at that period.

The separate `reference_unit_parameter_interval` divides the same measured
phase by the confidence box for the unknown reference gap. It measures
calibration uncertainty in reference units; it does not assert that the
source word evolved under the unknown exact generator. This local inverse
does not pay a full-state action error for the rare third sector, because
its enacted frequency is the recorded learned ratio and its scored response
isolates that frequency.

[`phase_clock.py`](phase_clock.py) continues from the training carrier
through both quadratures and both overlapping caps, retaining complete
records, terminal-state counts, stream hashes and every failure flag.
The numerical demonstration samples an explicitly quantized ideal-stage
population kernel with the finite native preparer. Its all-column exact
score audit supplies the conditional bias and variance contract.
Separate histories execute and certify the unmodified full 512-address
source arrays, including padding and all export branches. The population
stream is not presented as a microscopic raw-source trajectory, and a
seeded PRNG run does not prove independent physical randomness.

The retained protocol uses `2^24` training trials and `2^22` trials per
quadrature in each cap, with the fixed window `[-1,1]`. The decoder receives
neither the programmed parameter nor the analytical reference weights.
For the C alternative at parameter `3/4`, the two intervals are approximately
`[0.6771,0.7966]` and `[0.6792,0.8013]`; for E at `-3/4`, they are
`[-0.7918,-0.6714]` and `[-0.8148,-0.6952]`. These displayed endpoints are
outward rounded; the receipt retains exact rational bounds. Every interval
contains its validation target and passes the predeclared width cap `1/5`.
The shared-training, both-cap statistical failure bound is `1/128` for
each alternative, without conditioning on favorable records. Confidence
coverage comes from the theorem; a fixed seeded example does not measure
its repeated-run coverage. There is no new public-data fit.

Including failed trials, each alternative costs at most
`22,127,817,621,766,144` native attempts for training and readout,
`383,158,059,008` coherent source pulses and `766,316,118,016` proper-code
service events. Those are conservative costs of the represented protocol,
not counts of pulses physically executed in the accelerated simulation.
The receipt also bounds retained record/state counters and explicitly
excludes program storage, private tapes, stream hashes and arithmetic
workspace from that counter subtotal.

The completed ideal analyzer and native reads belong to the cap algebra;
locality of each intermediate encoded pulse is unproved. The result
is a blind inverse measurement of an actively programmed coherent source
word. It does not derive autonomous native time, strictly monotone
individual records, laboratory seconds, or physical source fidelity.

## The boundary includes noncommuting references

[`compatibility.py`](compatibility.py) removes a possible loophole in the
all-region test. Let a faithful state on the binary two-occupation shell over
`N>=3` ports preserve every bare two-port algebra under its modular group,
and put `H=-log(sigma)`. For each port `i`, its occupation `N_i` belongs to
every `A_{ij}`. Thus `[H,N_i]` lies in the intersection over `j != i` of
those algebras. Their off-diagonal supports are disjoint, so that intersection
is diagonal; `[H,N_i]` also has zero diagonal and hence vanishes. All joint
occupation eigenspaces on this shell are one-dimensional, forcing `H` and
`sigma` to be diagonal. This includes `N=3`, where the intersection need
not be only `span(I,N_i)`.

Spectator-independent modular ratios give exactly the positive product
weights `w_ij=lambda*r_i*r_j` characterized above. Consequently
**every faithful reference**, not only a preselected diagonal one, preserves
all bare two-port algebras precisely in that product family. The executable
certificate supplies forbidden-entry witnesses for every native
off-diagonal Hamiltonian entry and exact independent constraints for the
remaining diagonal energies. Connected symmetric repair selects only
the tracial member at equilibrium. Coherences cannot evade this specified
all-region boundary; protected regional families and changed observable
algebras are different contracts.

## Continuing the reference through recorded repair

[`history_transport.py`](history_transport.py) supplies the reference update
when a repair crosses the protected core. For either fixed C or fixed E,
a bounded controller chooses each native seam from its complete earlier
seam/coin/current-or-contrast tape. On binary66 every admitted history has
`K_h = sqrt(w_h) U_h P_h`, where `w_h=2^(-length(h))`, `P_h` is a diagonal
domain projector and `U_h` is injective on that domain. All terminal histories
are retained. The executable policy grammar admits at most 12 attempts and
65,536 reachable decision nodes; the proof applies to every finite bound.

For any faithful diagonal reference `sigma`, define the **actual** output
reference `tau_h=K_h sigma K_h*`. It is faithful on `range(K_h)`, although
that range can shrink. If `q_h=Tr(P_h sigma)`, then on this corner

```text
tau_h^(is) K_h = w_h^(is) K_h sigma^(is),
(tau_h/(w_h q_h))^(is) K_h = q_h^(-is) K_h sigma^(is).
```

Indeed `[P_h,sigma]=0` and the eigenvalue at `U_h e_i` is `w_h sigma_i`.
The branch scalar cancels in conjugation. Thus the complete recorded repair
channel intertwines the input modular flow and the direct sum of the output
corner flows, also with any untouched reference. Sequential updates compose
by restricting the injection and multiplying the coin weights; no fresh
reference preparation is necessary. Relative modular cocycles for two
diagonal references transport in the same way. Zero corners have no
logarithm or normalized state. The direct sum describes committed classical
records: this argument does not cancel relative phases between coherent
record blocks, nor identify a correlated probe's joint modular flow.

The binary carrier gives a canonical full-shell extension. A zero current
fixes the configuration; a nonzero current swaps the two seam bits.
Their recorded transpositions compose to a port permutation `pi_h`, and
`K_h=sqrt(w_h) U_pi_h P_h`. The faithful extension
`sigma_h=U_pi_h sigma U_pi_h*` obeys
`sigma_h^(is) K_h=K_h sigma^(is)` on the entire branch domain. Its compression
is `tau_h/w_h`. This does not fill zero posterior directions with a fictitious
prepared state: the extension specifies a source word whose restriction is
the proved corner action. It updates recursively from each current record,
without another preparation or an executed permutation service.

For the three-sector local reference only its two-port core moves. The
executed policy needs precisely four reference words, with cores `(0,1)`,
`(0,8)`, `(1,5)` and `(5,8)`. These transpositions are not graph symmetries,
and two reached cores are nonedges; adjacency and cap locality do not follow.
The extension fails for general full78 loads: repair sends `(2,0)` to
`(1,1)`, whereas a swap gives `(0,2)`. The adaptive full78 theorem below
uses partial injections, not this binary permutation shortcut.

The existing finite preparation supplies an operator sandwich
`(1-epsilon)sigma <= rho <= (1+epsilon)sigma`, `epsilon<1`. Every unnormalized
branch preserves it on the same support, even when `rho` is non-diagonal.
Operator monotonicity of the logarithm bounds its log error by
`epsilon/(1-epsilon)`. Duhamel's identity bounds the corresponding modular
channel half-diamond error by `min(1, |s|epsilon/(1-epsilon))`. Applying this
on both sides gives approximate intertwining error
`min(1, 2|s|epsilon/(1-epsilon))`, uniformly over rare records; there is no
division by their probability. This comparison is between reference-defined
linear channels, not a claim that posterior normalization is linear. The
sandwich concerns the native work support; compiled leakage into padding is
charged separately.

[`history_clock.py`](history_clock.py) executes the connection. The finite
local preparer retains its original success, herald-failure and two timeout
flags. On success a compiled balanced probe, two record-adaptive native
repairs crossing the core, the transported clock, a conditioned analyzer and
a native current read act on the same carrier. All 512 source coordinates
and native-export aborts are retained. C has 512 and E has 4,096 three-attempt
coin/raw-register words. A checked lazy zero extension represents impossible
slots without allocating thousands of empty quantum matrices; missing
attained records are errors. The preparer's fine auxiliary tapes remain
coarse-grained under its existing herald contract. Repairs execute through
the declared native Kraus interface; the coherent clock/analyzer source words
are replayed here, not the separately constructed native-repair pulse words.

At parameter `s=4`, the ideal clock can be moved across both repairs without
changing any terminal quantum block. Reusing the unchanged endpoint
reference instead changes the complete record law by more than `0.0003`
for C and `0.0006` for E. Removing the clock changes it by more than `0.004`
and `0.005`, respectively. These are conditional-on-success reference inputs;
the executed complete preparation keeps failures and charges its original
success weight `5/33`. Separate entangled-input controls check the full
quantum outputs, not only these classical marginals. Exact counterexamples
show why erasing the coin, reusing a fixed reference, or substituting a
general noncommuting reference cannot be silently admitted.

The reference is transported by the native repair channel. The probe rotation
disturbs the tested input, so its posterior is generally different from that
transported reference. An unperturbed reference is invariant under its own
modular flow; it cannot exhibit the nonzero phase signal. The continuing
reference theorem applies across the native repairs, not across arbitrary
coherent interventions as if they left the reference unchanged.

Every original compiled clock/analyzer array is checked against a rational
enclosure of its target, including off-code directions. The native Kraus
frame is checked against exact endpoint maps. For operator error `e`, use
gain `g=(1+e)^2` and half-diamond allowance `d=e(2+e)/2`; composition updates
`B <- g B+d`. With preparation error `B_p` and complete history error `B_h`,
the end-to-end bound is `G_h B_p+(5/33)B_h`. The receipt includes separate
attempt, source-pulse, service, tape, code and program-storage costs. The
five numerical error observations are freshly checked against their scientific
bounds and a fixed `2e-10` ceiling, then retained as exact `[0,2e-10]`
admission intervals. These summarize computed residuals; they are not rigorous
interval-arithmetic or implementation-error bounds. Replay also limits the
total variation of each complete retained probability distribution to `2e-10`.
The joined instrument bound compares against the ideal reference-defined
protocol. Comparing its action with the unknown finite-preparation
reference's own modular flow additionally uses the sandwich bound above.
The source bound certifies terminal logical words and native exports. Exact
proper-code intermediate exports have zero failure; floating per-pulse
apparatus errors and host workspace are not certified by this bound.

This closes finite record-conditioned clock continuation, not autonomous
physical time, universal bare-cap locality or protected defect dynamics.
The [adaptive coherence bound](../observer/RECORD_RECOVERY.md#persistent-exploration-with-adaptive-classical-control)
also prevents interpreting formal corner clocks as indefinitely recoverable
quantum memory under exploring record-only control. Active recovery,
restricted exploration and refinement require their own compatible joins.
Completion-program gate acceptance requires its declared predecessor
contracts; this component does not discharge CP01, CP02 or CP13 as a whole.

The [finite native-only construction](../completion/README.md) has a different
access contract. Its [stationary-reference obstruction](../completion/MODULAR_ACCESS.md)
excludes nontrivial modular implementation on the tested fixed full algebra
using its real operations alone. The compiled complex phase services here
are additional declared controls, so this continuation does not discharge
that obstruction or implement that construction's joint stationary reference.
The adaptive coherence bound applies to its actual repair events when the
exploration floor is retained; identity holds do not count as repair events.
The same repair count must be guaranteed on every compared branch. Halted
or unfinished histories cannot be discarded or replaced by an expected count.

## Replay

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference --check model/native_modular_reference/results.json
python -m model.native_modular_reference.preparation --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.full_flow --check model/native_modular_reference/full_flow_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.regional --check model/native_modular_reference/regional_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.completed_clock --check model/native_modular_reference/completed_clock_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.common_caps --check model/native_modular_reference/common_caps_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.local_preparation --check model/native_modular_reference/local_preparation_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.local_history --check model/native_modular_reference/local_history_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.cap_local_clock --check model/native_modular_reference/cap_local_clock_results.json
python -m model.native_modular_reference.reference_tomography --check model/native_modular_reference/reference_tomography_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.learned_clock --check model/native_modular_reference/learned_clock_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.phase_clock --check model/native_modular_reference/phase_clock_results.json
python -m model.native_modular_reference.native_census --check model/native_modular_reference/native_census_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.census_clock --check model/native_modular_reference/census_clock_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.history_clock --check model/native_modular_reference/history_clock_results.json.gz
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/native_modular_reference
```

Controls cover finite local-selector/wait timeouts, initially entangled
references, noncommuting compatibility controls, the faithful reference and
4+6 pair census, omitted native
coordinates, all padding columns, inverse/common phase, scalar amplification
that defeats an unweighted telescope, retained preparation failures, uniform
reference, native waiting, checkpoint continuation and hostile ndarray
subclasses. The retained receipt is software evidence, not measured data.
