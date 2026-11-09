# Complete population realization from short native records

The configured H66 source determines its finite stationary record law without
requiring a unique simulator architecture or a continuum trajectory. Its
complete one-record observable space has dimension 42. Two consecutive
records span all 66 population functions. Reflected moments of at most five
attempts therefore recover every one of the 721 selective population maps,
including all 421 zero maps. The recovered maps predict complete raw words
of arbitrary finite length under the same fixed scheduler. The executable
control predicts a complete 128-outcome seven-attempt sign law, with original
probabilities summing to one, and separate six-, seven- and eight-attempt
raw words.

This is a population realization of the actual weighted source in
[PROOF.md](PROOF.md), with its supplied stationary preparation and native
attempt clock. It supplies finite transfer energies. It does not identify
those energies with a coherent working Hamiltonian, the modular generator,
a cubic CAR dispersion or particle masses. The exact clock separation in
[TRANSFER.md](TRANSFER.md) applies. Coherent native effects, uncontrolled
hidden dimensions and arbitrary adaptive scheduling are outside this
reconstruction theorem.

## The source and a minimal observation depth

The native interpreter supplies 66 binary total-two configurations, port
weights `(25,25,8,...,8)`, all thirty seams, equal C/E law weights, and the
complete 721-label menu. Let `T_h` be its forward selective population map,
`P=sum_h T_h`, `pi` the stationary probability column, and `Pi=diag(pi)`.
All these objects are read from actual endpoint maps; the source is not
replaced by a fitted transition matrix. The verified identities are

    1ᵀP=1ᵀ,       P pi=pi,       P Pi=Pi Pᵀ.

The stationary reference has full support. For seam `e`, the score on an
original raw flag is

    Y_e(h)=60 1[law=E,seam=e] incoming_contrast(h).

It vanishes on every other flag, including the reserved abort slots. Its
conditional mean is the configuration function `d_e(x)=x_a-x_b`.
If `B_e=sum_h Y_e(h) T_h`, the conditional mean of the chronological
product of two actual records is

    f_(e,f)=B_eᵀ d_f.                                      (1)

The first measurement's disturbance is included in `B_e`. Replacing (1)
by the passive product `d_e d_f` changes this experiment.

From the constant, the thirty single scores, and all 900 chronological
pairs, exact elimination chooses one constant, eleven single scores and
54 pairs. Put these 66 conditional functions in the columns of `X`, in
that order. The code verifies both rational inverse identities and a
nonzero minor modulo 65521. Thus `X` is an exact invertible matrix. Each
feature has a complete two-slot interpretation: a single score ignores the
second slot and the constant ignores both. Feature magnitudes on records
are bounded by 3600.

Minimality is checked against the full raw menu, rather than just the
contrast scores. Stack the constant and all 721 functions `T_hᵀ1`.
Its exact rational rank is 42, also certified modulo 65521. Hence one
record cannot span all population functions, whereas two suffice. This
statement concerns the fixed native attempt, stationary source and its
unmodified raw grammar.

## Recovery of every selective map

Set `L_h=T_hᵀ`. The source reversal in TRANSFER.md reverses both record order
and each raw label. For feature cylinders `F_i,F_j`, marked detailed balance
gives the observed Gram and insertion moments

    G_ij=E[conj(F_i(reflected past)) F_j(future)],
    (M_h)_ij=E[conj(F_i(reflected past)) 1[middle=h] F_j(future)],

    G=Xᵀ Pi X,             M_h=Xᵀ Pi L_h X.                (2)

The first experiment uses four attempts and the second five. All raw
outcomes retain their original probability. Signed scores are expectation
statistics, not additional physical projective effects.

The inverse and the recovered coordinate insertion are

    G⁻¹=X⁻¹ Pi⁻¹ X⁻ᵀ,
    A_h=G⁻¹ M_h=X⁻¹ L_h X.                                (3)

Since the first column of `X` is the constant, its coordinate vector is
`e0`. For every chronological raw word, including an impossible word,

    Pr(h1,...,hm)=e0ᵀ G A_h1 ... A_hm e0.                  (4)

This follows by cancellation of adjacent `X X⁻¹`. The complete sum of the
maps recovers `Pᵀ`, so sums over any cylinder event or any unobserved slot
also agree. The fixed stationary Markov population source is the finite
realization premise. Five-event data alone do not identify unrestricted
hidden dimensions, non-Markov processes or quantum coherences.

`Lean/SourceLaw/RecordRealization.lean` proves the finite matrix identities,
uniqueness of each insertion from its Gram moment, and arbitrary finite
products. The exact H66 matrices and all 721 raw-map substitutions are
checked by `short_records.py`; their entries are not evaluated by the Lean
kernel. The execution reconstructs each `Pi L_h` from (2), including all
zero maps, before checking longer held-out laws. The seven-attempt coarse
law keeps the bit `1[current>0]` on every attempt. Coarsening is applied to
the complete maps, without conditioning on movement.

## Certified conditioning and transfer spectrum

Exact rational arithmetic gives

    ||X⁻¹||_1=7955/396,       ||X⁻¹||_infinity=464/11,
    G >= gamma I,            gamma=55756800/2129251819009.

For the last inequality the code computes the exact inverse in (3) and
uses `||G⁻¹||_2 <= ||G⁻¹||_infinity`; symmetry equates its one and infinity
norms. Thus `gamma` is a rigorous lower bound, approximately 0.0000261861.
It is not a floating singular-value admission.

The one-unread-attempt Gram is `G1=Xᵀ Pi Pᵀ X=sum_h M_h`.
The generalized eigenvalues of `(G1,G)` are exactly the native transfer
factors. The source bound is

    (296/495) G <= G1 <= G.

If separately estimated symmetric Grams have entry errors at most epsilon,
put `delta=66 epsilon`. If `delta<gamma`, their generalized Rayleigh
quotients differ from the exact ones by at most

    eta=2 delta/(gamma-delta).                              (5)

Indeed the denominator is at least `(gamma-delta)||v||²`, the numerator
perturbation is at most `delta||v||²`, and the exact quotient lies in
`[296/495,1]`. The min-max principle gives (5) for every ordered eigenvalue.
If `eta<296/495`, the logarithmic energies have absolute error at most

    eta/(296/495-eta).                                      (6)

These are native-attempt transfer logarithms. If the estimate gives a
factor above one, (6) remains an error bound; nonnegativity of every
estimated energy is not asserted. The exact transfer also admits the
rational degree-m logarithm approximation and tail in TRANSFER.md. No
real-time native implementation follows from taking this logarithm.

## Complete held-out prediction with finite observation errors

There is a direct finite inference procedure in the specified population
coordinates. The calibrated frame `X` and configured `pi` are supplied by
the source calculation. Let `J_h=Pi L_h`. Nonnegative candidate entries
for all `J_h`, with row sums of `sum_h J_h` equal to `pi`, give a complete
population instrument `T_h=J_hᵀ Pi⁻¹`. Interval constraints on
`M_h=XᵀJ_h X` are linear and rational. A feasible candidate can therefore
be found by a finite rational linear program. This is classical parameter
inference; it neither changes physical records nor renormalizes a selected
outcome. The exact source is feasible whenever the stated observation
intervals contain its moments.

Suppose both the source and a complete nonnegative candidate fit each
insertion entry within epsilon of the same observations. For
`A=X⁻¹`, each pair of moment entries differs by at most `2 epsilon`.
The source-coordinate reconstruction gives

    T_h=Aᵀ M_hᵀ A Pi⁻¹,
    ||T_h-T'_h||_1 <= 132 epsilon ||A||_infinity ||A Pi⁻¹||_1.

Summing all 721 outcomes and taking half the norm proves a uniform
one-attempt total-variation bound `C epsilon`, where the exact certificate
finds

    C=721*66 ||A||_infinity ||A Pi⁻¹||_1
     =1248314867975/264.

The complete Markov instruments contract population one-norm. Telescoping
therefore bounds the total variation of the entire original m-record law
by

    min(1, preparation_error + m C epsilon).                (7)

This includes all words simultaneously, including zero and failure labels.
It is stronger than a separate estimate for one selected cylinder. The
candidate's nonnegativity and completeness are necessary: (7) is not
asserted for an arbitrary signed inverse obtained from noisy moments.
The construction predicts the same fixed schedule; it does not infer a
quantum instrument or additional controller policies from population data.

For independent stationary repetitions, each scalar Gram/insertion score
has absolute magnitude at most `3600²`. Hoeffding and a union bound over
`66²(721+2)` quantities give the sufficient trial count per four-/five-event
menu

    N >= 2*3600⁴ log(2*66²(721+2)/failure_probability) / epsilon².

The producer replaces the logarithm by its integer binary-log upper bound.
This is a conservative planning bound. At epsilon `2^-40` and failure
probability `1/100`, it exceeds `10^40` trials per menu. No such sampling
experiment is executed or described as practical. The exact source replay
has no sampling error; finite mixing/preparation and empirical error are
separate. The result is an executable finite reconstruction with explicit
conditioning, not a claim of efficient tomography or unique physical
continuum selection.

## Calibration on a continuing native source and its resource boundary

A stationary independent reset is unnecessary for the finite inference
bound. The actual mixing theorem in `local_repair_family/mixing.py` is
uniform over the input source and an untouched reference, after restricting
away the preparation archive. Its exact twenty-attempt population
minorization exceeds `4/25`; the coherence contraction factor is `793/990`.
Thus n unread attempts give the reduced-source allowance

    delta_n <= (65/2)(793/990)^n + (21/25)^floor(n/20).       (8)

Before each four- or five-record probe, run 8192 attempts of the same
source. Keep their raw outputs in an archive which subsequent controls
and the calibration scores do not access. No erasure or source/archive
product reset is claimed. Conditioned on any earlier accessible probe
records, the input to this fixed burn-in is a valid source state. Uniformity
of (8) therefore bounds each conditional probe-score mean relative to its
stationary value by

    b <= 2*3600² delta_8192.

The factor two is the duality between half trace distance and bounded
signed scores. The code checks (8) with exact rationals and rounds its
reported upper bound upward on a `2^-128` grid. At entry tolerance
`epsilon=2^-40`, it verifies `b<=epsilon/2`.

For one scalar probe score Z_k, subtract its conditional mean given all
previous accessible probes. The resulting martingale difference has
conditional range length at most `2*3600²`. Conditional Hoeffding's lemma
and iterated conditional expectation give

    Pr[|N^-1 sum_k (Z_k-E[Z_k|past])| >= epsilon/2]
      <= 2 exp(-N epsilon²/(8*3600⁴)).                     (9)

The mean of the conditional drifts differs from the target stationary
moment by at most epsilon/2 on every past. Apply the same union bound as
above with statistical tolerance epsilon/2. This supplies the moment
intervals used in (7) without an IID assumption or a new prepared state.
It is a finite continued-source theorem. Preparation and probe records
retain their original weights; the comparison does not include the
unavailable preparation archive in a product-reset target.

The resource limit matters. The native interpreter admits at most `2^32`
raw attempts in one program. A worst-case trial uses 8197 attempts, so at
most `floor(2^32/8197)` complete trials fit in one admitted program. The
producer actually admits that maximal whole-trial prefix and reports its
resource count. The conservative two-menu trial requirement from (9)
exceeds this cap. A longer execution service or a composition of enough
admitted programs is not established here. The certificate therefore sets
`required_sampler_fits_one_admitted_program` to false and does not describe
the statistical sampler as executed.

The explicit interpreter also has a 64-attempt and 65536-history-block
evaluation cap. On the actual stationary input, the 8197-attempt request
executes one attempt, retains all 721 raw blocks with total probability one,
and then emits `RESOURCE_STOP` on every block because a second full
expansion exceeds its block limit. Zero calibration trials complete. This
stop is executed and verified. The rational short-record contractions,
all-map inverse calculation and held-out probabilities use finite matrix
identities, not a claim that this expanded sampling service completed.
