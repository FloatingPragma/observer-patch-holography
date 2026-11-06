# Thermodynamics of the actual binary native repair

## Exact source attachment and one common reference

Fix the existing twelve-port carrier and its thirty seams. Work on X, the
66 binary configurations with exactly two occupied ports. Total load two and
independent fair placement coins are supplied source data; collision states of
the 78-state carrier are transient under a covering scan. A seam (a,b) leaves the other ten bits and x_a+x_b unchanged.
If this sum is zero or two there is one possible endpoint configuration.
If it is one there are two, and the actual native rule chooses each with
probability one half. Thus the source kernel is

    P_s = (I + S_s)/2,

where S_s exchanges the two endpoint occupancies and is an involution.
Equivalently, P_s uniformly resamples the fibre of fixed outside bits and
seam total. There are 46 singleton fibres and ten two-point fibres for each
seam. The executable compares the actual source kernel, the independent
uniform-fibre formula and the hold/swap formula on every one of the 66^2
entries for all thirty seams. This is a source attachment, not replacement
of a desired dynamics by a heat-bath ansatz.

Each P_s is symmetric, bistochastic and idempotent. Consequently
u(x)=1/66 is stationary and satisfies detailed balance for every seam.
The token graph on X joins states exchanged by some source seam; exact
breadth-first enumeration verifies its connectedness. Hence functions fixed
by every S_s are constant. For a random scan P=sum_s w_s P_s with every
w_s>0 and sum_s w_s=1,

    <f,(I-P)f> = (1/4) sum_s w_s ||f-S_s f||_2^2.

Its stationary distributions are therefore constant, so u is the unique
normalized stationary law. This uniqueness is restricted to the declared
connected binary sector; neither the total occupancy nor a law between
sectors is selected. Zero scan weights are permitted by the numerical API
but need not give a connected scan and do not inherit this uniqueness claim.

A fixed exogenous word of seams preserves u without requiring a random
seam scheduler. A word visiting every seam also has u as its only stationary
law: equality in the norm contractions of its orthogonal projections forces
every intermediate projection to fix the vector, hence all seam projections
fix it. Its positive diagonal and positive probability of one allowed token
move per word make the covering-word chain irreducible and aperiodic. This
proves finite-state convergence, not a uniform bound under carrier refinement.
The existing `native_equilibration/sweeps.py` supplies separate quantitative
finite-source preparation bounds. A product or random mixture of different
P_s is generally not idempotent; the executable retains an exact P^2!=P
witness for the uniform random scan.

## State and transition entropy optimization agree

For any law p on X, using the usual zero-mass convention,

    H(p) = log(66) - D(p || u).

Nonnegativity of relative entropy, with equality only at equal distributions,
gives the unique maximum-entropy state u. For a probability row a supported
on a seam fibre F, the identical argument gives

    H(a) = log|F| - D(a || uniform_F).

Thus the actual native transition row is the unique entropy-maximizing row
on its allowed fibre. These two optimizations use the same uniform reference
that the source kernels preserve. They do not equate a state with a transition
kernel, infer that the aggregated random-scan row maximizes entropy on its
larger support, or assume an idempotent global transition. No numerical target
is used to choose a reference law.

## Entropy production under exogenous schedules

For one seam, q=pP_s is constant on each fibre and has the same fibre mass as
p. Therefore sum_x p_x log q_x=sum_x q_x log q_x, with zero-mass fibres omitted.
It follows exactly that

    H(q)-H(p) = D(p || q) >= 0.

The support inclusion p<<q is automatic because P_s has a positive diagonal.
For strictly positive p this identity is also the uniform-reference instance of the existing
`heatBath_kl_pythagorean` theorem in
`Lean/Thermodynamics/FluctuationTheorems.lean`; the zero-support extension is
the analytic fibre argument above, not a claim about that Lean theorem's
quantifiers. For any fixed exogenous seam
word, repeated application gives a telescoping sum of these nonnegative
relative entropies. Random exogenous schedules preserve the corresponding
conditional statement; arbitrary record-dependent feedback is excluded from
this unread-source assertion because it can correlate the kernel choice
with the state. Such feedback needs the complete memory/record ledger.

For one independent random scan let q_s=pP_s and q=sum_s w_s q_s. The entropy
of a mixture has the exact decomposition

    H(q)-H(p)
      = sum_s w_s D(p || q_s) + sum_s w_s D(q_s || q).

The second term is the mixing information, obtained from
sum_s w_s D(q_s||q)=H(q)-sum_s w_s H(q_s). Both sums are nonnegative.
Omitting this term, or applying the individual-seam projection identity to
the whole random scan, is incorrect in general. The executable retains a
case with nonzero second term. These checks use rational coefficients of
logarithms of primes: expanding every positive rational logarithm proves
the displayed equalities exactly. Positivity and uniqueness use the analytic
relative-entropy theorem, not a floating-point sign test.

## Fluctuation identity and its support guard

Let P be any bistochastic finite kernel and p a strictly positive initial
law, with q=pP. Column sums one make q strictly positive. For the uniform
reference, sigma(x,y)=log(p_x/q_y), and finite cancellation gives

    sum_(x,y) p_x P_xy exp(-sigma(x,y))
      = sum_(x,y) P_xy q_y = 1.

This is the source attachment of `integral_fluctuation` in the existing Lean
fluctuation module; it does not infer physical heat or temperature. Detailed
balance gives the separate pointwise reverse-weight identity if desired.
For support-deficient p, summing only its possible forward paths instead gives

    sum_y q_y sum_(x:p_x>0) P_xy <= 1.

Equality requires the corresponding reverse support to stay inside p's
support. A point mass on one moving configuration of a two-point seam fibre
gives exactly one half, so the strict-support hypothesis cannot be silently
discarded. No totalized logarithm of zero is used in this calculation.

## Complete tape accounting and the startup boundary

For n actual attempts with a fixed exogenous schedule, the map from initial
configuration and n coin bits to terminal configuration plus the complete
coin/current tape is injective: replay the recorded currents backwards.
Consequently

    H(X_n, tape) = H(X_0) + n log(2).

The same argument holds for a fixed-length public-prefix adaptive policy;
it does not imply that the public tape ALONE identifies the input. A retained
independently randomized seam word contributes its own Shannon entropy as
well. The unread marginal may lose distinguishability by discarding record
access, but the complete classical ledger has not destroyed that information.
Physical erasure and its energy cost are not inferred. This package checks
the exact entropy identity on the actual two-attempt E/F source tree.

The binary restriction is essential. The full total-two carrier has 78
configurations, including twelve double occupancies. Under a fixed seam the
states with endpoint loads (2,0), (1,1), (0,2) all become (1,1). The uniform
78-state law is not stationary: one output obtains mass 3/78 and two outputs
lose their mass. Its entropy change is exactly -log(3)/26. A covering scan
eliminates collision support, after which the binary heat-bath result applies;
the startup itself does not obey a blanket source-entropy-increase claim.

## Quantum and physical boundaries

The existing `model/native_equilibration/record_thermodynamics.py` supplies
the separate coherent C/E detailed-balance and full quantum record-entropy
argument under its specified cut, completed Palm reference, input-independent
seam choices and quantum semantics. Its source pins are rechecked here; its
quantum data-processing proof is not duplicated or relabelled as a scalar
consensus derivation. Its conditional quantum and reference scopes remain.

For the common uniform reference in this package, -log u=log(66) is a
constant modular energy. The proof therefore selects no nonconstant physical
Hamiltonian, energy unit, temperature, heat transport, physical clock, cooling
protocol or zero-temperature law. It proves an exact source-level finite
information-thermodynamic attachment, with sector, probability measure and
schedule premises explicit. Carrier refinement and physical identification
remain separate obligations.
