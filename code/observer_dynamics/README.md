# Observer dynamics: executable evidence interface

These tools verify bounded observer-like self-reading patches: local port
loads, repair attempts, boundary readback, accumulated records and public
receipt bundles. The immutable producers, specifications and compact numerical
inputs live in `../../evidence/observer_dynamics_20260925/`.

From the RER root:

```sh
python3 code/observer_dynamics/verify.py
python3 code/observer_dynamics/verify.py --replay-refinement --write-receipt
python3 evidence/exact_federation_tower_20260924/verify_tower.py
python3 evidence/source_net_causal_poset/verify_sampled_manifold.py
python3 evidence/observer_dynamics_20260925/lean/verify.py
```

The first command checks immutable file hashes, recomputes finite-chain and
finite-window receipts, checks the 48 terminal-state digest crosspins, verifies
the complete scale/angular population, estimator identities, null normalization
and joint schedule-moment consistency, and runs independent and mutation tests. The second
also reconstructs the fresh-noise streams for all 61 retained refinement chains
and writes a separate canonical verification receipt. Neither command modifies
historical evidence. The tower verifier streams load reconstruction in 8 MiB
chunks. The manifold command uses the original independently written verifier
and tests, recovered with hashes matching the sampled-run receipt.

`archive_adapter.py` answers only the historical producers' explicit `git show`
requests for the frozen simulator revision from the archived source blobs. No
sibling simulator checkout, network service or simulator Git history is needed.
It also serializes relative provenance paths with `/` on every operating system;
the frozen producer bytes and numerical comparisons remain unchanged. Install
the repository-root `requirements.txt` for the complete default controls,
including the pinned CAMB interpolation checks. Skipped controls fail verification.
The causal controls use a separate process because they vendor another source
revision. The tests can also be run through this entrypoint without introducing
those modules into an unrelated Python process.

The large terminal arrays and geometry caches remain in their original archive.
Digest crosspins and retained statistical identities do not constitute a fresh
microscopic replay of those arrays. Large transport spectra and boosted-clock
measurements are preserved with their original producers and checks; the
bounded verification reruns their small controls. The sampled manifold record
lacks the original large event samples and joint estimator covariance, so its
large-run uncertainty cannot be independently rebuilt here.

The conditional CMB comparison uses the retained standard-transfer calculation
and supplied source amplitude, tilt, cosmological background and temperature.
It is a reproduction of a conditional observer distribution and measurement
comparison, not a source-dynamics derivation of these inputs. Plot images are
not required by this numerical package. Exact finite-chain identities,
numerical receipts, retained large-run observations and compiled Lean witnesses
are distinct evidence types.

## Retained scale and angular consistency

The scale receipt contains 51 scale rows, each with 16 schedules, and 240
angular rows. The checker binds their census to the archived specification;
binds carrier count, port count, conserved mass and geometry digest to each
original settlement receipt; and checks their derived statistics independently
of the archived scale producer. It does not import that producer to calculate
the expected answers. The new controls run through `verify.py` in CI.

Two corruptions that pass the ratio checks alone show why the extra checks are needed.
Doubling the L8 `ancestor_1` equilibrium denominator while halving all 16
normalized variances changed their mean from **1.32427 to 0.66214**. Every
ratio identity passes, reversing enhancement versus suppression. Setting
schedule noise to zero and common power to total power also passes, although
the retained unequal responses alone require noise at least
**7.02851e-6**. The unchanged archived data pass the independent checks; these
are defects in the scientific cross-check, not discovered errors in the frozen
observations. The outer archive manifest rejects changed frozen bytes.

For `P` ports, `R` raised ports and disjoint groups with `c_i` cells, let
`a_i` be their normalized estimator weights. Exchangeable fixed-total binary
loads and iid initial loads give, respectively,

```text
C = sum_i (a_i-a_i^2)/(12 c_i),
V_exchangeable = R(P-R)/(P(P-1)) C,    V_iid = (35/12) C.
```

For `G` equally weighted groups, `C=(G-1)/(12 G h)` with harmonic cell mean
`h`. For cell-weighted groups, replace `h` with arithmetic mean `m`. These
identities hold even when the groups cover only part of the ports. Ancestral
group sizes follow from the specification and parent population. Cube checks
use the retained cell-count summaries and compare the weighted/unweighted
pair. They also enforce necessary integer-count bounds: the arithmetic mean
must admit an integer total with the declared minimum and maximum attained;
convexity of `1/c` then bounds the reciprocal sum and hence the harmonic mean.
Angular pixel extrema must likewise fit the total carrier population. These
checks do not reconstruct the absent geometry arrays or prove that every
interior summary corresponds to an integer histogram.

The schedule check uses a further connection: response, power and noise must
belong to **one vector population**. For the shared initial field `x`, set
`a_i=<x,y_i>/||x||^2`, `p_i=||y_i||^2/||x||^2` and let `n` be the unbiased
schedule noise in those same units. Writing
`y_i/||x||=a_i e_0+v_i` gives

```text
r_i = ||v_i||^2 = p_i-a_i^2,
R_sum = ||sum_i v_i||^2
      = S sum_i p_i - S(S-1)n - (sum_i a_i)^2,
max(0, 2 max_i sqrt(r_i) - sum_i sqrt(r_i))^2
    <= R_sum <= (sum_i sqrt(r_i))^2.
```

The polygon bounds are necessary and sufficient for some real vector
population when the residual space has dimension at least two: vectors of
specified lengths in a plane attain every resultant length between those
bounds. Centered scalar groups have residual dimension `G-2`, and real
degree-`ell` angular coefficients have residual dimension `2 ell`; all retained
rows meet this condition. This does not prove realizability by bounded integer
loads, reproduce an unseen covariance matrix, or authenticate the raw arrays.
It does reject impossible moments that separate row identities cannot detect.

The executable check is a numerical consistency diagnostic. It assembles the
supplied binary scalars as exact fractions and uses relative tolerance `1e-9`
on the relevant moment scale, without an absolute small-power floor. It
propagates residual subtraction uncertainty before taking square roots, with
outward rational root bounds. A tiny negative computed noise within roundoff
is allowed; reported statistics are never rewritten. This preserves valid
nearly collinear and identical schedules. For the response spread, the exact
sample variance `V` and variance about a rounded mean `m_hat` differ by
`S/(S-1) (m_hat-mu)^2`. The checker bounds the mean's binary64 summation and
division error, then permits that quadratic centering increment. It retains a
relative lower bound on `V` and binds the standard error to the reported
standard deviation; an erased or inflated resolvable spread still fails.
For example, 16 original vectors with responses `[1+ulp(1),1,...,1]` have a
NumPy variance about 6.7% above the exactly centered variance, although their
absolute discrepancy is only rounding of the common center. Controls exercise
both consumers at small, ordinary and large scales. Negative **unbiased common
power** and noise fractions greater than one are also possible and remain accepted.

The repair changes verifier acceptance and its live verification summary.
Frozen inputs, producer bytes, numerical observations, claim registrations,
papers and physical classifications are unchanged. It supports the standing
[evidence audit #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
without supplying a microscopic law, physical calibration or CMB prediction.
