# Global vacuum fidelity of local scalar read programs

The [derivation](DERIVATION.md) proves a global failure and a constructive
repair; the [contract](CONTRACT.md) fixes the exit requirements. Every energy
and overlap uses the original Hamiltonian and its original vacuum.

A stable second-order split can agree on smooth detectors while producing
many ultraviolet particles. For N sites and bulk frequency scale sigma,
uniform vacuum fidelity needs N*tau^4*sigma^4 to vanish; total excess energy
needs N*tau^4*sigma^5 to vanish. The standard symmetric fourth-order split
replaces these by N*tau^8*sigma^8 and N*tau^8*sigma^9. Matching time-average
lower bounds make these necessary as well as sufficient, with every sampled
time retained. A fixed-degree scalar polynomial symplectic update cannot
nontrivially preserve the original vacuum exactly on all refinements.

For the explicit repaired sparse stencil, sigma grows as q^(7/32), rather
than inverse read radius q^(1/8). The fourth-order schedule
tau=q^(-5/8)/64 gives total produced energy O(q^(-1/32)), global infidelity
O(q^(-1/4)), and Theta(q^(79/16)) ordered reference read incidences. It also
controls every bounded effect for the vacuum and bounded-energy coherent
inputs in a fixed frequency band. The same tick with the second-order method
has time-averaged energy Theta(q^(51/32)). The derivation gives matched
energy-tolerance cost comparisons with the dense family.

The support bound is sharp: an exact nonzero extreme signal travels (3j+1)
read radii in j fourth-order macros before wrapping. This proves that their
numerical durations cannot be an RG-compatible complete influence clock at
the same fixed speed. The quantum fidelity repair and this native-clock
obstruction are distinct proved results.

This work is standalone; the sparse construction of `code/m1_interfaces` is
fully defined again. It proves neither a native A1--A3 implementation nor a physical
clock. Inverse Hamiltonian gates are explicit operations. The classical
writer ledger is a reference computation, not a copy of unknown quantum
states; the quantum map factors into actual commuting two-mode phases and
onsite unitaries. Computational substeps are not faster physical transport.

Run from the repository root:

```sh
python code/m1_vacuum_fidelity/build.py
python code/m1_vacuum_fidelity/verify.py
python -m pytest -q code/m1_vacuum_fidelity
cd Lean
lake build Geometry.M1VacuumFidelityAxiomAudit
```

The compact receipt retains four complete finite periodic spectra (523 site
instances), both methods at both declared ticks, and all 32 observation times:
512 full-population observations, representing 66,944 mode-time evaluations.
The independent checker also executes 24 full spatial covariance comparisons;
its matrix construction and exact closed-walk moments do not use the producer's
Fourier implementation. The signed reference programs have 702 forward scalar
writes and 702 inverse writes, with every writer version included in the
trace commitment and 2,268 nonlocal read incidences. Independent exact cubic
arithmetic checks cancellation, stability and inverse execution.
Six additional 32-site controls execute the full non-wrapping signal rays
and independently verify nonzero extreme coefficients and zeros outside support.

The first q=65536 sparse level is checked by exact integer-ball column moments.
Its full q^3 population is not executed. The all-level fidelity thresholds and
cost exponents are analytic theorems; ten Lean results check finite oscillator,
composition, bulk-fraction and scale algebra. Their transitive axiom audit
allows only propext, Classical.choice and Quot.sound, with negative trust tests.

## Numerical contract

`model.quantum_case` and `check.quantum` accept a nonempty spectrum of positive
exact Python integers or `Fraction` eigenvalues, positive integer
multiplicities, a positive exact integer or `Fraction` tick, and order 2 or 4.
Booleans and inexact scalar types are rejected. The certified domain is
checked **before rounding**, as `0 < tick**2 * eigenvalue <= 1/100` for every
mode. All 32 times and the full multiplicities enter the reported particle
number, original-Hamiltonian excess energy (hbar=1), and negative log fidelity.

The producer uses the factored Chebyshev identity in the
[derivation](DERIVATION.md#evaluating-the-original-vacuum-signal-without-cancellation).
The checker constructs the original gate matrix from independently multiplied
exact cubic-field polynomials and evolves that matrix. They share the
input, interval arithmetic and output-formatting layer, not the mode formula.
Each uses a private outward-rounded interval context and recomputes from the
exact inputs at increasing precision. An observation or mean is emitted only
when both endpoints are positive and format to the same 23 significant
decimal digits. The positive series for `log(1+n)` includes an explicit
remainder bound; no absolute floor discards a small mode before weighting.

Precision starts at 90 decimal digits and doubles up to 5760. Unresolved
formatting or a zero-containing enclosure raises `ValueError` at that budget;
it is not an accepted zero or a mathematical domain exclusion. Neither path
changes the caller's ordinary or interval mpmath context. These numerical
enclosures rely on mpmath's directed arithmetic; they are not additional Lean
proofs or a certification of the physical source model.

The retained original-input controls reproduce failures on main commit
`0660c945`: both old paths erased second-order production at tick `10^-40`,
including an ordinary total of `1.5625` after multiplicity weighting. The
fourth-order producer reported a value about `2.56e23` times too large at
tick `10^-20`, while the checker reported zero. Near a step-32 revival both
were wrong by many orders of magnitude. The checker also accepted an exact
tick just above the certified upper boundary. Exact rational second-order
controls and independent elementary-gate covariance evolution cover these
cases, both orders, compensated units, mixed spectra and the retained catalog.
Canonical regeneration preserves the complete existing numerical evidence
and parent-claim hashes exactly; only source pins change. The four registered
vacuum, sparse-tick, polynomial-obstruction and signal-clock claims retain
their stated mathematical scope. No frozen acquisition, paper theorem,
physical classification or claim-registry payload changes.

The supplemental binary64 spatial covariance check uses relative
tolerance 0.002 and absolute tolerance 1e-24 for the very small produced
particle numbers; the exact algebra and high-precision receipt comparison
do not use this tolerance. Inverse reference execution is exact in the
checker and has residual below 1e-55 in the producer. No fit defines a graph,
accepted state, observed time, or refinement schedule.
