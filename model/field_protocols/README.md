# Complete adaptive field observations

This package executes a two-read, record-controlled protocol on the complete
one/two-site Weyl/contact Fock field. A first occupation record chooses the
sign of the next Hamiltonian evolution; a second occupation is then read.
The source implementation includes every six-dimensional primitive and its
sender-reset export. It retains all four joint records and every terminating
first-export failure, labeled by its actual record prefix. Separate
complement Kraus maps are added incoherently. No successful outcome is
renormalized. Initialized pointers, two durable record bits, the feedback
controller and access to signed evolution are declared services.

For each input basis column, the code calculates the complete source output.
Consequently the same maps act on any input entangled with an untouched
reference. If `A_r` and `B_r` are actual and ideal successful history maps,
and `F_j` the terminating abort maps, then

```
V_actual* V_actual = sum_r A_r* A_r + sum_j F_j* F_j = I,
||V_actual-V_ideal||² = ||sum_r (A_r-B_r)*(A_r-B_r) + sum_j F_j*F_j||.
```

The exact source isometries make the first identity telescope, even with
adaptive control. Contractivity after reading the orthogonal records gives
a complete-history half-diamond bound from the second identity. Numerical
evaluation checks trace preservation to its stated tolerance; its reported
dilation norm is not an outward floating-point or hardware certificate.
Ordinary execution also checks that this norm fits the uncapped ideal
word error plus its separate numerical allowance. A wrong but unitary
source therefore cannot pass merely by conserving probability. Deliberate
out-of-model leakage runs require `diagnostic_only=True`; their result
explicitly states whether the ideal source error budget was satisfied.

Independently, the existing rational overlapping-support product-formula
bound controls the ideal source words. If the first segment has error `e0`
and the two record-selected continuations have errors `e10,e11`, the complete
adaptive experiment has half-diamond error at most
`min(1,e0+max(e10,e11))`. This follows by telescoping channels on the direct
sum of retained record spaces: the norm of record-controlled continuation
is the maximum over its branches. No inverse outcome probability or number
of histories enters. This remains true for arbitrary untouched references.

Tests compare the complete instruments with independent matrix exponentials,
reject the wrong feedback sign, and inject genuine unitary source leakage.
The leakage control checks the reference marginal and demonstrates the
probability lost by deleting aborts. Resource counts give the largest cost
of one actual history, with two pointer initializations; evaluating every
counterfactual branch on the host is a separate computational cost.

This supplies a finite complete-observable-law attachment on the same action
as the particle and faithful-reference packages. The input state is supplied;
a preparation channel and its error can precede it by contractivity. It is
not a derivation of a chiral 3+1 continuum, a physical clock, or an empirical
particle dictionary. The modular and half-filled reference states keep
their distinct contracts.

Run `python -m model.field_protocols.protocol` and
`python -m pytest -q model/field_protocols` from the repository root.
