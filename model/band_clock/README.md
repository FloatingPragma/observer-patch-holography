# Moving the complete bound band with a record-derived rate

This component executes the **same periodic Weyl/contact action** on all 8,128
two-fermion states. The native count decoder supplies the rate
`beta_hat = log(theta_hat_1/theta_hat_7)`, and the declared transfer service sets
the Hamiltonian model duration to `t_hat = -beta_hat*s`, with `U=1`.
The hopping and contact words come from the existing source controls. No
band-only Hamiltonian, fitted spectral isometry, removed scattering state or
selected Weyl taste enters execution. A coherent 64-momentum pair packet is
entangled with an untouched reference qubit; every pointer output is retained.

This extends source execution from the earlier stationary eta-pair/dark M2
corner to a moving packet of the same isolated band. It supplies an operational
rate-transfer protocol **under an additional common-reference premise**. It
does not extend the native product reference by an intertwining isometry.

## Compatibility before the update

Fix the supplied periodic lattice, `c=1/8`, `U=1`, CAR representation, and
`H = dGamma(h) + U D`, where `h(k) = -c sigma.sin(k)`. Define the mathematical
reference on the complete N=2 algebra by

```
tau = exp(-beta H)/Z,       beta*U = log(25/8).
sigma_s(A) = tau^(is) A tau^(-is) = tau_H,-beta*s(A).
```

The Gibbs density is faithful on **all**8,128 states. Its spectrum is computed
through the 64 physical antisymmetric momentum fibers, of dimensions 120 or 128.
There is one global partition function. Giving each fiber its own normalized
state changes cross-momentum modular phases: the retained control has an
operator discrepancy around 0.071 on a single cross-fiber matrix unit.
The global Gibbs identity holds on the entire algebra, including these units.
The observed state is the entangled moving packet, never this Gibbs reference.

The reference is defined from the supplied action. Neither a local repair law
that selects it nor a source preparation of it is supplied here. Copying the
native reference cannot remove this premise: its modular gaps are integer
multiples of `log(25/8)`, whereas the moving band has noninteger energy gaps in
the eta-calibrated units. The [acceptance audit](../completion/acceptance_audit.py)
gives an exact spectral obstruction for every finite independent copy count. The present protocol uses the native records as a rate estimate
for external source control; it does not identify native repair with unitary
field motion.

## Actual source execution and uniform comparison

Each even periodic spatial matching is a disjoint matching of fermion modes.
Wrap bonds are reordered with conjugated coefficients. Independent sparse
position-space assembly checks the entire action before execution. The six
hopping groups have two-particle operator norms at most `c` each; the contact
projector has norm 1. A symmetric split with 32 steps uses the completed local
source hopping/contact words and counted FSWAP routes. Local words are
contracted on **every** vacuum, single-occupancy and double-occupancy input.
Fermionic route identities lift them to the full non-Gaussian N=2 state.

For `|s|<=1`, the admitted count-ratio window `[3,10/3]` gives
`|t_hat|<=5/4`. The simultaneous count confidence event gives
`|beta_hat-beta|<=3/80`; its failure allowance is
`381169/10^12`. Centering the full Hamiltonian gives

```
||H - I/2|| <= 1/2 + sqrt(3)/4 < 15/16,
epsilon_clock <= 9|s|/256,
epsilon_split <= (7/4)^3 (5|s|/4)^3 / (3 n^2).
```

The sum of these bounds and the statistical failure allowance bounds the
complete average half-diamond distance. This is uniform over all N=2 inputs
and arbitrary finite untouched references. Any common complete detector or
subsequent channel contracts the bound. Rejected training records cause an
identical stop in both comparisons, preserving the input and original tape.
No division by acceptance or pointer probability is made. The bound concerns
ideal completed source words; floating local replay errors and physical
hardware errors are distinct. Ideal proper-code export failures are zero by
the existing invariant-code identities, rather than discarded successes.

The retained numerical execution includes a two-bit source occupation pointer
at the origin, all four quantum outputs, full contact-site/complement
probabilities, sparse-action and Fourier-band comparisons, deleted-contact
and zero-clock controls. Costs include every evolution route and pointer,
with the 31,457,280 training attempts counted separately. Preparation of the
band packet, field-register initialization, coherent load/export services,
analogue control precision and the new Gibbs reference remain supplied.

This is a finite fixed-spacing interacting field calculation. It establishes
no physical seconds, electronvolt mass, nontrivial defect charge, Lorentz or
chiral continuum field theory. CP04 and CP06 still require the corresponding
native-reference, physical calibration and continuum joins.

Run:

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest model/band_clock -q
OPENBLAS_NUM_THREADS=1 python -m model.band_clock --check model/band_clock/results.json
```
