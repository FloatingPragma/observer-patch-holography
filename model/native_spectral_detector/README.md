# Native finite-clock readout of the complete ground-family spectrum

This package implements a complete four-outcome energy clock on the original
grounded cube. Its input is the invariant five-dimensional ground mixture,
prepared by finite native source words from blank registers. A native
charge instrument retains both addition and removal. The clock retains all
outcomes and executes coherent original-action histories through model time
three. Its probabilities include every one of the previously resolved
78 positive excitation lines. Four broad bins do not separately resolve
78 lines, and their centers are not particle masses.

The action remains the original zero-chemical-potential Weyl/contact cube
with `c=1/8`, `U=1` and exact ground energy `E_G=-sqrt(3)/4`.
It is not the separately considered centered action `H-N/2`. The existing
native six-level source controls, blank qubits, coherent proper-code
load/export, classical feedforward and computational readout are supplied
apparatus. Their literal pulse errors and all export failures are retained.
No claim about physical hardware accuracy is inferred from host arithmetic.

## A finite, exactly invariant ground preparation

The frozen source proof gives the unique spin-two ground multiplet with
highest-weight vector `G_0` and normalized ladder vectors `g_k`, `0<=k<=4`.
The existing native signed-tree word prepares `G_0` on hidden up modes;
the full local Fock inverse frame is applied after the following rotations.

Take five azimuths `phi_j=2*pi*j/5` and three polar cosines

```
x_l = -sqrt(3/5), 0, sqrt(3/5),
w_l = 5/18, 4/9, 5/18.
```

A four-bit native selector has squared amplitudes `[5,8,5]` repeated five
times, followed by zero, divided by 90. It is measured by the complete
product of the existing six-outcome qubit readout. All `6^4` slots are
retained; code record 15 and non-code outcomes stop the history. Their
ideal probabilities are zero, and any actual probability is included.
Conditional on record `(j,l)`, apply the same
`R_z(phi_j) R_y(acos(x_l))` to each of the eight hidden spin pairs.
Every rotation is compiled from two-body Pauli source words and fixes the
local vacuum and doublon, so this is the full Fock spin rotation.

The coefficient of `g_k` in the rotated vector is

```
sqrt(C(4,k)) ((1+x_l)/2)^((4-k)/2)
            ((1-x_l)/2)^(k/2) exp[-i phi_j (2-k)].
```

The five-azimuth sum kills every off-diagonal magnetic difference, since
`0<|k-k'|<=4`. The three-node polar quadrature integrates every polynomial
of degree at most five against `dx/2`: its even moments are `1,1/3,1/5`
and its odd moments vanish. Each remaining diagonal coefficient integrates
the degree-four binomial polynomial to `1/5`. The ideal finite mixture is
therefore exactly `P_ground/5`, with no group-average limit or projection.

The selector angles are only proposed by host trigonometry. Exact integer
weights, inverse Walsh coefficients and rational sine/cosine enclosures
independently certify the emitted tape. Polar rotation angles are similarly
checked against nested rational enclosures of their algebraic target
sines and cosines. Exact opposite Pauli coefficients retain the vacuum and
doublon phases. The already certified ground tree and physical frame are
reused without modification. The complete preparation error is the sum of
the tree/frame bound, selector bound and maximum conditional rotation
bound; it is below `10^-9`. Orthogonal selector records require the maximum
conditional error, not a sum over alternative histories.

## Both charge branches and a complete clock

At the origin the physical and hidden frames agree. Let
`gamma_1=Z_0 X_1=c_1+c_1^dagger`. Copy the input occupation of mode 1 into
a blank tag by a native CNOT, apply `gamma_1`, and read the tag. For every
full-Fock input, the two valid Kraus maps are exactly

```
tag 0: c_1^dagger,
tag 1: c_1.
```

Their effects sum to one by CAR. The remaining four native tag outcomes
are retained and have zero ideal weight. Measuring the input occupation
is part of this instrument; neither charge branch is normalized or
discarded. A missing Jordan-Wigner `Z_0` changes the instrument and is
tested as a negative control.

Initialize two clock qubits in the uniform superposition. On clock bit
`j=0,1`, apply controlled original-action evolution for time `2^j` and
the controlled scalar phase `exp(+i 2^j E_G)`. The scalar comes from the
independently proved finite ground energy; it sets the excitation-energy
origin and is not an empirical mass calibration. The final clock state
before its Fourier transform is

```
(1/2) sum_(n=0)^3 |n> exp[-in(H-E_G)] psi_tag.
```

Apply the positive-sign four-dimensional Fourier transform and measure
both clock qubits with the complete native readout. The effect at bin `k`
is exactly the previously proved positive kernel

```
K_k(E) = |(1/4) sum_(n=0)^3 exp[-in(E-2*pi*k/4)]|^2.
```

Finite Fourier orthogonality proves `sum_k K_k(E)=1` for every real energy.
This is one complete positive energy measurement, rather than a Fourier
reconstruction of separately noisy Ramsey expectations. The physical
excitation costs in both charge branches are nonnegative. Bin centers are
modulo `2*pi`; a negative principal center represents a clock outcome,
not a negative excitation cost. The signed retarded-response frequency
can be recovered by retaining the removal tag.
The zero-centered bin also receives positive-energy weight; observing it
does not certify a massless state.

The exact invariant source gives unnormalized addition measure
`mu_down/2+delta_lambda/4` and removal measure `delta_lambda/4`, where
`mu_down` has unit mass and `lambda=sqrt(3)/16`. The local spatial-projector
factor is `P_+(x,x)=1/2`. Thus the charge masses are `3/4` and `1/4`.
Every existing positive line contributes; at lambda the addition weight
is `3/10`, and removal weight is `1/4`. Exact interval kernel sums give
all eight joint charge/bin probabilities. They are compared to the
executed source measurement without fitting any line or probability.

## Controlled original-action words and the complete error budget

The controlled hopping word uses the original signed fermionic swap
routes. At their endpoints each `exp(-it alpha P)` becomes

```
exp(-it alpha P/2) exp(+it alpha P Z_clock/2),
```

which acts as the identity on clock zero and as the original hopping on
clock one. The swaps and their inverses are unconditional, so their
identity arm cancels exactly in the target circuit. The controlled
contact uses the eight commuting Pauli terms of
`n_up n_down |1><1|_clock`. All these words act on at most three local
work bits and are executed through the actual six-level source matrices.
Their templates, fixed phases, scalar terms and literal port errors are
independently checked. No controlled abstract Hamiltonian is a free
primitive.

Use 16 Strang steps for time one and 32 for time two. The controlled
ideal product has the same full-Fock operator error as the uncontrolled
product, by its two clock blocks. The two existing refined split bounds
sum exactly to `1887/65536`. Adding the complete preparation, charge,
clock, source-control and ground-phase bounds gives a complete output
state and total-variation error below `3/100` (approximately
`0.02879333504`). This applies before conditioning on any outcome and
preserves arbitrary untouched references during field evolution.

One longest successful history contains 612,666 source pulses: at most
5,902 preparation pulses, 13 charge-instrument pulses and 606,751 clock
pulses. The selector is read while the ground still occupies its eight
compact modes; the eight blank spin modes are embedded afterward. The
charge tag is measured before the two clock bits are initialized, so the
maximum simultaneous quantum register has eighteen bits. Pulse failures
stop that history and remain output events.
The final apparatus has `6^3=216` tag/clock slots for each valid selector
record. All 3,240 joint selector/tag/clock slots are recorded, including
zeros, in addition to selector stops and every first-export failure.
Failure probability tables marginalize the classical selector labels;
their complete vector family remains reconstructible from the frozen
initial columns and the literal prefix instruments.

## Execution and independent checks

The selector, all fifteen conditional ground words, charge instrument
and complete eighteen-bit clock are actually evaluated. The clock keeps
all 262,144 quantum rows and thirty classical selector/charge columns.
An auditable C++ contraction fuses the same local matrix multiplication
and reduced-density calculation used by the existing NumPy interpreter;
independent tests compare arbitrary complex inputs, nonordered targets,
all supported block sizes and up to thirty-two columns. It performs no
rank, number-sector, spectral or amplitude truncation. Five workers
propagate disjoint classical columns; linearity and orthogonal classical
labels make this an exact decomposition of the specified instrument.
All native export effects and failure slots are still evaluated.

The comparator alone uses the exact original Hamiltonian exponential in
the source-invariant number sectors, then assembles the complete field
and clock array. It supplies a diagnostic, not the literal-control error
proof. The full native output, complete failure tables, source hashes and
all probabilities are retained for cold replay. A records archive and
five ordered six-column state archives preserve every array entry while
keeping each repository file below the single-file size limit. Host roundoff receives
the declared numerical admission; it is not promoted into hardware
accuracy. No passive repair step is inserted during this active
measurement. Its initial reference and charge laws are the original
grounded-cube laws; active readout is not an autonomous vacuum-selection
or continuous-memory theorem.

For actual sampling, 100,000 independent complete histories suffice to
estimate all eight charge/bin probabilities and the coarse abort
probability within `1/100` simultaneously, with failure probability below
`1/1000`. This follows from Hoeffding and an exact positive-series bound
on the exponential. No random shot records are fabricated here.

This closes a finite native preparation-and-energy-POVM task. It does not
identify measured species, create additional resolution beyond its four
bins, derive relativistic dispersion or convert model energies to rest
masses. Infinite-ground preparation, bounded-frequency continuum
structure and the existing Einstein prerequisites remain separate.
