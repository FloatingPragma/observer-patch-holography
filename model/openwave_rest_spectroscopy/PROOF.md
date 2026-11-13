# Complete rest-energy records on one interacting source

## Native action and reference

The strong source is the five-input positive spin-gradient action in
`model/openwave_spin_gradient_nuclei`. The charge-changing contact, all of
its neutral-current terms and the reference are the native operators of
`model/openwave_dressed_baryon_lines/source.py`. The source is composed before
any spectral calculation. It retains 120 strong colour-singlet states and
all 70 four-fermion states of the eight light-family lepton modes. The
reference fills the four massless neutral modes; the charged rest energy
is one. The box kinetic coefficient is `(delta-1)/2`, with the exact same
delta as the strong source. The contact coefficient is 1/64.

The complete charge-zero and charge-one carriers have dimensions 1,400
and 360. Fixing twice the spin projection to one gives 398 and 102 states.
Exact cell reflection gives even blocks of dimensions 199 and 51. The
source checks every omitted block is reducing, retains the odd block and
all other charge, spin and Fock outcomes, and declares an idle compensating
charge reference. No equality with the massive three-family reference is
used.

A rational integer vector defines the actual preparation in the even block.
Its rank-one filter on the uniform fixed-spin carrier succeeds with
probability 1/398 or 1/102. The complementary filter is retained. The
nonselective preparation leaves this uniform state unchanged, so its mean
work is exactly zero. The full source, reference and input vector define
all conditional and original-input records.

## Exact finite spectral bounds

Every entry of the even Hamiltonian is an explicit algebraic number obtained
from integer native-current matrices and rational coefficients. The finite
spectral certificate constructs a dyadic frame O and eigenvalues E. If
eta bounds `||O*O-1||` and epsilon is its reported operator error, the polar
frame U obeys `||H-U diag(E) U*|| <= epsilon`. These bounds include the
algebraic-to-dyadic operator enclosure, the exact integer residual and Gram
error. Every eigenvalue and every cluster probability is retained. A
spectral gap is used only for a stated cluster probability enclosure,
never to remove a line. The contact shifts and splits the prepared energy
without changing the strong calibration coefficients.

The charged zero-mode addition in the quark-vacuum sector costs exactly
one source unit: the cross-current vanishes there. The common electron
calibration converts every energy into GeV. An expectation of the full
Hamiltonian is not identified with a single interacting particle mass.

## Finite register instrument

Let a=1792 and W=256 source units. The complete eigenvalue enclosure proves
that the whole admitted even block lies strictly inside `(a,a+W)`. For a
register of dimension M, define

`K_k = M^-1 sum_(n=0)^(M-1) exp(2pi i k n/M) exp(-2pi i n(H-a)/W)`.

The discrete Fourier identity gives `sum_k K_k* K_k=1` on the entire
admitted block. Each operator commutes with H, hence the nonselective
measurement changes no mean source energy. The outcome state is exactly
`K_k psi/sqrt(p_k)`. The complementary flag acts as the identity outside
the admitted block. This specifies the complete instrument, including
posterior states, even when a table contains only pointer probabilities.
The idle reference is unchanged.

The spectral probability kernel is
`[sinc(M delta)/sinc(delta)]^2`, with delta reduced modulo one. All M bins
are executed for M=1024, 4096 and 16384. They are neither clipped nor
renormalized. The interval certificate and the integer vector also define
the exact dyadic-frame approximation independently of floating evaluation.
Finite Duhamel bounds compare the controlled unitaries in the two sources:
with `tmax=(44/7)(M-1)/W`, an upper bound on the total-variation difference
between the exact instrument and the ideal dyadic-frame formula is
`4 eta + 4 tmax epsilon`. The inequality uses pi<22/7 and requires no
spectral gap. Reported floating sinc evaluations are numerical diagnostics;
the stated analytic error does not certify each machine transcendental
operation. Normalization checks and an independent direct controlled-unitary
regression verify their implementation.

The bin spacing W/M is an apparatus resolution. It is neither a decay
width nor a fundamental time scale. The register uses a declared controlled
source action, without assuming an external physical global clock.

## Complete integer proof witness

The pinned `spectral_witnesses.json` specifies each complete `2^40` dyadic
frame and the complete integer prepared state in the exact even basis.
The preparation is normalized from its exact integer squared norm and
independently checked against the rebuilt strong ground state. It is a
declared preparation, not an inferred exact eigenstate of the contact
action. The colour-string fixed-witness consumer rebuilds the exact source
and recomputes all integer Gram and residual bounds by the unchanged
six-quark complete-spectrum proof. It accepts neither stored intervals nor
stored error bounds. Every spectral column and every preparation coefficient
is retained, including zero coefficients and small-weight lines.

## Scope

The outputs are complete finite-box interacting spectral and pointer
records in electron-calibrated units. The contact source is not refitted to
physical interacting poles. The calculation supplies no infinite-volume
scattering pole, physical neutron lifetime, universal proper-time
calibration or transfer to a different neutral reference.
