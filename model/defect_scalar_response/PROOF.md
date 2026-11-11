# Native removal response to scalar occupation cost

The calculation uses the native thirty-mode local fermion source with its
five-state scalar register and all quark and lepton transition coefficients
one. The reference is its physical ground doublet, divided by two. The
literal annihilation maps and complete invariant targets are inherited from
`native_quark_response` and `chiral_joint/charged_ground_exchange`. The
calculation introduces no observed particle values.

## Scalar intervention and ground selection

Write the full local Hamiltonian as `H(mu)=V+mu N`, where
`N=diag(0,1,1,1,1)` acts on the scalar register and `mu=1` is the native
source. The interaction `V` connects scalar vacuum to scalar excitation,
so `U=2N-I` obeys `UV+VU=0`. The operator `N` commutes with every gauge
current and every fermion annihilator. Varying `mu` changes the scalar
occupation cost while preserving the native gauge constraints. This
intervention is distinct from a neutral Higgs field amplitude or a
Yukawa-channel coefficient.

On the actual sixteen-dimensional highest-spin source block, the complete
annihilator of `V^2` is

`x (x-105) (x^3-313x^2+31708x-1041468)`.

Let `e` be the least real root of

`[e(e-1)]^3-313[e(e-1)]^2+31708e(e-1)-1041468=0`,

and let `lambda=e(e-1)`. The local reference theorem selects this ground
on the full physical source. Every full local scalar block has energies
`(mu +/- sqrt(mu^2+4 lambda_j))/2`, together with scalar-dark energies
zero or `mu`. The negative branch decreases strictly with `lambda_j`.
Therefore the same maximal squared-transfer eigenspace gives the full
physical ground for every real `mu`; its spin multiplicity is two.
This uses the complete local ground theorem, without transferring an
eight-cell ground estimate to the varied action.

Let `Pi` be the complete rank-two spectral projector of `V^2` at `lambda`
in one magnetic block. The actual ground projector is

`P(mu)=Pi/2-Pi[2V+mu U]/(2 Delta(mu))`,
`Delta(mu)=sqrt(mu^2+4 lambda)`.

The code constructs `Pi` by exact polynomial division and verifies its
idempotence and its `V^2` eigenvalue equation in the degree-six source
number field. Its rank is two, and `tr(P(mu))=1`. At `mu=1`,
`Delta=1-2e`. No numerical eigenvector selects the reference.

## Complete energy-transfer moment

For a literal removal map `A_r`, let `H_t(mu)` be its complete target
Hamiltonian and let `epsilon_rj(mu)` run over every target energy, with
the complete spectral-projector weights `w_rj(mu)`. Define

`M2_r(mu)=sum_j w_rj(mu) [epsilon_rj(mu)-e(mu)]^2`.

This is an unconditioned second energy-transfer moment. Pump failure
retains probability one half and receives numerical readout zero. Pump
success also has probability one half; no division by this weight is
performed. A static charge-compensating reference contributes its stated
energy separately and is excluded from the source-only transfer. The
complement outside the complete invariant target has probability zero.
Degenerate and zero-residue target energies are included through their
full spectral projectors. No energy line or residue is selected.

The complete spectral theorem gives the exact sum rule

`M2_r(mu)=tr[P(mu) K_r^ddagger K_r]`,
`K_r=H_t(mu) A_r-A_r H(mu)`.

Here the adjoint uses the verified positive source and target Gram
metrics. Since `N_t A_r=A_r N`, the full operator `K_r` is independent
of `mu`. The code checks this equality on the entire admitted source
block, including reconstruction in the unreduced scalar and matter
coordinates. Its squared norm `O_r=K_r^ddagger K_r` commutes with `N`.
For quarks, exact separate-spin symmetry gives the physical mixed-ground
measure from either source magnetic block. For leptons, the calculation
averages the two positive spin-removal measures in source spin up; spin
covariance identifies this with a fixed spin probe on the mixed physical
ground. It does not use the signed spin difference employed for exchange.

Define `a_r=tr(Pi O_r)` and `b_r=tr(Pi U O_r)`. Scalar parity gives
`tr(Pi V O_r)=0`. The full response and its ground-following derivative
are therefore

`M2_r(mu)=a_r/2-mu b_r/(2 sqrt(mu^2+4 lambda))`,
`dM2_r/dmu=-2 lambda b_r/(mu^2+4 lambda)^(3/2)`.

The derivative includes the change of the actual ground projector.
Holding the reference fixed would omit this response. The complete
success effect has expectation one half and derivative zero in every
channel. All five moments at `mu=1` are independently checked against
the complete retained spectral moments `m2-2e m1+e^2 m0`, with the
declared factor one half for the lepton spin sum.

## Charge classes and interpretation

The removed fields `u^c,d^c,Q,e^c,L` carry centre classes
`2,2,1,0,3`, respectively; the created holes carry their negatives
modulo six. These are representation labels shared with the protected
flux classification. The response uses actual matrix elements, ground
polarization and complete transition images, which contain more
information than the centre label. No map creating a spatial flux or
winding defect is supplied by an elementary fermion removal.

The result is a native scalar-cost response in the declared source
clock. It does not identify a rest mass, infer a physical Yukawa
coupling, implement preparation of the varied ground, or establish
the varied eight-cell source response. The equal unit transition
coefficients are architectural inputs to the calculation.
