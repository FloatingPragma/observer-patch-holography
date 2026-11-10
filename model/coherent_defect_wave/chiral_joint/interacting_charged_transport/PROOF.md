# Charged propagation with a fixed compensator

## Source and complete preparation

The source is the actual eight-cell, thirty-mode-per-cell Peter–Weyl
Hamiltonian `H_lambda=H0+lambda T`, with `lambda=2^-34` and all onsite
Yukawa, scalar and link coefficients one. The static thirty-mode CAR
reference and its cost `4 N_ref` are those of
[interacting_charged_response](../interacting_charged_response/PROOF.md).
The reference is fixed at cell zero throughout the clock evolution.
No hopping term is added to H_lambda. The source clock is its actual
exponential, and the static reference clock is `exp(-4it N_ref)`.

Use the physical algebraic preparation from
[interacting_ground_preparation](../interacting_ground_preparation/PROOF.md).
Let q be its success probability, bounded below by
0.999999942661984023. Its conditional active source has fidelity at
least F with the actual unique interacting ground g_lambda, where
`F>=0.999999999999944496`. The second source and every reference and
ancilla are retained. Every preparation failure stays physical.

An explicit localization intervention follows. Let P0 be the tensor
product of all eight true local ground-doublet projectors and the
blank projector on every link. Measure `{P0,I-P0}`. The complete
failure is retained and subsequent controls on that branch are off.
Apply the same actual elementary removal pump B_r, for e^c or L, on
localization success. This extra intervention is part of the declared
preparation, rather than an implicit replacement of the interacting
ground by an atomic state.

The exact local occupation theorem gives
`P_ground N_r P_ground=d_r P_ground` on the full local doublet.
Consequently

`P0 B_r^dagger B_r P0=P0/2`.

Both pump success and pump failure have probability exactly one-half
conditional on localization, for every state or mixture in P0. The
successful packet is represented without renormalizing intermediate
branches by `B_r P0` acting on the inherited prepared source. All links
are blank at this point. The source has odd parity, and the static
reference has exactly one particle, so the joint parity is even.

## Arrival effects and charge accounting

There are six actual oriented source links incident to the origin. The
native compiler gives edges 0,1,2,5,7,12, leading to cells 4,2,1,1,2,4.
Their endpoint coordinates are `(1,0,0)`, `(0,1,0)`, `(0,0,1)`.
The two links to each cell are distinct periodic flux routes. Their
actual kinetic coefficients are respectively `+i/4` for outward source
terms and `-i/4` for inward source terms. Refilling the origin uses the
forward branch of an outward term or the adjoint branch of an inward
term. Hermiticity of each actual Pauli matrix gives the same spin
contraction, with the actual orientation and seam phases retained.

For each route e ending at y, define P_e to require the selected link
excited, every other link blank, and the true local matter ground at
every cell except y. It is identity on the entire thirty-mode/five-scalar
carrier at y and on the retained reference. This is a gauge-invariant
effect on the complete carrier. The six P_e are pairwise orthogonal
because their excited-link patterns differ. Their complement
`I-sum_e P_e` is retained in full, including states with other link
patterns, excited nonarrival cells and arbitrary arrival-cell matter.

The source at the origin returns to its local ground on arrival. The
reference does not move. For e^c removal its hypercharge label 6Y is +6;
the arrival matter has 6Y=-6. For L removal the corresponding values
are -3 and +3. The selected flux has the conjugate reference irrep at
its origin endpoint and the native removal-field irrep at its arrival
endpoint. Thus reference plus origin flux is a singlet, and arrival
matter plus arrival flux is a singlet. Every other site obeys its
uncompensated source Gauss constraint. The origin obeys the joint
source/reference constraint; source-only Gauss zero at that origin is
not asserted in the charged branch.

Raw lepton number is conserved by the actual clock. Every local source
ground has raw L=1. e^c removal raises total source raw L from eight
to nine, while the reference has raw L=-1. On arrival the seven other
cells have L=1 and the target has L=2. L removal lowers total source
raw L to seven and gives reference L=+1; the target then has L=0.
This is an actual relocation of the source defect with retained flux,
rather than changing the position label of a boundary reference.

All links are blank just after localization and pumping. H0 preserves
every link label. Therefore the exact operator identity

`P_e exp(-it H0) B_r P0=0`

holds on the entire input carrier for all real t. In particular every
arrival probability at t=0 is exactly zero, even with an imperfect
input preparation. The zero-hopping clock is a distinct control with
zero arrival at every time.

## Exact first-hop coefficient

Let u be the actual native-coordinate cube trial in P0. The source
ground theorem gives a phase with `||g_lambda-u||<=eta_g`, where
eta_g is approximately 0.00016931925. The local ground one-body tensor is

`P c_(s,a)^dagger c_(t,b) P = delta_ab[delta_st I/2 + b_r sum_i sigma_i,ts S_i]`,

with `b_ec=beta_ec` and `b_L=(1-beta_ec)/2`, evaluated from the exact
local ground algebraic field. No two native irreps have opposite
hypercharges. This excludes an origin pair-annihilation compression
to its neutral ground. The only hopping contribution surviving P_e
refills the removed family at the origin and removes that same family
at y. Every other native family and every other link has zero
contribution to this selected first-hop image.

The actual 119-state link coordinates satisfy

`<blank|U_ab^dagger U_cd|blank>=delta_ac delta_bd/d_r`.

Their forward and adjoint blank images are orthogonal. The code checks
this complete Gram matrix for both families. CAR grading gives only a
common phase in the selected refilling amplitude: all initial local
ground factors have fixed odd parity, and the compressed origin
one-body operator is even. No grading phase is omitted from a
probability interference term.

Contracting the actual spin matrix, retained reference components,
link Gram matrix and two local ground tensors gives the complete
two-spin first-hop effect

`M_axis=(1+3b^2)I/128 + b^2(1-b) S_0 dot R_axis S_y/16`,

where R_axis is +1 along the edge and -1 in its other two spin
directions. The implementation verifies all three 4-by-4 matrix
identities exactly. The factor 1/128 includes the actual hopping
coefficient squared 1/16 and the complete pump normalization.

In the unique cube ground, every one of the twelve undirected edge
correlations is epsilon_cube/12. The native bond-axis rotation converts
the twisted correlation above into that cube correlation. With
`x=4 epsilon_cube`, the exact route moment is therefore

`m_r=||P_e T B_r u||^2=(1+3b_r^2)/128+b_r^2(1-b_r)x/768`.

Exact root intervals give m_ec approximately 0.00774300793734 and
m_L approximately 0.01171374430756. Every route has this atomic
comparison moment. Equality of the full interacting route probabilities
is not required for the bounds.

## Full-H time bound

The full hopping norm is tau=360. Each selected actual oriented edge
has norm at most 15. From the full local norm 57, twenty-four unit link
penalties and `|E0|<12`,

`||H0-8E0||<=B0=576`,
`||H_lambda-8E0||<=B0+lambda tau`.

On the successful one-reference sector, subtracting the explicit
constant reference energy four in this norm estimate only removes a
common clock phase. That energy is retained in the energy records.

For any `psi=B_r P0 z` of norm at most one, the exact zero-kernel
identity and Duhamel's formula give

`||P_e exp(-it H_lambda)psi + i lambda t P_e T psi||`
` <= lambda tau t^2(2B0+lambda tau)/2`.

Both H and H0 in this argument act on the full source complement.
There is no truncated charged propagator. Only the selected edge can
contribute to `P_e T B_r P0`, so comparison with the cube trial costs
`15 lambda |t| eta_g`, rather than the full hopping norm times eta_g.
For an exact initial g_lambda the arrival amplitude is consequently
within `lambda |t| e_exact` of a vector of norm
`lambda |t| sqrt(m_r)`, where

`e_exact=15 eta_g+tau |t|(2B0+lambda tau)/2`.

## Physical preparation error before squaring

The preparation error must not be subtracted as a global probability
error, which would exceed this short-time signal. Retain the actual
preparation environment, including the second source. Its conditional
pure state has a purification comparison with `g_lambda tensor e`
of vector distance at most `epsilon_p=2 sqrt(1-F)`. The receipt uses
the certified upper value `epsilon_p<=4.711844974e-7`.

After P0 and B_r the distance to `B_r u tensor e` is at most
`eta_g+epsilon_p`. The exact no-arrival kernel holds for this actual
prepared state as well. Its amplitude error is bounded by

`e=15(eta_g+epsilon_p)+tau |t|(2B0+lambda tau)/2`.

Thus if m lies in [m_low,m_high], every selected route has conditional
probability, before normalization by localization and pumping, in

`[(lambda t)^2 max(0,sqrt(m_low)-e)^2,`
` (lambda t)^2 (sqrt(m_high)+e)^2]`.

Use the declared source time `t=2^-24`. Both lower amplitudes are
strictly positive. Multiplication of the lower probability by the
physical preparation success lower bound gives a strictly positive
unconditional probability from the original physical vacuum input:

| Field | One retained flux route | Union of all six routes |
| --- | --- | --- |
| e^c | [6.4300047848e-38,1.2745542710e-37] | [3.8580028709e-37,7.6473256254e-37] |
| L | [1.0483425938e-37,1.8251321419e-37] | [6.2900555633e-37,1.0950792851e-36] |

The probabilities are small because this certified source has tiny
hopping and the selected clock time is short. They are finite positive
same-H events. No physical speed or macroscopic transport law is
inferred from them.

## Complete retained law and the energy setting

Let a be the P0 success probability conditional on successful physical
preparation, and p_e each route probability conditional only on that
preparation. The graph norm r and the purification error give
`1-a<= (r+epsilon_p)^2`. The pump split is exactly a/2 on each branch.
There are ten coarse records, with exact probabilities

1. Physical preparation failure: `1-q`.
2. Atomic localization failure after preparation: `q(1-a)`.
3. Localized pump failure: `qa/2`.
4. Six separate flux-route arrivals: `q p_e` for each actual edge.
5. Complete pump-success nonarrival: `q(a/2-sum_e p_e)`.

Their sum is exactly one. The last probability has a strictly positive
lower bound. Every source copy, scalar/link memory, reference Fock
occupation, ancilla and failure output is retained. The preparation
controllers preserve source Gauss on all their branches. The charged
pump and subsequent detector preserve the joint source/static-reference
constraint and all other source constraints.

The same localized packet `B_r P0` also admits the actual energy-bin
measurement of interacting_charged_response. Since
`||B_r P0 g_lambda-B_r u||<=eta_g`, the same full-H window argument
applies without altering its bin endpoints. The inherited preparation
trace-distance bound expands the energy probabilities slightly and
preserves positive line and outside weights. Here
`p_line+p_outside=a/2`; localization and pump failures complete the law.
The static reference energy four is explicitly included in the total
energy and explicitly subtracted when reporting source energy.

Energy and arrival are separate observer settings on this same
preparation and clock. Their projectors need not commute. Separate
marginal bounds do not establish that the selected energy-bin component
has a positive arrival probability, so energy-resolved propagation is
not claimed. There is also no charged momentum band, free-particle
dispersion, rest-mass identification, larger-size transfer or continuum
limit in this finite side-two result. The neutral block-coupling
extension is a distinct source construction.
