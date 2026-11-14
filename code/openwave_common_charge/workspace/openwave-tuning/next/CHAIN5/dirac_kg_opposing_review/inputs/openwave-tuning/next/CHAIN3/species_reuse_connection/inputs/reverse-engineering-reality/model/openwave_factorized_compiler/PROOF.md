# Complete coherent compilation from native momentum factors

The source, controls and physical target are those of
[the kinetic compiler](../openwave_kinetic_phase_cycle/PROOF.md). Each cell
retains90 physical and324000 enlarged CAR modes. The bound applies to the
entire enlarged source with arbitrary external quantum references; no
ancillary measurement, reset or successful-outcome restriction is used.
The phase2/9, empirical electron anchor, uniform background and coherent
control services have the same ancestry as that construction.

## Exact even-factor decomposition

The invariant span of the physical flag and its image under the native
rank-one analyzer has dimension two for every character. Homogeneous
kinetics is identity on this flag. Fourier transformation pairs an l-leg
mode at momentum p with an r-leg mode at −p. Choose a common spin rotation
for that pair; invariance of the spin epsilon tensor preserves the pairing.
The eigenvalues of the Pauli kinetic matrix are ±|sin(ap)|/a. The paired
opposite-spin modes therefore have the same kinetic coefficient κ.
The two spin orientations, three characters and s momenta give6s disjoint
four-mode CAR factors. Self-inverse momenta also have two separate
opposite-spin leg pairs. All remaining CAR modes are free exact factors.

A factor has two flags on each of its two legs and Hamiltonian

```
H_r=κN+R_r* V R_r,
V=Σ_ij W_ij l_i* r_j* + adjoint,
W=|u><u|, ||u||=1.
```

Changing the one-particle flag coordinates reduces V to the creation and
annihilation of one normalized pair, so its full sixteen-row Fock norm is
one. The number spectrum is0,1,2,3,4. Subtracting the common scalar2κI
in both actual and target evolutions gives

```
||H_r−2κI|| ≤ 2|κ|+1 ≤ 4/a+1.
```

This scalar subtraction is only an error estimate: the original source
energy and reference are retained. Actual and ideal products have identical
total signed time and hence the same scalar phase. All even operators on
disjoint CAR factors commute. Each complete-source pulse, its four-phase
mean and the entire fourth-order word factor exactly into these even
blocks and the identical free-spectator evolution. Fermionic permutation
signs are included in the exterior coordinate change. No physical Fourier
or helicity operation is inserted into the apparatus; they prove an
identity for its existing homogeneous pulses.

## The complete-output bound

The exact fourth-order noncommuting-word certificate of the kinetic
compiler applies separately to every centered block. Its sum of absolute
segment weights is less than7/3. Let R=(7/3)(4/a+1). For R|t|/n≤1,
the block unitary error is at most R^5|t|^5/(20n^4). A tensor-product
unitary telescoping identity sums the6s errors. Multiplication by two
bounds the complete unhalved output trace distance, including any
reference, by

```
ε ≤ min(2, 6s R^5 |t|^5/(10n^4)).                (1)
```

The bound is linear in volume at fixed n and time. It requires disjoint
native factors; an interaction between factors invalidates the exact
decomposition. The executor adds such a coupling as a failing control.
Contractivity transfers(1) to every complete final instrument. A physical
energy marginal A obeys |Tr A(ρ_actual−ρ_target)|≤||A||ε without
conditioning on ancillary outputs.

## Localized energy and execution

The [compact packet](../openwave_localized_rest_packet/PROOF.md) has
physical excitation operator norm at most282s and mean excess bounded
by d. Choose ε≤d/(10^6·282s), using(1) to select an integer cycle count.
The energy transfer error is at most d/10^6. Its exact lower kinetic
moment, localization and phase records retain their original source
domain. The cycle counts are analytic resource bounds on the huge carrier,
not an execution of it. Coherent inversion, parity and phase services
retain their separate implementation costs and continuing controls.

The finite test evolves both sixteen-row factors on their full Fock spaces,
retains all256 rows and a sixteen-dimensional untouched reference, and
records all sixteen ancillary outcomes. Independent tests verify exact
rank-one norm identities, global-versus-factorized pulses, reference trace
bounds, volume scaling and rejection of inter-factor couplings. Numerical
matrix diagnostics are separate from the rational bound(1).

## State-sensitive physical energy estimate

The prepared packet permits a stronger energy estimate without changing its
complete-output instrument. Let A≥0 be its physical excitation energy,
tensored with identity on all ancillary/reference outputs, with ||A||≤Emax.
If the actual and ideal isometries differ by at mostδ and the ideal mean
is at mostμ, purification and Cauchy--Schwarz give

```
|Tr A(ρ_actual−ρ_ideal)| ≤ 2 sqrt(Emax μ) δ+Emax δ². (2)
```

Indeed, write the actual purified vector as v+w with ||w||≤δ, and expand
its A expectation. The cross term is bounded by
2||A^(1/2)v|| ||A^(1/2)w||. This includes arbitrary purifying references;
it bounds the physical energy marginal, not the energy stored in ancillary
hardware. For the conditional initial packet, μ≤m_upper+d, Emax=282s,
and δ≤3sR^5|t|^5/(10n^4) before saturation. Choose n so(2)≤d/10^6
and the complete trace bound is at most10^-12. Every original preparation
failure remains in the instrument. The explicit narrower state domain is
required for this energy improvement; the trace bound(1) remains uniform
on all physical states. The separate conserved phase of the target does
not establish an autonomous phase clock for continuing compiler controls.
