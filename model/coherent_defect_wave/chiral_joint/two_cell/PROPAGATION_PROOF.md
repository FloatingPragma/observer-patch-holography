# Controlled propagation on complete quantum two-cell sources

This construction measures actual finite motion under the uniform open
Hamiltonian of [the spectrum proof](SPECTRUM_PROOF.md) and, separately, the
two-oriented-link periodic Hamiltonian of [the symmetry proof](SYMMETRIES.md).
The two sources retain their distinct link memories and generators. Neither
uses a scalar expectation or a classical link as a replacement for a quantum
operator. The preparation below is physical but nonstationary; it is not
the neutral ground used by the charged spectral instrument.

## A prepared local pair with an empty receiver

Let Omega have empty lepton occupations, blank links and both scalars empty.
It lies in every local Gauss kernel, with L=0, 2Sz=0 and even matter parity.
Let `K_right=H_Y,right` be the actual quantum Yukawa operator at the right
cell, excluding its scalar-number term. Direct complete source action gives

`||K_right Omega||^2=4`, `K_right^2 Omega=4 Omega`.

Therefore the normalized state `psi_right=K_right Omega/2` can be prepared
by the complete unitary `exp[-i pi K_right/4]`, up to a global phase. The
controller acts through the same local even gauge-invariant pair operator.
It has no discarded branch. Preparation is a driven intervention; its work
supply is not calculated or confused with the subsequent idle clock.
The resulting state contains a right-cell matter pair correlated with its
quantum scalar, and the left cell is empty. All twenty-four actual local
Gauss generators annihilate it. Its L, 2Sz and matter parity remain zero,
zero and even. The periodic version uses the same local preparation with
both oriented quantum links blank.

The receiver detector asks whether left matter occupation is odd. Its two
projectors are R and `I-R`. They commute with all local Gauss currents:
source number commutes with its internal currents, and the scalar/link
currents act on independent factors. Local Yukawa terms change left matter
number by zero or two, so this detector cannot confuse a purely local pair
creation with transfer of one member of the prepared pair.

The exact initial responses are

`R psi_right=0`,

`||R H_open psi_right||^2=3/4`,

`||R H_periodic psi_right||^2=3/16`.

The receiver-odd term comes from actual link-assisted hopping of the prepared
right-cell matter. Its link occupation and the remaining source/scalar
records compensate the transported gauge charge. These derivatives alone
are not the finite-time conclusion; the remainder bound below controls the
entire declared unitary evolution.

## Complete finite-time probability bounds

The open full carrier has dimension 614400 and `||H_open||<=29`. The
periodic full carrier has dimension 3686400 and `||H_periodic||<=30`.
For the latter bound, forty quarter-unit hopping monomials including
adjoints contribute at most ten, the two local Yukawa operators at most
sixteen, and the two scalar plus two link projectors at most four.
The bound applies to every occupation, without restricting to a cyclic or
low-energy sector.

For either self-adjoint H, at `t=1/20` source time use the exact degree-eight
state polynomial

`psi_poly=sum_(n=0)^8 (-it)^n H^n psi_right/n!`.

The integral remainder and unitarity give the state-norm bound

`||exp(-itH)psi_right-psi_poly|| <= delta=(t B)^9/9!`,

where B is 29 or 30. The polynomial is calculated from the actual sparse
source columns, including every destination. It is not formed by projecting
onto a guessed moving-particle subspace. Its twenty-four Gauss actions are
also checked exactly, so the polynomial and exact evolution retain the same
physical constraint.

Write `p_poly=||R psi_poly||^2`, and let u be an exact rational upper bound
on `sqrt(p_poly)`. Projection is contractive, so the exact probability p obeys

`max(0,p_poly-2 delta u) <= p <= min(1,p_poly+2 delta u+delta^2)`.

The lower bound drops a nonnegative square term; it remains valid when the
triangle lower bound would be negative. The executable computes u by an
integer square root rounded upward. It then rounds the probability bounds
outward on the `2^-36` grid. Both final lower bounds are strictly positive.
The no-arrival outcome is retained with interval `[1-p_upper,1-p_lower]`.
The receipt reports the exact rational intervals and the full polynomial
support, rather than interpreting floating output as a proof.

## Controls that separate transport from local response

Remove only the open-source hopping terms and retain both local quantum
Yukawa operators, scalar energies and link energy. Each surviving operator
preserves left matter parity. Starting from psi_right, the receiver-odd
probability is therefore exactly zero at every time. The complete local
occupation/scalar domains verify the parity statement; the executed
finite-time polynomial is an independent zero-response control. The
positive response of the full H consequently requires transport through the
quantum link, rather than local pair creation at the receiver.

For the periodic source, its exact antiperiodic translation T sends
psi_right to the corresponding left preparation, exchanges the two receiver
projectors and commutes with the complete H. Thus the full opposite-direction
arrival laws agree at every time, with the translated preparation and
apparatus. The code verifies that relation on the complete finite-time
polynomial as well. The periodic ground and charged spectrum are not
borrowed from the open-source calculation.

These are controlled finite propagation records on explicitly prepared
physical states. They supply movement in the declared coordinate/register
readout. They do not establish an asymptotic scattering matrix, a
three-dimensional interacting continuum, a Lorentz light cone, a measured
particle identity or a relativistic rest mass. The finite physical clock
and coordinate assignment retain their specified source units.

## Executable sources and theorem interface

[propagation.py](propagation.py) computes the preparation, complete source
polynomial and outward probability intervals. It applies the actual open
columns in [hamiltonian.py](hamiltonian.py) and periodic columns in
[symmetries.py](symmetries.py), with the actual local currents from
[basis.py](basis.py). [test_propagation.py](test_propagation.py) checks both
sources, exact Gauss preservation, the local-only zero control and the
translated periodic experiment.

The finite momentum result in
[PeriodicMomentum.lean](../../../../Lean/QFT/PeriodicMomentum.lean)
consumes the supplied identity `HT=TH`. It then gives conserved polynomial
phase effects on the periodic source. The finite Taylor remainder above
is a direct self-adjoint spectral-calculus bound. Neither step transfers
the open ground to the periodic carrier or invokes a continuum locality
theorem.
