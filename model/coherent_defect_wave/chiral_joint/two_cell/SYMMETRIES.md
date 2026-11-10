# Spatial symmetry, charge conjugation and complete finite momentum

The open control is the full uniform two-cell source in `hamiltonian.py`,
with twelve lepton CAR modes, one six-state quantum link and two five-state
Higgs memories. Its dimension is 614400. A key is
`(mask,link,H_left,H_right)`. The actual charged-lepton Yukawa, scalar number
and link excitation coefficients are one at both cells. Its hopping is

`H_hop=(i/2) sum_(spin,a,b) sigma3_spin c_left,spin,a^dagger U_ab c_right,spin,b + h.c.`,

where the five admitted U_ab are the literal link raising matrices
`|a,b><0|`. All statements below retain every finite occupation. The
periodic completion is an explicitly different source, with two such links.

## 1. The single open link does not implement particle translation

A spatial operation implementing the same source fields must act on their
full CAR algebra, including occupied inputs. For example, the antiperiodic
two-cell shift is `T c_left T^-1=c_right`,
`T c_right T^-1=-c_left`. If it also swaps the two scalar memories, invariance
of the hopping requires `V U_ab V^-1=U_ba^dagger` on the link. An ordinary
site exchange has the same requirement with an additional phase. A spin
half-turn can change the sign of sigma3 but cannot remove the adjoint.

No unitary or antiunitary V on this six-state link has that property.
Indeed,

`sum_(a,b) U_ab^dagger U_ab=5 |0><0|`,

`sum_(a,b) U_ab U_ab^dagger=I-|0><0|`.

Their spectra are `(5,0,0,0,0,0)` and `(0,1,1,1,1,1)`; conjugation cannot
exchange them. Index permutations and phases do not affect the obstruction.
An operator implementing the stated complete CAR transformation differs
from its Fock implementer by an operator in the CAR commutant. Requiring
the actual scalar swap also fixes the scalar algebra, leaving precisely
the link operation tested above. Thus allowing an unspecified entangling
unitary does not evade this obstruction while retaining these field and
apparatus identifications. Arbitrary re-encoded observables have a separate
contract and are not excluded.

The executed candidate that merely transposes `|a,b>` fails directly.
On `|one right spin-up e^c, blank link, blank scalars>`, the coefficient of
`|one right spin-up e^c, occupied e^c link, blank scalars>` in `(HT-TH)`
is `-i/2`. Its conserved L and Sz sectors therefore do not supply momentum
through an unproved site-exchange assignment.

## 2. Axial spin labels and an exact Kramers symmetry

Define `J_+=sum_(site,a) c_site,up,a^dagger c_site,down,a` and
`J_z=(N_up-N_down)/2`. The real spin-singlet Yukawa interaction conserves
spin, while the one-axis hopping selects sigma3. Direct execution on one
right spin-down e^c gives coefficient `-i` in the left spin-up occupied-link
component of `[J_+,H]`. Thus `J_z` conservation alone does not give a
commuting full SU(2) spin action on this spatial control. Physical spin in a
relativistic theory requires the associated spatial transformations and
little-group action, rather than a count of internal labels.

There is nevertheless a useful exact antiunitary. Let Theta conjugate all
coefficients and act on source modes by

`Theta c_site,up,a Theta^-1=c_site,down,a`,

`Theta c_site,down,a Theta^-1=-c_site,up,a`.

It fixes the link and scalar configuration bases. In the hopping, complex
conjugation reverses i while the spin rotation reverses sigma3; the two
signs cancel. The real epsilon Yukawa tensors are invariant under this
spin rotation. Scalar and electric energies are fixed. Thus `Theta H=H Theta`
for this real-coupling source.

Actual current matrices show that Theta changes the sign of the imaginary
weak generator Q9 and fixes Q8, Q10 and hypercharge. The inactive color
currents vanish. Consequently every same-boundary Gauss kernel is retained.
Theta preserves L, reverses Sz and satisfies
`Theta^2=(-1)^N_matter=(-1)^L`. The code checks every one of the 49152
annihilator basis columns and all 4096 parity identities. Each odd-L energy
therefore has an orthogonal same-energy Theta partner. This is the finite
Kramers argument: an antiunitary squaring to minus identity cannot have an
eigenvector as its own partner. It proves degeneracy of the charged odd
sectors without treating a small numerical splitting as evidence. It does
not identify physical Lorentz spin or supply a dispersion relation.

## 3. An exact opposite-boundary particle/hole correspondence

Define another antiunitary A. It conjugates coefficients, exchanges the
cells, and acts by

`A c_site,spin,a A^-1=i c_other_site,spin,a^dagger`.

It transposes the link label `|a,b> -> |b,a>` and exchanges particle and
antiparticle states in each scalar memory while swapping the memories.
The link raising matrix stays raising, so this operation does not invoke
the impossible link-adjoint map of section 1.

For the hopping, coefficient conjugation changes i to minus i. Reordering
the transformed annihilator and creator gives another minus sign. The
transposed link indices restore the original term. For the scalar coupling,
the spin-singlet Weyl pair becomes its adjoint and phi becomes phi-dagger;
the complete Yukawa term is unchanged. All number energies are invariant.
Thus `AH=HA` on the full carrier, not only on a vacuum-cyclic subspace.

The implementation fixes `A|empty>=|all twelve occupied>` and constructs
all other columns by the transformed creators, retaining their complete
CAR permutation signs. It checks all 49152 annihilator columns and finds
`A^2=(-1)^N_matter`. The sector map is

`L -> 4-L`, `2Sz -> -2Sz`,

`(q_left,q_right) -> (-q_right,-q_left)`.

Here q is the actual compensating hypercharge entering each local Gauss
constraint. Each native matter, scalar and link endpoint current gives
`A G_left(q) A^-1=-G_right(q_transformed)` and conversely. The code verifies
these identities with the actual twelve currents and both boundary signs.
It also checks all twenty normally ordered hopping terms, the complete
local Yukawa input domain and occupied spectator controls.

Crucially for spectroscopy,

`A c_right,e^c,up^dagger A^-1=-i c_left,e^c,up`.

Right addition and left removal therefore have equal complete spectral
measures when the reference and boundary apparatus are transformed together.
For a selected stationary density rho the paired reference is `A rho A^-1`;
a unique invariant-sector ground makes this the same ground automatically.
Energy projectors, reservoir alternatives and failure effects must accompany
the transformation. Reality and positivity of the Born probabilities give
equality of their full retained weights. This equality does not compare
addition and removal at the same pinned right boundary. An unequal pair of
those same-side thresholds is consequently not, by itself, an intrinsic
particle/antiparticle mass difference. A is an exact symmetry of this finite
Hamiltonian; a relativistic CPT theorem is not inferred from it.

## 4. A periodic completion retaining the quantum links

Keep the six-state primitive link, and supply two independent oriented links
U01 and U10. A periodic key is `(mask,link01,link10,H_left,H_right)`, with
full dimension `2^12*6^2*5^2=3686400`. Declare

`H_periodic,hop=(i/4)c_0^dagger sigma3 U01 c_1
              -(i/4)c_1^dagger sigma3 U10 c_0 + h.c.`.

The two local Yukawa interactions and scalar energies are unchanged. Both
link excitation projectors contribute with unit coefficient. The relative
sign is the specified antiperiodic fermion boundary condition. As a separate
classical-background calculation, setting both link matrices to identity
gives coefficient `i sigma3/2` on the original side-two spectral edge. No
quantum vacuum or reference is asserted to implement that substitution.

The unitary T translates `c_0 -> c_1`, `c_1 -> -c_0`, swaps U01 with U10
and swaps the scalar memories. The first hopping term becomes the second
and conversely. A local pair acquires two wrap signs, so the Yukawa terms
also exchange without a sign. All scalar and electric terms are invariant.
The code expands every actual link matrix unit, including the untouched
second-link state, and verifies equality of the complete 240-term hopping
grammar. Full local occupation/scalar controls and occupied two-link
columns verify the same identity using the implemented Hamiltonian.

The local Gauss generator is the sum of matter, scalar and the appropriate
endpoint currents of both incident links. T exchanges the two complete
Gauss generators. This gives an exact translation on the zero-boundary
physical subspace. For charged apparatus, its boundary position must
translate too: `(q_left,q_right)` maps to `(q_right,q_left)`. A fixed charged
right boundary is a different physical sector. A retained apparatus-position
register or a joint source/apparatus translation supplies the complete
translated domain; omitting that register cannot define charged momentum
on the fixed-boundary kernel.

This construction adds one actual link and its electric term. It is a
concrete completion preserving the primitive link design, without a claim
that no smaller alternative design exists. Its spectrum is not obtained
by relabeling the open-edge energies.

## 5. Complete finite momentum and the remaining field test

On every source basis column, `T^2=(-1)^N_matter` and `T^4=I`. The complete
phase projectors are

`Q_k=(1/4) sum_(m=0)^3 exp(-i k m) T^m`,

for `k=0,pi/2,pi,3pi/2`. They are orthogonal, sum to identity and commute
with H. Odd sectors support phases `+/-pi/2`; even sectors support 0 and pi.
A branch convention gives momentum `k/a` modulo `2pi/a` for declared cell
spacing a. No laboratory scale is inferred by taking a logarithm.

An executed physical control prepares the normalized local quantum Yukawa
action on the empty Gauss-singlet state. Every one of the twenty-four local
Gauss actions annihilates the resulting pair/scalar state. The complete
momentum-phase probabilities are `(1/2,0,1/2,0)`; every zero outcome is
retained. This input is explicitly nonstationary and is not a ground state.
A separate odd one-mode algebraic control gives `(0,1/2,0,1/2)`; that naked
source input is not labeled a source-only Gauss state.

The supplied commutator `HT=TH` meets the premise of
`Lean/QFT/PeriodicMomentum.lean::polynomial_conserved`: finite Fourier
projectors and a specified spectral momentum logarithm commute with H.
It supplies an actual consumer of that finite theorem. It does not establish
three spatial translation directions, a common Lorentz shell, a controlled
interacting continuum or a particle mass. The mass-shell identities in
`Lean/Geometry/MassShellKinematics.lean` consume those physical assignments;
they cannot select them from L, Sz or a two-cell threshold. Likewise the
relative-evolution distinction in `Lean/QFT/FiniteUnitaryScatteringNoGo.lean`
permits finite scattering experiments but supplies no scattering theorem
for this Hamiltonian. Those boundaries keep the periodic completion,
charge-conjugate spectral correspondence and physical identification tests
on their respective source domains.
