# Coherent compilation with the complete native kinetic Hamiltonian

## Source, invariant coordinates and averaging

Use the full CAR source in
[lepton_mass_compiler](../lepton_mass_compiler/PROOF.md), with ninety
physical modes and 324,000 enlarged modes at each site of a finite periodic
cube. Its complete analyzer supplies W=U*QU. On character k, W has rank
one, and its physical-address compression is the exact native e_k.
The neutral scalar background norm is the declared value one. The phase
2/9 and three-character architecture retain their retrospective ancestry.

Let K be the homogeneous nearest-neighbor Weyl Hamiltonian on every
internal, character and address copy, with one-body symbol
`-sigma.sin(ap)/a`. Let V be the complete analyzer-conjugated pairing,
with full-Fock norm 6s on s sites. Every pulse uses this same K+V.
All fermion occupations, auxiliary addresses and untouched references
are part of the domain.

If N_A counts occupied ancillary addresses, R=exp(i pi N_A/2) commutes
with K. The average of `H_r=K+R^(-r)V R^r`, r=0,1,2,3, is
`Hbar=K+Vbar`. Pair monomials with one or two ancillary factors have
nonzero N_A degree modulo four and cancel. For the embedding J that
places every ancillary mode in vacuum,

```
Hbar J=J(K_physical+H_Y).                         (1)
```

This is a reducing-subspace identity on the entire physical Fock space,
including arbitrary entanglement with external references. It retains
the original analyzer coefficient Y; the kinetic action and clock do
not change at the embedding.

The span of the physical flag and W acting on that flag is two-dimensional
for each character. It is invariant under W and the phase rotations.
Homogeneous K is identity on flags and characters, so it preserves this
span at every site. The paired native lepton support has two internal
coordinates and two spins. Thus 24s CAR modes carry the noncommuting
calculation. Every other enlarged mode is an exact free kinetic
spectator. Its evolution commutes with the active generators and occurs
for the same total signed target time in the ideal and compiled laws.
This factorization is an operator identity, including states entangled
between active and spectator factors; it discards no mode.

Each one-body kinetic eigenvalue has magnitude at most3/a. There are
12s positive and12s negative rows counting multiplicity. Therefore
the full-Fock active kinetic norm is at most36s/a, and

```
||H_r,active|| <= B=36s/a+6s.                     (2)
```

## Separate sign controls for kinetic and pairing terms

Let I be the exact exterior lift of the site permutation x -> -x on the
periodic cube. It interchanges each forward hopping with its adjoint,
so `I K I*= -K`. Uniform onsite pairing is unchanged: `I V I*=V`.
The CAR exterior permutation includes its occupation-dependent signs.

Let S be parity of the complete native lepton-doublet leg. Every pairing
monomial contains one factor on that leg, whereas kinetic hopping
preserves it. Consequently `S V S=-V` and `S K S=K`. Both controls
commute with R and with each other. Their product C=IS obeys

```
C H_r C*= -H_r.                                 (3)
```

Thus a negative segment is implemented by conjugating a positive-duration
source segment by C. Leg parity alone does not reverse K+V. The exact
256-state control checks this failure, all four equations(3), every CAR
permutation identity, and a mutation that drops the fermionic signs.
An independent real-space check uses the complete periodic Weyl stencil.
These are declared coherent parity and spatial-permutation services.
No local primitive-repair implementation or free control overhead is
asserted. The construction does not require switching K off.

## Fourth-order product and complete error

Let S2(h) be the symmetric seven-factor product with letters
`(H0,H1,H2,H3,H2,H1,H0)` and weights
`(1/8,1/8,1/8,1/4,1/8,1/8,1/8)`. With
`p=1/(4-cuberoot(4))`, set

```
S4(h)=S2(ph)^2 S2((1-4p)h) S2(ph)^2.             (4)
```

The inherited exact noncommuting-word proof checks all 341 coefficients
through degree four. Its letters can be the full H_r, so noncommuting K
does not change the cancellation. The sum of absolute segment weights
is `8p-1<7/3`. Put `R0=(7/3)B`. For `R0|h|<=1`, the two Taylor tails
after degree four bound the local operator error by
`R0^5 |h|^5/20`: use `2 exp(1)/5!<1/20`. Unitary telescoping and the
operator-to-complete-trace inequality give

```
||actual output - ideal output||_1
    <= min(2, R0^5 |t|^5/(10 n^4)),   n>=R0|t|. (5)
```

Equation(5) includes all ancillary outputs and every untouched reference.
No reset, rare-success conditioning or measured intermediate flag is
used. The actual positive source-pulse duration is `(8p-1)|t|`, with
the parity/permutation/phase service durations accounted separately.
Target evolution time t is not thereby identified with laboratory time.

## Executed source blocks and resource meaning

At opposite periodic axial momenta, a spin-oriented pairing block has
four CAR modes: two physical charged legs and their two active ancillary
flags. Its kinetic term is `kappa N`, while its complete pair coefficient
is the actual rank-one W block. The opposite momentum has `-kappa N`.
The full eight-mode Fock carrier realizes the spatial inversion and
leg-parity controls in(3). The other spin-oriented blocks are related
by the fixed native epsilon signs and are retained in the all-site
norm and factorization.

For every native character, the executable evolves all sixteen Fock rows
and all four physical input columns at both axial momenta, including a
maximally entangled four-dimensional reference. It records all four
ancillary outcomes without postselection. High-precision matrix
evolution is a numerical check of the analytic bound, not an interval
rounding certificate. The target single-excitation gap is
`sqrt(kappa^2+e_k^2)`; its pair-sector gap is twice that value.

Reported large-volume round counts follow from(5); they are resource
bounds, not an execution of the enormous enlarged carrier. The declared
free source supplies controlled kinetic propagation with the same
analyzer-derived coefficients. Interacting charged minima, a selected
Higgs vacuum, mirror removal and observed-lepton identification are
separate assertions.
