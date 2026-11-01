# One spin factor survives both joined-field dynamics and complete repair

The unchanged six-site Weyl/contact field and the existing complete charged
block-repair operations preserve **one common full-Fock spin factor**. The
previous localized carrier still leaks, but its internal information need
not be lost: the two-hole spin-one sector has **105 multiplicity states
times a protected qutrit, dimension 315**. This is a finite same-process
protection result for CP-03/#1139, with spectrum work for CP-07/#1143.

The result concerns internal spin on a declared finite chain. It does not
make the old localized pair invariant, conserve the joined energy through
repair, or identify a physical particle. No primitive, field coupling,
reference column or repair outcome is changed.

## Exact common-factor theorem and its actual premises

Use three identical x-directed dimers, the original `c=1/8, U=1` action on
the joined six-site chain, and `K=H-N/2`. The global hidden spin is

```
2 S_a = sum_x (-1)^(sum_(j!=a) x_j) c_x* sigma_a c_x.
```

Along this chain its onsite frame repeats every two sites. Consequently
global spin is the sum of the same local dimer spins. The existing even
ground `g` is an exact spin singlet; both lower charged projectors are
doublets. Their projectors commute with every local spin generator, and
the unchanged graded whole-block swap commutes with the sum of the spins
on the swapped blocks. Every native C/E support projector, hold and
movement therefore commutes with global spin. So does the outside outcome
`Q=I-P`. The joined Hamiltonian commutes with that **same** spin, as verified
by exact Gaussian-integer full-Fock commutators.

These facts are proved locally without relying on a floating eigenvalue
test: the two analytic columns spanning the ground's singlet sector are
annihilated by all three integer spin generators; the six charged-projector
commutators and three complete 256-state swap commutators vanish exactly.
Tensor products and sums establish the full repair statement. Separate
numerical checks evaluate all **156 complete repair outcomes**, including
zero and outside slots, on all 4,096 states, and their **468** spin
commutators, under admission `2e-11`. Their zero numerical residuals are
checks of the exact proof, not outward floating-point certificates.

In each fixed number and total-spin sector, construct an orthonormal basis
of highest-weight vectors and their normalized lowering ladders. An
operator commuting with `S_z` and `S_+` maps the highest-weight space to
itself; commutation with `S_-` makes its multiplicity coefficients identical
on every ladder weight. Thus every admitted operator has the form

```
H = h_(N,S) tensor I_(2S+1),
L_r = l_(r,N,S) tensor I_(2S+1).
```

Completeness gives `sum_r l_r* l_r=I`. Any finite interleaving of this field
evolution and complete repair therefore leaves the reduced spin/reference
state unchanged, for arbitrary initial correlations with multiplicity and
an arbitrary untouched reference. Record-dependent choices of the admitted
operations also compose. Individual conditioned records can steer an
initially correlated state; the statement is about the complete output.
Merely commuting with `S²`, or covariance of an averaged channel, would
not suffice. A complete spin-z read keeps the spin-one sector but destroys
the qutrit/reference coherence, with trace distance exactly `2/3` on a
maximally entangled qutrit control.

The factor is exhibited, not claimed to be the largest possible joint
noiseless algebra of the field/repair process. Neither its global decoder
nor a new physical spin apparatus is implemented here.

## The leaky carrier is inside a larger protected sector

For six sites and fixed number `N`, the spin multiplicity is

```
m(N,S) = C(6,N/2+S) C(6,N/2-S)
       - C(6,N/2+S+1) C(6,N/2-S-1),
```

with inadmissible binomial indices zero. Highest-weight subtraction proves
this count; the full sum over number and spin is exactly 4,096. At `N=4`,
`S=1`, it gives 105 copies. The old three-position triplet carrier has
dimension nine and lies inside the 315-dimensional sector. It need not
remain in its initial three-dimensional multiplicity subspace.

The complete matrix experiment starts with the old localized hole-pair
triplet entangled with a three-state reference and executes

```
U(1) -> E repair on blocks (0,2) -> U(1/4) -> C repair on blocks (0,1).
```

The three block labels can be mapped to the literal native triangle
`(0,1,5)`. All **153 joint outcomes** remain, including both outside events
and every zero raw slot. After the first evolution the old small-carrier
outside probability is **0.005297146587459**. Nevertheless the complete
spin/reference marginal agrees with its initial value within `2e-10`.
Independent tests also use a general correlated 105-by-3-by-2 input.
This is evaluation of the same complete ideal operators, not native pulse
execution on a twelve-mode processor or a physical experiment.

The earlier two-dimer leakage identity
`L*L=(3 P_singlet+P_triplet)/512` remains correct. It excludes invariance of
the four-column local carrier, not the larger full-Fock spin factor.

## Complete finite spectrum and candidate probes

All number sectors and every highest-weight energy are retained. Their
spin-degeneracy-expanded dimension, trace and second trace moment match the
full action. Direct number-sector diagonalization and time-domain probes
independently check the spectrum. Numerical charts/energies use explicit
residual admissions; they are not certified interval eigenvalues.

The comparison reference is the **unique finite half-filled singlet ground
of the joined K**, with numerical energy `-3.057619842652825`. This supplied
pure reference is distinct from the product of dimer grounds and is not
faithful. It is not newly source-prepared in this packet. The local probes
are written in the exact Hubbard frame, with state and operators transformed
together. Complete response weights, including tiny and dark lines, remain.

| Joined-field candidate | Relative CAR number | Internal spin | Lowest model gap | Lowest multiplet probe weight |
| --- | ---: | ---: | ---: | ---: |
| End-site removal | -1 | 1/2 | 0.402792741312 | 0.043390177952 |
| End-site addition | +1 | 1/2 | 0.402792741312 | 0.043390177952 |
| Adjacent equal-spin pair removal | -2 | 1 | 0.847193432532 | 0.000402726621 |
| Adjacent equal-spin pair addition | +2 | 1 | 0.847193432532 | 0.000402726621 |

Each charge-one response keeps 210 highest-weight lines; each pair response
keeps all 105. These finite nonzero numerical residues attach actual probes
to the same action. The lowest neutral triplet gap is about
`0.007561476667`, below the charged gaps. This finite lattice spectrum is
neither a continuum mass gap nor a photon construction. Internal spin is
not Lorentz spin; its degeneracy is not a generation count. The fixed
repair projector acts on charge-two sectors: the charge-one candidates
take its **outside outcome only**. They have no nontrivial native repair
readout in this construction.

Spin protection also does **not** make these Hamiltonian poles stationary
under the complete repair process. For the retained E repair, the joined
energy pullback differs from K by Frobenius norm about **0.407239684974**.
The controller's energy/work account remains necessary. No new absolute
physical mass, species identification or empirical comparison is supplied.

## Exact spatial-frame criterion and the three-dimensional boundary

For a coordinate-identified block translation by `delta`, the hidden-spin
axis `a` changes by sign `(-1)^(sum_(j!=a) delta_j)`. Therefore all three
axes agree **if and only if**

```
delta_0 = delta_1 = delta_2 modulo 2.
```

Copies of a block with even side lengths, such as `2x2x2`, satisfy this
criterion under their tiling translations. Odd/odd/odd translations also
align. If the block reference/candidate projectors are spin invariant and
every complete repair uses their aligned permutations, the same proof
applies. This is a constructive criterion; it supplies no new cube
candidate, source run or replacement of the existing dimer law.

The unchanged transverse dimer placement on a `2x2x2` cube fails the
criterion. Assign its four x-directed dimers to native ports `(0,1,5,7)`.
The two actual E movement slots on seams `(0,1)` and `(0,5)` have exact
commutator Gram matrix `diag(64,32,32)` against the three global hidden-spin
generators, on the 24-dimensional two-doublet carrier. No nonzero linear
combination survives those slots. This defeats the automatic extension of
the chain proof to that three-dimensional placement; it does not classify
all possible noiseless algebras or exclude a properly specified future
aligned-block construction. No frame-dressed swap is silently substituted.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.common_spin_protection --check model/common_spin_protection/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/common_spin_protection
```
