# Same-candidate charged transport and even spectroscopy

This construction keeps the actual charged doublets of
[the interacting full-Fock dimer](../particle_response/README.md) through
localized transport, exchange and a neutral particle-hole response. Its
reference, charged projectors, energy gaps and local CAR residues are the
same objects as in that packet. It does not assign a different numerical
quantity to the electron-like and positron-like candidates.

The finite transport architecture consists of identical **decoupled** dimers
and supplied fermionic mode swaps. The ground in every unoccupied block is
the same half-filled singlet, rather than the empty Fock vector. The field
charge of a block means its excess occupation relative to that ground.
No transport by the uniform four-site interacting Hamiltonian is claimed.

## Analytic same-candidate columns

Write `t=1/16`, `U=1`, `D=sqrt(U²+16t²)=sqrt(17)/4` and `K=H-N/2`.
In the existing untwisted frame `d_0=c_0, d_1=i sigma_axis c_1`, the ground is

```text
g = a S + b D_pair,
a²=(1+U/D)/2, b²=(1-U/D)/2, a>0, b<0, ab=-2t/D,
E_g=-(U+D)/2.
```

Here `S` is the separated singlet and `D_pair` the symmetric doublon. The
hole's lower doublet has single-occupation columns
`(e_s, i sigma_axis e_s)/sqrt(2)`, with the left site's ordered internal
basis fixing the phase. The particle columns are the existing explicit
particle-hole unitary applied to these hole columns. Reversing the second
site's sign gives the upper doublet. Thus both lower projectors have the
unchanged gap `Delta=D/2-t` and local residue `r=(1+4t/D)/4`; the upper
doublet has gap `Delta+2t` and residue `s=(1-4t/D)/4`.

`algebra.py` constructs these columns analytically and compares them with
the predecessor's complete sixteen-state generator and projectors for all
three axes. No numerical degeneracy or visibility threshold defines the
candidate. Floating matrix residuals are checked against a stated numerical
admission; they are not outward interval proofs.

## Complete graded transport and route composition

For blocks of four CAR modes, let `pi` map each old block to its new place.
Permuting the occupied orbitals and sorting their exterior product gives

```text
Gamma(pi)|b_0,...,b_(L-1)>
  = (-1)^sum_(i<j,pi(i)>pi(j)) N_i N_j
    |b_(pi^-1(0)),...,b_(pi^-1(L-1))>.
```

The sign follows because crossing two blocks interchanges every occupied
orbital of one with every occupied orbital of the other. This proves the
formula on every full-Fock basis column. Linearity proves it for arbitrary
states, including superpositions of block numbers. Exterior lifting gives
`Gamma(pi)Gamma(sigma)=Gamma(pi o sigma)`; equivalently, neighboring swaps
obey the involution, braid and disjoint-commutation relations. Therefore
routes implementing the same block permutation agree. This is a statement
about the specified ideal FSWAP grammar, not invariance under arbitrary
continuous hopping schedules or their dynamical phases.

The two-block operation uses sixteen neighboring mode FSWAPs and keeps all
256 columns. Integer address checks compare it with independent orbital
permutations; tests also check every 4,096 address of three complete dimers.
The same formula defines finite larger routes without allocating their
exponential-size matrices. Such a route specification is not a claim of
larger source execution.

Since `g` has even occupation, the isometry `V_A` encoding an unknown charged
doublet at A and `g` at B satisfies exactly `Gamma(swap)V_A=V_B`. Its internal
spin is transported unchanged; the exterior block has its entire ground
density, not merely the same mean number. Tensoring this identity with any
untouched reference preserves arbitrary logical/reference entanglement.
The two-block generator `K_A+K_B` commutes with the block swap, so the same
gap follows the localized state. For two odd blocks,

```text
Gamma(swap)(v_A tensor w_B) = - w_A tensor v_B.
```

The charge-one hole and particle both satisfy this identity; the charge-two
pair states and the ground supply even controls. For arbitrary internal
states the result is **minus the internal swap**, not a universal scalar
`-1` on every state. Two equal pure internal states give a deterministic
minus sign. Orthogonal internal states have zero swap expectation, and two
independently maximally mixed doublets give expectation `-1/2`. These
distinctions matter for a complete Ramsey experiment.

Fermionic gates and swap networks are established constructions; see
[Bravyi and Kitaev](https://arxiv.org/abs/quant-ph/0003137) and
[Kivlichan et al.](https://arxiv.org/abs/1711.04789). The contribution here is
their explicit attachment to the same OPH candidate projectors, reference,
response and completed source interface, not a new general swap identity.

## Explicit dictionary to the native negative exchange lift

Choose one fixed normalized lower hole or particle column `v` and put twelve
identical dimer blocks at the twelve native port labels. Define the isometry
`W|{a,b}>` to have `v` in blocks a,b and the same even ground `g` everywhere
else. Its 66 columns are orthonormal: occupied and ground blocks have
different number. The complete field has 48 modes, and the embedded state
has relative charge `-2` for two holes or `+2` for two particles.

Writing `n_j` for the native occupied-defect projector and `Q=-1` for a
hole or `Q=+1` for a particle, `W†(N_j-2I)W=Q n_j`. A native defect current
`J` therefore transports CAR number `QJ`. This `Q` is excess field occupation,
not the separate source/core overlap label used by the native memory model.

Every reference block has even parity and each occupied block has odd parity.
Consequently the restriction of a block permutation to `im(W)` is exactly
the exterior-square port permutation used by the existing
[negative native exchange lift](../exchange_sector/STATISTICS.md). In
particular, moving one occupied block across a spectator gives precisely
that lift's ordering sign. The certificate checks all 600 directed movement
columns and 3,960 columns of the sixty proper rotations, using exact parity
addresses rather than allocating the 48-mode space.

This also intertwines complete native records. For each C/E raw flag r,
let `P_r` select its unchanged native input-support columns inside `im(W)`.
For a nonzero current use `L_r=Gamma(swap_ab)P_r/sqrt(2)`; for a hold use
`L_r=P_r/sqrt(2)`. It follows that `L_r W=W K_r^-`. All 720 raw slots,
including zero/abort slots, and their 7,920 supported columns are checked
against the actual native branches. Completeness holds on the specified
66-dimensional embedded carrier. Multiplication gives the identity for
every finite retained record history, and tensoring preserves arbitrary
untouched references. A full-Fock extension of these support projectors
and flags is not asserted or executed.

The dictionary attaches the native negative exchange character to the
**same charged column** used by the dimer spectrum. The larger field carrier
is not literally the native H66 source, and this is not a 48-mode source run.
The selected internal column and supplied permutation/projector services
are part of the construction. Allowing unknown internal doublets replaces
the fixed-column fiber by the corresponding internal-factor permutation;
it does not make every exchange a scalar. Neither the H66 recovery theorem
nor passive defect protection follows from the intertwiner.

## Neutral particle-hole response and its exchange sector

On two ground dimers use the neutral transfer
`T=c†_(A,j)c_(B,j)`, its real quadrature `X=T+T†`, and its normalized current
quadrature `Y=-i(T-T†)`. All commute with total number and parity. The name
current denotes the transfer quadrature of a supplied interblock hopping
control: the decoupled generator itself has no interblock current.

For the complete directed response about `g_A tensor g_B`, factorization of
the addition/removal spectral measures gives three lines:

| Frequency | Directed weight | Either Hermitian-quadrature weight |
| --- | --- | --- |
| `2 Delta` | `r²` | `2r²=(9+sqrt(17))/68` |
| `2 Delta+2t` | `2rs` | `4rs=4/17` |
| `2 Delta+4t` | `s²` | `2s²=(9-sqrt(17))/68` |

The directed weight sums to `1/4`; either Hermitian weight sums to `1/2`.
Every one of the 256 spectral states, including dark and wrong-charge
states, remains in the numerical check. Independent Jordan--Wigner matrices
and complete spectral moments verify the closed forms. The lowest line is
the same candidate particle-hole pair at gap `2 Delta`; the upper bands
are retained rather than folded into an assumed single pole.

Block exchange fixes the ground and sends `T` to `T†`. It fixes `X` and
negates `Y`. Thus each energy component of `X|g,g>` is exchange even, while
each component of `Y|g,g>` is exchange odd. Their identical spectral weights
do not make their exchange states identical. This provides a neutral,
number-preserving link between the same charged spectrum and the observable
exchange representation. A particle-hole pair consists of distinguishable
charge assignments; its symmetric or antisymmetric spatial superposition
does not alone measure identical-particle statistics.

`exchange_pair_column` exposes normalized analytic target columns for these
bands. It is a mathematical definition, not a postselected source protocol
or a claim that spectral filtering is free. An executed preparation must
account separately for its supplied state and controls.

## Source preparation and complete controlled exchange

`source.py` compiles the same analytic columns from fixed-occupation inputs.
It completes their basis only within the relevant number sector and uses
selective complex Givens rotations, expanded in commuting Pauli factors,
plus diagonal phases. The other number sectors remain present. Every
completed encoder preserves number; the expanded factors preserve parity.
The underlying encoded processor pulses are not identified with microscopic
CAR operations. Input particles and an unknown dual-rail internal state are
counted supplied resources, not generated from an empty vacuum for free.

Full sixteen-column source executions check each encoder with every export
complement. The ground encoder costs 4,168 pulses and each charged encoder
1,642. The block route costs 416 pulses, or 1,216 when controlled by a
neutral pointer. Its local complete instruments and all 256/512 successful
global columns are checked independently; complete preparation/transport
histories retain their export failures. Transport of an unknown charged
doublet, including ground preparation, costs 6,226 pulses. The corresponding
two-identical-candidate exchange read costs 4,510. Both pointer outcomes
remain: odd candidates give the minus port, the even ground gives the plus
port, replacing FSWAP by ordinary SWAP gives plus, and a which-path read
gives one half at each port. These are ideal predictions admitted against
full source output columns, not success-normalized probabilities.

`transport_readout` composes preparation, transport and a nondestructive
destination-presence read in one nine-bit source history. The inverse
charged encoder maps the lower doublet to its two occupation seeds; a
source-compiled controlled pointer marks their union, and re-encoding
preserves the unknown internal state. The complete 32-column local read
has blocks `[[I-P,P],[P,I-P]]`, with `P` the same lower charged projector.
It distinguishes that band from the upper doublet, ground and wrong-charge
sectors without distinguishing its internal columns. The joined history
costs 9,846 pulses and retains both pointer outputs and all 34,336
first-export Kraus slots. It uses eight initialized field-register bits
and one neutral pointer; no free spectral projector is inserted.

Each displayed complete history is checked with Stinespring Frobenius error
at most `1e-7`. A separate bound `16e-7` composes the sixteen local complete
FSWAP admissions for arbitrary route inputs and references; the complete
global success matrix is also checked directly. These are fixed numerical
implementation admissions, with no hardware noise or passive lifetime claim.

## Number-conserving charged spectroscopy

Let `a` be any dimer mode and `r` a supplied reservoir mode. With
`n=a†a` and `n_r=r†r`, define

```text
F_theta = I-n-n_r + exp(i theta) a†r + exp(-i theta) r†a.
```

The complete CAR identities give `F_theta†=F_theta`, `F_theta²=I` and
`[F_theta,N_field+n_r]=0`. In particular, its doubly occupied sector has
sign minus; an ordinary swap is not this operator. Evolve only the field
by the unchanged `K`. For a product reservoir number state, the correlator
`C_theta=<F_theta(t)F_0>` is

```text
reservoir 1: C_theta = <n(t)n> + exp(-i theta)<a(t)a†>,
reservoir 0: C_theta = <(I-n(t))(I-n)> + exp(i theta)<a†(t)a>.
```

These are operator identities on the complete field algebra, so the field
may be correlated with any untouched reference and need not be stationary.
The phase difference `(C_0-C_pi)/2` isolates the addition Green function;
the empty-reservoir version, conjugated, gives the existing removal
convention. This is linear record processing, without postselection.

`detector.py` executes both phases, both reservoir occupations and both
pointer quadratures: eight complete settings on six working bits. Every
setting includes a fresh source ground preparation, forward and reverse
split continuation, controlled **even** reflection probes, both pointer
outputs and every preparation/continuation export-failure branch. The
reservoir remains in each output state. All settings together cost 58,464
native pulses, with 58,464 loads and 58,464 exports; they require eight
field/reservoir/pointer initializations and four supplied reservoir particles.
At model duration `1/4` with four split steps, the ideal per-setting operator
error is at most `11/32768`. Full source vectors have a separate fixed
floating admission; these numbers are not hardware error certifications.

The measured addition correlation is approximately
`0.4960382759-0.0623579372 i`; the independent complete Lehmann value is
`0.4960382759-0.0623585494 i`. Removal gives its conjugate. All poles are
retained, including the same lower candidate residue `(1+1/sqrt(17))/4`.
Reservoir/bilinear Green measurements are established methods; see
[Bishop, Bagrets and Wilhelm](https://doi.org/10.1103/PhysRevA.111.062610).
The result here is their explicit number-preserving source realization on
the already specified OPH candidate/reference, rather than a new physical
particle identification or a general measurement-method claim.

## What is localized and what is protected

At a ready boundary, one charged dimer with ground blocks elsewhere has
strict finite support relative to the product reference. The supplied block
permutation transports that support exactly. This construction does not
show passive localization after turning on interblock couplings, a bulk
continuum limit, or protection under every admitted measurement and repair.

A direct counterexample uses only a complete, number-preserving local
occupation read. For the normalized lower-doublet projector `rho=P_lower/2`,
the channel `rho -> n rho n + (I-n)rho(I-n)` retains only `3/4` of its weight
in the lower doublet. Charge and parity survive, while the spectral candidate
does not. The package checks this for holes and particles on all axes.
Consequently neither CAR parity nor an exchange sign supplies passive
protection of this particular excitation. The H66 memory theorem has a
different encoding and explicit active recovery; it cannot be inherited
merely because both constructions use source operations.

Physical electromagnetic charge, Lorentz spin, autonomous selection of the
action, a protected defect sector and a common physical rest/propagation
clock remain additional requirements. The conditional electron/positron
calibration is unchanged; equal gaps are not new absolute mass predictions.

## Reproduction and evidence scope

```sh
python -m pytest -q -W error model/particle_transport
python -m model.particle_transport.packet --check model/particle_transport/results.json
```

The compact retained receipt binds the contract, implementation, tests,
import closure and actual source data. Verification executes the algebra,
source histories and detectors afresh and compares every retained key and
value; it rejects nonfinite numbers, altered types, missing/extra fields,
stale loaded code and changed source inputs. Exact integer identities,
analytic finite formulas and floating implementation admissions are labeled
separately. This is CP-05/CP-13 component progress. No physical comparison
or integrated completion gate is added.
