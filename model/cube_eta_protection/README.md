# Three-dimensional eta protection and local interacting response

This packet gives an exact protected sector of the **original open2x2x2
Weyl/contact field**, including all three hopping axes, with a finite
four-dimer native repair service. It also prepares a fixed-number reference
with the actual source and gives a one-site disturbance with eight
interacting spectral lines. No hidden-spin alignment of the dimers is used
for protection.

Run from the repository root:

```sh
.venv/bin/python -W error -m pytest model/cube_eta_protection -q
.venv/bin/python -W error -m model.cube_eta_protection --check model/cube_eta_protection/results.json
```

## Exact full-field protection

The field has eight sites, sixteen CAR modes, c=1/8, U=1, and the full
65536-dimensional Fock space. `cube.py` consumes the integer Gaussian
matrices of `field_symmetry`, with integer H16=16H. In the **physical**
mode basis define

`eta+ = sum_x (-1)^(x0+x1+x2) c_up(x)^dagger c_down(x)^dagger`,
`eta-=(eta+)^dagger`, and `eta_z=(N-8)/2`.

For the supplied i sigma_i links, onsite singlets are invariant under the
SU2 link rotation; opposite sublattices cancel each hopping contribution.
Consequently `[H,eta+]=U eta+` and `[H-(U/2)N,eta+]=0`. The packet checks
the latter on the **complete**65536-dimensional action, not a compressed
candidate or a spin-aligned replacement.

Let X=eta+ from N4 to N6. The exact Gram identity

`X†X = (eta-|N4)†(eta-|N4) + 4I`

proves injectivity and rank C(16,4)=1820. Also
`H_N6 X = X(H_N4+U)`. Thus the polar isometry
`V=X(X†X)^(-1/2)` gives an exact normalized intertwiner. Its range is the
entire eta>=2 part of N6; eta_z=-1 there. This contains the full interacting
four-particle parent action, not only free Slater energies.

The repair domain consists of four x-directed dimers `(0,4),(1,5),(2,6),(3,7)`
in canonical site order. Blocks use existing native ports `(0,1,5,7)`.
All five inherited seams `(0,1),(0,5),(0,7),(1,5),(1,7)` are admitted.
Exactly two blocks carry the existing lower-hole doublet and two carry
the unchanged dimer ground, giving a24-dimensional N6 support P. The
unnormalized ground is

`-i|up,up> + i|down,down> + (4-sqrt17)(|double,empty>+|empty,double>)`,

and the two hole columns are `|up,empty>+i|empty,down>` and
`|down,empty>+i|empty,up>`. `repair.py` expands the four-block columns
exactly as A+sqrt17 B, with integer Gaussian A,B, and applies the actual
fermionic permutation to the physical cube ordering. Both `eta- A=0`
and `eta- B=0` hold exactly; hence P has eta1. A numerical check binds
these algebraic columns to the **consumed analytic normalized** ground
and hole columns; this separately admitted check is not the exact proof.
The normalization itself is exact: with r=4-sqrt17, the consumed analytic
ground formula has b/a=r, a>0 and a²(1+r²)=1. Two hole and two ground
factors give the common column normalization n=4(1+r²)=136-32sqrt17.
The integer Gram identities are `A†A+17B†B=35904I` and
`A†B+B†A=-8704I`, hence `(A+sqrt17B)†(A+sqrt17B)=n²I`.
The floating binding residual is therefore a diagnostic of the code's
evaluation, not an unproved mathematical normalization premise.

Every coin/raw slot of C and E is read from the original twelve-port
native rule, restricted to these four blocks. The resulting matrices are
checked exactly against the physical graded whole-block permutation.
For each seam/law they obey completeness and unitality, including zero
slots. Embedded native Kraus maps are `V_P M_r V_P†/sqrt2`; the outside
outcome is `Q=I-P` on the full Fock space. There are130 outcome slots
over all ten services:120 native and10 outside. Since `P X=0`, every
native outcome annihilates range X and Q acts identically. The complete
repair instrument is the identity there with a deterministic outside
record, including all operators entangled with an arbitrary reference.
Interleaved field evolution is therefore the unitary interacting field
evolution on this sector, with every repair record retained.

The orthogonality argument applies to any specified dimer covering whose
local ground is eta0 and whose two chosen hole columns are eta-lowest;
this packet executes the stated x-oriented covering and all its five
seams only. It does not assert arbitrary repair invariance, maximality,
or a65536-state representation of all twelve original blocks.

## Exact protected spectrum and one-site disturbance

On N2 let `v=eta+|vac>`, so `v†v=8`, and `Q3=I-vv†/8` has rank119.
For `T=(eta+)²:N2->N6`, the exact identities are

`T†T=60I+24vv†`, and `H_N6 T=T(H_N2+2U)`.

Thus `T/sqrt60` is isometric on Q3. Its eta3 image has35 hidden-spin
singlets and84 triplet states. Hidden spin is used to calculate this
spectrum, not to infer repair protection or Lorentz spin.
In the actual Hubbard onsite frame the singlet basis is the eight
doublons plus28 unnormalized antisymmetric spin pairs; its metric is
diag(1 on doublons,2 on separated pairs). The integer action preserves
this basis. Exact characteristic polynomials give, for `z=16(E_N6-2)`,

`z^57 (z²-12)^21 (z²-16z-12)
 (z³-16z²-12z+64)^3 (z³-16z²-12z+128)^3`.

This is119 states at11 distinct energies:57 at z0,21 at each of ±sqrt12,
one at each quadratic root, and three at each of six cubic roots.
The last20 dimensions have contact-shifted energies. Exact coprimality,
squarefreeness and rational real-root isolation keep this conclusion
separate from the displayed approximate weights.

Use the fixedN6 reference `|R>=(eta+)³|vac>/sqrt2016`, with H energy3.
It has56 signed equal occupation amplitudes, one for each three-doublon
configuration. At one physical site0 let n=n_up+n_down. Then

`<n>=3/4`, `Var(n)=15/16`, and `|q>=(4n-3)|R>/sqrt15`

is a normalized eta3 singlet. Its parent is `(8d_0-v)/sqrt56` and the
packet checks the exact lifting identity. Its mean H is3 and variance
is3/112. The eight-column parent Krylov Gram determinant is exactly

`1488395072912184839188800532906554156908544 > 0`.

Its annihilator is exactly the product of the quadratic and both cubics.
Therefore this particular local probe has **positive spectral weight
at all eight interacting energies**, including the degenerate roots,
and zero weight at the three free-energy values. The rational resolvent,
root-isolating intervals, and normalized weights are retained.

The one-site number-preserving phase `F=exp(i*pi*n/2)` obeys
`F|R>=(|R>-sqrt15|q>)/4`. Its elastic/inelastic weights are1/16 and15/16.
The full seven-site exterior reduced state is exactly unchanged by F;
the integer coefficient-matrix Gram equality verifies this on the
entire support. Globally, the pure-state trace distance is sqrt15/4.
This is a **finite local unitary disturbance of an entangled reference**,
not a claim of boundary charge or a localized asymptotic particle.
The disturbed state remains in the protected sector and has H variance
45/1792, so its time response is nontrivial. A fixedN6 Ramsey construction
with a pointer, controlled F, common field evolution, and controlled F†
has correlation

`C(t)=1/16+(15/16) sum_j w_j exp[-it(z_j/16-1)]`.

H=3 for the reference fixes the clock: replacing H by K without shifting
the reference changes the answer. No vacuum/charged-state coherence is
needed during this spectroscopy. No physical vacuum or mass is selected.

## Actual source preparation; explicit detector boundary

`source.py` compiles the signed56-amplitude Dicke state on eight even
occupation bits, while all eight odd modes remain blank. All5379 native
tree pulses and all their exports are evaluated. Relabelling compact
bit j as physical bit2j with blank spectators is an exact isometry on
every retained input/output, not a success projection. Eight literal
native CNOT words then copy into the odd modes in24 pulses, executed on
all65536 rows. The one-site F word has another3 pulses, including its
scalar phase. No charge postselection occurs.

The independent exact trig/port certificate checks all255 tree splits,
all fixed copy/phase angles and the actual represented source controls.
The complete source-state error is at most the retained rational bound,
less than8.48e-13, including all export failures. The numerical full-source
diagnostic is separately admitted. All first-failure records are kept,
including true zero outcomes and explicitly identified table padding.
All130 following repair outcomes are also retained; outside records
continue to the next service without being counted twice as terminal
probability. Source failures and native repair branches stop; the final
outside continuation completes the probability sum.

The complete native seventeen-bit time-dependent detector and a native
field-continuation approximation are **not executed in this packet**.
The local probe spectrum and Ramsey formula are exact analytic results;
the source execution covers the reference preparation and local kick.
There is no hardware-accuracy claim. The old small-carrier leakage and
transverse hidden-spin frame mismatch remain valid: eta protection is a
different larger-sector mechanism. Growing-volume control, transported
deformations, calibrated species/masses and interacting continuum QFT
remain open. The receipt binds all consumed laws, focused tests and the
paper fragment and replays in a fresh interpreter.
