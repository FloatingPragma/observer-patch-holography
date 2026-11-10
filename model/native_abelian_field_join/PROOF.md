# Native odd charge coupled to the finite transverse source

This construction is a declared extension of thirteen independent copies
of the actual eight-cell native Hamiltonian. It retains every native
source state and adds an explicitly enlarged reference, a finite outer
abelian field and spectral matrix-unit hopping. It is not the uniform
larger native lattice, and does not inherit its Hamiltonian by a change
of notation. The main result below concerns a selected reducing charged
sector of the declared extension. A distinct penalty alternative is
specified separately.

## Native spectral input and the complete local reference

Let h be the actual full eight-cell native H_lambda at lambda=2^-34,
with physical neutral ground g and energy E_g. Its kinematic dimension is
`D=2^240 5^8 119^24`. Its two native e^c annihilators at the origin have
indices9 and24. The actual e^c internal representation has twelve
current values `(0,...,0,1)`, and its one-reference spin space has
exact dimension two. On that reference sector every gauge generator is
the charge value times I2 and the source Hamiltonian is h tensor I2.
The old reference energy4 is therefore an exact spectator shift.

The original complete thirty-mode CAR reference is denoted F_R. Retain
all its states and its reference Hamiltonian4N. Add the conjugate of
each nonvacuum number sector and two even boundary ports. Share the
single vacuum of the original and dual Fock carriers:

`R_loc=C|vac> direct_sum (sum_(n=1)^30 [F_R,n plus conjugate(F_R,n)])`
`      direct_sum C|port_+> direct_sum C|port_->`.

Its dimension is `2^31+1`. On the original occupied sector use its actual
gauge, spin and raw-lepton matrices and energy4n. On the dual sector use
the negative complex conjugates of each Hermitian current and the same
energy4n. Fermion parity is (-1)^n on both sectors. The ports have
native charges +1 and -1 respectively, spin zero, raw lepton number zero,
even parity and reference energy zero. This is a direct-sum reference
architecture; no global CAR algebra is asserted on the enlarged carrier.
The original full CAR carrier is an unchanged invariant subspace of the
reference Hamiltonian and currents. The joined transport below explicitly
couples its vacuum to the added ports.

An antiunitary A_R fixes the common vacuum, maps an original occupied
basis vector to its dual, and maps that dual back with phase(-1)^n.
It swaps the two ports. It satisfies, on the complete declared carrier,

`A_R H_ref A_R^-1=H_ref`,
`A_R G_ref A_R^-1=-G_ref`,
`A_R S_ref A_R^-1=-S_ref`,
`A_R L_ref A_R^-1=-L_ref`,
`A_R^2=(-1)^N_ref`.

The port parity is even in the last identity. The source antiunitary
Theta from the public charge/reflection package, together with A_R,
preserves the complete local Gauss kernel and h+H_ref. This is a new
full-reference proof. Applying particle-hole to the old reference would
instead change4N to120-4N and is not used here.

The original occupied-reference physical sector has the frozen gap
`12069208469156468531631400952645460636513962326353 /
 73287825008760668595250553296982590732615706214400 >0.16468`.
The dual sector has exactly the same spectrum by the just-defined
antiunitary. Both port sectors have gap at least
`a=144855983017581/3518437208883200 >0.0411705465` above E_g.
For port+ this is the complete old one-e^c source sector after removing
the two-state spectator and its energy4. The source charge/reflection
identity proves the port- bound. The common-vacuum sector is the original
physical neutral source. Hence the complete physical local block has
unique ground `g tensor vac`; neither charged complements nor duplicated
reference vacua are discarded.

## Actual odd G1 spectral matrices

The frozen full-H e^c bin has source excitation endpoints

`L=6370399497542706231/9223372036854775808`,
`U=12776828565198489711/18446744073709551616`.

Its residue is strictly positive. The native source ground has trivial
proper-cubic character, so at every energy the two origin annihilator
vectors have scalar Gram matrix. The actual antiunitary pairs their
addition and removal measures. Choose E_* in this bin with maximum
per-spin residue, breaking ties by the lower energy. The selection is a
finite spectral service on the actual source. Its weight is positive,
at least the frozen bin lower weight divided by D.

Normalize the two removal vectors and their two addition partners at
E_g+E_*. Attach removal to port+ and addition to port-. Denote the resulting
physical block vectors by e_(q,s), where q=-1 is removal and q=+1 is
addition, and s has two values. They are odd, orthonormal and have local
energy E_g+E_*. Their finite proper-rotation type is G1. After a fixed
pseudoreal basis identification, each pair has the same two-dimensional
G1 matrices. The physical source spin at arbitrary continuous angles is
not proved by this identification.

Define `b_(q,s)^dagger=|e_(q,s)><g,vac|`, extended by zero on the complete
local complement. These are actual native spectral matrix units. Their
odd grading is inherited from the full native Fock space. Replicated
blocks use its graded tensor product. On one block they obey the
exclusion algebra

`b_a b_b^dagger=delta_ab P_g`,
`b_b^dagger b_a=|e_b><e_a|`.

They are not canonical multi-flavor CAR operators on that five-state
code. At most one selected charged state occupies a block. The full
native complement remains present and reducing. No claim that the entire
native bin has rank two is required.

## Outer field, compensator and complete Gauss constraints

Use the public cone's exact13 vertices,42 edges,50 faces, D,C,G and
thirty transverse modes, with counting metrics and total photon cutoff
three. The field dimension is5456. Let the block's outer charge be
`Q_y=P_(port-)-P_(port+)`, on all states of those port sectors, and zero
on its other reference sectors. The reference port joins native charge
to this declared outer abelian charge; it does not identify all native
hypercharge with an empirical electromagnetic field.

At the fixed outer reference vertex r=1 retain a full52-mode CAR
compensator: thirteen flavors of a two-component spinor for each of the
two opposite charges. Its charge is Q_R=N_+-N_-, its raw lepton number
is -Q_R, and its energy is4(N_++N_-). Opposite charges carry conjugate
spinor matrices. The number-preserving antiunitary exchanges these
opposite-charge mode spaces, with square minus one on one-particle
states. Its exterior lift preserves4N_R and reverses Q_R and raw L on
the entire Fock space. This is a separate reference from R_loc.

For each y other than r, retain a ternary longitudinal register
`ell_y in {-1,0,+1}`. Define

`rho_ell=sum_(y!=r) ell_y (delta_y-delta_r)`,
`E_L=D G rho_ell`, `E=E_L-E_T`, `B=C A`,
`Gauss_v=(D^T E)_v-Q_v-Q_R delta_(v,r)`.

Every native block also retains all of its original local Gauss
constraints, now including its complete R_loc currents at the origin.
The outer physical kernel is exactly
`ell_y=Q_y (y!=r), Q_R=-sum_y Q_y`.
All transverse coordinates are divergence free on the complete cutoff
carrier. The finite Peierls matrices are literal exponentials
`W_e(q)=exp(i q A_e)`. No infinite-oscillator Weyl identity is used.

The complete physical space includes every matched old or dual local
reference occupation and every global compensator occupation satisfying
these equations. Extra Fock and longitudinal states remain in the
kinematic carrier; truncated shift operators vanish at both boundaries.
A compensated pump raises Q_y by q, raises the longitudinal register by
q when y!=r, and creates reference charge -q at r. A number-conserving
hop lowers the source and longitudinal labels at u and raises both at v;
it leaves the fixed compensator untouched. Thus every operation commutes
with every Gauss constraint, including on kinematic boundary inputs.

## Declared Hamiltonian and actual neutral ground

Subtract the common constant13E_g+E_vac when reporting excitation energies.
The unpenalized joined Hamiltonian is

`H_join=sum_y(h_y+H_ref,loc,y)+4N_R+E_C(ell)+E_vac+H_ph`
`       -eta sum_(e=(u,v),q,s) [b_(v,q,s)^dagger b_(u,q,s) W_e(q)+h.c.]`,

with eta=1/24 and the just-described longitudinal shifts implicit in each
hop. Each matrix unit is a finite spectral service; an efficient native
implementation is not asserted. This Hamiltonian explicitly adds outer
transport to the replicated native sources. The original intra-block
Hamiltonian and all its energies are retained.

Let N_active be the number of selected charged blocks. The graded
Cauchy-Schwarz form bound gives
`V>=-eta sum_y degree(y)P_active,y>=-(1/2)N_active`.
Since E_*>=L>1/2, the native block lower bounds, nonnegative reference
cost, Coulomb energy and photon energy prove that the complete physical
neutral ground is the product of the thirteen actual native grounds,
the common local reference vacua, the global compensator vacuum, zero
longitudinal registers and photon vacuum. It is unique for this declared
source/reference Hamiltonian. Additional zero-Hamiltonian preparation
apparatus, if retained, produces ground degeneracy in that enlarged
apparatus and must not be included in this uniqueness statement.

The product actual native ground is the same one used by the spectral
operators. It is not replaced by an atomic trial when defining H_join.
A finite preparation can use the inherited physical native preparation,
with its source copies, references, failures and fidelity error retained;
the quantitative service-level join is specified below. An efficient
native implementation remains a separate problem.

## Antiunitary symmetry on the entire declared source

Use the native charge/reflection antiunitary on each block, the complete
A_R above on each local reference, the number-preserving opposite-charge
exchange on the full global compensator, ell->-ell on every ternary
register, and ordinary coefficient conjugation on the real normal-mode
photon basis. This is one antiunitary on the complete tensor product,
formed from the product of its unitary parts and one overall conjugation.
It sends Q_y and Q_R to their negatives, A to A, transverse momentum to
its negative and hence E to -E. Every outer Gauss generator is negated;
the native Gauss generators have the already-proved reflected action.

All reference costs are invariant on every sector, including old occupied,
dual occupied and both port sectors. Coulomb and photon energies are
preserved. Native spectral matrix units are exchanged between opposite
charges with matching spin-unitary factors. Those factors cancel in the
spin-contracted hop, while antiunitarity sends W_e(q) to W_e(-q).
Thus the complete H_join is invariant, including its full complements.
The evolution parameter is reversed. No unitary physical-CP or physical
CPT inference is made. The optional penalty below is also invariant,
because its full port and selected projectors are exchanged in pairs.

## Exact charged-sector embedding

Fix charge q and one compensator flavor of charge -q, with precisely one
reference fermion. Restrict the source to exactly one selected e_(q,s),
all other blocks g, matched longitudinal registers and the full photon
carrier. This is a reducing sector. Its isometry records the charge
position, the native G1 index s and every one of the5456 photon states.
All original native complements remain elsewhere in the complete source.

In these coordinates the actual joined Hamiltonian is exactly

`H_join-(13E_g+E_vac) = [H_affine_mobile,N=1 +(E_*+2)I] tensor I_G1`,

with the charge-sign Peierls orientation. The shift is explicit: native
excitation E_* plus new reference energy4 replaces the affine source
cost1 plus reference cost1. No static charged line is transferred after
changing the Hamiltonian. The equality proves the full operator transfer
on the declared reducing sector, including the complete photon cutoff.
The negative-charge Peierls orientation is unitarily conjugate to the
positive one by photon parity(-1)^Nph; the initial photon vacuum and
position/occupation projectors are unchanged. Therefore their complete
position and photon-number records agree at the same positive time.
This last statement uses that explicit unitary equivalence, not the
antiunitary's reversal of the time parameter.

Every certified affine one-charge spectral interval shifts by E_*+2;
the photon and position projectors transfer exactly. Subtracting the
retained compensator energy4 leaves source energy E_*+E_affine-2.
The full70928 residual result therefore supplies the charged orbital and
photon ground within this selected sector, with its G1 degeneracy and
retained-reference multiplicities explicit. The spin-independent hopping
supports an internal finite-group doublet. It is not a Dirac dispersion,
a continuous physical spin representation or a relativistic mass shell.
The fixed-reference cone has finite spatial levels rather than a proven
translation-invariant momentum band.

This sector's lowest energy is not automatically the lowest energy in the
complete outer-charge sector. Inactive native charged-port states have
only the lower bound a, which is below E_*. That obstruction remains in
the unpenalized construction.

## Whole-bin coupling with an unnormalized native CAR signal

A second explicit coupling choice avoids the small residue of a single
chosen eigenvalue. Let J be the actual inherited native spectral bin in
each compensated charge sector. Its source excitation energies lie in
[L,U]. At the origin its Gram matrix is

`<g,c_s^dagger P_J c_t g>=w_J delta_st`,

where w_J is the unnormalized squared norm per spin and has the frozen
lower bound above0.41970415. The total squared norm of c_s g is1/2.
Therefore the conditional bin probability after the native family pump
is2w_J. It would be incorrect to multiply w_J by another1/2.

Define `psi_s=P_J c_s g/sqrt(w_J)` for removal and its exactly paired
addition vectors. Use `|psi_s,port><g,vac|` as the transport matrix units.
This changes the inter-block coupling from the eigenchannel choice above;
it leaves the actual native h unchanged. Every state in J, including all
states orthogonal to these two vectors, remains present. The subspace
with one charged block in J is reducing, although the two-vector span
alone need not be invariant under h. All charged states outside J remain
immobile under the added hopping and retain their full native evolution.
The same complete references, Gauss proof and antiunitary proof apply.

Put `E_c=(L+U)/2` and `epsilon=(U-L)/2`. For comparison only, replace h-E_g
on J by E_c times the identity. The difference of actual and comparison
Hamiltonians has norm at most epsilon in the whole one-bin-charge
sector. The comparison psi sector is exactly the full affine70928
operator shifted by E_c+2 and tensored with its G1 index. The other bin
vectors are retained as immobile states. Duhamel's unitary identity gives

`||(U_actual(t)-U_comparison(t)) P_one_bin||<=epsilon |t|`.

This bound includes every photon state and every bin state. It is not a
truncated perturbative evolution or a replacement of the actual h.

A literal native CAR pump at outer vertex y is the contraction
`B_q=(1/sqrt2)sum_s r_(-q,s)^dagger c_(q,s)` with the appropriate local
port creation and longitudinal shift. Here c_(-,s)=c_s and c_(+,s)=c_s^dagger;
the vacuum projections on the local and global references are retained.
Its full-carrier Halmos completion exists because its reference-vacuum
effect is N_ec/2 or(2-N_ec)/2. On the actual native ground its success
probability is exactly1/2. Its bin component has squared norm w_J and
normalized source spin vectors psi_s. The outside-bin component has
squared norm1/2-w_J and cannot move to a different outer vertex.

If A_xy(t) is a certified lower arrival amplitude in the affine source
conditional on its successful pump, then the actual native CAR arrival
probability, before native preparation factors, is bounded below by

`w_J [max(0,A_xy(t)-epsilon |t|)]^2`.

The affine universal value is A_xy(1/4)>=47/9216. Its sharper full finite
first-integral controls give larger positive amplitudes. They transfer
by this uniform comparison, with the actual native bin width retained.
The same statement applies to arrival jointly with nonzero transverse
occupation, using the corresponding affine amplitude. This gives a
large positive finite signal without a1/D eigenvalue-weight selection.

The full unpenalized charged-sector ground remains unidentified because
out-of-bin native charged states can have lower energy. A positive
finite motion record does not remove that complement. The whole-bin
coupling is spin independent and supplies repeated internal G1 content;
it does not supply a Dirac Hamiltonian or continuous physical spin.

## Physical native preparation and the complete arrival law

Prepare thirteen native source blocks using the frozen physical two-copy
protocol, retaining all thirteen partner sources, every spin reference,
all clock/controller registers and every filter outcome. The new local
charge references and global compensator begin in their common vacua,
with zero longitudinal registers and photon vacuum. No port is populated
on any native-preparation branch, so the new hopping is exactly zero
throughout this preparation. The original native Hamiltonians and clocks
are therefore the actual ones used by the inherited preparation.

Let Q be the probability that all thirteen original filters succeed.
Its lower bound is the original filter success lower bound to power13.
Their conditional two-source states have at least the inherited
ground-pair fidelity. Independence gives the tensor-product overlap
with all twenty-six actual native grounds at least that fidelity to
power13.

Add a separate binary spectral apparatus measuring
`P0=product_(all26 native copies) |g_y><g_y|`, including the thirteen
retained partner sources. The active blocks also have their neutral
local reference vacua. Let p be
its success probability conditional on the preceding successes. Keep
its failure as a separate record. Conditional success is exactly g^26,
and p has the product-fidelity lower bound just stated. The failed
branch is also Gauss physical. This is not an assertion that the earlier
approximate preparation already made an exact ground.

P0 is an additional finite spectral service with its own pointer and
control cost. For example, its controlled-pointer unitary is
`P0 tensor I+(I-P0) tensor X`; an interaction
`(pi/2)(I-P0) tensor(I-X)` produces it in one declared source-time unit.
It commutes with the free source Hamiltonian and all Gauss constraints.
The source also undergoes its ordinary free evolution during that unit.
The complete free Hamiltonian is H_join on the thirteen active blocks
plus the thirteen unchanged native partner Hamiltonians, with the
original zero Hamiltonians on spin references, counters and pointers.
Conditional success leaves the partner energy13E_g in the common
ground baseline. It is not silently removed or assigned zero energy.
No efficient native synthesis or zero-cost projection is claimed. A
finite approximate spectral filter can instead use the actual native
gap and controlled source evolution, with its residual preparation error
retained in every subsequent bound. That variant is a separate apparatus
calculation.

The inherited neutral spin-reference Hamiltonians are zero and their
full Fock spaces remain retained. Their native particle-hole action,
with the reflected site labels, preserves those Hamiltonians. Together
with conjugation on the clock/pointer registers and the source-copy
antiunitaries, it extends the free-clock symmetry to the full retained
preparation apparatus. It does not make the preparation schedule or its
initial state invariant, and its zero-energy multiplicities do not alter
the narrower source-ground uniqueness statement above.

After the literal native CAR pump and bin flag, let alpha be the actual
conditional arrival probability from the normalized bin state. The
complete six-record law of the whole-bin construction is

```
native preparation failure          1-Q
native preparation pass, P0 failure Q(1-p)
P0 pass, native CAR pump failure    Q p/2
CAR success, outside native bin     Q p(1/2-w_J)
bin success, arrival at x           Q p w_J alpha
bin success, nonarrival at x        Q p w_J(1-alpha)
```

The six probabilities sum exactly to one. The bin flag tests membership in J across all possible charge positions,
rather than membership at the original vertex alone. It commutes with the
joined evolution in the one-charge sector; outside-bin failure keeps its
entire source state, compensator and clock. The same time t is used on
arrival and nonarrival branches. Preparation and pump failure branches
also retain their autonomous free source evolution. The pointer flags
make these branches orthogonal rather than discarding their states.
Every branch satisfies all native and outer Gauss constraints.

The strictly positive unconditional lower bound is
`Q_lower p_lower w_J,lower [A_xy(t)-epsilon|t|]^2` whenever the bracket
is positive. The residue factor w_J is unnormalized; no additional1/2
or1/4 multiplies this expression. The affine amplitude already removes
its own separate1/4 pump probability. `structure.py` checks this
normalization and evaluates the exact rational lower bounds, including
arrival jointly with a nonvacuum transverse occupation. Energy and
arrival are separate detector settings. The exact eigenchannel
construction above additionally permits direct transfer of a selected
band energy-window law; the whole-bin source requires its stated
comparison error for energy or ground-state conclusions.

## Retained finite storage

Let `D=2^240 5^8 119^24`, `R=2^31+1`, and let M be the inherited exact
native agreement-filter counter dimension. One complete kinematic
storage bound for the six-record protocol is

`D^26 R^13 2^52 3^12 5456 2^416 M^13 2^17`.

It includes all twenty-six native copies, thirteen complete old/dual/port
references, the full global compensator, every ternary longitudinal
register, the full finite photon carrier, sixteen native spin-reference
modes per source copy, all thirteen filter counters, thirteen retained
filter flags, and the P0, CAR, bin and arrival flags. A separate binary
photon-occupation output multiplies this bound by two. This is a finite
carrier bound, not an efficient synthesized gate count. The algebraic
native preparation and spectral services remain declared apparatus
operations with their explicit source/reference domains.

## Explicit alternative with a charged-complement penalty

This alternative concerns the exact eigenchannel coupling, whose
selected eigenspaces reduce the native Hamiltonian. It does not assert
a full charged-ground theorem for the whole-bin coupling.

A different declared source may add
`2 sum_y(P_port,y-P_active,y)` to that eigenchannel H_join. This positive spectral penalty
acts on every inactive charged-port state, including any omitted partner
at the same native energy. It changes the Hamiltonian on that complement;
it preserves the neutral ground, all selected-band operators and the
operator embedding above. It is not included in H_join.

For net outer source charge+1 or-1, the global compensator costs at least4.
An inactive charged port then costs at least2+a, so such sectors lie above
6+a before any positive active-charge remainder. If all charged ports are
active, their number is odd. Three or more have lower energy
`4+3(L-1/2)`, which exceeds the selected-sector upper bound
`U+2+1.87398449270` by a positive margin. More active charges only raise
this form lower bound. Neutral block complements add positive energy
and remove their vertices from the active hopping graph; this is a
Dirichlet restriction of the full selected one-charge operator, so it
cannot lower that operator's ground energy. Extra compensator fermions
also add positive cost.

These comparisons identify the lowest net-charge energy with the selected
sector for the explicitly penalized alternative. They do not make it a
unique joint charged vector: the native G1 doublet and the retained static
reference flavors remain degenerate. A full charged-ground statement
must retain those multiplicities. The unpenalized construction retains
only the reducing-sector result.
