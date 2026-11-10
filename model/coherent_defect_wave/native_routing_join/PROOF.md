# One finite clock for routed flux and native interaction

This is the specified architecture `coherent-native-routed-paired-flux-v1`.
It joins the [complete flux route](../mobile_flux/ROUTING_COMPLETION.md) to
the actual [mobile-source occupation unitary](../mobile_flux/PROOF.md).
The native update acts on the same link word that the route transports.
The construction changes the clock explicitly. It proves a finite joint
interaction and retains the free one-object invariant sector; it makes no
interacting relativistic limit or Standard Model identification.

## Source, code, reference and observation domain

Each spatial cell has two link lanes. Each lane contains the twelve literal
link words, neutral routing memory and an independent native CAR register
with the same 7,200 modes as the occupation source. Independence matters:
using one shared native register in both lane Gauss operators would require
a further compatibility proof for each native update.

Write a link word as `r xor sigma*STAR`, with r the canonical representative.
For native occupation n let `chi(n)=(-1)^(number of port-5 doublet modes)`.
The exact logical encoding of each lane is

\[
 E|r,b;n\rangle=
 \frac{|r,b;n\rangle+\chi(n)|r\mathbin{\mathrm{xor}}STAR,b;n\rangle}{\sqrt2}.
\]

For the flat orbit `r=0`, routing memory b is empty. For the five odd orbits
b has two values. The logical space is therefore `C11 tensor F_native`.
This is an isometry into the positive interior Gauss sector. It retains
both gauge partners and every native occupation, with exactly two link
registers and two native registers per cell. The same link register occurs
in E, in routing and in the native update; no frozen or cloned field is used.

The common reference is the flat-orbit, empty-routing-memory, native-empty
product. Every gate below fixes it. The normalization is defined in the
logical code and preserved by E. A specified preparation can implement E
deterministically by a Hadamard on the gauge-pair bit followed by the
occupation-dependent sign. This is a coherent preparation service, not a
projection with its failure removed. Alternatively, starting from definite
literal words, all products of Gauss projectors form a complete instrument;
on L lanes its `2^L` signs each have probability `2^-L`. The all-positive
branch has that success probability and all other signs must be retained.

The admitted detector includes all alternatives of the face-syndrome
readout and the regional pair-position readout, together with specified
native occupation projectors. These commute with the interior constraints.
No missing boundary Gauss law or continuous gauge field is supplied by
these interior constraints.

## Exact local native action and complete routing

The literal native occupation unitary `U_F` preserves regional charge Q,
native number and the interior Gauss constraint. It therefore preserves E's
range: flat words retain empty routing memory, and odd words retain their
occupied routing memory while visiting any of the five allowed face orbits.
The compressed unitary `V=E* U_F E` is an exact restriction. In
`native_column`, both source gauge partners execute the actual U or U*;
decoding is accepted only when re-encoding reproduces the entire output.
Agreement of a projected matrix alone would be insufficient.

The routing rule acts on the two C11 factors while leaving both native
registers at their cells. Its twenty-dimensional one-object coin, vacuum
identity and one-hundred-dimensional coincident-object identity are all
retained. The existing two local swap layers implement each countertranslation.
Thus routing is a complete finite unitary for every object occupation, with
the native register as a retained spectator during each routing operation.
Lifting a logical gate W as `E W E*+(I-E E*)` gives a unitary on the full
physical support that commutes with each admitted Gauss constraint. This
dressing is necessary when the fixed native occupations change the Gauss sign.

Let `R_a` denote the complete three-direction routing tick, including its
face-count mass coin `exp(-2ia beta)` in the one-object sector. Define

\[
 K_\ell=I-\frac{V_\ell+V_\ell^*}{2},\qquad
 K=\sum_\ell K_\ell,\qquad
 J_{a,\lambda}=\exp(-ia\lambda K)R_a.
\]

This order, spacing/tick duration a and real coupling lambda define the
common clock. Distinct K_l act on independent lanes and commute. Each K_l
is positive with norm at most two and vanishes on the native-empty sector.
Since `K_l<=2 P_(native occupied,l)<=2 N_l`, the fixed total native-number
sector satisfies `0<=K<=2N`. Each individual native number N_l is conserved
by the full clock; routing moves the flux while native registers stay local.
Total paired-flux number and all interior constraints are also conserved.

The weak exponential is a declared local Hamiltonian service. It is neither
the literal U_F tick nor an energy logarithm inferred from its phase bins.
The contact spectrum of the separate positive-Hamiltonian construction is
not transferred to this routed clock by renaming an energy.

## Controlled finite calculation

For `epsilon=a lambda`, compute

\[
 T_r=\sum_{j=0}^r\frac{(-i\epsilon K)^j}{j!}.
\]

Taylor's integral remainder uses the exact unitary `exp(-itK)` and gives

\[
 \|\exp(-i\epsilon K)-T_r\|
 \le \frac{(2N|\epsilon|)^{r+1}}{(r+1)!}=\delta_r.
\]

This is an operator-norm bound on the fixed-N sector. On normalized input
it is also a vector error; no exponential prefactor is required because K
is self-adjoint. For arbitrary input its vector bound is multiplied by
the input norm. The exact joined tick remains unitary. The Taylor vector
is neither renormalized nor passed off as its exact state.

For any detector projection P, if `b=||P T_r R_a psi||`, its exact probability
lies in `[(max(0,b-delta_r))^2,(b+delta_r)^2]`, intersected with `[0,1]`.
For a binary detector the complementary interval is obtained by subtracting
the first interval from one. Exact complete probabilities sum to one;
approximate squared norms need not do so.

The executed side-three cube places one paired flux at `(1,1,1)` and one
native wave, in mode `(port7,internal3,repair0,address0)`, at lane0 of
`(2,2,2)`. Its spacing is `a=pi/12` and coupling `lambda=3/(2500pi)`, hence
`epsilon=1/10000`. The normalized routed state has sixteen amplitudes.
The actual encoded generator produces 1,440 amplitudes and satisfies

\[
 \|K R_a\psi\|^2=\frac{283}{72},\qquad
 \|P_{\rm changed\ face}K R_a\psi\|^2=\frac1{19200}.
\]

P detects any of the four odd face orbits different from the prepared
orbit, retaining its complement. Routing alone leaves that face label
unchanged. The native gate acts on the routed object's actual link word.
At order one, `delta_1=1/50000000`. Since `138<sqrt(19200)<139`, its exact
yes probability lies between

\[
 \left(\frac1{1390000}-\frac1{50000000}\right)^2
 \quad\hbox{and}\quad
 \left(\frac1{1380000}+\frac1{50000000}\right)^2.
\]

The lower endpoint is strictly positive. Both records are bounded and
retained. The first-order vector has squared norm
`7200000283/7200000000`, and its approximate yes weight is
`1/1920000000000`; the receipt labels these approximations explicitly.
No normalization step converts them into claimed exact output probabilities.
With the native register empty, or with lambda zero, the exact face-change
probability is zero. This demonstrates native influence on an actual
spatially routed protected object under the declared common clock.

## Exact free restriction and physical boundaries

If every native register is empty, K vanishes identically. That sector is
invariant under both factors of J, so the full finite joined clock equals
R_a there, for every object occupation and every lambda. Its one-object
restriction has the proved free window
`omega=+/-sqrt(|p|^2+4)` with error `ta(3P^2+6P)` on `|p_i|<=P`.
The whole many-object completion is not thereby a free fermionic field.

The same invariance has a direct mass consequence: this native coupling
does not generate a vacuum self-energy or Yukawa mass correction in the
native-empty free sector. Its rest-frequency magnitude two persists
exactly there. A mass-renormalization claim requires a different populated
reference or a coupling that can access additional sectors.

Total paired-flux number conservation also prevents object production from
an initially empty routing vacuum. The demonstrated object is inserted by
the specified preparation. Its finite response is not a calculation of
early-universe particle production.

The joined construction supplies no interacting relativistic limit.
The native registers have no covariant spatial propagation; their numbers
are fixed at their cells. The coupling's Lorentz transformation law is
unproved. Rotationally scalar density coupling is not automatically a
Lorentz scalar, and beta insertion alone does not provide a covariant field
for its coefficient. Native action can mix momentum windows, so the free
fixed-window estimate does not bound full interacting histories without a
window/leakage control. Positive many-body energy, physical statistics,
chiral currents, a scattering interpretation and an independent conversion
to GeV require their own constructions. The five face locations remain
configurations of the same source object, without five-species identification.
Nontrivial continuum response also requires a specified native-state density
and coupling scaling; a localized finite native preparation can lose its
physical support under refinement.
