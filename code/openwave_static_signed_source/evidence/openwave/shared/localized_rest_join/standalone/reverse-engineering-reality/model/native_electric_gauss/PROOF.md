# Integer charge and electric records on native observer links

## Carrier and architectural resources

The source is the native two-load space on twelve ports, with the thirty
physical seams and the active six-state memory from
`model/openwave_active_native_charge`. Each unordered occupied pair receives a
transported tag. An integer register (q\in\mathbb Z) assigns charge (q) to
the tagged load and (-q) to the other load. This retains the compensating
reference explicitly. The active memory character is the declared quotient
(q\bmod6) of the distinguished-load/reference label. It is not the total
electric charge, which is zero in this preparation. The real-valued charge
is (q_*q), with supplied unit (q_*). No nonzero additive homomorphism from
$\mathbb Z_6$ to $\mathbb R$ exists: (6f(1)=0) implies (f(1)=0).

Each oriented physical seam (a<b) carries an integer electric flux
(E_{ab}\in\mathbb Z). The source, memory and these thirty links form one
extended carrier. The tag and integer links are added architectural registers;
they are not inferred from the source's six memory characters. Controlled
integer addition, a neutral-pair preparation swap, an electric diagonal
phase and complete orthogonal pointer readouts are declared services. The
integer translations and their finite-window permutations are executed on
complete finite local bases. This does not supply a gate decomposition of
unbounded arithmetic into the source's smaller native alphabet.

Let (B) be the oriented incidence matrix, with (+1) at an edge's origin.
For a basis state (s), define

\[
 \rho_v(s)=q\,1_{v=\mathrm{tag}(s)}-q\,1_{v=\mathrm{other}(s)},\qquad
 G_v(s)=(BE)_v-\rho_v(s).
\]

The integer spectrum defines the continuous local action
$\exp(i\sum_v\theta_vG_v)$. The physical Gauss sector is (G=0).
The operations below preserve every Gauss sector, including nonzero ones.
There is no physical electron assignment, SI unit, continuum field or
independently calibrated physical clock in this theorem.

## Exact native transport and complete records

Consider an actual native collision-free hop along (a<b). Its number
current is (j=+1) for (a\to b), and (-1) for the reverse. Let (q_h)
be the charge on the moving load. Apply the existing native source
permutation and its active-memory voltage, transport the tag with its load,
and set

\[
 J=q_hj,\qquad E_{ab}'=E_{ab}-J.
\]

All other links are unchanged. Since $\Delta\rho=-JB_{ab}$ and
(B\Delta E=-JB_{ab}), every (G_v) is preserved. Collision and coherent-event-bit-zero
columns are the identity. The coherent bit is the recoded identity/swap
selector. The original C/E coin chooses the final endpoint; it is not this
permutation selector. A committed raw branch selects the swap from the
absolute value of its decoded number current. The reverse native hop has opposite current and
inverse memory voltage. Each fixed-seam event is therefore a full basis
permutation, and consequently a unitary on all matrix units and all quantum
references. The tag is essential: reassigning signs by sorted particle order
would fail to transport charge through an ordering exchange.

For a finite register (-L\le E_e\le L), execute an admitted transposition
only when its target is inside the window. An outward proposal is the
identity on source, memory and flux. The inverse partner of an admitted
column is inside the same window, so this complete map is an involution.
No flux is wrapped modulo the register dimension. Its proposed current and
realized current are separate records. This boundary policy differs from
the unbounded native extension exactly on flagged proposals. A word of
(n) events agrees with the unbounded extension when every input link has
margin at least (n|q|).

For committed records, retain every original C/E raw slot, including zero
slots, and split each slot into admitted and overflow outcomes. Write
(D_r) for its original diagonal source-domain projector and (P_b) for
its boundary flag. The extended Kraus operator is

\[
 K_{r,b}=\sqrt{w_r}\,S_{|j(r)|}P_bD_r.
\]

Each (S) is the complete permutation proved above. Thus
(K_{r,b}^*K_{r,b}=w_rP_bD_r), and summing both flags and all twenty-five
raw slots gives the identity. The original source-domain weights sum to
one; no successful output is normalized. These partial injections specify
the entire matrix-unit instrument, not only its populations. All other
links tensor the identity. Nonadmitted encoded source/memory addresses are
the identity complement of the program's explicit admission flag; they
are not accepted as physical pair states.

The finite certificate executes every pair, tag, active-memory coordinate,
link value and coherent event bit for (q=-1,0,1), (L=1), on all thirty seams. It also
executes (q=-2,-1,0,1,2), (L=2), on seam ((0,1)). It checks full
permutation uniqueness, involution, Gauss preservation, all native raw
slots and the exact electric-work identity. Arbitrary integer charges,
spectator links and wider windows follow from the algebra above. Removing
the flux translation supplies a nonzero Gauss-defect control.

## Neutral preparation and electrostatic work

Choose a simple native path (x=v_0,\ldots,v_d=y). Prepare (+q) at (x),
(-q) at (y), and electric flux (q) along the oriented path. Its
divergence is exactly the charge density. A source/reference-controlled
swap between this vector and its retained vacuum address, with the identity
on every other address, is a complete Gauss-preserving preparation. The
memory passes through; a specified character preparation includes all six
memory coordinates. A fixed preparation target is part of the program,
so distinct paths or (q)'s are distinct preparation settings.

For supplied (kappa\ge0), define

\[
 H_E=\frac\kappa2\sum_eE_e^2.
\]

The preparation's electric work is (kappa q^2d/2). Its source-sector and
apparatus work require their Hamiltonians and are not set to zero. A native
hop has the exact electric-energy change

\[
 \Delta H_E=\kappa(J^2/2-E_{ab}J).
\]

This holds for every admitted column; a boundary identity has realized
(J=0) and zero change.

For a unit neutral pair, every integer solution (BE=\rho) satisfies
(sum E_e^2\ge\sum|E_e|\ge d(x,y)). The second inequality follows by
pairing (BE) with a unit-Lipschitz graph-distance function. A shortest
integer string attains the bound, proving the exact integer electric ground
energy (kappa d(x,y)/2). On the same graph, the *real-valued relaxation*
has minimizer (E=B^TL^+\rho), (L=BB^T), and energy
(kappa\rho^TL^+\rho/2). These are distinct feasible sets. At
(kappa=1/8), an adjacent pair has integer minimum (1/16), while the
Poisson lower bound is (11/480), with strict excess (19/480). The receipt
solves all eleven pairs with origin zero exactly. Neither graph energy
proves a Euclidean (1/r) law.
For a normalized quantum state supported on this fixed integer-charge sector,
the diagonal energy operator obeys the same lower bound: its expectation is
an average of integer-string energies and is at least kappa d(x,y)/2.
The energy of the mean flux is a different nonlinear quantity and may lie
below that value. A real Poisson relaxation cannot replace the quantum
expectation of the integer quadratic energy.


## Operational integer boundary charge and the half-split control

For every native region (R), measure the oriented boundary flux

\[
 Q_R(E)=\sum_{e=(a,b)}E_e(1_{a\in R}-1_{b\in R}).
\]

The Gauss constraint gives (Q_R=\sum_{v\in R}\rho_v). This readout
uses the actual link records rather than the integer label as its input.
It has an integer spectrum. Adding an integer cycle (c), (Bc=0), leaves
every (Q_R) unchanged. Deforming a region without crossing a charged
endpoint also preserves its result. This is an integer relative-chain
boundary index, with localized endpoints transported by the native events.
It does not assert a continuum texture or winding classification.

The executable checks all 4,096 regions on eight actual prepared, looped or
transported source states. Four successive native loop deformations have
identical full cut tables. For the unit neutral pair each table contains
1,024 outcomes with charge (+1), 1,024 with (-1), and 2,048 with zero.
No pair of disjoint-region sharp pointer records can equal
((1/2,1/2)), because both pointer spectra are integer. This is a statement
about the complete readout; it does not forbid fractional expectation
values. A supplied Hadamard reference controls the actual native swap on
((0,5)). The resulting two regional expectations are both (1/2), while
the complete joint charge records are exactly ((1,0)) and ((0,1)), each
with probability (1/2). All twelve reference/memory outcomes are retained.
The half-half sharp record has probability zero.

This proves the finite integer-flux and unsplittable-unit portions of the
charge criterion on the declared extension. Architectural uniqueness is
not a premise. Equivalence to the pinned platform construction, or explicit
review of the declared extension, and identification of its endpoint defect
with the physical electron are separate requirements. No platform icon
change follows from this finite proof alone.

## Native loop, field quadratures and complete observations

Prepare the pair ((0,1)), tagged port zero, (q=1), the direct electric
string on ((0,1)), and the full active character (k=1). The three literal
recoded coherent native swap events

\[
 (0,5),\quad(5,11),\quad(0,11)
\]

return the source pair and its tag. They change memory by the inherited
voltage (2\bmod6), and change flux by the divergence-free cycle (c),
with components (-1) on ((0,5),(5,11)), and (+1) on ((0,11)).
Let (U) denote this actual three-event permutation. Its inverse is the
reversed event word. On the unbounded fixed-matter orbit,

\[
 U|n\rangle=|n+1\rangle\quad\text{with the active-memory action retained},
 \qquad H_E(n)=\frac{1+3n^2}{16}.
\]

The finite executor uses (L=2) for the two-block experiment, retains all
six memory coordinates, and executes every branch of both blocks. A
Hadamard reference, controlled (U), sixth-root phase (zeta^r), and
Hadamard readout give the complete Kraus pair

\[
 M_\pm^{(r)}=(I\pm\zeta^rU)/2,
 \qquad \sum_\pm M_\pm^{(r)*}M_\pm^{(r)}=I.
\]

A reference is supplied for each block. These blocks use the declared
coherent native-record policy; committing which-path source records is a
different instrument. The first phase is (r=1). The second phases (r=0)
and (r=1) produce, respectively, the exact joint probabilities

\[
 (5/16,3/16,3/16,5/16),\qquad(3/8,1/8,1/8,3/8).
\]

Every terminal source, link and memory probability appears in the receipt.
For each second-reference setting their sum is one. The nonselective
electric energies are (1/16,5/32,11/32), giving electric work (3/32)
and (3/16). The first outcomes each have probability (1/2) and
conditional loop expectation $\langle U\rangle_\pm=\pm\bar\zeta/2$.
The two measured quadratures determine its real and imaginary parts;
these expectations are not inferred from a selected-only experiment.

For a supplied coefficient $\beta\ge0$, define the finite loop generator
(H_B=\beta[I-(U+U^*)/2]) and its current
(J_B=i\beta(U-U^*)/2). This generator is a declared finite magnetic
interaction, without an autonomous source implementation claim. Its
conditional current is $\pm\beta\sqrt{3}/4$, with nonselective current
zero. The chosen loop uses seams (0,5), (5,11), (0,11); it does not
contain the prepared unit-string seam (0,1). Their electric vectors are
orthogonal. On this loop orbit the string adds the same constant electric
energy to every vacuum-orbit basis state. Consequently its vacuum-subtracted
pair ground energy is independent of beta, including with the declared loop
generator. The quadrature/current identities do not establish a beta-dependent
pair force or interpolation to the real Poisson minimum. Such a comparison
requires a specified cycle intersecting the string. The same source loop obeys

\[
 [E_e,U]=c_eU,\qquad
 [H_E,U]=U\kappa(c\cdot E+|c|^2/2),\qquad
 i[H_B,E_e]=c_eJ_B.
\]

These are finite Maxwell operator identities on the unbounded integer
loop orbit, and on finite input supports whose required words do not meet
a boundary. All generators preserve Gauss's law. At a boundary, the actual
complete permutation and its inverse are executed; a constant-cycle
identity is not imposed on a changed source orbit. Boundary controls
exercise and retain such modifications.

An independent electric phase readout uses a controlled forward loop,
the declared diagonal (e^{-i\tau H_E}), the inverse loop and complete
reference readout. At $\tau=16\pi/9$ source-time units, initial cycle
sectors (n=0,1,2) give energy differences (3/16,9/16,15/16) and pointer
probabilities ((3/4,1/4),(0,1),(3/4,1/4)). The calculation executes all
native forward/inverse columns in the prepared state and uses exact
relative sixth-root phases. The source and links return to their input
state; the field work of the forward and inverse controlled loops cancels.
The reference and phase service are consumed apparatus resources, without
an assertion that their Hamiltonian work vanishes. At $\kappa=0$ the
pointer is deterministically plus. This tests the supplied electric
coupling rather than deriving its normalization.

## Attachment boundary

This construction combines the actual native source/memory process with an
explicit integer electric extension and finite complete observations. It
has exact additive charge, local gauge invariance, charge-current continuity
and a measured native-loop response. The signed carrier, coupling, phase
clock and pointer are declared. A physical electric unit, source-derived
magnetic dynamics, continuum Maxwell/Coulomb limit, photon identification
and any OpenWave physical row validation are unproved.
