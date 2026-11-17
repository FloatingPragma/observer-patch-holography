# A stationary dressed charge kernel with finite current memory

Platform row11 remains untested. This packet establishes an exact finite
observer algebra and source-record kernel. The physical charged-lepton minima
and mass ratios are absent. The two registered N+ memory policies retain the
same bounded native/Gauss motion and differ in whether the ACTUAL overflow
workspace becomes an output record. Their distinction is essential.

## 1. The full finite state domain is six paths

Keep the original triangle seams (0,5),(5,11),(0,11), a negative compensator at1,
and a positive tag on0,5,11. The integer charge is q=1. Every active memory
m=0,...,5 is included. Fix E_(0,1)=1, all links outside the triangle and(0,1)
to zero, and each triangle link to[-L,L]. Impose the actual source Gauss law
G=BE-rho=0. This is a complete finite domain under the specified triangle
source, not a list of selected successful paths.

Put k=E_(0,11). The Gauss equations give precisely these states:

    position 0:  (E05,E5,11,E0,11)=(-k,-k,k),       -L<=k<=L;
    position 5:  (-k-1,-k,k),                      -L<=k<L;
    position11:  (-k-1,-k-1,k),                    -L<=k<L.    (1)

Let d be0,1,2 at the respective positions and define

    r=3k+d,
    c=m-2k-h(d) mod6,      h(0)=0, h(1)=h(2)=5.              (2)

The native active-memory voltage is5 on0->5,0 on5->11 and3 on11->0.
Reverse moves have opposite voltage modulo6. Thus every admitted move
preserves c and changes r by one. The six c values give six disjoint paths,
each with vertices r=-3L,...,3L and length n=6L+1. At the two endpoints the
outward move is the source's whole-event overflow hold. Neither the electric
link nor the source/tag/memory is wrapped or partly updated. Each of the
three seam operations is an involution on the complete finite domain.

The finite checks enumerate all Gauss-zero states independently of(1), then
test(2) against every actual source event. They find42 states at L=1 and78 at
L=2, including all six memory components. There are12 overflowing seam/selector
columns per window, two in each component. No boundary state is removed.

## 2. Source, reference and all transported flavors

For each tagged load use the already declared flavor a=3s+f with spin s=0,1
and family f=0,1,2. Both loads' flavors are retained, so there are36 flavor
pairs. In canonical CAR order the tagged basis vector is

    |S_(i,a)> = epsilon(i)|CAR mask(i,a),m_i,E_i>,
    epsilon(i)=+1 if tag<compensator, otherwise -1.            (3)

The actual CAR sign obeys epsilon(i) sign(hop)=epsilon(j). With the finite
whole-event hold included, the complete source unitary U_e is consequently
the unsigned path permutation in the tagged basis(3), tensor identity on
all36 flavor pairs. The code checks the literal pinned endpoint proposal,
both exact C/E CAR columns and realized bounded source on every state,
flavor pair, seam and selector. A native proposal at an overflow is retained
separately from the realized identity hold. Counting that proposal as the
realized hop would violate the source contract.

For the reference let R_a be the tagged CAR state of the original
pair_string((7,0,1),q=1,m=0), with the same flavor pair a. Each scheduled seam
fixes R_a, its phase, memory and all links. It never overflows. The admitted
source space is the direct sum of these reference states and every S_(i,a).
Any other tagged/full-Fock input has a separately recorded identity complement;
no signed-pair electric interpretation is assigned to it.

The relational matrix unit |S_(i,a)><R_a| preserves total charge and Gauss:
both electric strings have the same G=0. It raises the charge in the triangle
by one and lowers the reference-region charge by one. This is operational
integer charge with the source's supplied unit. It is not a physical electron
identification. All flavor blocks use the same source law and clock.

## 3. A finite sufficient ghost and a common stationary reference

The observer's ghost stores the CURRENT index i, hence its position, full
integer electric vector and active memory. It updates by the same disclosed
seam and recoded selector as the source. It does not measure which charge arm
is occupied. Source and ghost may be different in a control; the ghost rule
and comparator remain frozen to the native design.

The original CAR sign telescopes to epsilon(initial)epsilon(final), while
all active-memory and electric path information needed by this source is
already in i. Two histories with the same i therefore need no extra phase
tape for the comparator in the basis(3). This sufficiency uses the stated
phase-exact source contract. It does not say every possible coherent lift or
other memory law has no additional holonomy.

Use the finite quantum/classical working algebra

    A = direct_sum_i B(H_source),

with ghost i labelling its classical blocks. An emitted seam/selector/overflow
letter is an external outcome of the complete instrument. All such outcomes
are retained in the finite execution. The stationarity below is a property
of the WORKING algebra. An ever-growing one-sided external archive does not
have a stationary distribution on its finite lengths; that earlier statement
is unchanged. The complete stationary instrument instead supplies consistent
time-translation-invariant probabilities for finite output windows.

A single common positive reference is

    rho_0 = (1/36) sum_a |R_a><R_a| tensor (I_ghost/N),
    N=6n.                                                    (4)

Uniform preparation and the source/reference dressings are explicit apparatus
services. This is one reference across every flavor and memory block, not
separate choices of spectral origin for different families. For

    O_(i,a)=|S_(i,a)><R_a| tensor |i><i|_ghost,

its Gram form is

    omega(O_(i,a)† O_(j,b)) = delta_ij delta_ab/(36N).         (5)

The full source operation fixes each R_a, and the ghost operation is a
permutation of its finite state set. Hence rho_0 is stationary after summing
all instrument outcomes, under BOTH memory policies. Actual source overflow
is always0 on this reference; the ghost overflow can still be1.

## 4. Actual overflow and ghost overflow are different records

Choose a seam uniformly and a coherent recoded selector fairly. Selector0
is identity and has no overflow. Selector1 executes the full source U_e and
the ghost permutation S_e. Let B_e be the projector on SOURCE states whose
proposed event overflows, zero on every R_a. Such states are fixed by U_e, so
B_e commutes with U_e. The ghost analog is a classical predicate b_e(i).

For A-stationary-dressed-read, the source overflow workspace is read. At a
fixed ghost input i the selector1 source Kraus maps are

    (1/sqrt6) U_e(I-B_e),   actual overflow0,
    (1/sqrt6) U_e B_e,      actual overflow1,                  (6)

and the new ghost is S_e i. Its own current/overflow is emitted in both
source branches. Selector0 supplies I/sqrt6. Equations(6) are complete on
the entire source, including all references, boundary states and arbitrary
external references. They are the complete source operation with its actual
orthogonal overflow flag, not a postselected survivor.

If ghost i is an overflow state, the two arms of O_(i,a) have different actual
source flags: S_(i,a) flags1 and R_a flags0. Summing the flags removes that
cross block even though both source arms undergo identity holds. If i is not
an overflow state, both flags are0 and the cross block is transported.
Consequently the unconditioned instrument on the charged operator span is

    O_(i,a) -> sum_j K_ji O_(j,a),
    K=(1/6) sum_e [I+S_e(I-B_e^ghost)].                      (7)

Here B_e^ghost is the diagonal classical predicate. Each summand is positive:
on an admitted swapped pair it is I+Swap; on a nonoverflow fixed point it
is2I; on an overflow fixed point it is I. Thus K>=0 independently of the path
calculation. It is symmetric but substochastic; its lost cross-block amplitude
is actual recorded which-arm information, not lost probability.

For A-stationary-dressed-erase, the actual overflow workspace is coherently
uncomputed BEFORE it is read. This is possible as a declared reversible
service because B_e is invariant under U_e. The whole source then uses
U_e/sqrt6 at selector1; the same ghost overflow/current is still reported.
There is no classical actual-overflow flag in this different instrument.
Its charged kernel is

    P=(1/6) sum_e(I+S_e).                                    (8)

Coherent uncomputation is not deletion of a committed record and is not
asserted to be available without the declared service. The comparison packet
calculates both instruments and retains their distinct output alphabets.
Reading a source-derived current or another which-arm flag would likewise
change these charged kernels. Original native raw C/E labels in the column
certificate are not simultaneous classical outcomes of(6)-(8).

After summing labels, these finite working channels are self-adjoint in the
uniform block Hilbert-Schmidt inner product. Indeed each S_e is an involution,
and each U_e(I-B_e), U_e B_e is Hermitian. The stationary GNS two-point kernel
of(5), for either transfer T=K or P, is therefore

    C_((i,a),(j,b))(t)=omega(O_(i,a)† Phi^t(O_(j,b)))
                     =delta_ab (T^t)_ij/(36N).              (9)

The same result follows directly from(7)-(8), without an assumption that a
classical population kernel always equals a quantum charged kernel.

## 5. Exact spectra, reflection positivity and the cutoff gap

On each of the six native paths, an interior vertex has four held selector
columns and two neighboring columns. An endpoint has five held columns in P,
one of which reports actual overflow and is absent from the charged K.
Thus the complete matrices are exactly

    P=I-L_path/6,
    K=(4I+A_path)/6.                                        (10)

The implementation tests(10) on all entries of both full42/78-state matrices.
The standard discrete cosine and sine vectors give the exact eigenvalues

    P: 2/3+(1/3)cos(j pi/n),       j=0,...,n-1;
    K: 2/3+(1/3)cos(j pi/(n+1)),   j=1,...,n.                 (11)

Each list is repeated in all six memory components and all36 flavor blocks.
In particular both transfers are strictly positive and at least I/3. P has
one constant eigenvector in each memory component. Its charged zero-response
mode is accessible through the corresponding sum of the O_(i,a). The ghost
and full electric records have not eliminated it.

K has no eigenvalue one at a finite window. Its lowest attempt-log energy is

    Delta_L=-log[2/3+(1/3)cos(pi/(6L+2))]
            ~ pi²/[6(6L+2)²] as L->infinity.                (12)

It comes from the ACTUAL overflow readout and disappears as the finite window
is enlarged. Coherent overflow uncomputation restores the zero mode at the
same finite L. No extra Hamiltonian was evolved to obtain either transfer;
the energy here means the specified logarithm of these measured record
kernels. It has not been attached to a physical particle rest energy.

For finite vectors v_t, site and link reflection forms are positive:

    sum_st <v_s,T^(s+t)v_t> = ||sum_t T^t v_t||²,
    sum_st <v_s,T^(s+t+1)v_t>
      = <sum_s T^s v_s, T sum_t T^t v_t> >=0.                (13)

Equation(13) gives temporal reflection positivity on this finite charged
operator span with the common reference(4). It does not prove all OS axioms,
a relativistic field hierarchy, or a spatial mass dispersion. Products of two
O operators vanish in this finite source/reference sector, so a many-particle
CAR or Gaussian hierarchy cannot be borrowed from its two-point kernel.

## 6. Complete pointer experiment and controls

For a public initial ghost j and flavor pair a, prepare
`(|S_(j,a)>+|R_a>)/sqrt2` with ghost j. The terminal comparator is the SINGLE
fixed current-controlled observable

    X=sum_i (O_(i,a)+O_(i,a)†)

on the selected flavor block, zero on its orthogonal source complement.
The complete effects are(I±X)/2. A finite implementation service can use the
same full Naimark construction as the dressed-history packet: for J=X,
P0=I-J², the Hermitian unitary [[J,P0],[P0,-J]] and a retained Hadamard pointer
give maps [(I±J)psi|0> ± P0 psi|1>]/2. The query and source remain outputs.
No new M6 compilation of that readout is asserted here.

The certificate retains the FOUR fine pointer/query outcomes. For a branch
of original weight w, c=<J>, t=<J²>, their weights are

    w(1+t+2c)/4, w(1-t)/4, w(1+t-2c)/4, w(1-t)/4.           (14)

On a matched coherent pair c=t=1. If an actual source overflow has been
recorded, the source and reference belong to separately labelled pure-arm
branches, each with c=0; their probabilities are still retained. A HOLD
control with a mismatched ghost has c=0, t=1/2. It would be incorrect to use
a complement-identity swap in the Hadamard test instead of(14).

For L=1,2 the program uses the c=0 component's left endpoint and center and
enumerates every word through horizons0,...,4. Every seam/selector word has
its original weight6^(-t). Source versus reference overflow flag words can
split it into two original half-weight outcomes. The full6,220 word rows
retain79,656 fine pointer slots across the two native policies and HOLD.
Their endpoint pointer differences reproduce the corresponding columns of
K^t and P^t. All original probabilities sum to one at each horizon.

The one-step moved-neighbor contrast is1/6 on the native source and0 on
HOLD, with the same frozen ghost and comparator. HOLD's TOTAL contrast at
the left endpoint is5/6, and is retained rather than falsely called zero.
At the same boundary the self-endpoint contrast is2/3 when actual overflow
is read and5/6 when it is coherently uncomputed. Their difference1/6 tests
the stated memory policy on the exact same source object. Correctly dressed
native transport is the positive plant.

The original run makes26,352 literal pinned endpoint calls and52,704 C/E
coherent proposal-column evaluations, including all36 flavor pairs and all
reference columns. Later word propagation reuses this exhaustively executed
one-step table: its native path-execution share is explicitly zero and its
table-compression share100%. Markov composition of the exact full instrument
proves the compression; it is not a fresh call count. Normalization, matched
comparator identities and coordinate changes are useful checks but carry no
physical row clause by themselves.

## 7. What this adds and what remains

A finite sufficient CURRENT ghost supports one stationary positive dressed
charge reference and actual complete record instruments. No infinite path tape
is required by this particular readout. The retained actual-overflow and
coherently uncomputed policies are distinct; their charged spectra differ
although their bounded native population motion is the same. Every electric
state, active-memory value, flavor and overflow outcome is accounted for.

The source/Gauss/recoding services, finite windows, uniform preparation and
coherent source/reference/readout services remain declared N+ apparatus.
The family factor is supplied and the unchanged source is identical on it.
The zero modes or cutoff gaps in(11)-(12) are not identified charged-lepton
minima, and no mass ratios are inferred. Physical clock/energy, electric unit,
family selection and all other accessible particle branches remain obligations.
A full Lorentz continuum is not imposed as an extra prerequisite on static
row11. Same-group proof and exact replay confer no G6 or platform award.
