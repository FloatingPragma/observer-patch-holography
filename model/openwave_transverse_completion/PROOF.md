# Native matter and a compact transverse plaquette

Each external vertex carries the complete 192-CAR longitudinal member of
`openwave_common_electromagnetism`. Its internal two-cell Gauss rotors and
quantum reference are retained. The external spatial graph is separate.
Every positive Dirac excitation has its actual native electric charge;
colourful annihilators carry the native fractional charges. The formal
quartic Ward certificate proves that the complete internal source commutes
with its vertex charge, including the semileptonic terms.

For an oriented external edge use L2(U(1)), integer electric operator E and
unitary multiplication U. Charges are measured in thirds: q=3Q is integer.
The complete constraints are div E=q+R delta_0, with a quantum compensating
reference R and the full constraint-failure outcome. Add

    K=t sum_edges,j (a_xj-U_edge^qj a_yj)* (a_xj-U_edge^qj a_yj),
    H_field=(a/2) sum_edges E_edge^2+b sum_plaquettes(1-Re U_loop),
    t=1/16, a=1/16, b=1/8.

Here a_j denotes the positive Dirac annihilator or native colourful
annihilator specified above. The internal Gauss isometries dress those
operators too. All terms preserve every external and internal constraint.
The kinetic form is positive on the full Fock space, annihilates its actual
Dirac/colourful vacuum, and does not alter the inherited lower bounds.
The coefficients are declared architecture, without electromagnetic target
selection. This external finite-range kinetic is distinct from the
square-root kinetic in `openwave_boosted_records`.

The closed nonnegative form defines the full rotor source. On a finite
compact graph the field Hamiltonian is an elliptic Schrödinger operator
with a real bounded potential. Its unique strictly positive ground is
invariant under gauge transformations, which commute with the operator.
Thus its global ground energy E_gamma is attained in the zero-divergence
sector. Subtracting this one field vacuum gives the inherited matter floor
on every physical charge sector. The vacuum is attained by the matter
vacuum and that field ground. This finite-graph statement is not a positive
vacuum gap for a physical massless-photon limit.

## Exact one-excitation reduction

The colourful vacuum with one positive Dirac excitation is an exact reducing
matter sector: its native semileptonic action vanishes and the normalized
internal electric term is zero. Every other Fock sector remains in the
source. For the four-cycle and a charge q at x with reference -q at zero,
all integer solutions of Gauss law are

    E_e=n-q 1_(e<x), n in Z, x in {0,1,2,3}.

The executed carrier is consequently C4 tensor ell2(Z), with no declared
flux cutoff. A forward hop changes x to x+1 and leaves n unchanged except
for 3 to 0, where n changes to n-q. The reverse is its adjoint. The magnetic
plaquette changes n by plus or minus one. These actions give the literal
integer matrix H-m=A/32 in the code. Every charge in the native list is
checked; an undressed charged hop fails Gauss law. The reference-bound
charged object is not identified as an isolated free electron.

The compression to |n|<=C retains the complete diagonal and every allowed
matrix element. It is an approximation to this infinite rotor source.
Outside the compression, electric energy is at least
2a(C+1-|q|)^2. Positivity of K and of the magnetic potential preserves this
bound. The compression-to-complement norm is at most 2t+b/2. If [l,u]
encloses the entire compression's lowest eigenvalue and D>u is that outside
floor, the Schur complement gives

    l-(2t+b/2)^2/(D-u) <= inf(H-m) <= u.

The field-vacuum calculation uses D=2a(C+1)^2 and boundary norm b/2.
Subtracting its enclosed vacuum gives the one-electron rest-energy interval.
This is a minimum on that complete reducing sector, not a minimization over
every charged Fock sector of the source.
No counterterm retunes that external dressing to the electron anchor.

## Exact arithmetic for every compressed eigenvalue

Let Z/d_v be the dyadic complete eigenvector matrix and Lambda the dyadic
reported diagonal. The code forms both Gram error and matrix residual with
integer arithmetic, checking every int64 product bound before use. Their
Frobenius norm upper bounds g,r use integer square-root enclosures. If g<1,
the complete change of basis is invertible. With K>=||Lambda||,

    epsilon=[2K g+(1+g)r]/(1-g)

bounds every ordered eigenvalue error. Indeed the polar unitary of the
approximate basis differs from it by the Gram error; equivalently, expand
the Rayleigh quotient U*HU=Lambda+(U*U-I)Lambda+U*R and use
(1-g)I<=U*U<=(1+g)I. The same quotient bound on every ordered subspace
and min-max gives the stated enclosure. No eigenvector, spin copy or
complement is selected away.

## Infinite-source time records

Separate the diagonal electric energy from bounded off-diagonal terms B,
whose norm is at most 2t+b. An interaction-picture Dyson term of order k
changes n by at most k max(1,|q|). From n=0, leaving |n|<=C requires
r=floor(C/max(1,|q|))+1 transitions. For z=(2t+b)|time|<r+1, the omitted
path sum has norm at most z^r/r! divided by 1-z/(r+1). Comparing infinite
and compressed evolutions gives twice this tail. Any complete terminal
instrument has total-variation error at most that vector bound; the
retained preparation success of one half multiplies the original-input
bound by one half. This state-dependent path bound includes the full
unbounded diagonal exactly and does not infer a bound from convergence
of two finite cutoffs.

This exact tail is the Galerkin truncation component. It does not by itself
bound floating evaluation of the finite propagator. For each displayed
position/flux table the code separately certifies its numerical error.
Write Q=Z/d_v, D=diag(E/d_e) for the complete dyadic witness, with
||Q*Q-I||<=g and ||HQ-QD||<=r. Duhamel gives

    ||exp(-itH)-Q exp(-itD) Q*|| <= g+|t| r sqrt(1+g).

Each scalar phase is enclosed by rational Taylor sums after range reduction
with a fifty-decimal rational enclosure of pi. All finite matrix-vector
products and differences from the displayed complex amplitudes are then
evaluated as exact rationals. Their vector error delta gives half-L1 error
at most delta(2+delta)/4 for the original-input position/flux table, including
its success weight one half. Each live calculation must fit the fixed
outward vector budget dimension/2000000 and probability-rounding half-L1
budget 1/10^14. Publishing these verified budgets avoids changing exact
receipt strings when a LAPACK implementation chooses a different complete
eigenframe. A calculation that exceeds either budget fails; no numerical
error is assumed from a floating residual alone. Adding the tail gives the reported infinite
position/flux table bound. Every row and every eigenvector participates.
The separate displayed native-probe weights remain floating diagnostics;
the exact Galerkin-tail theorem applies to the fixed complete native
instrument, while this additional numerical enclosure concerns the
position/flux table rather than certifying the probe-weight evaluation.

The declared preparation input is the matter vacuum with reference charge
zero and the bare electric loop state n=0. Its mean field energy is b;
relative to the true field vacuum it is b-E_gamma. This initial field is
externally prepared, not identified with the field ground. The success
Kraus operator is 1/sqrt(2) times the dressed native creation at x=0 times
the projector on this input; sqrt(I-K*K) is the complete failure operator.
The failure branch initially retains the input. Its subsequent full
vacuum-sector unitary conserves mean field energy b. On success the initial
mean energy is m+2t+b. Thus the actual preparation work, including both
branches, is (m+2t)/2, rather than half the success state's total energy.

The instrument records x, n and all 96 native particle/hole probe rows.
Their weights factor through the literal isometry (A; conjugate B), not
through substituted particle labels. The displayed factorization specifies
every outcome, including exact zeros. Preparation failure has probability
one half; Fock and constraint complements are retained. Absorption transfers
charge q at the measured x to a retained local detector of declared zero
Hamiltonian. The constraint during absorption is
div E=q_matter+q_detector+R delta_0; the compensating reference and flux
remain. This is an explicit apparatus extension of the same constraint.
Projecting
onto flux destroys magnetic coherence: conditional terminal work is
-m-<K>+b-<V_mag>. Including success probability one half and the untouched
failure record, the final mean energy is b+<E_electric>/2. Consequently
W_prepare+W_readout=<E_electric>/2. The constant E_gamma subtraction cancels
from every work difference. The receipt verifies this total mean-energy
identity and reports the initial external preparation energy separately.
The finite source time is declared; no autonomous observer-clock or
relativistic continuum equivalence is supplied.


The retained rest certificates use complete dyadic proof frames for the
charged cap-18 compression (148 coordinates) and neutral vacuum caps
18, 24, and 48 (37, 49, and 97 coordinates). Every replay checks their
integer matrix identity, complete square shape, bounded integer entries,
sorted centres and exact binary round trips before applying the unchanged
`exact_cap_enclosure` consumer. That consumer recomputes the full integer
Gram and residual matrices, verifies its integer-overflow bounds and
proves the frame invertible from Gram error below 1/100. No stored error
bound or spectral interval is trusted. The original infinite-flux Schur
bounds and exact physical intervals are retained. Other exploratory caps
continue to construct a fresh frame and apply the same proof. Finite-time
probability records continue to use live propagation and independently
verified outward numerical budgets. These proof coordinates carry no
additional empirical calibration.
