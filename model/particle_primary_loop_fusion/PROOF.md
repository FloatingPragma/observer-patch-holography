# Complete source loop fusion and its energy records

Two adjacent native charge-six Wilson loops fuse by canceling their
common edge. The resulting loop has a lower source energy than two
disjoint face loops. Full-complement estimates place their actual energy
records in disjoint windows at interaction coefficient one. This is a
finite static pair-creation experiment in the empty-matter reducing
sector. No moving incoming packets, physical binding mass or scattering
cross section is asserted.

## 1. Source, orientation and complete creation instruments

Keep member zero of the [primary contract](../particle_common_contract/PROOF.md),
including all480 native CAR modes, the full gauge group, degree-one
product link cutoff, physical Gauss carrier, selected polynomial density,
retained reference and Hamiltonian clock. The
[isolated gauge-sector proof](../particle_primary_reduced_sector/PROOF.md)
establishes the exact empty-Fock reducing sector and

    H_empty = H_E+B,              ||B||<=6/7.                  (1)

On each outward oriented square f, let chi_f be the native one-dimensional
charge-six character of its loop holonomy. Define the actual compressed
contraction C_f=P_word chi_f P_word, restricted to Gauss support. Closed
Wilson loops preserve Gauss support. Projection makes C_f a contraction;
it is not assumed unitary on the finite cutoff.

For distinct faces f,g, execute the complete sequential instrument

    C_g C_f P_empty,
    sqrt(I-C_g* C_g) C_f P_empty,
    sqrt(I-P_empty C_f* C_f P_empty).                       (2)

The last branch stops after the first failure. Squaring and summing the
three Kraus operators gives I on the complete physical source. Every
branch and the untouched reference are retained. Selecting each of the
30 ordered distinct face pairs with probability1/30 gives90 complete
probe records before the energy-clock labels.

Let Omega be constant links and empty matter. C_f Omega has four signed
unit charge-six fluxes and is admitted at degree one. For adjacent
outward faces the common edge has opposite orientations. Their product
cancels that edge, leaving a simple six-edge loop. Thus

    C_g C_f Omega=chi_hexagon Omega,
    H_E C_g C_f Omega=(189/2) C_g C_f Omega.                 (3)

For opposite faces their edge supports are disjoint and

    H_E C_g C_f Omega=126 C_g C_f Omega.                    (4)

Both vectors have norm one by independent Haar integration. Every
intermediate and final word is in the actual cutoff, so (3)-(4) are
identities for the compressed operators. They require no claim that a
smaller subspace is invariant under H. The code enumerates all30 words
on the literal cube incidence matrix:24 adjacent pairs and6 opposite
pairs, with zero divergence and only flux0,+1,-1 on every link.

## 2. A full-source isolated fusion cluster

The degree-one physical electric spectrum is classified by the same
native irreps and cube support argument as the isolated gauge result.
Every nonconstant link has electric energy at least63/4. At a vertex
with just one nonconstant incident link, Gauss invariance vanishes.
All nonempty cube supports with at most five links are the six squares;
all sixteen six-link supports without a degree-one vertex are simple
hexagons. A simple cycle must carry the same irrep around the cycle,
and its invariant contraction has multiplicity one.

The actual native one-link energies are

    42, 189/8, 441/8, 147/4, 63/4,                         (5)

with the corresponding five duals. Therefore H_E has exactly the
following levels below441/4:

    0       multiplicity1,
    63      multiplicity12,
    189/2   multiplicity44=32+12.                          (6)

The32 states are the two charge-six orientations of each of16 hexagons.
The12 states are the two native weak-charge orientations on each of6
squares. The next irrep on a square and the next irrep on a hexagon
exceed441/4; every support with at least seven links has energy at
least7*(63/4)=441/4. Some such supports may have no invariant state;
the lower bound does not require their existence.

The full perturbation norm in (1) is smaller than half every separation
in (6). Finite spectral perturbation and continuous variation of the
spectral projection along H_E+tB, 0<=t<=1, therefore give a rank44
cluster in

    I_f=[189/2-6/7,189/2+6/7]=[1311/14,1335/14].           (7)

The spectrum outside this cluster is at most63+6/7 or at least441/4-6/7.
Its distance from189/2 is at least417/28. For the normalized hexagon
vector phi in (3), (H-189/2)phi=Bphi. The spectral theorem gives

    ||(I-1_(I_f)(H))phi||² <= ( (6/7)/(417/28) )²
                            =64/19321.                    (8)

No truncated Hamiltonian supplies either the cluster rank or (8).
For the normalized disjoint-loop vector psi in (4), the same argument
without a cluster claim gives

    ||1_([124,128])(H)psi||² >=1-(6/7)²/4=40/49.           (9)

## 3. Actual gauge ground and positive original records

Let g be the exact empty-sector ground, phase aligned to Omega. The
isolated gauge result supplies

    -1/45675<=E_g<=0,
    ||g-Omega|| <(10/7)/(27*62)<1/1000.                    (10)

Every product C_g C_f is a contraction. Equations(8)-(10), and
`sqrt(1-64/19321)>998/1000`, imply for every adjacent pair

    ||1_(I_f)(H) C_g C_f g||² >(997/1000)²>99/100.          (11)

For every opposite pair, `sqrt(40/49)>9/10` gives

    ||1_([124,128])(H) C_g C_f g||² >(899/1000)²>4/5.       (12)

These inequalities refer to the exact interacting ground in the rare
reducing sector. The original selected density is not replaced by g.
Its authenticated energy-local population bound is

    p(E_g)>=c_g=2^(-32406376).                             (13)

An initial complete source energy clock, the instrument(2) and a final
complete source energy clock give the actual original experiment. Use
the unchanged clock plan from the isolated gauge result: resolution
eta=1/64, one-clock error delta<1/512, with the implementation and tail
errors combined. Retain initial labels in
`[-1/45675-eta,eta]`, and final labels in I_f expanded by eta or
`[124-eta,128+eta]`, according to the recorded pair geometry.

For the source eigenvector g, the initial good-label probability is at
least1-delta. The final positive clock effect is at least1-delta on
the respective spectral window. Thus the original joint weights exceed

    adjacent: c_g*(99/100)*(1-delta)² > (49/50)c_g,
    opposite: c_g*(4/5)*(1-delta)²   > (3/4)c_g.            (14)

In the complete30-pair menu each individual pair row therefore has
original probability greater than c_g/40. This statement includes both
creation successes and the specified energy labels; the other records
are retained. No tiny sector probability is divided out. Clock errors
multiply the eigenvector probability before the original factor c_g;
no fixed absolute clock error is subtracted from c_g. The bound uses the
source marginal, not an operator domination of its retained pure
source-reference state. Additional initial-state clock leakage is not
excluded.

## 4. The energy deficit and its physical boundary

The isolated A1u gauge transition has energy E_1-E_g with
`63-59/1351350<=E_1<=63`. Every fused-cluster energy E_f in (7)
therefore satisfies the exact source comparison

    429/14-118/1351350 <= 2(E_1-E_g)-(E_f-E_g)
                       <=453/14+1/45675.                  (15)

This is a positive energy deficit relative to twice that isolated
gauge transition, approximately30.64277 to32.35716 source units.
The adjacent-loop and disjoint-loop energy windows are disjoint, and
the full source measures can distinguish their geometry. Equation(15)
is not a two-asymptotic-particle scattering threshold: the intermediate
loop was not energy projected and no incoming moving packets were
constructed. The44-level cluster may contain several energies; neither
the individual poles nor a momentum dispersion are inferred from its
finite enclosure.

The empty-sector ground near zero is not the ambient primary vacuum,
whose energy is below-723. The inherited original extreme-sector
probability bound vanishes on the selected refinement. Every bounded
probe in (2) remains supported there, so this finite source response
has no demonstrated surviving particle weight, GeV mass or physical
binding-energy interpretation.
