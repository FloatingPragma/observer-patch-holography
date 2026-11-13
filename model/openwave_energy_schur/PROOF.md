# Energy-dependent full-Haar Schur enclosure

## The actual source and its full complement

Use the coarse and refined operators, clocks and complete Haar frames of
`openwave_refined_haar_spectrum`. The refined projection contains all 201
states of electric energy below 5, including every intertwiner at that energy.
Write the unmodified full source in blocks as H=[[A,T*],[T,C]]. The complete
complement satisfies C>=d with d=208 for the coarse member and d=80 for the
refined member. These are source-energy bounds, without a GeV calibration.

Let D=alpha E_Q+F beta and V=beta Q(sum_f W_f)Q. Then C=D-V,
D>=d+M, ||V||<=M=F beta. The existing complete off-block row/column proof
supplies the exact rational b2>=||T||^2. Every direct omitted action row is
retained. No bound treats the omitted representations as absent.

For real z<d define Sigma(z)=T*(C-z)^(-1)T and
G(z)=T*(D-z)^(-1)T. The program evaluates G(z) at each nonzero endpoint
from the actual complete Haar action. A row of electric energy E contributes
its exact outer product divided by alpha E+F beta-z. This positive Gram
construction includes off-diagonal matrix entries and the entire directly
reachable complement.

The resolvent identity gives

    ||Sigma(z)-G(z)|| <= b2 M / ((d+M-z)(d-z)).

This holds on the whole real half-line z<d, including negative z. It does not
use a Taylor approximation in z. With K=(D-z)^(-1/2)V(D-z)^(-1/2), its norm
is at most r=M/(d+M-z)<1. Retaining words 0 through N in the Neumann series
gives the all-representation remainder

    b2 M^(N+1) / ((d+M-z)^(N+1)(d-z)).

The executed matrix is N=0. Higher-order numbers in the receipt are remainder
bounds for that formula, not claims that the additional word matrices ran.

The endpoint certificate also computes J(z)=T*(D-z)^(-2)T from exactly the
same rows. Its largest absolute row sum, enclosed with rational algebraic
bounds, supplies j(z)>=||J(z)||. Expanding the resolvent between two factors
of (D-z)^(-1)T gives the sharper bound

    ||Sigma(z)-G(z)|| <= M j(z)/(1-M/(d+M-z)).

Indeed the middle series is V+V(D-z)^(-1)V+..., with norm at most M/(1-r),
and the two outer factors have squared norm ||J(z)||. The smaller of this
bound and the unweighted bound is used in the actual inertia check. This
weights the full directly reachable complement by its individual electric
energies; it does not introduce an omitted-state cutoff.

## Full-source eigenvalue brackets from finite inertia

For x in P and y in Q, completing the quadratic form gives

    <(x,y),(H-z)(x,y)>
      = <x,(A-z-Sigma(z))x>
        + ||(C-z)^(1/2)(y+(C-z)^(-1)Tx)||^2.

Hence the number of eigenvalues of H below z equals the negative inertia of
the finite Schur matrix A-z-Sigma(z). The bounded T and strictly positive
C-z justify this identity on the form domain; the second square adds no
negative directions. The finite frame checker authenticates every eigenvalue
of A-G(z), using a complete fixed integer frame, exact Gram and residual
bounds, and outward-rounded algebraic source coefficients. Proposed LAPACK
coordinates carry no authority until this check passes.

At the lower endpoint exactly one eigenvalue remains negative after the
full self-energy error is included. At the upper endpoint at least two do.
The first excited full-source eigenvalue lies between these endpoints, with
multiplicity retained. Subtracting the authenticated parent ground interval
and intersecting the parent gap interval gives the reported sharper gap.
The parent ground bound is reused as a theorem, not relabelled as a separate
calculation. Neither bracket follows from diagonalizing a truncated Hamiltonian
alone.

## A third patch and the global-norm obstruction

The incidence control evaluates cubic patches with L=1,2,3 subdivisions,
F=3L^2(L+1) faces, spacing 1/L, alpha=32/L and beta=L^3/32. These declared
couplings have g^2=64/L^2 and L alpha=32 in the same source clock. They are
not a derived running coupling. The third patch calculation is an exact
incidence and Haar lower-bound control, not a spectrum calculation.

Choose the normalized fundamental character state |p>=chi_p Omega of one
elementary plaquette. For every elementary face f vertex-disjoint from p,
W_f|p>=(1/2)chi_f chi_p Omega is a normalized two-loop state times 1/2,
with electric energy 6. It lies in Q=1-1_(E<5). Distinct f give orthogonal
Peter-Weyl supports, even when the f's intersect each other. For any other
face g, an edge of f absent from g and p makes the matrix element vanish
by Haar orthogonality. Thus the contributions cannot cancel. Bessel gives

    ||QHP||^2 >= beta^2 N_far(p)/4.

At most 12 elementary faces meet one cubic-lattice vertex. The four vertices
of p therefore give N_far(p)>=F-48. The executable also counts N_far exactly
on each patch. The global off-block norm for this fixed total-electric-energy
projection grows with volume at fixed beta. Local incidence cannot make it
uniform, and the global Neumann ratio approaches one as the face count grows.
This obstructs that approximation strategy, not the existence of a physical
mass gap.

The local refinement route instead reuses the interaction-picture theorem
in `openwave_su2_propagation`: each link meets at most four bounded plaquette
terms, with word rate 32 beta independent of box volume at fixed beta. The
receipt exposes its coupling dependence on all three patches. That theorem
controls local observable dynamics; it neither substitutes a product vacuum
for the interacting vacuum nor proves a uniform ultraviolet spectral gap.

The result is a finite interacting Haar-source spectral improvement. It is
not a four-dimensional continuum construction, a Clay solution or a Nature
validation.
