# An isolated interacting gauge transition in the actual source

Member zero has a simple A1u gauge level within that symmetry sector
of its empty-matter reducing carrier. Its energy difference from the
empty-matter gauge ground is
contained in

\[
 [\,63-59/1351350,\ 63+1/45675\,].
\]

A specified Wilson-loop contraction couples those two exact eigenstates
with squared matrix element greater than \(4/49\). The unchanged selected
density therefore has a positive transfer atom, with its original
rare-sector probability retained. This is a finite source-unit gauge
transition, not a charged fermion or a GeV rest mass.

## 1. The actual reducing carrier

Use [primary member zero](../particle_common_contract/PROOF.md), with the
same physical Gauss projection, 480 matter modes, native gauge group,
degree-one link word spaces, selected polynomial density, idle reference,
Hamiltonian and clock. The
[extreme-sector proof](../particle_primary_survival/PROOF.md) establishes

\[
 [H,N_F]=0,\qquad
 H|_{N_F=0}=H_E+B.                                      \tag{1}
\]

Every matter bilinear vanishes on the empty Fock vector. The fully filled
sector is a second, conjugate copy, but the detector below labels the
empty sector explicitly. Zero external ghost charges reduce every
dangling link to its constant vector. The remaining graph is the actual
eight-vertex, twelve-edge cube.

No source action is removed by hand in (1). It is a reducing sector of
the full specified operator. Its density is the restriction of the
original \(\rho=p(H)\), with **no conditional normalization**. It carries
the full native product group, not a selected compact-simple factor.

The five native irreducible representations have
\[
 (d,C,Q_6)=(3,4/5,-4),(2,9/20,-3),(6,21/20,1),
          (3,7/10,2),(1,3/10,6).
\]
Each actual degree-one link space is the constant plus all coefficients
of these representations and their duals. Exact Weyl integration checks
their ten characters are orthonormal. Thus they are distinct irreducibles,
and the link dimension is
\(1+2(3^2+2^2+6^2+3^2+1^2)=119\).
The electric energy on an occupied link is \(105C/2\). Its two smallest
values are \(63/4\) and \(189/8\). This uses the actual native metric.

## 2. Gauss support and an isolated electric cluster

Decompose the link space into its representation blocks. With no matter
charge, a vertex with exactly one nontrivial incident representation has
no invariant vector. A nonconstant support of a physical state therefore
has no degree-one vertex. Every finite nonempty such graph contains a
cycle. The cube is simple and bipartite. Among its edge subsets with at
most five edges, exactly its six square faces have no degree-one vertex.
The code independently enumerates all4096 subsets.

The electric ground is the unique constant-link vector \(\Omega\).
The first nonzero physical electric energy is
\[
 E_1=4(63/4)=63.                                         \tag{2}
\]
For a four-edge square every vertex has degree two. Schur's lemma forces
a single representation, with its appropriate dual orientation, around
the cycle. The invariant contraction is its character, once. At energy63
that representation is the one-dimensional charge-six character or its
dual. Hence this level has exactly twelve orthonormal vectors: both
orientations of each of the six faces.

Any support with at least six edges costs at least \(6(63/4)=189/2\).
A four-edge face in any other representation costs at least
\(4(189/8)=189/2\). There is no five-edge physical support. Consequently
the next electric eigenvalue is at least \(189/2\).

The actual magnetic multiplication is
\[
 B=-\frac1{210}\sum_{p=1}^6 S(U_p),\qquad
 S=\sum_r(\chi_r+\overline{\chi_r}).
\]
Since \(\sum_r d_r=15\), its full parent norm is at most \(6/7\).
Compression to the source and restriction to Gauss support contract
this norm. Finite spectral perturbation therefore gives one eigenvalue
in \([-6/7,6/7]\), twelve in \([63-6/7,63+6/7]\), and all other
eigenvalues at least \(189/2-6/7\). These separated clusters hold at the
actual interaction coefficient one. The interpolation \(H_E+tB\),
\(0\le t\le1\), proves constancy of each cluster's rank; it is not an
extrapolation beyond a small-coupling theorem.

## 3. Exact cube labels and residual moments

Orient each face outward. Proper rotations permute its six outward
normals; inversion exchanges opposite faces and reverses loop
orientation. On the real face characters inversion therefore acts by
face exchange, and on the imaginary characters by its negative. The
literal 48 proper source lifts and inversion give

\[
 \begin{array}{c|c}
 \text{real characters}&A_{1g}\oplus E_g\oplus T_{1u}\\
 \text{imaginary characters}&A_{1u}\oplus E_u\oplus T_{1g}.
 \end{array}                                                \tag{3}
\]

Every summand occurs once. The Hamiltonian commutes with these finite
symmetries. Its isolated cluster carries the same representation along
the interpolation, because the character traces of its continuous
spectral projector are integer combinations of a fixed finite character
table. The ground remains \(A_{1g}\). The first \(A_{1u}\) level is simple
within its symmetry sector, and the next \(A_{1u}\) energy is at least
\(189/2-6/7\). No nondegeneracy across other cube representations is inferred.
These are CAR-even gauge states; \(A_{1u}\) denotes spatial inversion,
not odd fermion number or Lorentz spin.

Put \(y_p=(\chi_6(U_p)-\overline{\chi_6(U_p)})/i\) and
\[
 \phi=\frac1{\sqrt{12}}\sum_p y_p.
\]
Distinct face characters are Haar orthogonal, so phi is normalized and
has type \(A_{1u}\). The actual character calculations give
\[
 \int S^2=10,\quad \int y^2=2,\quad
 \int Sy^2=0,\quad \int S^2y^2=18.                       \tag{4}
\]
For the last identity, \(y^2=2-\chi_{12}-\chi_{-12}\).
The only two native charges summing to12 are6+6, and their character
product is \(\chi_{12}\); the dual statement is identical. Thus (4)
is \(2(10)-1-1=18\). The code also obtains it from the complete exterior
weight polynomial and exact Weyl constant term, retaining all native
representations.

Products involving distinct face holonomies are evaluated using private
edges. Every fourth-order multiplicity pattern other than four equal
faces or two doubled faces has an edge with a single nontrivial word,
so its integral vanishes. Two different doubled faces have independent
Haar holonomies, even when they share an edge, by integrating one private
link of each. Terms containing \(Sy\) have zero integral by conjugation.
Consequently
\[
 \langle\Omega,B\Omega\rangle=\langle\phi,B\phi\rangle=0,
 \qquad \|B\Omega\|^2=\frac1{735},\qquad
 \|B\phi\|^2\le\frac{59}{44100}.                         \tag{5}
\]
The second bound is the exact **parent** integral
\((18+5\cdot20)/(2\cdot210^2)\). Final product-cutoff projection contracts
it. Products of two characters on the same link need not lie in the
degree-one space, so equality is not asserted for the compressed second
moment. The first vacuum action is admitted exactly.

For a normalized vector with Rayleigh mean mu, variance at most v and
a unique lowest eigenvalue below a certified second level b>mu, the
finite spectral inequality
\[
 E_{\min}\ge\mu-\frac{v}{b-\mu}                           \tag{6}
\]
follows by averaging \((H-E_{\min})(H-b)\ge0\). Apply it to Omega in
the full gauge sector and phi in the \(A_{1u}\) sector. Rayleigh upper
bounds and (5) give
\[
 -\frac1{45675}\le E_g\le0,\qquad
 63-\frac{59}{1351350}\le E_u\le63.                       \tag{7}
\]
This is a spectral enclosure of actual eigenvalues, not a Ritz label.
The same-symmetry complement of the \(A_{1u}\) level lies at least
\(189/2-6/7-63=429/14\) above its Rayleigh center.

## 4. A nonzero exact transition residue

Use the source contraction
\[
 A=\frac1{12i}\sum_p(\chi_6(U_p)-\overline{\chi_6(U_p)}).
 \qquad \|A\|\le1,\quad A\Omega=\phi/\sqrt{12}.           \tag{8}
\]
The source compression preserves these statements. Let g,u be normalized
exact eigenvectors for (7), with phases chosen to give positive overlaps
with Omega,phi. For either anchor the spectral component orthogonal to
its lowest eigenvector has norm at most its residual norm divided by
the distance from its mean to the rest of its spectrum. Both residual
norms in (5) are below1/27, and the respective distances exceed62 and30.
Thus
\[
 \|g-\Omega\|+\|u-\phi\|
 \le\sqrt2(1/1674+1/810)
 <\frac1{364}.
\]
Since \(1/\sqrt{12}>15/52\) and \(15/52-2/7=1/364\),
\[
 |\langle u,A g\rangle|>2/7.                            \tag{9}
\]

Within the physical empty sector define the exact spectral projectors
\(P_g=1_{H<1}\) and
\(P_u=P_{A_{1u}}1_{62<H<64}\). They have rank one by the proved
separations. Their extensions vanish on all other native occupation
sectors. The resolved source operator \(B_{\rm tr}=P_u A P_g\) is a
contraction and has exactly
\[
 \mu_{B_{\rm tr}}=p(E_g)|\langle u,A g\rangle|^2
                        \delta_{E_u-E_g}.               \tag{10}
\]
This is an atom in the original source probability law. The raw loop
probe also contains that positive transition, but may have other
frequencies. Equation(10) does not claim that two bare loop states form
an invariant dynamical subspace. It uses the actual isolated spectral
projectors, defined by fixed rational separating intervals.

For \(B_{\rm tr}+B_{\rm tr}^*\), the positive-frequency commutator
coefficient is \([p(E_g)-p(E_u)]|\langle u,A g\rangle|^2\).
The exact selected modular slope gives
\(p(E_u)/p(E_g)\le\exp[-\beta_{\rm low}(E_u-E_g)]<1/2\).
It is a nonzero source response pole as well. Neither (10) nor the
commutator identifies spatial momentum, relativistic dispersion or a
physical species.

## 5. Original weights and complete clock records

The [energy-local density bound](../particle_primary_branch_certificate/PROOF.md)
applies to this same selected polynomial. With the actual global ground
lower bound L0, ambient dimension bound \(2^{b_D}\), and \(E_g\le0\), put
\[
 c_g=2^{-b_D-\lceil2\beta(0-L_0)\rceil}.
\]
Then \(p(E_g)\ge c_g\). The original positive transfer atom is larger than
\(4c_g/49\), and its positive commutator coefficient larger than
\(2c_g/49\). This is a source-marginal inequality; no Loewner bound on the
pure source-reference density is asserted.

A complete finite experiment uses the initial source-energy clock, the
two Kraus operators \(AP_{\rm empty}\) and
\(\sqrt{I-P_{\rm empty}A^2P_{\rm empty}}\), and complete final cube and
source-energy measurements. All native-sector, symmetry, clock, failure
and reference outputs are retained. The clock uses the unchanged H,
resolution eta=1/64, and authenticated one-clock error
delta<1/512. Its implemented kernels are functions of H. Starting in the
actual eigenvector g, the accepted initial window
\([-1/45675-\eta,\eta]\) therefore retains at least \(1-\delta\).
The u contribution to the final \(A_{1u}\) window
\([63-59/1351350-\eta,63+\eta]\) retains at least
\((4/49)(1-\delta)\). Positivity gives original joint probability
\[
 p_{\rm event}>
 c_g\,\frac4{49}(1-\delta)^2>c_g/13.                    \tag{11}
\]
Clock errors are included in each eigenvector's positive effect before
multiplying by c_g. A fixed absolute error is not subtracted from a tiny
original probability. Other initial states can produce clock leakage;
(11) is a guaranteed contribution, not a perfect finite detector for an
unknown exact eigenvalue. The reported energy difference is in the
interval (7) for \(E_u-E_g\), widened by2eta. No extreme-sector success
probability is divided out.

## 6. The physical boundary

The empty-sector ground has energy near zero in the source convention,
whereas the global source ground is below -723. This sector's energy
origin cannot be substituted for the physical global vacuum. The source
proves
\[
 \operatorname{Tr}\rho(P_{\rm empty}+P_{\rm full})
 \le \frac{\epsilon a}{90V-P/7}
 \le \frac7{627\,2^{8n+7}}.
\]
Thus every contraction supported in these extreme sectors has total
original spectral weight tending to zero on the selected refinement
trajectory. Conditioning on the sector, changing its preparation or
injecting it with an apparatus has a different probability/work account.

The isolated finite gauge transition is compatible with that obstruction.
It supplies neither a surviving continuum particle nor a pure compact-
simple Yang-Mills construction. The CAR-odd trion's multiplicity-space
complement requires a separate bound. The known 48 bare minimal trion
states and the scalar retained-penalty estimates do not isolate that
branch at interaction one. The GeV mass table is empty.
