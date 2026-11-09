# A common gauge-coordinate field, current and recorded experiment

The construction chooses the gauge-coordinate completion of the native
link algebra. It then fixes a Hamiltonian, a state and a common model clock
and computes complete observer records. None of these software choices
requires a uniqueness theorem. The distinction from a particle prediction
concerns the outputs: physical chirality, the interacting continuum and
the observed spectrum require the corresponding dynamical and measurement
proofs.

[`field_route.py`](field_route.py) supplies a finite-cutoff bridge for this
chosen field program and executes an interacting spinor/quantum-flux
control. The control uses the actual native singlet charge, four Dirac CAR
modes and every admitted flux record. It has a nonzero exact current and
retains both chiralities. Its one spatial cell is not a demonstration of
three-dimensional propagation.

## The constructive route and its inherited domains

The [retained coordinate completion](../interacting_extension/chiral_lattice/PROOF.md)
embeds the native 60-level link into matrix coefficients in
\(L^2(G,\mathrm{Haar})\), with \(G=S(U(3)\times U(2))\). Its multiplication
operators include transitions outside the native finite representation
menu. Choosing this completion specifies a different field program from
the raising-only link dynamics; it does not assert equality of their
readouts. Haar integration, electric Casimirs, source charge
representations and the spatial-overlap matrix conventions have explicit
definitions in that package.

The [free spatial-overlap census](../interacting_extension/chiral_lattice/light_branch.py)
retains every momentum, spin and internal component. It has one free zero
momentum, with the other corners at cutoff energy, and a controlled free
continuum symbol. The [Euclidean overlap implementation](../common_source/overlap_regulator.py)
supplies the modified projectors, gauge covariance and a uniform polynomial
approximation on its admissible gap domain. These results are reused;
none needs a second pole census here.

The [electric connection](../interacting_extension/chiral_electric/PROOF.md)
preserves a fixed-rank Euclidean chiral bundle and includes its geometric
potential. It does not by itself conserve the spatial Hamiltonian's
modified axial charge. The [canonical quantum-link dressing](../particle_calibration_selection/PROOF_FAMILIES_CHIRALITY.md)
supplies another finite route to exact CAR while retaining its dark
complement. It does not remove mirror modes. Architecture choice removes
no obligation to check those outputs on one Hamiltonian.

## A bounded overlap polynomial on the full coordinate Hilbert space

On a fixed finite cubic lattice, let \(H_W=\beta(aD_W-I)\) be the Hermitian
Wilson kernel in the retained fifteen-state representation. Its norm is
at most 8 for every unitary coordinate-link configuration. In three
directions, the Wilson operator is \(aD_W=3I-\sum_\mu W_\mu\), where
each \(W_\mu\) is a unitary sum of the two orthogonal spin projections
times forward and backward covariant shifts. Hence
\(\|aD_W-I\|\le5<8\), independently of plaquettes or a spectral gap.
Use the same Legendre polynomial as the existing overlap implementation:

\[
 q_m(H_W)=\frac4{17}H_W\sum_{j=0}^{m}
       (15/17)^jP_j\!\left(\frac{257I-8H_W^2}{255}\right),
 \qquad h_m=\frac{\beta+q_m(H_W)}a.                     \tag{1}
\]

It is a gauge-covariant Hermitian polynomial of coordinate-word degree at
most \(d=2m+1\). On the admitted gap \(H_W^2\succeq I/4\), the retained
locality proof gives

\[
 \|q_m(H_W)-\operatorname{sign}H_W\|
 \le\epsilon_m=16(15/17)^{m+1}.                         \tag{2}
\]

There is also a useful bound on the whole coordinate domain, including
small Wilson gaps:

\[
                  \|q_m(H_W)\|\le1+\epsilon_m.          \tag{3}
\]

For a real spectral value \(|x|\ge1/2\), (3) follows from (2). For
\(0<|x|<1/2\), the Legendre argument lies between 1 and \(257/255\).
All \(P_j(z)\) are positive there, as follows from
\(P_j(1+u)=\sum_{k=0}^j\binom jk\binom{j+k}k(u/2)^k\) for \(u\ge0\).
Their infinite generating series
equals \(1/|x|\) after multiplication by \(4/17\). A finite partial sum
therefore obeys \(0\le|q_m(x)|\le1\); at \(x=0\), the polynomial is
zero. This proves (3) without a gap assumption. Thus
\(\|h_m\|\le(2+\epsilon_m)/a\) globally. Equation (3) does not extend
the accuracy statement (2) across gap closures.

Use the full fermionic Fock space of the finite Dirac/M15 lattice. With
\(n_f=4\cdot15\,|\Lambda|\) CAR modes,
\(\|d\Gamma(h_m)\|\le n_f\|h_m\|\). Positive electric Casimirs plus
this bounded interaction and bounded coordinate-independent gauge-invariant
local terms define a
self-adjoint Hamiltonian on the electric domain. A chosen plaquette term
\(-\kappa\sum_p\operatorname{Re}\operatorname{Tr}U_p/15\) has norm at
most \(|\kappa|\,|\mathcal P|\) and coordinate degree four. These
coefficients define the program; they are not obtained from particle data.

## A controlled representation cutoff for complete observer histories

Let \(\mathcal V_K\) be the finite span of products of at most \(K\)
native M15 matrix coefficients or their conjugates across all links. This
is total representation-word degree, not a coordinate-grid cutoff.
The spaces are nested and invariant under endpoint gauge transformations.
Electric Casimirs preserve each space. Their Haar-orthogonal projectors
therefore preserve the Gauss constraints. Full Fock occupations are not
truncated.

A compactly supported nonempty Gauss-zero preparation exists without
particle masses or a vacuum fit: fill all 60 spin/internal CAR modes at
one site, leave the other sites empty and use the constant Haar wavefunction
on every link. The native M15 representation has zero trace for every
generator, so its determinant character is trivial; the full local Slater
state is invariant. This preparation has coordinate degree zero and
nontrivial allowed hopping. It is a chosen finite initial state, rather
than an interacting Lorentz-invariant vacuum.

A degree-\(d\) interaction maps \(\mathcal V_K\) into
\(\mathcal V_{K+d}\). From initial support in \(\mathcal V_{K_0}\),
reaching the omitted complement of \(\mathcal V_K\) needs at least

\[
 r=\left\lfloor\frac{K-K_0}{d}\right\rfloor+1
\]

bounded-interaction insertions. The electric interaction picture does not
change that support count. If \(X\) bounds the integrated interaction
norm, the full and compressed evolutions agree through order \(r-1\),
giving

\[
 \|\psi(T)-\psi_K(T)\|
 \le2\sum_{n\ge r}\frac{X^n}{n!}.                    \tag{4}
\]

Occupation measurements and neutral pointer readouts preserve the degree
filtration. Their complete record dilations obey the same bound, including
arbitrary untouched references. Bounded polynomial current-pointer
couplings are included as additional interaction intervals in \(X\) and
in the maximum degree \(d\). Consequently the total variation of the
entire recorded distribution is bounded by the right side of (4), capped
at one. No particular outcome is postselected. The code encloses the
exponential tail with rational summation and a geometric remainder; no
floating exponential supplies that bound.

The compressed Hamiltonian is \(\Pi_KH_m\Pi_K\). Computing the sign
polynomial after compressing its Wilson kernel generally gives a different
operator: \(q_m(\Pi_KH_W\Pi_K)\ne\Pi_Kq_m(H_W)\Pi_K\). The independent
Fourier-mode control detects that boundary error. A halo of the polynomial's
full word degree preserves all intermediate paths before compression.

For the same polynomial field, local number and its hopping current obey

\[
 N_x=\sum_\alpha c_{x\alpha}^\dagger c_{x\alpha},\quad
 J_{x\to y}=i(c_x^\dagger h_{m,xy}c_y-
                     c_y^\dagger h_{m,yx}c_x),\quad
 i[H_m,N_x]=-\sum_yJ_{x\to y}.                        \tag{5}
\]

Equation (5) holds for the kinetic/electric/plaquette program and local
terms commuting with \(N_x\). Pair-creating Higgs vertices require their
explicit local number source term. The complete bilinear current is gauge
invariant. A neutral pointer coupled to it measures a function of the
same Hamiltonian, rather than an independently assigned particle current.
Its polynomial coordinate degree is at most \(2m+1\).

## Executed common-state spinor and quantum-flux control

Choose the native \(q=6Y=6\) singlet. On one periodic spatial cell, give
one direction a quantum rotor \(U=e^{iq\theta}\) and set the other two
directional links to identity. This is a declared one-coordinate control,
not an asserted invariant reduction of the complete three-dimensional
gauge field. The source gamma matrices give

\[
 H_W=-\beta\cos(q\theta)+i\beta\gamma_1\sin(q\theta),
 \quad H_W^2=I,\quad
 h=\beta+A U+A^\dagger U^\dagger,
 \quad A=(-\beta+\beta\gamma_1)/2.                    \tag{6}
\]

Thus the exact spatial-overlap sign is available here. Retain all four CAR
modes and all integer fluxes, with

\[
 H=\frac{E^2}{2}+d\Gamma(\beta)
       +d\Gamma(A)T_q+d\Gamma(A)^\dagger T_q^\dagger,
 \qquad T_q|\ell\rangle=|\ell+q\rangle.              \tag{7}
\]

The electric term and the spinor/flux hopping do not commute. This is a
quantum gauge-coordinate interaction, rather than a classical background
average. All terms share the declared duration parameter. There are no
GeV or laboratory-time inputs.

The one-cell charge constraint with a static neutralizing charge \(-2q\)
is \(G=q(N_F-2)\). It commutes with (7). Its complete physical occupation
block has six states; the parent operator retains all sixteen occupations.
Restricting the computation to this invariant Gauss-zero block is exact.
At cutoff \(|\ell|\le K\), every boundary vector is retained and the
compressed shift's loss of unitarity is explicit.

The same Hamiltonian defines an electric-flux current

\[
 \dot E=i[H,E]=iq\big(d\Gamma(A)^\dagger T_q^\dagger
                              -d\Gamma(A)T_q\big).       \tag{8}
\]

In a source-derived spin basis with \(\beta=\operatorname{diag}(1,1,-1,-1)\),
prepare the two positive-beta occupations and the rotor state
\((|0\rangle+i|q\rangle)/\sqrt2\). It lies in the complete Gauss-zero
space and gives the exact initial values

\[
             \langle E\rangle=3,\qquad
             \langle\dot E\rangle=6,\qquad
             \langle H\rangle=11.                      \tag{9}
\]

The initial state, relative phase, electric coefficient and clock are
specified apparatus/program data. Equation (9) follows from them; it is
not a measured particle-mass identification.

The implementation evolves the state for \(1/64\), measures the complete
signed flux, evolves for a further \(1/64\), and measures again. At
\(K=18,30\) it retains all 1,369 and 3,721 joint outcomes, including zeros.
Flux modulo six is an exact invariant of this specified Hamiltonian.
The prepared residue is evolved as a complete reducing block; all detector
labels in the other residues retain their exactly zero probabilities.
Only a hopping insertion changes flux, by at most six. The conservative
full-Fock interaction norm bound is eight, so (4) applies with
\(K_0=6\) and \(X=8/32=1/4\). The code compares the entire distributions,
padding absent cutoff labels by zero, against the sum of their rational
infinite-cutoff error bounds.

The finite evolution has a separate rational error enclosure. In the
chosen spin basis the doubled Hamiltonian \(2H_K\) has Gaussian-integer
entries. The implementation applies the order-64 Taylor polynomial
\(P_{64}(-itH_K)\) using Gaussian-integer Horner recursion and a common
integer denominator. It retains every first-readout branch without
renormalization before applying the second pulse. Squared amplitudes give
nonnegative rational record weights; their sum differs from one only by
the explicitly bounded Taylor error. The displayed decimal entries are
conversions of these exact weights, rather than floating matrix
exponentials.

The onsite term has norm two: \(\beta\) has two eigenvalues \(+1\) and
two \(-1\), and full-Fock occupation sums range from \(-2\) to \(2\).
Each hopping lift has norm at most four. Therefore

\[
 \|H_K\|\le K^2/2+10=:B_K,\qquad
 \|e^{-itH_K}-P_{64}(-itH_K)\|
 \le \frac{(tB_K)^{65}}{65!}=:\varepsilon_K(t).       \tag{10}
\]

The second inequality follows from the integral Taylor remainder and the
unitarity of the Hermitian evolution; it needs no exponential norm factor.
If the two pulse errors are \(\varepsilon_1,\varepsilon_2\), the complete
record-isometry error is bounded by
\(\delta=\varepsilon_1+\varepsilon_2+\varepsilon_1\varepsilon_2\).
For exact normalized record vector \(\psi\) and computed vector
\(\widetilde\psi\), the trace norm of their rank-one density difference
is at most \(2\delta+\delta^2\). Dephasing into all record labels gives

\[
 \frac12\sum_{a,b}|p_{ab}-\widetilde p_{ab}|
 \le \delta+\delta^2/2=:\eta_K.                      \tag{11}
\]

Equation (11) compares a normalized unitary distribution with the
unnormalized rational weights. It bounds their normalization error by
\(2\eta_K\), every event error by \(2\eta_K\), and each signed-flux
mean error by \(2K\eta_K\). At \(K=30\), the rational enclosure gives
\(\eta_K<1.155\times10^{-35}\). The cutoff comparison uses exact
rational differences and adds the two bounds (11) to the cutoff
allowance (4); it does not treat decimal agreement as an error proof.

The field Hamiltonian is executed directly as the specified finite
coordinate program. No compilation and execution of its complete native
C/E pulse tape is claimed by this calculation.

## Output obligations beyond this finite construction

The configured polynomial gauge-coordinate field has controlled finite
representation cutoffs and complete current/occupation records. Its
Dirac module is vectorlike. A physical chiral theory requires an
anomaly-compatible measure or canonical regulator with controlled
complementary modes. A Ginsparg-Wilson identity or omission of detector
weight on those modes does not supply that result.

The exact-overlap limit also needs gauge-state control near small Wilson
gaps. If a normalized evolving state has probability at most \(\eta\)
on configurations violating the gap in (2), a direct good/bad
configuration split gives the one-body-to-full-Fock comparison

\[
 \|(d\Gamma(h)-d\Gamma(h_m))\psi\|
 \le\frac{n_f}{a}\big[\epsilon_m+
                           (2+\epsilon_m)\sqrt\eta\big]. \tag{12}
\]

A common-clock Duhamel bound integrates (12) along the actual state and
instrument history. A physical application must supply \(\eta\) for
that experiment; the one-cell identity \(H_W^2=I\) does not establish it
for all interacting nonabelian configurations.

Uniform interacting volume and refinement limits, the appropriate positive
state, and the same detector's particle poles, widths and scattering
response are separate output requirements. Calibration to physical units
and comparison with measured species follow those definitions. No
architecture-uniqueness requirement is imposed on any step here.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_observer_outputs.field_route
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/particle_observer_outputs/test_field_route.py
```
