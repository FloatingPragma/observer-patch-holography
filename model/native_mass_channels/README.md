# Onsite masses of the native Weyl hopping

This package classifies every Hermitian onsite one-body perturbation that
anticommutes with the existing three-dimensional Weyl hopping. The answer is
four real coefficients, with all eight Weyl valleys retained. It gives the exact
uniform-mass spectrum, three connected gapped coefficient regions, a constructive
spin-wall unwinding path, and a finite exact degree witness for a continuous
spin-direction texture. These are statements about a supplied quadratic mass
family. They do not derive contact-induced order, an interacting mass shell,
physical particle species, a defect zero mode, or native repair protection.

Run the exact checks from the repository root:

```sh
python -W error -m pytest model/native_mass_channels -q
python -W error -m model.centered_source_ground --check
python -W error -m model.native_mass_channels --check
```

`algebra.certificate()` returns a JSON-compatible report. Public numeric inputs
are Python integers or canonical `fractions.Fraction` values; booleans, floats,
complex numbers, malformed fractions, and incomplete coordinate lists are
rejected. Matrices are immutable SymPy matrices. Reports are reconstructed on
each call. No floating eigensolver, fitted tolerance, or rounded winding number
is used.

Receipt validation has two stages. The centered-ground packet independently
replays its full integer Gram and charged-sector proofs through its own
`model.verify` entry. This packet validates that parent's declared conclusions
and complete current source manifest, then replays the entire mass catalogue.
It does not rerun the parent's large Gram calculation during each mass check.
The mass validator compares every catalogue field with a fresh algebra report,
including nested scope flags, and requires the complete current source manifest.

## Exhaustive real-space classification

For the actual source link $h_{x,x+e_j}=ic\sigma_j/2$, with $c>0$, an onsite
Hermitian $B_x$ satisfies $\{h,B\}=0$ precisely when

\[
B_{x+e_j}=-\sigma_j B_x\sigma_j.
\]

The three real-linear maps on the right commute and square to the identity.
Consequently an arbitrary Hermitian $B_0$, with four real coordinates,
extends uniquely to the entire connected cubic lattice and has period two.
Writing $\eta_x=(-1)^{x_1+x_2+x_3}$, the complete answer is

\[
B_x(m_0,\mathbf m)=m_0\eta_x I+
\sum_{j=1}^3m_j(-1)^{x_j}\sigma_j.
\]

The finite certificate retains all 96 real edge equations on the eight-site
period-two cell: their rank is 28 on 32 real onsite unknowns. The four displayed
columns span its kernel. The recursion proves extension and exhaustiveness
without presupposing periodicity. The nonzero-speed assumption matters: with
zero hopping every onsite field would anticommute.

`source_anchor_certificate()` separately reconstructs the actual open-cube
`FieldGeometry.single_h()` at $c=1/8$. Exact rational conversion of its binary
matrix entries verifies all four anticommutators and computes rank 28 for the
complete real onsite linear system formed directly from that source matrix.
The tests reject a Hermitian source link with a wrong Pauli axis. The source
anchor also verifies the existing hidden-frame matrices, rather than importing
the mass basis as a premise about those matrices.

For the source frame $R_x=(i\sigma_1)^{x_1}(i\sigma_2)^{x_2}
(i\sigma_3)^{x_3}$, with $d_x=R_xc_x$,

\[
(-1)^{x_j}\sigma_j=\eta_xR_x^*\sigma_jR_x.
\]

Thus the three spin masses are staggered hidden-spin (Néel) fields. The centered
local quadratic observable is

\[
c_x^*B_xc_x-\tfrac12\operatorname{Tr}B_x
=m_0\eta_x(N_x-1)+\eta_x\mathbf m\cdot(d_x^*\boldsymbol\sigma d_x).
\]

This centering is essential for the scalar channel: at a half-filled state with
local density one, the uncentered scalar term has expectation $m_0\eta_x$.
The centered-ground symmetry consequences are owned by
[`centered_source_ground/particles.py`](../centered_source_ground/particles.py).
Zero order expectation in a finite symmetric ground does not prove an absence
of an interaction gap or exclude infinite-volume symmetry breaking.

## Complete folded spectrum

Fold $p+\pi b$, $b\in\{0,1\}^3$, into one reduced Brillouin zone, retaining
the two spin rows of every valley. Let

\[
\gamma_j|b\rangle=-(-1)^{b_j}|b\rangle\sigma_j,
\quad\beta_0|b\rangle=|b\oplus(1,1,1)\rangle I,
\quad\beta_j|b\rangle=|b\oplus e_j\rangle\sigma_j.
\]

The resulting matrix is $16\times16$:

\[
H(p)=c\sum_j\sin(p_j)\gamma_j+m_0\beta_0+\mathbf m\cdot\boldsymbol\beta.
\]

Each generator is a Hermitian involution. The three kinetic generators
anticommute pairwise and with every mass. The three spin masses anticommute
pairwise; the scalar mass commutes with every spin mass. For $r=|\mathbf m|>0$,

\[
Q=\beta_0\mathbf m\cdot\boldsymbol\beta/r,
\qquad Q^2=I,\quad[Q,H]=0,
\qquad H^2=c^2\sum_j\sin^2(p_j)+m_0^2+r^2+2m_0rQ.
\]

The exact trace identities $\operatorname{Tr}Q=\operatorname{Tr}H=
\operatorname{Tr}(HQ)=0$ give two eight-dimensional $Q$ sectors with balanced
positive and negative levels. Therefore the complete spectrum is

\[
E=\pm\sqrt{c^2\sum_j\sin^2(p_j)+(m_0\pm r)^2},
\]

with multiplicity four for each of the four generic signed levels. Equal levels
have their multiplicities added, including zero levels; `spectrum()` also covers
$r=0$. The gap from zero in the infinite uniform bulk is exactly

\[
\Delta=\bigl||m_0|-r\bigr|.
\]

The tests compare a direct full-matrix characteristic polynomial with this
formula and check its degenerate cases. This is a bulk gap, not a claim about
arbitrary open boundaries or spatially varying mass fields.

For a pure third-axis spin mass, fixing the first two valley bits gives four
invariant four-row blocks. The 16 products of the four Clifford generators in
each block span its full matrix algebra. `dirac_decomposition_certificate()`
checks this irreducibility and the complete 16-row partition. Rotations generated
by $\beta_i\beta_j$ commute with every kinetic generator and rotate the spin
mass direction, so every nonzero pure spin mass has four massive Dirac blocks.
These are algebraic copies of a supplied free massive operator; no Weyl valleys
are discarded and no physical species assignment follows.

## Restricted coefficient-space topology

The gapped coefficients have exactly three connected components. This is the
topology inside the specified onsite anticommuting family, not a classification
under arbitrary local Hamiltonian perturbations or added degrees of freedom.

| Coefficients | Deformation retract | Fundamental group | Second homotopy group |
| --- | --- | --- | --- |
| $m_0>r$ | $(1,0,0,0)$ | $0$ | $0$ |
| $m_0<-r$ | $(-1,0,0,0)$ | $0$ | $0$ |
| $\lvert m_0\rvert<r$ | $S^2$ | $0$ | $\mathbb Z$ |

The first two cones are convex. For the central component set
$s=|m_1|+|m_2|+|m_3|>0$. The rational homotopy implemented by
`component_retraction(m,t)` is

\[
(m_0,\mathbf m)\longmapsto
\bigl((1-t)m_0,[(1-t)+t/s]\mathbf m\bigr),\qquad0\leq t\leq1.
\]

It preserves $r>|m_0|$ and retracts onto the octahedral sphere
$m_0=0,|\mathbf m|_1=1$, which is radially homeomorphic to $S^2$.
The imported sphere facts are $\pi_1(S^2)=0$, and classification of maps
$S^2\to S^2$ by their integer degree; see Hatcher,
[Proposition 1.14, Proposition 2.30, and Corollary 4.25](https://pi.math.cornell.edu/~hatcher/AT/AT.pdf).
Their use is limited to these continuous coefficient-space statements.

There are three component labels, so walls between different components cannot
be joined by a gapped path within this family. Within the central component a
spin sign wall has no such component obstruction. `wall_rotation(t,M)` gives an
exact rational path from $(0,M,0,0)$ through $(0,0,M,0)$ to $(0,-M,0,0)$,
with $|\mathbf m|=M>0$ throughout. All components have trivial fundamental
group; a loop of masses in this family carries no nontrivial winding class.
The central component admits integer point-texture degree on a surrounding
two-sphere.

## Finite texture degree and its limits

`texture_degree()` uses the six vertices and eight outward-oriented faces of an
octahedron. On each domain face it linearly interpolates three rational target
vectors and normalizes radially. It requires an invertible target-column matrix
$A$ on every face, which guarantees that the interpolant never reaches zero.
For a nonzero rational ray $q$, compute $\lambda=A^{-1}q$. If every component
of $\lambda$ is positive, that face contains one preimage of the ray direction,
with local degree $\operatorname{sign}\det A$. Summing these signs computes the
integer degree. The evaluator conservatively rejects a zero component of
$A^{-1}q$ on any face, including faces not intersecting the ray; callers may
choose another generic ray. Degenerate image faces require a different
triangulation or evaluator and are not silently assigned a degree.

The identity and antipodal vertex assignments have degrees $+1$ and $-1$.
A separate hemisphere control has degree zero. `hedgehog_stability()` supplies
an exact perturbation certificate: if every target vertex differs by at most
$\delta<1/3$ in each coordinate from the identity (or antipodal) assignment,
each face has a signed separating normal along the entire straight homotopy.
Its dot product with the interpolated target is at least $1-3\delta>0$.
Radial normalization therefore remains defined and the degree is preserved.
This homotopy certificate also applies when the ray evaluator's narrower
nondegeneracy domain does not apply.

A nonzero degree obstructs filling the ball by a continuous mass direction that
stays in the central component. It does not by itself imply a fermionic index,
bound state, energy minimum, dynamical stability, or survival under a native
complete instrument. An arbitrary discrete vertex assignment has no such
continuous invariant without the declared interpolation and its nonvanishing
condition. In particular, position-dependent coefficients add gradient terms
to $\{h,B\}$; the tests explicitly reject extending the uniform bulk-gap
identity to a spatially varying rotation. Additional restrictions, boundary
conditions, and a separate full-operator theorem would be needed to identify
protected physical defects.
