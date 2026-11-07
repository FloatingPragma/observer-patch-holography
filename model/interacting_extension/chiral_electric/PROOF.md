# Electric motion on a changing chiral subspace

## What is constructed

The source provides a Euclidean overlap projector \(P(U)\). This note supplies
a local kinetic completion when a gauge configuration coordinate itself is
quantum mechanical. It uses one native charge \(q=6\), the complete source
\(2^4\) site stencil (64 spin/site components), and a fixed-rank, admissible
holonomy patch. The retained chiral rank is 32.

The coordinate interval, Dirichlet boundary, kinetic coefficient and penalty
are declared controls. Neither their source selection nor a many-body
Lorentzian chiral gauge theory is established. In particular, the projector
here is not the modified axial charge of a vectorlike spatial-overlap
Hamiltonian.

## 1. A first-order electric derivative needs a connection

Differentiate \(P^2=P\). Then

\[
P P' P=0,\qquad P'=PP'+P'P.
\]

For a section \(\psi=P\psi\), the ordinary derivative has a complementary
component \((1-P)\psi'=P'\psi\). Define the anti-Hermitian connection

\[
C=[P,P'],\qquad \nabla=\partial_\theta+C.
\]

Since \([[P,P'],P]=-P'\),

\[
\nabla P=P'+[C,P]=0.
\]

Thus this covariant derivative preserves the chosen bundle. The parallel
transport equation is \(V'=[P',P]V=-CV\) in a horizontal frame: the generator
has the opposite sign from the connection. The implementation checks this
sign independently using centered derivatives of the source overlap sign
kernel and canonical transport.

This is a statement about consistent local bundle transport. It does not
assert conservation of the physical axial charge in an interacting gauge
theory, cancellation of anomalies, or existence of a global fermion measure.

## 2. Compressing the parent second-order electric kinetic energy

Let \(V^\dagger V=I\), \(VV^\dagger=P\), and write a physical section
\(\psi(\theta)=V(\theta)\chi(\theta)\). Set

\[
A=V^\dagger V',\qquad A^\dagger=-A.
\]

Applying the full parent kinetic energy to a section and compressing gives

\[
V^\dagger(-\partial_\theta^2)V
=-(\partial_\theta+A)^2+\Phi,
\qquad
\Phi=V'^\dagger(1-P)V'\succeq0.
\]

To check the sign, use
\(A'=V'^\dagger V'+V^\dagger V''\) and
\(V'^\dagger P V'=A^\dagger A=-A^2\). Therefore

\[
\Phi=V'^\dagger V'+A^2,
\]

with a **plus** sign before \(A^2\). Discarding this term changes the projected
parent model. A connection-only kinetic law can instead be postulated, but
that is a different microscopic completion; the source has not selected it
over the projected parent law.

The source patch evaluated at \(\theta=0.037\) has derivative leakage norm
approximately \(0.713671\), and the eigenvalues of \(\Phi\) range from
\(0.0122187\) to \(0.509326\). This is a nonzero effect, not a rounding residue.
The covariant projector derivative has Frobenius residual below
\(5.2\times10^{-15}\). Halving the centered-difference step reduces the
independent derivative error by approximately four, as expected.

## 3. Complete finite grid and frame covariance

On an interior Dirichlet grid, retain every 64-component parent fiber. Let
\(\mathcal V=\operatorname{diag}(V_1,\ldots,V_m)\), and let \(K\) be the full
positive nearest-neighbor parent Laplacian. The exact compressed matrix is

\[
K_{\rm comp}=\mathcal V^\dagger K\mathcal V.
\]

Its forward off-diagonal blocks are proportional to
\(-V_j^\dagger V_{j+1}\). Under arbitrary independent frame changes
\(V_j\mapsto V_j Q_j\), the whole matrix transforms by unitary conjugation.
An identity matrix on each forward edge fails this covariance.

Replacing the overlaps by their unitary polar factors gives a separate
connection-only discretization. **The finite-grid difference of these two
matrices is not asserted positive semidefinite.** Positivity is the
differential statement about \(\Phi\).

This frame covariance concerns redundant coordinates within the same
projector range; it is not a proof of the full gauge Ward identity.

## 4. Retaining the complementary space

Rather than deleting the complement, take the finite parent

\[
H_\Lambda=K+\Lambda(1-\mathcal P),\qquad
\mathcal P=\mathcal V\mathcal V^\dagger.
\]

All parent states remain present for every finite \(\Lambda\).
In a complete \(P\oplus Q\) basis, let
\(A=PKP\), \(B=PKQ\), \(D=QKQ+\Lambda Q\).
Since \(K\succeq0\), \(D\succeq\Lambda Q\). For \(E<\Lambda\), the Schur
operator is

\[
A-B(D-E)^{-1}B^\dagger.
\]

If \(a_k\) is the \(k\)-th eigenvalue of \(A\) and
\(\Lambda>a_k\), the first \(k\) low eigenvalues obey the common conservative
bracket

\[
0\le a_j-\lambda_j(H_\Lambda)
\le\frac{\|B\|^2}{\Lambda-a_k},
\qquad j\le k.
\]

This follows from the min-max principle and the Schur complement below the
complement threshold. The numerical audit evaluates this analytic bound in
float64; the reported values are not interval-arithmetic certificates.

For the nine-point grid, the full dimension is 576 and compressed dimension
288. Increasing \(\Lambda/\|K\|\) from 10 to 40 to 160 reduces the maximum
first-eight energy discrepancy from approximately \(2.95\times10^{-6}\) to
\(7.44\times10^{-7}\) to \(1.86\times10^{-7}\). The complete complement is
included in every run.

## 5. What this resolves and leaves open

It gives an explicit, checked way to join gauge-coordinate electric motion
to a changing finite chiral projector, and exposes a term that a naive
join omits. It does not select the parent kinetic law, its coefficient, the
global admissible gauge measure, the anomaly-compatible many-body lift,
reflection positivity, or an interacting continuum limit. It also does not
derive a measured particle mass.
