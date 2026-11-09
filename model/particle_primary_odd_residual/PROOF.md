# Scalar-odd spectral enclosure and coherent original response

The actual eight-cell primary source has an excited physical
\(A_{1u}\) eigenvalue with excitation energy in
\([0.54794769,15.00864937]\) source units. This enclosure follows from
the full residual of a specified twelve-vector trial, bounded by
\(5.365163597\). A separate calculation on the unchanged selected
polynomial state gives an original coherent scalar-operator response of
weight greater than \(5/1024\) in transfer interval \([3/4,45/4]\).
Neither statement isolates a spectral line. Their overlapping energy
windows do not identify the selected response with the trial eigenvalue.

## 1. Fixed source and inherited trial

Use member zero of the [primary contract](../particle_common_contract/PROOF.md):
spacing one, eight cells, twelve interior bonds, six plaquettes, all sixty
CAR modes per cell, degree-one product link cutoff, physical Gauss carrier,
the selected polynomial density, its idle conjugate reference and source
clock. The Hamiltonian has its original scalar energy zero and interaction
coefficient one. The [branch calculation](../particle_primary_branch_resolution/PROOF.md)
supplies the exact twelve-vector \(A_{1u}\) compression and rational trial
coefficients. Its trial \(\psi\) is an auxiliary vector; it does not
replace the selected preparation \(\rho\).

Let \(\Omega\) fill all negative-beta rows with constant links. Write
\(\phi=T\Omega/4\), where \(T=\sum_{x,\alpha}B_{x\alpha}\) is the sum
of the sixteen matched charge-six pair-raising operators. The inherited
physical trial is

\[
 \psi=(\phi-v_t)/\sqrt{1+s},\qquad
 v_t=X_t\phi+W\phi+M_t\phi .                         \tag{1}
\]

Here \(X\) is hopping in the four native blocks other than charge six,
\(X_t\) has their four inherited weights, and \(Y\) is charge-six
hopping. Put
\(W\phi=(t_2P_{D_6=2}+t_4P_{D_6=4})Y\phi\).
The coefficients \(t_2,t_4\) retain the different electric denominators
\(63/4\) and \(87/4\). Magnetic multiplication and its five weighted
characters are \(M,M_t\). Thus \(V=X+Y+M\).
No commutation between \(Y\) and these defect projections is assumed.

The exact local ground-trial moments are reused from
[the residual calculation](../particle_primary_residual_refinement/PROOF.md).
Their full-Fock Haar identities, and the finite cooling transport from
`Lean/QFT/GroundRecordTail.lean` through
[mass readout](../particle_primary_mass_readout/PROOF.md), retain their
source domains. A projected Ritz residual alone would omit the required
outside space; every residual below includes that space.

## 2. Exact charge-six second action

The charge-six native block is one dimensional. Its independent link
characters have modulus one and Haar-orthogonal integer powers. Use all
32 charge-six CAR modes on the actual cube. Apply the literal Wilson
matrix entries to \(\phi\), keeping every resulting occupation and
every oriented link exponent. The first-action norms are pointwise
constant in the link coordinates:

\[
 \|P_2Y\phi\|^2=3,\qquad
 \|P_4Y\phi\|^2=21/2,\qquad
 \langle\phi,Y\phi\rangle=\langle\phi,W\phi\rangle=0. \tag{2}
\]

The two sectors are orthogonal. Their second-action Haar Gram is

\[
 G_{kj}=\langle YP_kY\phi,YP_jY\phi\rangle
   =\begin{pmatrix}93/2&123/2\\123/2&348\end{pmatrix}_{kj},
       \qquad k,j\in\{2,4\}.                            \tag{3}
\]

The code executes (2)--(3) with integer Gaussian amplitudes. Multiplying
\(\phi\) by four and each hopping action by two removes every
denominator. A basis key contains the full 32-bit occupation and the
twelve link exponents. After two actions each exponent is in
\(\{-2,-1,0,1,2\}\); balanced base five is injective on this set.
Haar integration pairs equal exponent vectors. Pointwise first-action
norms are checked as complete Laurent coefficient tables, including
cancellation of every nonconstant coefficient. The retained vectors have
48 and 672 first-action terms, then 1928 and 16456 second-action terms.
This is an exact action on these vectors, not a dense simulation of the
primary physical Hilbert space. Independent tests use tuple exponents
and literal CAR action to check the same Gram.

## 3. The full parent norm and the physical residual

The other native blocks factor from the charge-six matter state. Write
their pointwise hopping variance, weighted variance and covariance as
\(z_X,s_X,q_X\), and their exact fourth mixed moment as \(M_X\).
With native dimensions \(d_j=(3,2,6,3)\), twelve bonds and the inherited
open-cube connected coefficient 18,

\[
 \begin{split}
 z_X&=168,\quad s_X=12\sum_jd_jt_j^2,\quad
 q_X=12\sum_jd_jt_j,\\
 M_X&=z_Xs_X+2q_X^2-18\sum_jd_jt_j^2.             \tag{4}
 \end{split}
\]

Equation (2) gives
\(z_6=27/2\), \(s_6=3t_2^2+(21/2)t_4^2\) and
\(q_6=3t_2+(21/2)t_4\). Equation (3) gives
\(M_6=(t_2,t_4)G(t_2,t_4)^T\). Expanding the four vectors
\(XX_t\phi,XW\phi,YX_t\phi,YW\phi\), factorization and (2) give

\[
 M_J=M_X+M_6+z_Xs_6+s_Xz_6+4q_Xq_6.              \tag{5}
\]

The two covariance cross terms are each \(2q_Xq_6\); the remaining
cross terms vanish because a first hopping expectation on either
matter factor vanishes pointwise. No statistical independence of native
representations on a common link is assumed.

Set \(z=z_X+z_6\), \(s_h=s_X+s_6\), \(q_h=q_X+q_6\).
The inherited exact character calculation supplies
\(z_l=\langle M^2\rangle\), \(s_l=\langle M_t^2\rangle\),
\(q_l=\langle MM_t\rangle\) and \(M_l=\langle M^2M_t^2\rangle\).
Magnetic multiplication commutes with matter projections. The bipartite
sign reverses every hopping and preserves \(\phi\), so odd total
hopping terms vanish. Thus the complete parent norm is

\[
 \|Vv_t\|^2=M_J+zs_l+s_hz_l+4q_hq_l+M_l.         \tag{6}
\]

Both \(\phi\) and \(v_t\) lie in the actual product cutoff. Its final
compression contracts the norm in (6). Let \(q=q_h+q_l\),
\(s=s_h+s_l\), \(g=6/22050\), and \(C_l=-N/420\) be the inherited
loop correction. The exact projection of \(Vv_t\) into the trial space
has norm squared \(q^2+g\|C_lt_l\|^2\). The free onsite/electric
operator preserves the trial space. For its rational root endpoint
\(d_{\rm lo}\) and Rayleigh shift \(e\),

\[
 \|(H+714-e)\psi\|^2
 \le\frac{(q+e)^2+(d_{\rm lo}+e)^2s+
       \|Vv_t\|^2-q^2-g\|C_lt_l\|^2}{1+s}.       \tag{7}
\]

Every omitted physical state is included. At degree one the parent
norm is an upper bound, not an asserted equality after compression.
The authenticated trial has center approximately
\(-717.8547086510794\) and residual at most
\(5.365163596221\). Finite self-adjoint spectral theory on the reducing
\(A_{1u}\) space therefore puts some eigenvalue in

\[
 [-723.219872248,-712.489545054].                  \tag{8}
\]

The same member's authenticated full ground upper is approximately
\(-723.7678199426574\), below the lower endpoint in (8). Its ground
lower is approximately \(-727.4981944205706\). Using the exact rational
intervals gives some strictly excited \(A_{1u}\) energy in
\([0.54794769,15.00864937]\). This is not a positive lower bound on
every excitation, ground uniqueness, an isolated pole, or a rest mass.
No subtraction of two trial upper bounds is used as an excitation bound.

## 4. A selected-state coherent moment

Use the existing constant-spatial, scalar channel of the complete
[pair apparatus](../particle_primary_pair_channels/PROOF.md),

\[
 K=K_{00}=\frac{T}{16\sqrt2},\qquad
 T=\sum_{x,\alpha}c^\dagger(p_{x\alpha})c(m_{x\alpha}). \tag{9}
\]

Its normalization is unchanged when read as the selected outcome of that
apparatus or in the complete two-outcome instrument
\(\{K,\sqrt{I-K^\dagger K}\}\). Sixteen independent quasispin pairs
give \(\|T\|^2=72\), hence \(\|K\|^2=9/64\). The full-Fock identity is

\[
 [K^\dagger,K]=(16-D_6)/512.                      \tag{10}
\]

The source polynomial state is stationary and the authenticated native
defect bound obeys \(\langle D_6\rangle<46/25\). In particular its
original success mass is at least \((16-\langle D_6\rangle)/512\).
No division by the success probability or by a spatial channel rank is
made.

Only charge-six hopping contributes to \(R_T=[H,T]-6T\). In the actual
beta chart its one-body blocks are

\[
 h=\begin{pmatrix}A&B\\B^\dagger&D\end{pmatrix},\quad
 S=\begin{pmatrix}0&I\\0&0\end{pmatrix},\qquad
 [h,S]=\begin{pmatrix}-B^\dagger&A-D\\0&B^\dagger\end{pmatrix}. \tag{11}
\]

Here \(A-D\) is minus the link adjacency with a two-spinor identity,
so its Hilbert--Schmidt norm squared is exactly 48 at every coordinate.
Each off-diagonal spin block of \(B\) has norm \(1/2\), and each cube
vertex has degree three, giving \(\|B\|\le3/2\).

In beta particle/hole variables, (11) splits into pair creation \(L\)
with coefficient \(A-D\), and a number-preserving \(N=d\Gamma(n)\).
The 32-dimensional one-body matrix is
\(n=\operatorname{diag}(-B^\dagger,-(B^\dagger)^T)\), has trace zero,
and has norm at most \(3/2\). The defect operator is precisely the
number of these quasiparticles. Singular-value decomposition of the
pair-creation coefficient reduces \(L\) to commuting hard-core pair
creators. On a sector with \(q\) quasiparticles, at most \(q/2\) pairs
are occupied. Embedding their zero/one occupations in bosonic Fock
space and compressing its creation operator gives

\[
 L^\dagger L\le48(1+D_6/2).                      \tag{12}
\]

Blocked single occupations only remove available pair modes. For \(N\),
the ordinary particle bound and the trace-zero hole bound give
\(\|N|_{D_6=q}\|\le(3/2)\min(q,32-q)\). Since
\(\min(q,32-q)^2\le16q\),

\[
 N^\dagger N\le36D_6.                            \tag{13}
\]

These are full-Fock coordinatewise form inequalities, hence hold on the
Haar carrier. The product link projection commutes with \(T,D_6\), so
compression and then Gauss restriction preserve the resulting bound.
With Young parameter \(6/7\), divide the bound for \(L+N\) by 512:

\[
 ([H,K]-6K)^\dagger([H,K]-6K)
 \le\frac{624+858D_6}{3584}.                     \tag{14}
\]

It follows on the selected state that
\(\int(\omega-6)^2\,d\mu_K<13767/22400\). This uses its actual
stationarity and occupation bound. The auxiliary trial is absent from
the expectation.

## 5. Original coherent windows and complete records

Equations (10) and (14), with the symmetric interval \(6\pm21/4\), give

\[
 \mu_K([3/4,45/4])>
 \frac{177}{6400}-\frac{13767/22400}{(21/4)^2}
 =\frac{7057}{1317120}>\frac5{1024}.              \tag{15}
\]

Both negative and positive out-of-window transfers are charged to the
same second moment. No polynomial-slope small-frequency approximation
is used. The same-source cooling transport with initial cut two costs
at most \(\|K\|^2\epsilon/2=9/2048\), so the original final-excitation
law has weight

\[
 \nu_K([3/4,53/4])>
 \frac{10151}{10536960}>\frac1{2048}.             \tag{16}
\]

For the two complete energy clocks, choose per-clock resolution \(1/64\),
tail allowance \(1/32768\), and complete implementation allowance
\(1/32768\). The inherited source-bound clock plan then has total event
loss at most \(1/8192\). All \(2J^2\) probe and clock records retain
the reference. The reported transfer interval
\([23/32,361/32]\) has original probability at least \(39/8192\).
An extra preparation, compiler or sampling error is deducted separately.
The exact finite services are specified; a large primary histogram is not
executed.

The final interval in (16) lies within the broad excitation enclosure
from (8), but this overlap supplies no spectral identification. The
coherent probe has cubic operator type \(A_{1u}\); a mixed input state
does not force its output to have that state irrep. The trial spectrum
and the selected response therefore remain separate consequences of one
Hamiltonian. An isolated component and its original residue require a
spectral projection or excluded-complement estimate joining them. Rest
energy also requires physical momentum and controlled dispersion of that
same component. The GeV rest-mass table is empty; no continuum or
refinement-uniform version of these member-zero estimates is asserted.
