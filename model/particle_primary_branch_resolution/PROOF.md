# Concentrated pair response and the odd variational sector

On primary member zero, the original charge-six pair response has weight
greater than \(3/20\) in the transfer interval \([4,8]\). Its original
autocorrelation has real part less than \(-3/8\) at source time \(1/2\).
A complete finite apparatus retains these signs and probabilities with
explicit errors. An independent twelve-vector physical \(A_{1u}\) trial
resolves the different electric denominators dressing a scalar pair.
Its Rayleigh bound concerns the bottom of an odd sector; neither that
bound nor the concentrated response isolates a pole or a rest mass.

## 1. One source and two distinct uses of vectors

Fix member zero of the [primary contract](../particle_common_contract/PROOF.md):
eight cells, spacing one, all sixty matter modes per cell, the complete
degree-one product link cutoff, physical Gauss carrier, the selected
polynomial density, idle conjugate reference and source clock. The
interaction coefficient is one throughout. Operators are compressed from
complete parent products; no multiplication identity is assumed for
arbitrary compressed link operators.

The [occupation proof](../particle_primary_occupation_bound/PROOF.md)
and [refinement response](../particle_primary_refinement_response/PROOF.md)
supply the full-occupation electric Schur inequality, including inactive
Dirac rows and all number sectors. The [pair response](../particle_primary_response_moments/PROOF.md)
supplies the actual matched beta chart and exact local residual norm.
The [cubic/parity construction](../particle_primary_pair_channels/PROOF.md)
supplies the physical cube action and its inversion. Only local character
integrals and rational solves are reused from the
[ground trial](../particle_primary_variational_readout/PROOF.md).

Selected probabilities below always use \(\rho\). The beta-filled product
vector \(\Omega\), and its scalar pair excitation \(\phi\), are auxiliary
physical trial vectors. Their Rayleigh quotients give full-sector energy
upper bounds. They supply no lower bound on a selected-state residue.
The finite cooling transport is the existing
`Lean/QFT/GroundRecordTail.lean` argument, consumed through
[mass readout](../particle_primary_mass_readout/PROOF.md).
`Lean/QFT/PeriodicMomentum.lean` requires a supplied conserved translation;
it supplies no momentum of the open cube. None of these finite results
identifies a source eigenvalue with a physical rest energy.

## 2. Retaining the native charge-six defect penalty

Let \(D_r\) count positive-beta particles and negative-beta holes in native
block \(r\). The native dimensions, Casimirs and charges are exactly

\[
 (d_r,C_r,Q_{6,r})=(3,4/5,-4),(2,9/20,-3),(6,21/20,1),
                 (3,7/10,2),(1,3/10,6).
\]

Retain onsite fractions \(r_r=(0,0,0,0,4/5)\). The cube has interior
degree three, so each incident bond receives endpoint coefficient

\[
 m_r=1-r_r,\qquad
 \Delta=4\sum_r d_r(\sqrt{m_r^2+1/4}-m_r),\qquad
 w_r=(105C_r/2-\Delta)^{-1}.                     \tag{1}
\]

Every radical in the certificate is rounded upward. The full nonconstant
link block has positive margin \(63/4-\Delta>0\). The inherited inverse
order argument reduces its constant-link Schur loss to

\[
 m_r(D_{r,x}+D_{r,y})-w_r Q_r\ge-c_rI,\qquad
 Q_r=N_x+N_y-N_xN_y/d_r-4S_x\cdot S_y/d_r.        \tag{2}
\]

For the four larger blocks use the proved four-defect vacuum crossblock
bound. For charge six use the exact rational positive-pivot test on all
256 two-site Fock states, in its 25 conserved-number blocks. This is the
same quartic as (2), with two active and two inactive rows per endpoint.
It is not a truncation to the product vacuum or the two-particle sector.
Summing the operator inequalities, then restricting to the physical
carrier, gives

\[
 H\ge-720-12\sum_r c_r-6/7+(12/5)D_6.           \tag{3}
\]

The magnetic term is bounded on the full source by \(-6/7\), and the
boundary has twelve interior hopping bonds. The exact source ground
upper \(U_0\) and selected excess \(1/16\) give

\[
 \operatorname{Tr}\rho D_6
 \le\frac{U_0+1/16+720+12\sum_r c_r+6/7}{12/5}
 <1.840<15/8.                                  \tag{4}
\]

The upper \(U_0\) is regenerated from its actual variational producer.
Its trial is used only in this energy inequality. Equation (4) holds
for the unchanged selected polynomial state, not for a preparation
chosen to maximize pair success. This allocation uses degree three;
no transfer to another primary member is asserted.

## 3. A narrower original response and a fixed-time sign

For the sixteen matched charge-six pairs

\[
 B_{x\alpha}=c^\dagger(p_{x\alpha})c(m_{x\alpha}),\qquad
 \mu=\frac1{16}\sum_{x,\alpha,i,j}p_i
           |(B_{x\alpha})_{ji}|^2\delta_{E_j-E_i},
\]

the full-Fock effect inequality gives

\[
 w:=1-\frac{\langle D_6\rangle_{\rm upper}}{16}
 \le\mu(\mathbb R)\le1,\qquad w>113/128.          \tag{5}
\]

The actual corner residual in the inherited chart obeys

\[
 \|[H,B_{x\alpha}]-6B_{x\alpha}\|^2
 =\hbox{at most }3/2+\sqrt2<102/35=:M.
\]

The last strict inequality follows from \(\sqrt2<99/70\). Compression
and physical restriction contract the parent norm. Thus

\[
 \int(\omega-6)^2d\mu\le M,\qquad
 \mu([4,8])\ge w-M/4>3/20.                      \tag{6}
\]

The receipt keeps the stronger rational value, above 0.156. These are
original weights for a uniform local-address experiment. They are not
the spectral measure of the individual coherent \(A_{1u}\) component.

There is a source-time consequence without a shrinking filter. Put

\[
 F(t)=\int e^{it\omega}d\mu(\omega)
     =\frac1{16}\sum_{x,\alpha}
       \operatorname{Tr}\rho B_{x\alpha}^\dagger\alpha_t(B_{x\alpha}).
\]

At \(t=1/2\), write \(u=(\omega-6)/2+3-\pi\). Since

\[
 \cos(\omega/2)=-\cos u\le-1+u^2/2,\qquad
 0<\pi-3<1/7,\qquad\sqrt M<12/7,
\]

Cauchy--Schwarz and (5) give

\[
 \operatorname{Re}F(1/2)
 \le-\frac{97}{98}w+\frac{M}{8}+\frac{12}{98}
 <-3/8.                                        \tag{7}
\]

The total mass is at most one in bounding the mixed first moment. In the
negative coefficient of the total mass its lower bound \(w\) is used.
Equation (7) is a fixed finite-time oscillatory response on member zero.
It establishes no periodicity of the complete response and no isolated
line. The original real contrast obeys

\[
 F(0)-\operatorname{Re}F(1/2)
 \ge w-\left(-97w/98+M/8+12/98\right).
\]

The cooling transport with initial-excitation cut one gives final
excitation weight at least \(7/80\) in \([4,9]\). Initial momentum is
uncontrolled, so this is not a rest-energy interval.

## 4. Complete measurements of the band and the sign

Reuse the seventeen-outcome probe instrument

\[
 K_j=B_j/4\quad(1\le j\le16),\qquad
 K_f=\sqrt{I-\sum_jK_j^\dagger K_j}.
\]

The complete two-energy instrument from
[joint spectroscopy](../particle_primary_joint_spectroscopy/PROOF.md)
has per-clock resolution \(1/32\), tail allowance \(1/512\) and complete
implementation allowance \(1/512\). Its total event loss is at most

\[
 2/512+2/512=1/128.
\]

All \(17J^2\) outcomes and the idle reference are retained. The reported
transfer interval \([63/16,129/16]\) therefore has original probability
at least \(9/64\). Any extra compiler, preparation or sampling event
error subtracts separately. The initial energy measurement preserves the
stationary source marginal and changes its purification coherences.

The sign (7) has a separate complete Ramsey experiment with 34 outcomes:

\[
 L_{b,\pm}=\tfrac12(K_b\pm\alpha_{1/2}(K_b)).      \tag{8}
\]

The two signs cancel the cross terms, and completeness of the two
conjugate probe instruments proves \(\sum L^\dagger L=I\). Score the
successful branches by \(+1,-1\), with score zero on both failure
branches. Its unconditioned score expectation is exactly
\(\operatorname{Re}F(1/2)\). It retains the original source state as
input and the complete quantum reference. This experiment does not use
an initial energy record.

The existing [Cayley service](../particle_filtered_probe/PROOF.md) acts on
the original \(H\), including its scalar energy zero. Choose a finite
number of subdivisions with operator evolution error at most \(1/1024\).
Conjugation changes the complete Stinespring map by at most twice that
error. The same factor bounds the signed correlation score directly by
operator Cauchy--Schwarz, since the aggregate success map is contractive.
Thus the complete-output half-trace error is at most \(1/512\), including
an idle reference, and the specified finite rational service has score expectation
at most \(-191/512\). This is a declared finite algebraic service; its
large primary tape and histogram are not executed. A further complete
compiler error \(\sigma\) costs at most \(2\sigma\) in the signed score.

## 5. An actual scalar-odd electric packet calculation

Let \(\Omega\) fill both negative-beta rows of every native block with
constant links, and set

\[
 \phi=\frac14\sum_{x,\alpha}B_{x\alpha}\Omega.
\]

The sixteen occupied particle-hole vectors are orthogonal, so \(\phi\)
is normalized. It is physical, has the same five native occupations as
\(\Omega\), and has operator/state type \(A_{1u}\) under the physical
cube and inversion. The product anchor is inversion even: each cell has
an even number of occupied negative-beta modes, cell permutations carry
no odd Fock exchange sign, and its links are constant. The one-link ghost
constraints are the inherited physical restriction. No rotation of the
asymmetric ambient dangling-link allocation is imposed.

Write \(H=H_0+V\), with the frozen dimensionless split, and
\(H_0\phi=(-720+6)\phi\). The magnetic and hopping expectations of
\(\phi\) vanish by Haar integration. For a fixed bond and one-dimensional
charge-six block, use its full eight-mode Fock algebra. The local vector
\(\Phi=(B_{x1}+B_{x2}+B_{y1}+B_{y2})\Omega_{xy}\) has norm squared four.
With the actual Wilson matrices, the two oriented hopping images have
combined squared norms

\[
 \sum_{\pm}\|P_{D=2}T_\pm\Phi\|^2=4,\qquad
 \sum_{\pm}\|P_{D=4}T_\pm\Phi\|^2=2.             \tag{9}
\]

The executed exact CAR calculation verifies (9) for all three axes.
The pair lies outside this bond with probability \(3/4\); there the
vacuum hopping variance is one and it produces a second pair. Inside
the bond its reduced vector contribution is \(|\Phi\rangle\langle\Phi|/16\).
Haar integration removes opposite link characters. Thus the global
one-pair and two-pair packets per bond have squared norms \(1/4\) and
\(3/4+2/16=7/8\). Distinct bonds have distinct link characters, hence
orthogonal Haar images even if their fermion occupations overlap.

For twelve bonds the charge-six packets have norms squared and energies
relative to \(-714\)

\[
 (3,63/4),\qquad(21/2,87/4).                    \tag{10}
\]

The first moves the existing pair and costs only one electric excitation;
the second creates an additional pair and includes six onsite units.
For another native block \(r\), the independent anchor factor gives the
ground packet norm \(12d_r\) and relative energy \(6+105C_r/2\).
The five loop packets have norm squared \(g=6/22050\) and relative
diagonal \(210C_r\).

The resulting trial space has dimension twelve. The bipartite sign
operation preserves \(\phi\) and reverses hopping, so every hopping
triple expectation vanishes, including the one/two-pair projections.
Two hopping factors cover at most two links; a plaquette has four.
A private link removes every hopping-hopping-loop matrix element.
Loop triples therefore retain exactly the inherited native character
matrix \(N\). In an orthonormal packet basis the endpoint compression is

\[
 \begin{pmatrix}
 0&u_h^\dagger&\sqrt g\,\mathbf1^\dagger\\
 u_h&\operatorname{diag}(\Delta_h)&0\\
 \sqrt g\,\mathbf1&0&A
 \end{pmatrix},\qquad
 |u_{h,j}|^2=n_j,\quad A_{rs}=210C_r\delta_{rs}-N_{rs}/420.   \tag{11}
\]

All packet vectors have word degree at most one and are retained by the
actual member-zero cutoff. Matrix elements of \(V\) between them equal
the parent matrix elements. The full trial space is not invariant.

Matrix (11) has one negative Ritz root \(-d\), where

\[
 d=\sum_j\frac{n_j}{\Delta_{h,j}+d}
     +g\mathbf1^\dagger(A+dI)^{-1}\mathbf1.       \tag{12}
\]

Positive rational pivots and dyadic bisection enclose it to width
\(2^{-60}\). The code evaluates the rational Rayleigh vector at the
lower endpoint. This proves a full physical \(A_{1u}\)-sector energy
upper bound at coupling one, with no perturbative extrapolation.
The difference between this upper bound and the ground trial's upper
bound is not a bound on an excitation energy. The receipt prints that
difference only with this distinction; the full ground lower bound
instead supplies an upper bound on the odd-sector bottom above the
full ground. The full ground's own parity is not established.

## 6. The excluded complement and the resolution boundary

Let the rational trial be
\(\psi=(\phi-v_t)/\sqrt{1+s}\), with packet coefficients obtained from
(12), squared packet norm \(s\), and linear overlap \(q\). Write its
relative Rayleigh energy as \(e\), and its rational root endpoint as
\(d_{\rm lo}\). The full perturbation obeys
\(\|V\|\le J=360+6/7\). Its projection from \(v_t\) back to the trial
space has squared norm \(q^2+g\|(-N/420)t_l\|^2\). Hence

\[
 \|(H+714-e)\psi\|^2
 \le\frac{(q+e)^2+(d_{\rm lo}+e)^2s
       +J^2s-q^2-g\|(-N/420)t_l\|^2}{1+s}.       \tag{13}
\]

This includes every omitted finite state and requires no commutation of
the projected one/two-pair dressing with \(V\). Its evaluated error is
too large to isolate a spectral component. The bare pair trial's exact
variance is \(363/2+6/4410\), which also fails to isolate a line. These
are failures of the displayed complement estimates, not a statement
that the source has no particle branch.

Equations (5)--(7) also permit several frequencies: a measure of mass
\(113/128\) equally split at \(9/2\) and \(15/2\) satisfies their
moment cap and fixed-time sign, as does a single line at six. Such
finite measures are realized by explicit ground-to-level contraction
operators in the tests. This comparison diagnoses the information in
the moment inequalities; it is not substituted for the primary source.

A resolved source branch requires a controlled omitted-space resolvent
or an equivalent spectral exclusion and selected original residue. An
experimentally resolved threshold requires an admitted frequency
resolution and a bound on the probability below its claimed onset; (6)
provides neither an onset nor exclusion of lower transitions. Rest mass
further requires that same channel's physical momentum, preparation
momentum control and finite-domain dispersion errors. A declared SI tick
can immediately convert energy units, but does not supply those spectral
premises. The physical-mass table in GeV is empty. Uniform refinement is
an additional requirement for a continuum claim, without being imposed
as a prerequisite on a scoped finite-capacity comparison.

## 7. An exact obstruction to retained-penalty isolation

Consider the spatially uniform per-native allocations underlying (1)--(3),
with arbitrary retained onsite fractions \(0\le r_j<1\) and a common
retained electric fraction \(0\le s<1\). At member zero each endpoint has
degree three. Thus \(m_j=1-r_j\) and

\[
 \Delta=4\sum_jd_j(\sqrt{m_j^2+1/4}-m_j),\qquad
 w_j=(105(1-s)C_j/2-\Delta)^{-1}.                 \tag{14}
\]

Only allocations satisfying the inherited positive nonconstant-link
margin are admitted. Suppose the full-occupation inequalities are proved
with any scalar losses \(c_j\):

\[
 m_jD_j-w_j\mathcal Q_j\ge-c_jI.
\]

The two-site beta-filled vector has \(D_j=0\) and
\(\langle\mathcal Q_j\rangle=d_j\), by the exact active/inactive-row
covariance calculation in the occupation proof. Taking this expectation
forces \(c_j\ge w_jd_j\). This necessary condition also applies to an
exact solution of every occupation block. It does not replace
\(\mathcal Q_j\) by its expectation in an operator inequality: its vacuum
crossblock is nonzero. The resulting comparison has the form

\[
 H\ge-720-\mathcal L+K,\quad
 \mathcal L=12\sum_jc_j+6/7,\quad
 K=3\sum_jr_jD_j+\frac{105s}{2}\sum_\ell C_\ell. \tag{15}
\]

The physical penalty has an explicit minimum on the anchor complement. A single
matter defect changes total \(Q_6\) by a nonzero native charge, whereas
the full Gauss carrier has total matter \(Q_6=0\). Hence \(D=1\) is
absent. At \(D=0\) the matter vector is its unique local product singlet.
Ghost Gauss constraints force dangling links to their constant words.
In an edge-irrep decomposition, the support of nonconstant interior links
has no degree-one vertex: one nontrivial irrep cannot couple to a singlet
at a matter-singlet vertex. Any nonempty such support contains a cycle.
The interior cube has girth four, and each nonconstant link has Casimir
at least \(3/10\). Its full electric energy is therefore at least
\(4(105/2)(3/10)=63\). This argument uses physical Gauss invariance;
an ambient one-link excitation costs only \(63/4\).

A same-native particle-hole singlet with constant links attains the matter
penalty \(6r_j\). Its spatial and Dirac scalar sum is a nonzero physical
\(A_{1u}\) vector. A charge-six character around one plaquette attains
electric penalty \(63s\), with \(D=0\). Both vectors lie in the actual
product cutoff. Consequently, writing \(r_* = \min_jr_j\),

\[
 \inf\operatorname{spec}(K|_{\Omega^\perp})
   =\min(6r_*,63s),\qquad
 \inf\operatorname{spec}(K|_{A_{1u}})\le6r_* .    \tag{16}
\]

Zero retentions are included: the corresponding infimum can be zero.
The argument does not identify the physical ground of \(H\) with the
product anchor. Min-max applied to (15) supplies a second-eigenvalue floor,
and restriction to \(A_{1u}\) supplies an odd-sector floor. For either
isolation test the obtainable floor is at most
\(-720-\mathcal L+6r_*\).

There is an exact two-case ceiling over the complete parameter domain.
The function \(\sqrt{m^2+1/4}-m\) decreases for \(m\ge0\), so

\[
 \Delta\ge60(\sqrt{5/4}-1)>7,
 \qquad
 r_*\ge1/2\ \Longrightarrow\
 \Delta\ge60(\sqrt{1/2}-1/2)>12.                 \tag{17}
\]

The strict rational checks are \((67/60)^2<5/4\) and
\((7/10)^2<1/2\). Setting \(s=0\) enlarges the positive denominators and
decreases \(w_j\). Equations (14)--(17) imply

\[
\begin{array}{c|c|c}
 r_*\text{ range}&\mathcal L\text{ strict lower bound}
     &\text{comparison-floor strict upper bound}\\ \hline
 [0,1/2]&12\sum_j\dfrac{d_j}{105C_j/2-7}+6/7
       =\dfrac{184218}{24871}
       &-\dfrac{18016725}{24871}<-724.4069\\[6pt]
 [1/2,1)&12\sum_j\dfrac{d_j}{105C_j/2-12}+6/7
       =\dfrac{2867416}{274505}
       &-\dfrac{198863986}{274505}<-724.4457 .
\end{array}                                                       \tag{18}
\]

The authenticated full-source ground upper is
\(U_0\simeq-723.7678199426574>-724\). Thus neither floor reaches
\(U_0\), for any admitted parameters in this allocation family, even
with exact local occupation losses. It cannot certify an isolated ground
or place the scalar-odd sector strictly above the full ground by this
comparison. These ceilings concern the lower-bound method, not the actual
excited energies. A sharper ground upper could change the comparison.
A sharper electric inverse, coupled-cluster estimate or controlled
omitted-space resolvent is outside this obstruction. The
selected state, source action, reference and clock are unchanged.
