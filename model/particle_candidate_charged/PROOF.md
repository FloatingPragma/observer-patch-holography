# Charged representations, mass coefficients and composite dressing

The [primary construction](../particle_common_contract/PROOF.md) has a
sixty-dimensional one-body fiber at each cell, four source spin rows times
fifteen internal rows. Its selected density, canonical idle reference and
Hamiltonian clock are held fixed. The calculation classifies exact local
coefficient constraints and constructs a complete gauge-neutral apparatus
with a charged odd regional factor. It supplies no physical charged-lepton
or quark mass. Native chirality and source cube spin retain their finite
meaning until an effective field/current dictionary is established.

## 1. The actual local coefficient commutant

In the integer exterior-mask chart, the twelve anti-Hermitian gauge
generators have central operator \(Q_6=-2T_{12}/i\). Its five eigenspaces,
in native occupation order, are

| Conditional all-left name | Internal representation | Q6 | Dimension |
|---|---|---:|---:|
| u conjugate | (anti-3,1), Y=−2/3 | −4 | 3 |
| lepton doublet | (1,2), Y=−1/2 | −3 | 2 |
| quark doublet | (3,2), Y=1/6 | 1 | 6 |
| d conjugate | (anti-3,1), Y=1/3 | 2 | 3 |
| e conjugate | (1,1), Y=1 | 6 | 1 |

Let P_r be the corresponding spectral projectors. They are computed as
the exact Lagrange polynomials

\[
 P_r=\prod_{q\ne q_r}{Q_6-qI\over q_r-q}.                 \tag{1}
\]

Any gauge-commuting matrix preserves these distinct central eigenspaces.
Within a block of dimension d, the code solves the complete exact system
\(XT_a-T_aX=0\), for all twelve literal generators. Row-major
vectorization gives the matrix
\(I\otimes T_a^T-T_a\otimes I\). Its ranks are respectively
8,3,35,8,0; every block nullspace is the span of I_d. This is an exact
algebraic rank calculation, not a tolerance-based numerical test.

The actual onsite and oriented hopping spin coefficients are

\[
 {3\beta\over a},\qquad {-\beta+\beta\gamma_i\over2a},
 \qquad {(-\beta+\beta\gamma_i)^\dagger\over2a}.           \tag{2}
\]

Their sixteen-variable commutant system has rank15 and scalar nullspace.
Equivalently beta and the differences of the two hopping orientations
generate all four Clifford matrices and hence M4. It follows that

\[
 \{X:[X,I_4\otimes T_a]=[X,S_b\otimes I_{15}]=0
                  \text{ for every coefficient }a,b\}
 =\{I_4\otimes\sum_r z_rP_r:z_r\in\mathbb C\}.           \tag{3}
\]

The local gauge commutant before imposing (2) has dimension80; the
simultaneous commutant in (3) has dimension5. The generated complex
coefficient algebra is the direct sum of M_(4d_r), of dimension944.
The central projectors (1) separate its summands. Thus the declared
elementary fiber has no additional three-dimensional multiplicity factor
commuting with these coefficients. The four source spin rows do not form
three flavor copies. Equation (3) does **not** compute the commutant of
the global Hamiltonian alone, exclude composite or spatial excitations,
or exclude an explicitly enlarged source.

The [separate chiral attachment](../standard_model_source_map/KINETIC_ATTACHMENT.md)
declares a three-family tensor factor in a different carrier. Its factor
does not transfer to this sixty-mode source. The finite Schur mechanism in
[`GaugeIrreducibleBorn`](../../Lean/QFT/GaugeIrreducibleBorn.lean) concerns
commutants under its explicit irreducibility hypotheses. Here the literal
matrix ranks supply the finite hypotheses directly.

## 2. Local Dirac mass coefficients and native chirality

Set \(\alpha_i=i\beta\gamma_i\) and
\(\gamma_5=\gamma_1\gamma_2\gamma_3\beta\). All are Hermitian;
\(\beta\gamma_5\beta=-\gamma_5\). Solving the exact linear
equations \(\{M,\alpha_i\}=0\) gives a two-dimensional complex
space. Its Hermitian part is

\[
 M=m\beta+ni\beta\gamma_5,\qquad m,n\in\mathbb R.       \tag{4}
\]

Imposing the actual inversion's spin action, \([M,\beta]=0\), leaves
the single beta direction. Gauge preservation and (1) consequently give

\[
 M_{\rm even}=\beta\otimes\sum_r m_rP_r.                 \tag{5}
\]

These are algebraically allowed parity-even local mass **coefficients**.
The frozen source uses the single onsite coefficient3/a in (2); its
hopping and gauge dynamics determine measured transfer laws. Reading3/a
as a physical rest mass is unjustified. Adding independent m_r to the
action changes that action. Bare off-diagonal quadratic terms between
the inequivalent gauge blocks fail the gauge commutant condition.

The [actual parity construction](../particle_primary_odd_sum_rules/PROOF.md)
proves that beta exchanges the rank-two native chirality projectors
\(\Gamma_\pm=(I\pm\gamma_5)/2\), while preserving each native gauge
representation. The source and polynomial state are invariant under the
full finite inversion. Complete parity-conjugate apparatuses therefore
have equal original classical energy/readout laws after reflecting
positions. The retained reference output is conjugated, not assumed equal.
The onsite term makes native gamma5 nonconserved. The finite result
constrains a parity-equivariant one-handed field assignment; a selected
broken state, different action or proved different field/current dictionary
requires its own construction. It is not a general impossibility theorem
for physical chirality.

## 3. Five native Ward identities and a common scalar constraint

The [source number calculation](../particle_primary_channel_labels/PROOF.md)
shows that the five global block counts N_r commute with all terms of H,
the word cutoff and Gauss projection. Hence \([N_r,\rho]=0\), and the
same is true of all source-energy clock factors. If
\([N_r,A]=\Delta_{A,r}A\), trace invariance under
\(\exp(i\theta\cdot N)\) implies

\[
 \operatorname{Tr}\rho A^\dagger f(H)B=0
            \quad\hbox{unless }\Delta_A=\Delta_B,          \tag{6}
\]
\[
 \operatorname{Tr}\rho A f(H)B=0
            \quad\hbox{unless }\Delta_A+\Delta_B=0.        \tag{7}
\]

This is a scalar source-marginal statement. A charge-coherent reference
or compensated apparatus can carry the opposite increment; its
reference-valued correlation need not vanish. A Wilson link carries
gauge flux but has zero native CAR occupation increment, so a link alone
cannot cancel a nonzero vector in (6) or (7).

Write the unit native vectors in order \((u^c,L,Q,d^c,e^c)\).
For the three conventional bare Yukawa monomials to preserve these
source charges, one conjugate-scalar H-dagger would require

\[
 \begin{array}{c|c|c}
 \text{bare vertex}&\Delta_{H^\dagger}&Q_6(H^\dagger)\\
 Le^cH^\dagger&-e_L-e_{e^c}&-3\\
 Qd^cH^\dagger&-e_Q-e_{d^c}&-3\\
 Qu^cH&+e_Q+e_{u^c}&-3.
 \end{array}                                               \tag{8}
\]

All three have the same required hypercharge, but distinct native
occupation vectors. A single scalar of **sharp native occupation
increment** cannot make all three bare vertices neutral. The statement
does not exclude differently dressed fermions, multiple scalar sectors,
a state breaking the number symmetry, or a different source action.
It diagnoses the literal unbroken source and the stated bare dictionary.
The finite one-Yukawa construction in
[`ElectroweakBreakingComposition`](../../Lean/Screen/ElectroweakBreakingComposition.lean)
uses its declared coefficient, scalar and vacuum premises. Equations
(6)–(8) test whether that premise can be attached to this source without
silently adding a symmetry-breaking interaction.

## 4. An actual conjugate-doublet scalar and dressed odd factor

Take four negative-beta lepton-doublet modes l_(a,s), with weak index
a=0,1 and source-spin index s=0,1, and two negative-beta charged-singlet
modes e_s. These are annihilators in the actual native blocks. All64
occupations of these six modes are retained. With epsilon_01=1, put

\[
 \Phi_a={l_{a0}e_1-l_{a1}e_0\over\sqrt2},\qquad
 D_s=\sum_a l_{as}^\dagger\Phi_a.                         \tag{9}
\]

The code lifts the twelve actual native generators to the full Fock
space and verifies their commutators with every operator in (9). Phi
is a weak anti-doublet, color singlet, Q6=−3, even CAR polynomial with
native increment \((0,-1,0,0,-1)\). Its epsilon contraction is invariant
under every actual48 binary-cube lift and under the three infinitesimal
spin matrices. It has source scalar coefficients and the internal
comparison weights of H-dagger. This supplies the first required vector
in (8) as an actual source operator. It does not supply a Higgs vacuum
or a scalar pole. Its one-point function in the invariant original
source vanishes by (7) and by its gauge charge. A vanishing charged
one-point function in an invariant finite density does not exclude a
Higgs phase defined through gauge-invariant dressed, vector or radial
correlations. The required attachment is such an operational response;
it does not require replacing this exact zero by a nonzero expectation
of the same forbidden local charged operator.

D is odd, weak/color neutral, with Q6=−6 and native increment
\((0,0,0,0,-1)\). It transforms in G1, equivalently the dual spinor
related by epsilon. For an anti-Hermitian spin matrix K the exact
commutator is \([\widehat K,D_s]=\sum_tK_{ts}D_t\).
The epsilon identity \(g^T\epsilon g=\epsilon\) is checked on all48
source lifts. These are finite source rotation statements.

The full-Fock operator identity is

\[
 D_s={1\over\sqrt2}\sum_{a,t,u}\epsilon_{tu}
                           l_{as}^\dagger l_{at}e_u.      \tag{10}
\]

On the filled negative-beta anchor Omega, it reduces to

\[
 D_s\Omega=\sqrt2\sum_u\epsilon_{su}e_u\Omega.            \tag{11}
\]

The exact full64-dimensional positive-effect spectra give

\[
 \|\Phi_a\|=1,\quad\|D_s\|^2=3,\quad
       \left\|\sum_sD_s^\dagger D_s\right\|=9/2.          \tag{12}
\]

Equation (11) is an anchor identity, not replacement of D by a constant
times a bare hole in the original state. On an input with the lepton
modes empty, D vanishes while a bare e hole can be nonzero. Equation
(10) explicitly supplies the density and spin dependence discarded by
such a replacement. D has the same native charge vector as the e^c
hole, so (6) allows their normal correlation. The following exact positive
identity supplies a nonzero original correlation, beyond that selection
rule. Set \(h_s=\sum_u\epsilon_{su}e_u\). Direct expansion, with the
same CAR order on both sides, gives

\[
 \sum_s h_s^\dagger D_s
       =\sqrt2\sum_a\Phi_a^\dagger\Phi_a\ge0.            \tag{12a}
\]

The right side has expectation2 square-root2 on the filled anchor.
Source-marginal faithfulness therefore gives original expectation at
least2 square-root2 times c. This is the gauge-invariant product before
Gauss compression. Inserting P_G between its two charged factors would
erase the product; it is not an undressed charged physical spectral
kernel. Section5 gives a physical compensated comparison.

## 5. Complete compensated readout on the original physical source

A single charged D does not map the physical Gauss carrier into itself.
Choose an interior edge xy, retain a positive-beta e^c creator at x and
the D factor at y, and use the actual charge-six character link and
constant-link projector. Define four rows

\[
 K_{ps}={1\over3}\chi_6(U_{xy})P_{0,xy}
                           c_{x,p}^\dagger D_{y,s},
       \qquad p,s=0,1.                                   \tag{13}
\]

The inside and outside Q6 increments are−6 and+6. Their full five native
increments cancel. The weak currents annihilate D and the e^c factor;
the link character cancels the Abelian endpoint transformations. Thus
each K preserves all Gauss constraints. On the three-word link control,
\(E=\operatorname{diag}(-1,0,1)\), \(L=|1\rangle\langle0|\),
\([E,L]=L\) and \(L^\dagger L=P_0\). The explicit endpoint generators
are \(6N_p-6E\) and \(Q_y+6E\), which commute with (13).
The action is admitted in every primary word cutoff. Gauge-covariant
link compression preserves these identities; its first constant-link
action is exact. Ghost links are not used.

On the complete eight-mode CAR space let
\(E_D=\sum_sD_s^\dagger D_s\). Disjoint even factors commute, and
the graded CAR signs give the exact effect

\[
 \sum_{p,s}K_{ps}^\dagger K_{ps}
       =P_{0,xy}{(2-N_p)E_D\over9}\le I.                  \tag{14}
\]

The two factors in (14) have spectra in[0,2] and[0,9/2]. The complete
instrument is (13) and
\(K_f=\sqrt{I-\sum K_{ps}^\dagger K_{ps}}\).
Every success is globally even; its two regional factors are odd. The
outside charge and link are retained, so this is not an isolated charged
particle creation or half of a justified rest-energy measurement.

With all negative-beta rows filled, positive-beta rows empty and links
constant, Omega is a normalized physical source vector. The four
vectors in (13) obey

\[
 \langle K_{ps}\Omega,K_{p's'}\Omega\rangle
              ={2\over9}\delta_{pp'}\delta_{ss'}.        \tag{15}
\]

The actual selected marginal has \(\rho\ge cP_G\), with the
authenticated finite dyadic c from the
[support bound](../particle_spatial_modular_response/PROOF.md). Therefore

\[
 \Pr(p,s)\ge2c/9,\qquad\Pr(\text{any success})\ge8c/9.  \tag{16}
\]

These probabilities use the original density. No anchor preparation
substitutes for it. The anchor is only a test vector in the positive
trace inequality. The idle reference is kept and no joint pure-state
Loewner inequality is asserted. The complete source-energy instruments
can be appended before and after (13). Their operators are functions of
the same H; the unread first clock fixes the stationary rho and the
second clock is trace preserving. Summing all complete clock outcomes
therefore leaves (16) unchanged. This faithfulness argument alone gives
no energy-window location. The floor c deteriorates with capacity
and is not a useful uniform continuum residue.

The [regional spinor calculation](../particle_primary_regional_fermion/PROOF.md)
supplies an original energy-window bound for the bare compensated e^c
spinor apparatus at member zero: every raw joint spin row has weight
greater than1/16 in its declared pair-transfer window. Section6 attaches
a dressed complete-window bound through a full original-state residual
inequality. Common anchor overlap alone does not transfer that law.

There is a positive original connection to a normalized bare charged
probe. Write A_ps=chi P0 c_p-dagger h_s/2 and B_ps=K_ps from (13).
The bare summed effect is P0(2−N_p)N_e/4≤I. Their exact mixed product is

\[
 \sum_{p,s}A_{ps}^\dagger B_{ps}
  ={\sqrt2\over6}P_0(2-N_p)\sum_a\Phi_a^\dagger\Phi_a
  \ge0.                                                   \tag{17}
\]

All factors on the right commute and are positive. Its anchor value is
2 square-root2/3. Define the eight rows and one failure

\[
 L_{\pm,ps}=(A_{ps}\pm B_{ps})/2,\qquad
 L_f=\left(I-\sum_{\pm,p,s}L_{\pm,ps}^\dagger
                                      L_{\pm,ps}\right)^{1/2}. \tag{18}
\]

The summed success effect is half the sum of the two complete
subinstrument effects, hence at most I. Every row is Gauss invariant,
globally even, native-number neutral, and retains the compensating
outside endpoint. Its complete original records obey

\[
 \sum_{p,s}\Pr(+,p,s)-\sum_{p,s}\Pr(-,p,s)
 =\operatorname{Re}\operatorname{Tr}\rho\sum_{p,s}A_{ps}^\dagger B_{ps}
 \ge{2\sqrt2\over3}c>0.                                  \tag{19}
\]

This is unconditioned contrast on the original selected state, with all
failures and the reference kept. Summing complete initial and final H
clock outputs preserves (19). In an ideal energy eigenbasis the trace
in (19) is a sum of \(p_i\overline{(A_{ps})_{ji}}(B_{ps})_{ji}\).
Thus the compensated bare and dressed probes share at least one nonzero
energy-transition matrix element. Its transfer frequency may be zero.
Its energy location and narrow-band weight are not obtained by this
positive sum rule.

The onsite part of (2) gives [H_on,c_p-dagger D_s]=6c_p-dagger D_s/a.
The constant-to-character link adds63/(4a), so (13) has exact onsite and
electric center87/(4a), the same center as the bare compensated pair.
The hopping commutator is absent from that center. Its spectral response
cannot be replaced by the bare-pair response or by a rest-mass value.

## 6. A positive dressed energy window in the original state

Define the full local residual and negative weak-hole count

\[
 E_s=D_s-\sqrt2h_s,\qquad n_{L,-}^{\rm hole}=4-N_L.
\]

The exact positive matrix \(G=\sum_sE_s^\dagger E_s\) commutes
with N_L. On the full blocks with weak-hole counts0,1,2,3,4 its maximum
eigenvalues are0,1,3,3,4, respectively; both e-mode occupations and all
weak spin states are included. Therefore

\[
 G\le{3\over2}n_{L,-}^{\rm hole}.                         \tag{20}
\]

The coefficient3/2 is attained in the two-hole block. This exact
full-occupation inequality replaces the invalid shortcut of using the
filled anchor alone.

At member zero retain all24 directed cube edges and all four endpoint
spin pairs. Put

\[
 A_{xyps}={\chi P_0c_{x,p}^\dagger h_{y,s}\over2\sqrt{24}},
 \qquad
 B_{xyps}={\chi P_0c_{x,p}^\dagger D_{y,s}\over3\sqrt{24}},
 \qquad r={2\sqrt2\over3}.                                \tag{21}
\]

Each menu has96 successes and its own square-root failure. A is the
actual bare raw-spin apparatus up to epsilon's signed permutation of
the inside spin labels. Its complete unconditioned event is the recorded
transfer interval \(J=[77/8,271/8]\), with original weight
\(p>43/160\) supplied by the regional receipt. The recorded endpoints
include the error allowances of its actual source-energy clocks.

For each edge the residual effect has the same link projection and
outside vacancy factor as (14). Using (20), and the fact that each
origin has exactly three directed neighbors, gives

\[
 \sum_{xyps}(B_{xyps}-rA_{xyps})^\dagger(B_{xyps}-rA_{xyps})
 \le {1\over72}\sum_{x\sim y}n_{L,-}^{\rm hole}(y)
 = {1\over24}\sum_y n_{L,-}^{\rm hole}(y).                \tag{22}
\]

The actual selected density is invariant under source charge conjugation.
This exchanges positive-beta occupations and negative-beta holes within
each native block. The
[selected odd bound](../particle_primary_selected_odd_response/PROOF.md)
gives the full weak-plus-e defect count \(\delta<4.99118\) at this
member. Consequently the expectation of the last sum in (22) is at
most delta/2 and

\[
 \epsilon:=\operatorname{Tr}\rho\sum(B-rA)^\dagger(B-rA)
       \le\delta/48<0.103983<(129/400)^2.                 \tag{23}
\]

Use the **same** complete initial and final clock maps in both menus.
On the actual source-reference purification, write the vector of
success histories as \(f_j(H)A_bf_i(H)|\Psi\rangle\) or
\(f_j(H)B_bf_i(H)|\Psi\rangle\), retaining b,i,j and the reference.
Completeness of the final clock and stationarity of rho under the unread
initial clock imply that the squared norm of the full difference vector
is the left side of (23). Projecting onto the same recorded event J can
only reduce it. The Hilbert-space triangle inequality therefore gives

\[
 \sqrt{\Pr_B(J)}\ge r\sqrt{\Pr_A(J)}-\sqrt\epsilon
 >{61\over125}-{129\over400}={331\over2000},              \tag{24}
\]
\[
 \Pr_B(J)>{109561\over4000000}>{1\over40}.                \tag{25}
\]

The rational margins follow from
\((61/125)^2<(8/9)(43/160)\); the receipt checks them exactly against
the original rational p and delta. There is no extra clock subtraction:
p is the recorded-event probability and (24) compares
its full implemented histories. Different failure operators remain in
the two complete instruments; neither is inserted into the matching
success-history norm. No success or link conditioning occurs.

This is a positive original finite energy response of an actual
scalar-dressed charged regional spinor with its compensator, on the
source whose bare response was computed. It identifies a populated
pair-transfer window, with no rest-frame limit, isolated pole or
single-endpoint energy division. The statement is member-zero specific;
the fixed-edge faithfulness result of Section5 has its separate
all-member domain.

## 7. Candidate conclusions

Every charged-lepton or quark target has an exact internal comparison
row in Section1. The source also supplies actual finite spinor factors,
Gauss-complete compensated probes and the constructive scalar (9).
None of these statements selects electron versus muon versus tau, up
versus charm versus top, or down versus strange versus bottom. The local
coefficient algebra supplies one copy of each internal irrep, while
effective composite or spatial families require their own spectral
identification. The direct bare chirality and Yukawa attachments are
constrained by Sections2–3. A charged species confirmation requires its
original energy/spatial law, effective transformation and interaction
dictionary and a justified scale conversion. No GeV value is inserted
or inferred by this calculation.
