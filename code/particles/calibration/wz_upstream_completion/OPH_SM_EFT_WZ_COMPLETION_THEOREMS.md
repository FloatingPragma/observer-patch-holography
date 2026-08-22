# OPH → renormalized Standard Model → Physical W/Z complex poles

## Mathematical completion theorem stack and computational receipt contract

**Version:** 4.2 (integrated scientific specification)
**Date:** 2026-07-20
**Status:** corrected *draft sufficiency specification*; production receipt schemas and a production aggregate verifier remain open
**Permitted present claim:** `DRAFT_SUFFICIENCY_STACK_DEFINED__SIMULATION_RECEIPTS_OPEN__NO_OPH_NATIVE_POLE_PROMOTION`

---

## 0. Executive result

The downstream strict-one-loop map already constructed is mathematically adequate once it is fed a complete renormalized two-point packet. The missing bridge is upstream. It consists of:

1. a source-emitted local gauge-invariant Standard-Model effective action and active-field census;
2. a source-to-renormalized-parameter and EFT matching map;
3. a precise one-point/tadpole coordinate, preferably the Fleischer--Jegerlehner (FJ) coordinate;
4. complete Yukawa matrices or a certified approximation with a remainder;
5. a fully specified renormalization and chiral-regularization prescription;
6. two independent one-loop realizations of the W and neutral two-point sectors;
7. executable BRST, Slavnov--Taylor, Ward and Nielsen receipts;
8. an analytic-continuation and gauge-invariant amplitude receipt identifying the complex resonance poles; and
9. a source law and, for GeV output, a source-derived operational clock.

The present OPH structural theorem fixes at most a branch of geometry, compact gauge group, representations and a one-Higgs matter pattern. Those data do not determine the coefficients or quantum generating functional. Consequently, the current theorem cannot entail a unique W/Z pole pair. This document proves that non-entailment and states one explicit sufficient augmented branch. No minimality claim is made.

A crucial type correction is:

> **OPH must emit the gauge-invariant EFT action. It need not and should not emit a unique gauge-fixing function.** Gauge fixing, ghosts and gauge-parameter values are calculational coordinates. The simulator chooses a declared BRST family and proves that the complex pole is invariant within it.

The same distinction applies to the counterterm scheme: a scheme is a coordinate chart on the renormalized theory. It must be frozen and transformed consistently; it is not itself an observable selected by OPH.

---

# Part I: Exact non-entailment results

## Theorem 1: Structural gauge data do not determine a local action

Let `T_struct` be any theory whose relevant outputs are:

- a four-dimensional oriented Lorentzian event manifold or an imported local Minkowski chart;
- the global gauge group
  \[
  G_{\rm SM}=\frac{SU(3)_c\times SU(2)_L\times U(1)_Y}{\mathbb Z_6};
  \]
- three copies of the Standard-Model chiral representations;
- one Higgs doublet \(H\sim(1,2)_{1/2}\); and
- the statement that a local gauge-covariant action is admissible.

Then `T_struct` does not determine the numerical local action, the renormalized one-point coordinate, or the W/Z complex poles.

### Proof

For every parameter point

\[
p=(g_s,g,g',m^2,\lambda,Y_u,Y_d,Y_e)
\]

in an open domain with \(g_s,g,g'>0\), \(\lambda>0\), and \(m^2<0\), the same fields and representations admit the power-counting-renormalizable gauge-invariant action written in Theorem 4 below. Two choices \(p_1,p_2\) have identical manifold, gauge group, quotient, representations, generation count and Higgs multiplicity, so they satisfy the same structural sentences. Yet already at tree level

\[
w(p)=\frac{g^2}{4}\frac{-m^2}{\lambda},\qquad
z(p)=\frac{g^2+g'^2}{4}\frac{-m^2}{\lambda}
\]

vary continuously with \(p\). Choosing, for example, the same \(g,g'\) and different ratios \(-m^2/\lambda\) gives different W and Z poles. At one loop the self-energies also depend on \(\lambda\) and on the Yukawa matrices. Hence there are at least two models of `T_struct` with different pole outputs. By elementary model-theoretic non-definability, no unique pole pair is entailed. ∎

## Corollary 1.1: A representation witness is not a generating functional

An exterior-algebra decomposition, anomaly cancellation, a compact-group reconstruction, or a finite response rank can constrain which local operators are allowed. None assigns coefficients to all allowed operators or fixes a quantum state/measure. Therefore none alone supplies the renormalized 1PI functional \(\Gamma\).

## Theorem 2: Normal-form closure does not select a source law

Let \(N:Q\to Q_{\rm nf}\) be an idempotent repair/normal-form map. The pushforward operator on laws, \(\mu\mapsto N_\#\mu\), is idempotent and fixes every law supported on \(Q_{\rm nf}\). Therefore repair confluence cannot select a unique probability law or Euclidean action on the normal-form set.

### Proof

Since \(N\circ N=N\), one has \(N_\#N_\#\mu=(N\circ N)_\#\mu=N_\#\mu\). If \(\mu(Q_{\rm nf})=1\), then \(N(q)=q\) \(\mu\)-almost surely, hence \(N_\#\mu=\mu\). Distinct laws supported on \(Q_{\rm nf}\) are therefore all fixed. ∎

## Theorem 2.1: A MaxEnt state does not determine the Lorentzian dynamics

Let \(\rho\) be a faithful finite-dimensional density matrix. Knowledge of \(\rho\), even together with its MaxEnt representation

\[
\rho=Z^{-1}e^{-K},
\]

does not determine a unique physical Hamiltonian, clock normalization, Lorentzian action or 1PI functional.

### Proof

Every Hermitian \(H\) satisfying \([H,\rho]=0\) makes \(\rho\) stationary under \(e^{-itH}\), and there are infinitely many such \(H\) whenever the Hilbert space is nontrivial. The logarithm \(K=-\log\rho+\log Z\) is a modular/statistical generator, but replacing \(K\) by \(aK+bI\), with a corresponding inverse-temperature/clock rescaling, leaves the same normalized Gibbs state after the thermodynamic parameter is changed. More generally, a single equal-time state does not specify time-ordered correlation functions or the analytic continuation needed for a 1PI action. Therefore an additional transfer/clock/dynamical receipt is necessary. ∎

### OPH consequence

The local-MaxEnt/refinement clause can constrain a finite operator family and produce a Gibbs state once its moments are supplied. It does not, without a transfer/reflection-positivity and continuum packet, identify \(-\log\rho\) with the physical Lorentzian Standard-Model action. A simulator must reject any source bridge that performs this relabeling without the missing theorem.

## Theorem 3: Gauge fixing is not source-identifiable and need not be

Suppose the gauge-invariant action \(S_{\rm inv}\) is fixed. Let \(\Psi_1\) and \(\Psi_2\) be two admissible gauge-fixing fermions connected within a nonsingular BRST family. Then

\[
S_i=S_{\rm inv}+s\Psi_i
\]

are different gauge-fixed actions with the same BRST cohomology and, when the extended Slavnov--Taylor/Nielsen identities hold and the pole is simple, the same physical complex pole.

### Proof

The difference \(S_2-S_1=s(\Psi_2-\Psi_1)\) is BRST exact. BRST-closed observables are unchanged under such a deformation up to the standard source terms required by the Nielsen identity. The determinant identity proved in Theorem 14.2 and the simple-root argument of Theorem 14.3 show that a simple zero of the physical inverse two-point block is independent of every gauge parameter. Thus no unique gauge-fixing action is physically identifiable from the source, and none is required. What is required is an executable proof that the chosen quantization family preserves the physical pole. ∎

---

# Part II: One explicit sufficient augmented source branch

## Definition 4: `OPH+SM-EFT+FJ` source packet

A sufficient source-side augmentation is an equivalence class

\[
\mathfrak G_{\rm src}
=[M,g,\mathcal P_{G_{\rm SM}},\mathcal F,S_{\rm inv},p_F(Q),\mathfrak R_{\rm EFT},\mathfrak L_{\rm src}]
\]

with:

1. **geometry:** an oriented time-oriented four-dimensional Lorentzian spin manifold \((M,g)\), or an explicitly imported Minkowski validation chart;
2. **bundle:** a principal \(G_{\rm SM}\)-bundle compatible with the selected global quotient;
3. **field census:** exactly the declared active gauge, Higgs and chiral matter fields on each EFT interval;
4. **invariant action:** the local action of Theorem 4, including all source-retained operators and explicit exclusions;
5. **renormalized FJ basis:**
   \[
   p_F(Q)=\{Q,g_s,g,g',m^2,\lambda,Y_u,Y_d,Y_e\}_{\overline{\rm MS},\,FJ};
   \]
6. **EFT transport:** thresholds, schemes, beta functions, finite decoupling maps, Jacobians and remainder enclosures;
7. **source law:** either a proved deterministic point law or a target-independent joint law; and
8. **source ancestry:** a hash-pinned DAG with no path from W, Z, Higgs, top, \(G_F\), measured weak mixing, measured VEV, or calibrated proxies to any claimed source output.

Two packets are physically equivalent when related by a bundle isomorphism, a
BRST-exact gauge deformation, or an order-by-order local invertible field
redefinition that stays inside the frozen EFT accuracy.  For a field
redefinition the receipt must include the transformed sources and parameters,
the functional Jacobian and any local anomaly term, every induced operator in
the common term mask, and a bound on discarded higher-order operators.  An
unqualified analytic field redefinition is not an equivalence at finite EFT
order.

This definition is sufficient input structure. It is **not** claimed to follow from current OPH.

## Proposition 4A: Conditional finite-source-to-1PI reconstruction route

A genuine nonperturbative OPH derivation of \(\mathfrak G_{\rm src}\), rather
than an explicit augmentation, would be sufficient if it supplied a cofinal
regulator family and a gauge-invariant Euclidean observable algebra with all of
the following:

1. quotient-intrinsic finite actions or weights \((m_r,S_r)\), not only normal-form maps;
2. exact gauge invariance and a consistent anomaly-free chiral measure at the regulator level;
3. Schwinger functions satisfying the relevant Osterwalder--Schrader hypotheses: Euclidean covariance, permutation symmetry, reflection positivity on the physical observable algebra, tempered regularity, clustering/vacuum uniqueness, and compatibility with the continuum limit;
4. local source couplings sufficient to recover every physical field/operator insertion used in the W/Z current amplitudes, plus a declared BRST gauge-fixed extension for off-shell 1PI coordinates;
5. refinement-compatible convergence of all required Schwinger distributions, with uniform bounds and a nontrivial limit;
6. a renormalization map fixing field normalization, operator basis, scale \(Q\), and the finite-order equivalence class defined above;
7. differentiability of the connected functional and a Legendre construction on the gauge quotient or on a fixed BRST slice, including treatment of zero modes; and
8. a reconstruction and uniqueness theorem for the resulting physical Lorentzian theory and local renormalized jet, modulo precisely the restricted equivalences above.

Then functional differentiation of the resulting renormalized \(\Gamma\) emits D10, D11, RG/matching and two-point pole packets as views of one source object.

### Proof

By the assumed Osterwalder--Schrader reconstruction statement, items 1--5
produce the physical Lorentzian observable theory. Item 6 fixes the otherwise
free renormalized coordinate chart. Item 7 supplies the additional gauge-fixed
coordinate construction needed for an off-shell 1PI functional; its local
derivatives give the potential/Yukawa jet, inverse propagators and vertices.
Item 8 makes those outputs unique in physical content within the declared
finite-order equivalence. Therefore the downstream packets are views of one
reconstructed source object. ∎

### Scope boundary

This proposition assumes the hard reconstruction and uniqueness statements; it
does not prove that the present OPH finite carriers satisfy them. A full
nonperturbative chiral Standard Model construction is much stronger than the
bounded strict-one-loop validation lane. For the latter, a source-emitted
formal perturbative action plus the algebraic-renormalization receipts of Part
VI is sufficient; for a foundational nonperturbative OPH claim, the stronger
finite-measure/continuum hypotheses remain necessary.

---

# Part III: The gauge-invariant local action

## Theorem 4: Canonical local Standard-Model action

Choose the perturbative chart with Minkowski signature \((+---)\), \(Q=T_3+Y\), standard hypercharges, and \(g'\equiv g_Y\) (not the GUT-normalized \(g_1=\sqrt{5/3}\,g'\)). Freeze

\[
D_\mu=\partial_\mu-i g_sG_\mu^AT^A-i gW_\mu^at^a-i g'B_\mu Y,
\qquad e=g s_W=g'c_W,
\]

with the corresponding plus-sign nonabelian commutator term in the field
strengths. Let

\[
\widetilde H=i\sigma^2H^*.
\]

The source action is

\[
S_{\rm inv}=\int d^4x\,\mathcal L_{\rm inv},
\]

\[
\begin{aligned}
\mathcal L_{\rm inv}={}&
-\frac14G^A_{\mu\nu}G^{A\mu\nu}
-\frac14W^a_{\mu\nu}W^{a\mu\nu}
-\frac14B_{\mu\nu}B^{\mu\nu}
+(D_\mu H)^\dagger D^\mu H
-V(H)\\
&+\sum_{\psi=Q_L,u_R,d_R,L_L,e_R} i\bar\psi\gamma^\mu D_\mu\psi\\
&-\Bigl[
\bar Q_LY_dHd_R+
\bar Q_LY_u\widetilde H u_R+
\bar L_LY_eHe_R+\mathrm{h.c.}
\Bigr],
\end{aligned}
\]

with

\[
V(H)=m^2H^\dagger H+\lambda(H^\dagger H)^2.
\]

The packet must explicitly state whether it excludes or retains the QCD and electroweak theta terms, the Weinberg operator, right-handed neutrinos, higher-dimensional operators, additional scalars and vectorlike matter. An omitted operator is not zero by silence: it needs a symmetry/selection receipt or an explicit branch axiom.

### Proof of gauge invariance

The kinetic terms are constructed from curvatures and covariant derivatives. Each Yukawa monomial has vanishing total hypercharge and contracts the color and weak indices to singlets:

\[
\bar Q_LHd_R,
\quad \bar Q_L\widetilde H u_R,
\quad \bar L_LHe_R.
\]

The potential depends only on \(H^\dagger H\). Thus every term is invariant under the local Lie algebra, and the selected hypercharge lattice makes it well defined under the global \(\mathbb Z_6\) quotient. ∎

## Theorem 5: Perturbative anomaly cancellation and conventional SU(2) Witten check

Using left-handed Weyl fields

\[
Q:(3,2)_{1/6},\quad
u^c:(\bar3,1)_{-2/3},\quad
d^c:(\bar3,1)_{1/3},\quad
L:(1,2)_{-1/2},\quad
e^c:(1,1)_1,
\]

the perturbative gauge and mixed gravitational anomalies vanish generation by generation, and the global \(SU(2)\) anomaly vanishes for three generations.

### Proof

Up to common positive trace normalizations,

\[
\begin{aligned}
SU(3)^3:&\quad 2-1-1=0,\\
SU(3)^2U(1):&\quad 2\left(\frac16\right)-\frac23+\frac13=0,\\
SU(2)^2U(1):&\quad 3\left(\frac16\right)-\frac12=0,\\
\mathrm{grav}^2U(1):&\quad
6\left(\frac16\right)+3\left(-\frac23\right)+3\left(\frac13\right)
+2\left(-\frac12\right)+1=0,\\
U(1)^3:&\quad
6\left(\frac16\right)^3+3\left(-\frac23\right)^3+3\left(\frac13\right)^3
+2\left(-\frac12\right)^3+1^3=0.
\end{aligned}
\]

There are four left-handed \(SU(2)\) doublets per generation after color multiplicity, hence twelve for three generations, an even number. The conventional Witten \(SU(2)\) anomaly therefore vanishes. ∎

This calculation does not classify every possible global anomaly of the chosen
\(\mathbb Z_6\) quotient. A stronger global claim requires a separate
global-anomaly or bordism receipt for the precise spacetime and bundle class.

**Receipt consequence.** The anomaly arithmetic is a necessary parent of the algebraic-renormalization receipt. It is not by itself a quantum measure or continuum existence theorem.

---

# Part IV: Electroweak vacuum and full Yukawa packet

## Theorem 6: FJ tree coordinate and electroweak mass matrix

Assume \(m^2<0\) and \(\lambda>0\). Define

\[
v_F=\sqrt{-\frac{m^2}{\lambda}}>0,
\qquad
H=\begin{pmatrix}\chi^+\\(v_F+h+i\chi_3)/\sqrt2\end{pmatrix}.
\]

Let

\[
g_Z=\sqrt{g^2+g'^2},\quad
s_W=\frac{g'}{g_Z},\quad c_W=\frac g{g_Z},
\]

\[
W^\pm_\mu=\frac{W^1_\mu\mp iW^2_\mu}{\sqrt2},\quad
A_\mu=s_WW^3_\mu+c_WB_\mu,\quad
Z_\mu=c_WW^3_\mu-s_WB_\mu.
\]

Then the tree quadratic action has

\[
w=\frac{g^2v_F^2}{4},\qquad
z=\frac{g_Z^2v_F^2}{4},\qquad
m_A^2=0,\qquad
m_h^2=2\lambda v_F^2.
\]

### Proof

The stationarity condition is

\[
\frac{\partial V}{\partial h}\Big|_{h=0}=v_F(m^2+\lambda v_F^2)=0,
\]

and the positive broken solution is the displayed \(v_F\). Expanding \((D_\mu H)^\dagger D^\mu H\) gives a charged mass term \(g^2v_F^2W^+W^-/4\) and the neutral matrix

\[
\frac{v_F^2}{8}(W^3_\mu,B_\mu)
\begin{pmatrix}g^2&-gg'\\-gg'&g'^2\end{pmatrix}
\binom{W^{3\mu}}{B^\mu}.
\]

Its eigenvalues are \(0\) and \(g_Z^2v_F^2/4\) with eigenvectors \(A\) and \(Z\). The second derivative of the potential gives \(2\lambda v_F^2\). ∎

### Mandatory OPH correction

The existing transmutation quantity

\[
v_{\rm chart}(P)=E_{\rm cell}(P)
\exp\!\left[-\frac{2\pi}{\beta_{\rm EW}\alpha_U(P)}\right]
\]

must remain typed `v_chart` until a source theorem identifies it with the positive root \(v_F\) of a canonically normalized Higgs action. Numerical proximity or a tree mass formula is not such a theorem.

## Definition 7: Full Yukawa receipt

At scale \(Q\), a complete receipt contains three complex \(3\times3\) matrices and unitary diagonalizations

\[
U_{uL}^\dagger Y_uU_{uR}=D_u,\qquad
U_{dL}^\dagger Y_dU_{dR}=D_d,\qquad
U_{eL}^\dagger Y_eU_{eR}=D_e,
\]

with nonnegative diagonal entries and

\[
V_{\rm CKM}=U_{uL}^\dagger U_{dL}.
\]

A convenient frozen weak basis is

\[
Y_u=D_u,\qquad Y_e=D_e,\qquad Y_d=V_{\rm CKM}D_d,
\]

with the right-handed down basis fixed by declaration. Running masses are

\[
m_{f}(Q)=\frac{v_F(Q)}{\sqrt2}y_f(Q).
\]

## Theorem 7.1: Basis covariance of the pole

Unitary flavor-basis changes conjugate the fermion kinetic/Yukawa blocks and leave all closed-loop traces and hence W/Z pole positions invariant, provided the CKM matrix and counterterms are transformed consistently.

### Proof

A flavor basis change is an invertible unitary field redefinition with unit Jacobian in perturbation theory. Every fermion loop contains products contracted into traces or the invariant charged-current combination \(V_{\rm CKM}\). Conjugation cancels under the trace, and the determinant zeros of the bosonic Hessian are unchanged under an invertible field redefinition. ∎

### Completeness rule

A physical one-loop pole/width receipt must include all active fermion eigenvalues and CKM data. A `top-only` packet is not complete. Light fermions may be set massless only under an explicit approximation theorem that:

1. names the affected diagrams and thresholds;
2. bounds the omitted real and absorptive contribution;
3. preserves the correct count of open decay channels; and
4. places the bound in a nonstochastic error object.

---

# Part V: BRST quantization and nonlinear gauge family

## Definition 8: BRST differential in the gauge basis

Introduce ghosts \(c_s^A,c^a,c^B\), antighosts \(\bar c^I\), and Nakanishi--Lautrup fields \(b^I\). The classical BRST operator is

\[
\begin{aligned}
sG_\mu^A&=\partial_\mu c_s^A+g_sf^{ABC}G_\mu^Bc_s^C,
&sc_s^A&=-\frac{g_s}{2}f^{ABC}c_s^Bc_s^C,\\
sW_\mu^a&=\partial_\mu c^a+g\epsilon^{abc}W_\mu^bc^c,
&sc^a&=-\frac g2\epsilon^{abc}c^bc^c,\\
sB_\mu&=\partial_\mu c^B,
&sc^B&=0,\\
sH&=i\left(g c^a\frac{\sigma^a}{2}+g'c^B\frac12\right)H,\\
s\psi&=i(c_s^Ag_sT^A+c^agT^a+c^Bg'Y_\psi)\psi,\\
s\bar c^I&=b^I,
&sb^I&=0.
\end{aligned}
\]

It is extended as a graded derivation. Jacobi identities and the representation property imply \(s^2=0\).

## Theorem 8.1: BRST invariance of the invariant action

\[
sS_{\rm inv}=0.
\]

### Proof

BRST acts as an infinitesimal gauge transformation with the commuting gauge parameter replaced by a ghost. Gauge invariance of every term in Theorem 4 therefore gives zero. Nilpotence follows from closure of the gauge algebra and the matter representations. ∎

## Definition 9: Frozen nonlinear gauge family

Use

\[
\begin{aligned}
F^+={}&(\partial_\mu-ie\widetilde\alpha A_\mu
-igc_W\widetilde\beta Z_\mu)W^{+\mu}\\
&+\xi_W\frac g2(v_F+\widetilde\delta h+i\widetilde\kappa\chi_3)\chi^+,\\
F^-={}&(F^+)^\dagger,\\
F^Z={}&\partial\!\cdot Z+
\xi_Z\frac{g}{2c_W}(v_F+\widetilde\epsilon h)\chi_3,\\
F^A={}&\partial\!\cdot A,\\
F_s^A={}&\partial^\mu G_\mu^A.
\end{aligned}
\]

With gauge-fixing fermion

\[
\Psi=\int d^4x\left[
\bar c^-\left(F^++\frac{\xi_W}{2}b^+\right)
+\bar c^+\left(F^-+\frac{\xi_W}{2}b^-\right)
+\bar c^Z\left(F^Z+\frac{\xi_Z}{2}b^Z\right)
+\bar c^A\left(F^A+\frac{\xi_A}{2}b^A\right)
+\sum_{A=1}^{8}\bar c_s^A
 \left(F_s^A+\frac{\xi_s}{2}b_s^A\right)
\right],
\]

define

\[
S_{\rm gf+gh}=s\Psi.
\]

After eliminating the \(b\)-fields, the gauge-fixing terms are

\[
\mathcal L_{\rm gf}=-\frac{F^-F^+}{\xi_W}
-\frac{(F^Z)^2}{2\xi_Z}-\frac{(F^A)^2}{2\xi_A}
-\sum_{A=1}^{8}\frac{(F_s^A)^2}{2\xi_s},
\]

and the ghost action is the corresponding \(-\bar c\,sF\) expression in the frozen sign convention.

## Theorem 9.1: Gauge-fixed BRST invariance

\[
s(S_{\rm inv}+S_{\rm gf+gh})=0.
\]

### Proof

The first term vanishes by Theorem 8.1. The second is \(s^2\Psi=0\). ∎

## Extended gauge-parameter doublets

For every

\[
\eta_i\in\{\xi_A,\xi_Z,\xi_W,
\xi_s,
\widetilde\alpha,\widetilde\beta,
\widetilde\delta,\widetilde\kappa,\widetilde\epsilon\},
\]

introduce a Grassmann source \(\chi_i\) and set

\[
s\eta_i=\chi_i,\qquad s\chi_i=0.
\]

The 45-point sidecar is an eight-parameter **electroweak** stress grid and
holds \(\xi_s=1\).  That is sufficient only for the strict one-loop W/Z lane,
where no internal gluon or QCD ghost occurs.  Any retained mixed QCD correction
must also verify the \(\xi_s\) Nielsen identity (and should vary \(\xi_s\) in a
separately frozen grid).  The QCD ghost vertices must be generated from
\(s\Psi_s\), not supplied by a handwritten table.

These doublets generate the Nielsen identities. They are not new physical fields.

### Simulator rule

The ghost action and all ghost vertices must be generated by symbolic application of \(s\) to the exact hashed \(F^I\). A handwritten ghost vertex list cannot prove completeness.

---

# Part VI: Renormalization and counterterm theorem

## Definition 10: Frozen one-loop renormalization chart

Use dimensional regularization in \(d=4-2\epsilon\) with a declared chiral \(\gamma_5\) treatment (the reference contract uses BMHV) and \(\overline{\rm MS}\) subtraction for the independent invariant parameters. Bare objects are written

\[
\begin{aligned}
X_0&=Z_X^{1/2}X,\\
g_{s,0}&=\mu^\epsilon(g_s+\delta g_s),
&g_0&=\mu^\epsilon(g+\delta g),
&g'_0&=\mu^\epsilon(g'+\delta g'),\\
m_0^2&=m^2+\delta m^2,
&\lambda_0&=\mu^{2\epsilon}(\lambda+\delta\lambda),\\
Y_{f,0}&=\mu^\epsilon(Y_f+\delta Y_f).
\end{aligned}
\]

Gauge-fixing parameters, BRST sources and nonlinear-gauge parameters must also have an explicit bare/renormalized declaration. No universal choice of their finite counterterms is assumed; the generated Slavnov--Taylor restoration fixes the chosen chart.

Because \(v_F\) is derived,

\[
\boxed{
\frac{\delta v_F^{\rm par}}{v_F}
=\frac12\left(\frac{\delta m^2}{m^2}
-\frac{\delta\lambda}{\lambda}\right)
}
\]

at first order. The finite FJ tadpole shift of Part VII is a separate object and must not be merged silently into this parameter counterterm.

## Theorem 10.1: One-loop algebraic renormalizability, conditional form

Assume:

1. the local action of Theorem 4;
2. the anomaly cancellation of Theorem 5;
3. a regularization satisfying the quantum action principle;
4. a complete local counterterm basis of ghost number zero and dimension at most four; and
5. the BMHV (or equivalently explicit) chiral prescription and normalization conditions in the receipt.

Then all one-loop ultraviolet divergences can be removed and the renormalized 1PI functional can be made to satisfy the Slavnov--Taylor identity by local counterterms.

### Proof

Let \(\Gamma=S_{\rm cl}+\hbar\Gamma^{(1)}+O(\hbar^2)\). A regulator may break the Slavnov--Taylor identity by

\[
\mathcal S(\Gamma)=\hbar\Delta^{(1)}+O(\hbar^2).
\]

The quantum action principle makes \(\Delta^{(1)}\) an integrated local polynomial of ghost number one and power-counting dimension at most four. Nilpotence of the linearized Slavnov operator gives the Wess--Zumino consistency condition

\[
\mathcal S_{S_{\rm cl}}\Delta^{(1)}=0.
\]

The cohomology consists of gauge-anomaly classes plus exact terms. The anomaly coefficients vanish by Theorem 5, so

\[
\Delta^{(1)}=\mathcal S_{S_{\rm cl}}\widehat\Delta^{(1)}.
\]

Adding the local finite counterterm \(-\hbar\widehat\Delta^{(1)}\) restores the identity. Ghost-number-zero divergences decompose into renormalizations of the invariant parameters and fields plus BRST-exact gauge-fixing terms, all contained in the declared basis. ∎

### Computational consequence

The counterterm action must be generated by substituting the bare maps into

\[
S_{\rm inv}+S_{\rm gf+gh}+S_{\rm ext}
\]

and expanding to first order. `S_ext` couples antifield/external sources to nonlinear BRST variations. A manually curated counterterm list is not a completeness proof.

The independent checker must verify:

- exact cancellation of every \(1/\epsilon\) coefficient;
- the one-loop linearized Slavnov--Taylor identity;
- the photon Ward conditions;
- the declared finite symmetry-restoring counterterms; and
- equality of the action, term-mask and convention hashes consumed by both pole engines.

---

# Part VII: FJ tadpoles and the two-engine equivalence theorem

## Definition 11: Direct FJ route

The FJ coordinate is the positive renormalized tree-parameter root

\[
m^2(Q)=-\lambda(Q)v_F(Q)^2.
\]

The direct route expands around the bare minimum and retains all explicit tadpole diagrams and their induced contributions to every Green function. Define the finite shift \(\Delta v^{(1)}\) by the sign-safe equation

\[
\boxed{
T_h^{(1)}+\Gamma_{hh}^{(0)}(0)\,\Delta v^{(1)}=0
}
\]

in the declared 1PI sign convention.

## Definition 12: Converted tadpole-free route

Let \(p_L\) denote a separately calculated tadpole-free parameter chart, for example a Landau effective-potential-minimum chart. The complete finite map is

\[
p_L=p_F+\hbar\Delta p^{(1)}+O(\hbar^2),
\]

including field shifts, parameter shifts, mass arguments, counterterms, mixing angles and every induced Taylor term. The map orientation is part of the receipt.

## Theorem 12.1: Exact field-redefinition invariance of pole locations

A local invertible translation of the Higgs integration variable, with the action and sources transformed consistently, does not change the zeros of the physical propagator determinant.

### Proof

The path integral is unchanged by an invertible change of integration variable up to a field-independent Jacobian for a translation. Connected and 1PI functionals are related by the corresponding source/field reparametrization. Their Hessians are conjugate by the invertible Jacobian of the field map, up to terms proportional to the transformed equations of motion that are included by the complete source transformation. Hence the determinant is multiplied by a nonvanishing analytic factor and has the same zero set. ∎

## Theorem 12.2: Strict finite-order reparametrization rule

Let a pole coordinate have the expansion

\[
s(p)=s_0(p)+\hbar s_1(p)+O(\hbar^2)
\]

and let \(p_L=p_F+\hbar\Delta p^{(1)}\). Then, expressed at the same physical point,

\[
\boxed{
 s_{1,F}(p_F)=s_{1,L}(p_F)
 +\Delta p^{a(1)}\frac{\partial s_0}{\partial p^a}(p_F)
}
\]

for this declared map orientation.

### Proof

Substitute \(p_L=p_F+\hbar\Delta p\) into

\[
s_L(p_L)=s_0(p_L)+\hbar s_{1,L}(p_L)+O(\hbar^2)
\]

and Taylor expand:

\[
s_L=s_0(p_F)+\hbar\left[s_{1,L}(p_F)+\Delta p^a\partial_as_0(p_F)\right]+O(\hbar^2).
\]

Equality of the exact physical pole gives the result. ∎

### FJ receipt pass condition

The direct and converted engines must agree, as complex balls, for

\[
s_{W,1},\quad s_{Z,1},
\]

and separately for their real and imaginary parts. Replacing only `v` in tree masses is insufficient; the complete \(\Delta p\cdot\nabla s_0\) and all transformed one-loop terms are mandatory.

---

# Part VIII: Charged and neutral complex-pole mathematics

## Definition 13: Frozen inverse-propagator convention

Let \(\kappa=(16\pi^2)^{-1}\) be the loop-counting factor (distinct from the
Higgs fluctuation \(h\) and from Planck's constant). Use

\[
\Gamma_W^T(s)=s-w+\kappa\Pi_{WW}^{(1)}(s)+\kappa^2\Pi_{WW}^{(2)}(s)+O(\kappa^3),
\]

and

\[
\Gamma_N^T(s)=
\begin{pmatrix}
 s+\kappa\Pi_{AA}^{(1)}(s)+\kappa^2\Pi_{AA}^{(2)}(s)&
 \kappa\Pi_{AZ}^{(1)}(s)+\kappa^2\Pi_{AZ}^{(2)}(s)\\
 \kappa\Pi_{ZA}^{(1)}(s)+\kappa^2\Pi_{ZA}^{(2)}(s)&
 s-z+\kappa\Pi_{ZZ}^{(1)}(s)+\kappa^2\Pi_{ZZ}^{(2)}(s)
\end{pmatrix}+O(\kappa^3).
\]

All \(\Pi\)'s are renormalized, analytically continued complex functions on the declared sheet. Define

\[
D_N(s)=\det\Gamma_N^T(s).
\]

## Theorem 13.1: Strict charged pole coefficients

If \(w\ne0\) and the self-energies are analytic near the reference root on the chosen sheet, the root

\[
s_W=w+\kappa s_{W,1}+\kappa^2s_{W,2}+O(\kappa^3)
\]

has

\[
\boxed{s_{W,1}=-\Pi_{WW}^{(1)}(w)}
\]

and

\[
\boxed{s_{W,2}=-\Pi_{WW}^{(2)}(w)
-s_{W,1}\Pi_{WW}^{(1)\prime}(w)}.
\]

### Proof

Insert the series into \(\Gamma_W^T(s_W)=0\), Taylor expand, and equate powers of \(\kappa\). The coefficient of \(\kappa\) gives the first formula; the coefficient of \(\kappa^2\) gives the second. ∎

## Theorem 13.2: Strict neutral pole coefficients and mixing order

The massive neutral root

\[
s_Z=z+\kappa s_{Z,1}+\kappa^2s_{Z,2}+O(\kappa^3)
\]

has

\[
\boxed{s_{Z,1}=-\Pi_{ZZ}^{(1)}(z)}
\]

and

\[
\boxed{
 s_{Z,2}=-\Pi_{ZZ}^{(2)}(z)
-s_{Z,1}\Pi_{ZZ}^{(1)\prime}(z)
+\frac{\Pi_{ZA}^{(1)}(z)\Pi_{AZ}^{(1)}(z)}{z}
}.
\]

### Proof

Expand

\[
D_N=(s+\kappa\Pi_{AA}+\cdots)(s-z+\kappa\Pi_{ZZ}+\cdots)
-(\kappa\Pi_{AZ}+\cdots)(\kappa\Pi_{ZA}+\cdots).
\]

At order \(\kappa\), the first factor is \(z+O(\kappa)\) and the second gives \(s_{Z,1}+\Pi_{ZZ}^{(1)}(z)=0\). At order \(\kappa^2\), the off-diagonal product first contributes and yields the displayed term after division by \(z\). ∎

**Order firewall.** The product \(\Pi_{ZA}^{(1)}\Pi_{AZ}^{(1)}\) must be computed for matrix identities but must not enter a strict-one-loop root. It is mandatory at strict two loops.

## Theorem 13.3: Strict mass/width series

Let

\[
s_V=m_{V,0}^2+\kappa s_{V,1}+\kappa^2s_{V,2}+O(\kappa^3)
\]

and choose the square-root branch with positive real part and negative imaginary part. Then

\[
\sqrt{s_V}=m_{V,0}
+\kappa\frac{s_{V,1}}{2m_{V,0}}
+\kappa^2\left(\frac{s_{V,2}}{2m_{V,0}}
-\frac{s_{V,1}^2}{8m_{V,0}^3}\right)+O(\kappa^3).
\]

Writing \(\sqrt{s_V}=M_V-i\Gamma_V/2\),

\[
\delta M_V^{(1)}=\Re\frac{s_{V,1}}{2m_{V,0}},\qquad
\Gamma_V^{(1)}=-\frac{\Im s_{V,1}}{m_{V,0}}.
\]

### Proof

Taylor expand the analytic square root around \(m_{V,0}^2\). ∎

The exact square root of a truncated \(s\) is a useful display coordinate but contains kinematic terms beyond the retained loop order. Gauge and engine-equivalence tests must compare strict coefficients.

---

# Part IX: Ward, Slavnov--Taylor and Nielsen theorems

## Definition 14: Extended Slavnov--Taylor identity

Let \(K_\Phi\) be the external source coupled to \(s\Phi\). The renormalized 1PI functional obeys

\[
\boxed{
\mathcal S(\Gamma)+\sum_i\chi_i\frac{\partial\Gamma}{\partial\eta_i}=0
}
\]

where \(\mathcal S\) is the usual antibracket/BRST Slavnov functional including the \(b\)-field terms.

The simulator must derive all two-point ST residuals by differentiating this exact master identity in the declared field/source basis. This is safer than maintaining convention-dependent component formulas by hand.

## Theorem 14.1: Photon Ward protection

On the unbroken electromagnetic branch, the renormalized transverse neutral block satisfies

\[
\Gamma_{AA}^T(0)=0,\qquad \Gamma_{AZ}^T(0)=0.
\]

### Proof sketch

Differentiate the electromagnetic Ward identity with respect to the photon and neutral fields, evaluate at vanishing fields, and use the unbroken \(U(1)_Q\) vacuum. The first identity protects the massless photon; the second prevents a photon--Z mass mixing at zero momentum. In the executable receipt these are recomputed from the full counterterm-completed 1PI functional. ∎

## Theorem 14.2: Matrix Nielsen identity

Differentiating the extended Slavnov--Taylor identity with respect to \(\chi_i\) and the neutral fields gives matrices \(\Lambda_i,\widetilde\Lambda_i\), regular at a simple physical pole, such that

\[
\boxed{
\partial_{\eta_i}\Gamma_N^T
=\Lambda_i\Gamma_N^T+\Gamma_N^T\widetilde\Lambda_i
}.
\]

Consequently

\[
\boxed{
\partial_{\eta_i}D_N
=\operatorname{tr}(\Lambda_i+\widetilde\Lambda_i)D_N
}.
\]

### Proof

The first relation is the two-point projection of the differentiated extended Slavnov--Taylor identity. Where \(\Gamma_N^T\) is invertible, Jacobi's determinant identity gives

\[
\partial_\eta\det\Gamma
=\det\Gamma\,\operatorname{tr}(\Gamma^{-1}\partial_\eta\Gamma)
=\det\Gamma\,\operatorname{tr}(\Gamma^{-1}\Lambda\Gamma+\widetilde\Lambda)
\]

and cyclicity of the trace gives the formula. Analytic continuation extends it to a neighborhood of a simple zero. ∎

## Theorem 14.3: Gauge independence of a simple complex pole

If \(D_N(s_Z,\eta)=0\) and \(\partial_sD_N(s_Z,\eta)\ne0\), then

\[
\frac{\partial s_Z}{\partial\eta_i}=0.
\]

The analogous statement holds for the charged transverse pole.

### Proof

Implicit differentiation gives

\[
\frac{\partial s_Z}{\partial\eta_i}
=-\frac{\partial_{\eta_i}D_N}{\partial_sD_N}\Big|_{s_Z}.
\]

The determinant Nielsen identity makes the numerator proportional to \(D_N(s_Z)=0\). ∎

### Finite-order receipt boundary

This theorem proves gauge independence for the correctly constructed theory. It does not prove that software included every diagram, tadpole, counterterm, ghost, Goldstone, longitudinal block, analytic branch or strict re-expansion. The receipt must verify the finite implementation directly.

---

# Part X: Physical resonance interpretation

## Definition 15: Second-sheet continuation

The scalar integrals and self-energies are first defined by the Feynman \(+i0\) prescription in the Euclidean/upper-half-plane domain. The receipt names every physical cut and defines the continuation path to the resonance sheet. A sheet is identified by its cut-crossing vector, not by an informal string such as `physical`.

## Theorem 15.1: Laurent pole of a mixed propagator

Let \(\Gamma(s)\) be an analytic \(n\times n\) matrix near \(s_p\). Assume

\[
\operatorname{rank}\Gamma(s_p)=n-1,
\]

so its right and left kernels are one-dimensional, spanned respectively by
nonzero vectors \(r\) and \(\ell\):

\[
\Gamma(s_p)r=0,\qquad \ell^\dagger\Gamma(s_p)=0.
\]

Assume also

\[
\ell^\dagger\Gamma'(s_p)r\ne0.
\]

Then

\[
\Gamma(s)^{-1}
=\frac{r\ell^\dagger}
{(s-s_p)\,\ell^\dagger\Gamma'(s_p)r}+O(1).
\]

### Proof

The rank hypothesis makes the complementary block nonsingular after choosing
bases with \(r\) and \(\ell\) as the first right/left directions. Its Schur
complement is scalar and vanishes at \(s_p\). The nonzero derivative hypothesis
makes that zero simple, with derivative proportional to
\(\ell^\dagger\Gamma'(s_p)r\). Inverting the block matrix gives the displayed
rank-one principal part; rescaling either null vector cancels between numerator
and denominator. ∎

## Theorem 15.2: Physical pole criterion for unstable W/Z bosons

A two-point determinant zero is promoted to a physical W or Z resonance pole only if:

1. it is simple and lies on the declared continuation sheet;
2. the Laurent coefficient in Theorem 15.1 is nonzero;
3. at least one BRST-invariant physical scattering amplitude in a charged-current or neutral-current channel has nonzero coupling to that pole; and
4. the same pole appears in that gauge-invariant amplitude after the nonresonant part is separated analytically.

No positivity condition is imposed on the propagator residue at a complex unstable-particle pole.

### Proof

Insert the Laurent expansion between amputated BRST-closed production and decay vertices. If their contractions with \(r\) and \(\ell\) are nonzero, the amplitude has a nonzero simple pole at \(s_p\). Slavnov--Taylor identities cancel unphysical-state contributions, and Theorem 14.3 makes the pole position gauge independent. Conversely, a determinant zero decoupled from every physical channel is not an observed resonance. ∎

This corrects an overly strong “positive-residue pole” requirement sometimes appropriate for stable particles: W and Z are unstable resonances, not asymptotic Hilbert-space states.

## Theorem 15.3: Certified pole uniqueness by Rouché

Let \(C\) be a closed contour on the declared sheet. If a reference analytic determinant \(D_0\) has exactly one zero inside \(C\), no zero on \(C\), and

\[
|D-D_0|<|D_0|\quad\text{on }C,
\]

then \(D\) has exactly one zero inside \(C\), counted with multiplicity.

### Proof

Rouché's theorem. ∎

The contour receipt must use complex balls and prove the strict inequality pointwise on an interval/mesh with an analytic interpolation bound.

---

# Part XI: Matching, source law, covariance and units

## Theorem 16.0: One-loop pure-SM gauge beta coefficients from the census

For three generations, one complex Higgs doublet and the standard hypercharge convention \(g'=g_Y\), define

\[
\frac{dg_i}{d\log\mu}=\frac{b_i}{16\pi^2}g_i^3.
\]

Then

\[
\boxed{b_{g'}=\frac{41}{6},\qquad b_g=-\frac{19}{6},\qquad b_{g_s}=-7.}
\]

### Proof

For a nonabelian group,

\[
b=-\frac{11}{3}C_2(G)+\frac23\sum_{\rm Weyl}T(R_f)+\frac13\sum_{\rm complex}T(R_s).
\]

For \(SU(3)\), the Weyl index sum is \(2\) per generation and therefore \(6\), giving \(-11+4=-7\). For \(SU(2)\), the Weyl index sum is \(2\) per generation and therefore \(6\), while the Higgs contributes \(T=1/2\), giving

\[
-\frac{22}{3}+4+\frac16=-\frac{19}{6}.
\]

For \(U(1)_Y\), include component multiplicities. Per generation,

\[
6(1/6)^2+3(2/3)^2+3(1/3)^2+2(1/2)^2+1^2=\frac{10}{3}.
\]

Three generations give \(10\), and the Higgs contributes \(2(1/2)^2=1/2\). Therefore

\[
\frac23(10)+\frac13\left(\frac12\right)=\frac{41}{6}.
\]

∎

In the GUT-normalized convention \(g_1=\sqrt{5/3}\,g'\), the pure-SM coefficient is \(41/10\). The MSSM triple \((33/5,1,-3)\) is not a pure-SM one-loop packet.

## Reference one-loop matrix RG system on a pure-SM interval

With

\[
T=\operatorname{tr}(3Y_u^\dagger Y_u+3Y_d^\dagger Y_d+Y_e^\dagger Y_e),
\]

\[
H=\operatorname{tr}\left(3(Y_u^\dagger Y_u)^2+3(Y_d^\dagger Y_d)^2+(Y_e^\dagger Y_e)^2\right),
\]

a frozen reference implementation may use

\[
\begin{aligned}
16\pi^2\beta_{Y_u}&=\frac32\left(Y_uY_u^\dagger-Y_dY_d^\dagger\right)Y_u
+Y_u\left[T-\frac{17}{12}g'^2-\frac94g^2-8g_s^2\right],\\
16\pi^2\beta_{Y_d}&=\frac32\left(Y_dY_d^\dagger-Y_uY_u^\dagger\right)Y_d
+Y_d\left[T-\frac{5}{12}g'^2-\frac94g^2-8g_s^2\right],\\
16\pi^2\beta_{Y_e}&=\frac32Y_eY_e^\dagger Y_e
+Y_e\left[T-\frac{15}{4}g'^2-\frac94g^2\right],\\
16\pi^2\beta_\lambda&=24\lambda^2+4\lambda T
-(9g^2+3g'^2)\lambda
-2H+\frac38\left[2g^4+(g^2+g'^2)^2\right],\\
16\pi^2\beta_{m^2}&=m^2\left[6\lambda+2T-\frac92g^2-\frac32g'^2\right].
\end{aligned}
\]

The production checker must nevertheless derive or independently verify the equations from the exact interval census and convention; this displayed system is not a license to reuse it under a different EFT or normalization.

## Definition 16: `EFT-1` matching packet

For every scale interval \([Q_j,Q_{j+1}]\), emit:

- active fields, spins/statistics, representations and multiplicities;
- the complete local action/operator basis;
- scheme and gauge-independent invariant parameter basis;
- beta functions with a contribution-level monomial mask;
- source-derived threshold intervals;
- finite decoupling and scheme-conversion maps \(D_j\);
- Jacobians \(J_j\);
- interval RG/matching/truncation remainders; and
- a final map into `SM_MSbar_FJ(Q)`.

The checker independently derives the one-loop beta coefficients from the census. An MSSM coefficient triple may not be consumed by a pure-SM pole engine without an actual MSSM/equivalent interval and finite matching map.

## Theorem 16.1: Deterministic matching uniqueness and stability

If the beta vector field is locally Lipschitz on each interval, threshold maps are deterministic and Lipschitz, and the ordered threshold list, schemes and term mask are fixed, then one source point determines one matched parameter point. Perturbations are bounded by the product of Grönwall and threshold-map Lipschitz factors plus the declared remainders.

### Proof

Picard--Lindelöf gives uniqueness between thresholds; Grönwall bounds the flow difference. Apply each finite threshold map and iterate. ∎

## Theorem 16.2: Source covariance is not an interval

A non-singleton source enclosure does not determine a covariance matrix.

### Proof

For distinct \(z_1,z_2\) in the enclosure, \(\delta_{z_1}\) and \(\tfrac12(\delta_{z_1}+\delta_{z_2})\) have different covariances but the same support enclosure. ∎

## Deterministic source mode

Only after global root uniqueness, selector uniqueness, exact primitive selection and deterministic matching are proved may the source law be

\[
\mathrm{Law}(z)=\delta(z-z_0),\qquad C_z=0.
\]

Then

\[
C_{x,\rm source}=J_xJ_{\rm match}J_{\rm source}
C_zJ_{\rm source}^{T}J_{\rm match}^{T}J_x^{T}=0,
\]

where the pole output is explicitly real-vectorized as

\[
x=(\Re s_W,\Im s_W,\Re s_Z,\Im s_Z).
\]

This transpose formula is for real covariance. A native complex formulation
must emit both covariance and pseudocovariance.

Root enclosures, branch alternatives, perturbative truncation, scale variation, continuation uncertainty and floating-point rounding remain separately named nonstochastic error objects.

## Stochastic source mode

A target-independent law or hash-pinned ensemble with justified weights is required. Equal branch weights cannot be assigned merely because branches are listed.

If \(F\) is the verified source-to-pole map, the stochastic output is the joint
pushforward law

\[
\mu_{WZ}=F_\#\mu_{\rm src},
\]

not one pole pair. Its covariance, when finite, is a second-moment summary and
is not a confidence region or a support enclosure. Linearized covariance uses
the real Jacobian of \(x\); an exact nonlinear pushforward requires direct
law/ensemble propagation. Certified truncation, continuation, and numerical
remainders remain separate deterministic set-valued objects.

## Theorem 16.3: Operational-unit attachment

Before a source clock closes, the only licensed outputs are

\[
\frac{s_W}{E_\star^2},\qquad \frac{s_Z}{E_\star^2}.
\]

If a source-derived dimensionless clock gap \(\varepsilon_{\rm clk}\) corresponds to the defining physical frequency \(\nu_{\rm clk}\), then

\[
E_\star=\frac{h_{\rm P}\nu_{\rm clk}}{\varepsilon_{\rm clk}}
=\frac{\hbar_{\rm P}\omega_{\rm clk}}{\varepsilon_{\rm clk}}
\]

and

\[
M_V=E_\star\Re\sqrt{s_V/E_\star^2},\qquad
\Gamma_V=-2E_\star\Im\sqrt{s_V/E_\star^2}.
\]

### Proof

This is dimensional factorization by the operational reference energy, using
\(h_{\rm P}\) for Planck's constant and
\(\omega_{\rm clk}=2\pi\nu_{\rm clk}\). ∎

The clock receipt must bind the adopted exact SI value of \(h_{\rm P}\), the
joule-to-eV conversion, the frequency convention, and every rounding step.

The current OPH clock checksum does not instantiate this theorem until its electromagnetic, absolute-electron, QCD/nuclear and atomic parents are source-produced on one no-target DAG.

---

# Part XII: Exact simulator/producer receipt obligations

## 17.1 Claim lanes

Every run must select exactly one lane:

1. `OPH_CHART_ONLY`: current D10/D11 coordinates; no pole interpretation;
2. `EXTERNAL_SM_EFT_VALIDATION`: imported complete SM parameter packet; validates field-theory machinery but is not OPH-native;
3. `OPH_NATIVE_DIMENSIONLESS`: source/matching/FJ/BRST gates pass, clock open; emits \(s/E_\star^2\);
4. `OPH_NATIVE_PHYSICAL`: all prior gates plus clock and physical-amplitude pole interpretation;
5. `TARGET_COMPARISON_ONLY`: immutable post-processing with no source write access.

The lanes are exclusive. A target-comparison process must not share a writable directory or environment variables with source producers.

## 17.2 Raw action receipt

The action producer emits:

- exact field census and global-representation labels;
- canonical kinetic normalization;
- full invariant operator list and exclusions;
- parameter basis and scale;
- metric, Fourier, \(+i0\), hypercharge and mixing conventions;
- action AST and canonical hash;
- source ancestry DAG and blacklist result.

## 17.3 Feynman-rule and diagram-universe receipt

Two independent rule generators are required. At least one must derive vertices by functional differentiation of the action AST. The other must not import its generated rule tables.

Each diagram record contains:

```text
external_block
loop_order
topology_id
internal_fields
vertex_ids
symmetry_factor
fermion_loop_sign
coupling_monomial
nonlinear_gauge_monomial
integral_family
counterterm_or_loop
source_action_hash
```

Completeness is proved against an independently enumerated universe of one-loop 1PI two-point topologies:

- two three-point vertices;
- one four-point/seagull vertex;
- fermion bubble;
- ghost bubble;
- scalar/Goldstone/Higgs bubble;
- gauge and mixed gauge-scalar bubble;
- tadpole-induced self-energy contributions in the direct FJ route; and
- all one-loop counterterm insertions.

A unique tag only prevents duplicate IDs; it does not prove that a topology is present.

## 17.4 Integral receipt

For every scalar integral:

- analytic definition and normalization;
- \(d=4-2\epsilon\) UV pole coefficient;
- cut and sheet label;
- 128-, 192- and 256-bit complex-ball enclosure;
- nesting and radius-shrink proof;
- threshold-distance diagnostic.

Exact algebraic/UV identities must reduce to exact zero. Numerical identities pass by ball containment of zero, not by a decimal coincidence.

A numerical residual \(R=f(x)\) must use a residual-specific validated
extension. In normalized real coordinates, a mean-value enclosure may take the
form

\[
R(B)\subseteq f(x_0)+J_f(B)(B-x_0),
\]

where the interval/ball Jacobian is computed for that residual and all terms
carry the same declared units. The receipt stores the derived input-propagation
radius, rounding radius, and a pre-frozen modeling tolerance separately. Exact
algebraic identities must reduce to exact zero; for a numerical comparison the
residual ball must contain zero and satisfy its derived, dimensionally typed
bound at every requested precision. There is no universal factor such as
`8 * sum(input radii)`. Near a threshold, increase precision or return
`UNCERTIFIED_THRESHOLD_POINT`; do not relax the bound silently.

## 17.5 Direct-FJ engine

Must generate from its own action/rule objects:

- one-point function;
- explicit tadpoles;
- WW, AA, AZ, ZA and ZZ transverse blocks;
- charged/neutral longitudinal, Goldstone and mixed blocks;
- ghosts;
- all counterterms;
- complex absorptive parts.

## 17.6 Converted engine

Must be separately compiled and use a separate integral implementation. It emits the tadpole-free coefficients, the full finite map \(\Delta p\), every tree derivative, and the re-expanded FJ coefficients.

It must not call the direct engine or load its generated expressions.

## 17.7 BRST checker

A small independent checker consumes the canonical action, BRST rules and coefficient records, but never imports the diagram generator. It recomputes:

1. \(s^2=0\) on every field;
2. \(sS_{\rm inv}=0\);
3. \(S_{\rm gf+gh}=s\Psi\);
4. UV-pole cancellation;
5. the linearized one-loop ST identity;
6. the charged and neutral two-point ST projections;
7. \(\Gamma^T_{AA}(0)=\Gamma^T_{AZ}(0)=0\);
8. charged Nielsen identity;
9. neutral matrix and determinant Nielsen identities;
10. direct/converted FJ equality;
11. strict pole series and order masks; and
12. gauge-grid constancy of strict \(s_{W,1},s_{Z,1}\).

The predeclared grid contains the \(3^3\) combinations

\[
\xi_A,\xi_Z,\xi_W\in\{1/2,1,2\}
\]

plus one-at-a-time and mixed nonlinear-gauge stress points. The grid is secondary evidence; the Nielsen identity is proof-bearing.

## 17.8 Physical-pole checker

The physical interpretation receipt must include:

- sheet/cut vector;
- certified contour and argument count;
- simple-zero derivative ball excluding zero;
- left/right null-vector residuals;
- a BRST-invariant physical amplitude channel;
- nonzero production/decay vertex contractions;
- equality of the amplitude and two-point pole balls; and
- the mass/width square-root branch.

Do not require a positive propagator residue for W/Z.

## 17.9 Mandatory mutation suite

Every mutation below must fail for a distinct reason:

1. omit one ghost loop;
2. omit one Goldstone loop;
3. flip one fermion-loop sign;
4. flip the tadpole-shift orientation;
5. replace only `v` without the full \(\Delta p\cdot\nabla s_0\);
6. drop \(AZ\) or \(ZA\) from the neutral matrix;
7. leak \(\Pi_{ZA}^{(1)}\Pi_{AZ}^{(1)}\) into a one-loop root;
8. omit that product from a two-loop root;
9. partially resum the root in a strict-coefficient comparison;
10. switch \(g'\) to GUT-normalized \(g_1\) without conversion;
11. feed MSSM beta coefficients to a pure-SM interval;
12. omit a fermion species/open channel;
13. replace full Yukawas by a top-only packet without a remainder;
14. change the analytic sheet;
15. remove the BMHV symmetry-restoring finite counterterm;
16. attach a one-loop BRST receipt to an order-2.5 result;
17. label `v_chart` as `v_F`;
18. set \(C_z=0\) without source uniqueness;
19. introduce a W/Z/H/top/\(G_F\) ancestor;
20. enable GeV output before the clock gate;
21. set `positive_residue_required=true` for an unstable-boson pole;
22. allow target comparison to write source inputs.

## 17.10 Promotion conjunction

`OPH_NATIVE_PHYSICAL` is true only when all of the following receipt hashes are present and mutually consistent:

```text
EVENT_GEOMETRY_OR_IMPORTED_CHART
SM_EFT_ACTION_1
ANOMALY_AND_PERTURBATIVE_BRST_1
EFT_MATCHING_1
FULL_YUKAWA_1
FJ_DIRECT_1
FJ_CONVERTED_1
FJ_EQUIVALENCE_1
RENORMALIZATION_ST_1
GENERAL_GAUGE_BRST_1
WARD_ST_NIELSEN_1
POLE_SERIES_1
ANALYTIC_CONTINUATION_1
PHYSICAL_CURRENT_POLE_1
SOURCE_LAW_1
SOURCE_COVARIANCE_1
SOURCE_CLOCK_1
NO_TARGET_ANCESTRY_1
```

Every receipt must carry the same action hash, field-census hash, scheme hash, FJ convention hash, perturbative term-mask hash, analytic-sheet hash and source-root hash.

## 17.11 Draft machine-contract inventory

The package supplies ten fail-closed Draft 2020-12 validation instances backed
by nine distinct schema documents (W and Z reuse the pole schema):

```text
SM-EFT-ACTION-1
EFT-MATCHING-1
FULL-YUKAWA-1
FJ-EQUIVALENCE-1
RENORMALIZATION-ST-1
GENERAL-GAUGE-BRST-1
PHYSICAL-CURRENT-POLE-W-1
PHYSICAL-CURRENT-POLE-Z-1
SOURCE-LAW/COVARIANCE-1
SOURCE-CLOCK-1
```

The bundled checker is a **fixed-template specification linter**. It recomputes
the pure-SM one-loop gauge coefficients on the named template interval,
enforces basic interval/source-mode conditions, requires two pole-instance
declarations, compares selected digest strings, and rejects placeholder hashes.
It does not resolve artifacts, recompute their digests or derive the evidence
booleans. The integrated checker therefore refuses scientific promotion
unconditionally, even if every declaration is forged true. The shipped
templates are deliberately valid at the schema level while failing their
candidate conjunction. A future production verifier must resolve immutable
artifacts, recompute their digests and equations, and derive rather than trust
every evidence result, as required by
`../simulator/UPSTREAM_MATH_SIMULATOR_REQUIREMENTS.md`.

---

# Part XIII: OPH paper-side implications

## 18.1 What must be changed in the papers

1. Rename the current electroweak transmutation output `v_chart` everywhere it is upstream of a pole claim.
2. Treat all D10 W/Z pairs as chart coordinates until `SM-EFT-ACTION-1`, `EFT-1`, FJ, BRST and pole-interpretation receipts pass.
3. State that the structural Standard-Model theorem selects a gauge/matter **type**, not the renormalized action coefficients.
4. Move gauge fixing and counterterm choices into a quantization/receipt section rather than presenting them as recovered observables.
5. Quarantine the target-ancestral D11 core.
6. Add the complete-Yukawa/open-channel requirement to every width claim.
7. Replace any “positive residue” requirement for unstable W/Z by the gauge-invariant current-amplitude criterion of Theorem 15.2.
8. Keep the H3/event-base, event-chart, cone, mixed-GNS and clock repairs as parents of an OPH-native spacetime/unit claim. An imported Minkowski validation lane may bypass those paper-side geometry parents but must be labeled imported.

## 18.2 What simulation can and cannot prove

A symbolic/numerical simulator can instantiate:

- a declared action and gauge family;
- perturbative Feynman rules and diagram completeness;
- renormalization and finite identities;
- certified integral enclosures;
- strict complex pole coefficients;
- gauge and FJ equivalence; and
- source-DAG separation.

It cannot by itself prove that OPH selects the action, the Yukawa matrices, the source law or the operational clock. Those require source theorems whose outputs become simulator inputs. A successful external-SM validation proves the computational QFT bridge, not OPH-native source closure.

---

# Final theorem

## Theorem 19: Sufficient OPH-native Physical W/Z complex-pole theorem

Assume:

1. the augmented source packet of Definition 4 is uniquely selected modulo the stated physical equivalences;
2. an artifact-resolving verifier, independent of the producers, recomputes the complete action, Yukawa, matching, FJ, counterterm and BRST receipts rather than trusting producer flags or digest strings;
3. both independent one-loop engines agree on the complex strict coefficients;
4. all Ward/ST/Nielsen and UV identities pass;
5. each charged/neutral determinant has exactly one certified simple resonance-sheet zero inside its declared contour, satisfies the rank-\(n-1\) Laurent hypothesis, and has a nonzero residue in a BRST-invariant physical amplitude;
6. no-target ancestry passes, and the source is typed either as a unique deterministic point or as a target-independent joint probability law with any claimed moments and tail/confidence objects explicitly established; and
7. the source clock passes.

Then the verified deterministic source-to-pole map \(F\) assigns to every
admissible source point \(z\) one physical pair, unique inside the declared
contours,

\[
F(z)=(s_W(z),s_Z(z)),
\]

and the corresponding operational coordinate quadruple

\[
(M_W,\Gamma_W,M_Z,\Gamma_Z),
\qquad s_V=(M_V-i\Gamma_V/2)^2,
\]

through the certified perturbative order and deterministic error ledger. If the
source law is \(\delta_{z_0}\), this gives one source-defined pair \(F(z_0)\).
If the source law is stochastic, the scientific output is instead the joint
pushforward law

\[
\mu_{WZ}=F_\#\mu_{\rm src};
\]

it is not one deterministic pair. Gauge-parameter and direct-versus-converted
FJ independence hold pointwise on the verified source domain.

### Proof

The source packet and law determine either one renormalized invariant parameter
point or a joint parameter law. Matching transports each point uniquely into
the common FJ basis. The algebraic-renormalization and BRST receipts define
finite renormalized two-point blocks. The FJ theorem identifies the two
tadpole realizations. The implicit-function series fixes the strict pole
coefficients, while Ward/ST identities separate the massless photon and
unphysical sectors. The Nielsen determinant identity gives pointwise
gauge-parameter independence. Rouché/argument-principle receipts isolate one
simple zero per declared contour, and the strengthened Laurent hypothesis plus
BRST-invariant physical amplitudes identify those zeros as physical
resonances. The operational clock fixes the scale orbit. Validated interval/
ball maps propagate deterministic error sets. In stochastic mode, measurability
of the verified map gives the pushforward law; covariance is derived from that
law when finite and is not used as an enclosure. ∎

## Current status against Theorem 19

```text
draft_sufficiency_stack_defined = true
blocking_mathematical_corrections_applied_in_integrated_copy = true
production_receipt_schemas_complete = false
production_aggregate_verifier_complete = false
current_OPH_entails_stack = false
external_SM_validation_engine_complete = false
OPH_native_dimensionless_pole = false
OPH_native_physical_GeV_pole = false
```

The next efficient milestone is the `EXTERNAL_SM_EFT_VALIDATION` lane: it can close the direct/converted FJ, counterterm, BRST, Nielsen, analytic-continuation and physical-current pole machinery without waiting for OPH to emit the source action. Once validated, the imported packet can be replaced by a genuinely source-produced packet without changing the downstream map.


## Frozen computational sidecars

- nonlinear gauge grid: `data/nonlinear_gauge_grid_v1.json`
- count: 45
- canonical SHA-256: `6e0265eda5d55e5def430548441777f8809671fbbddeeb1ea4c5f0ace588abfd`
- receipt dependency DAG: `data/receipt_dependency_dag_v4.json`
- DAG canonical SHA-256: `79a8ec98d77286e8063e57b2395ddc36f9507cdc4e0721eb7995b0b8f519dfca`
