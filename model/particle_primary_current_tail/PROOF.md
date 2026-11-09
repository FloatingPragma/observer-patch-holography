# Common local-phase histories and energy-form control

Use the unchanged open Wilson/Gauss family, polynomial density, idle
conjugate reference and source Hamiltonian clock of the
[primary contract](../particle_common_contract/PROOF.md). Write
\(a=2^{-n}\), \(R=1/a\), \(L=2R^2\), and
\(K=H-E_0\ge0\) on the actual physical carrier. Subtraction of the true
finite ground energy occurs only in inequalities; no source clock or
recorded energy zero is changed. The selected state is stationary and
\(\operatorname{Tr}\rho K\le\epsilon=a/16\).

The [local phase construction](../particle_primary_local_fluctuation/PROOF.md)
supplies physical unitaries
\[
 F_X=\pi a\sum_y f_R(y-X)(N_y-30I),\qquad U_{X,\sigma}=e^{i\sigma F_X},
 \quad\sigma=\pm1.                                           \tag{1}
\]
Only the apparatus tent coefficients wrap on the finite address list;
there is no wrapping source bond. The same retained origin controls the
whole history. All estimates below also hold on each individual origin
block, including boundary profiles. Number phases commute with the product
link cutoff and Gauss projection, so parent identities descend to the
actual carrier.

## 1. An actual uniform double commutator

Onsite, electric and magnetic terms commute with \(F_X\). For each
Hermitian hopping bond \(T_b\), including its \(1/a\), put
\(\delta_{Xb}=\pi a\Delta_b f_R\). Then
\[
 [F_X,[H,F_X]]=-\sum_b\delta_{Xb}^2T_b.                       \tag{2}
\]
The full-Fock norm of a bond is at most \(30/a\). The exact tent sum in
one direction is
\[
 D_R=\sum_y(\Delta_1 f_R(y))^2
       ={2(2R^2+1)^2\over9R^3}.                            \tag{3}
\]
Discarding the absent wrapping bonds only decreases the sum of squared
coefficients. For \(R\ge2\), using \(\pi^2<10\),
\[
 \|[F_X,[H,F_X]]\|\le90\pi^2 aD_R
       <900aD_R\le{2025\over2}.                            \tag{4}
\]
Indeed \(aD_R=(2/9)(2+R^{-2})^2\le9/8\). At \(R=1\) the tent
occupies a single cube corner, with exactly three source bonds. Its
bound is \(90\pi^2<900\), also below (4). This estimate is independent
of volume and link cutoff and applies to every physical input vector.
It is not a bound on \(\|[H,U_X]\|\).

## 2. Semibounded energy controls the unitary form domain

For finite self-adjoint \(K\ge0,F\) assume
\(\|[F,[K,F]]\|\le C\). For a normalized vector \(\psi\), let
\(e(t)=\langle e^{itF}\psi,K e^{itF}\psi\rangle\). Then
\(e(t)\ge0\) and \(|e''(t)|\le C\). Taylor's inequality gives
\[
 0\le e(t+s)\le e(t)+s e'(t)+Cs^2/2.
\]
The quadratic on the right is nonnegative for every real \(s\), hence
\(|e'(t)|^2\le2Ce(t)\). Integrating the corresponding inequality for
\(\sqrt{e(t)+\eta}\) and taking \(\eta\downarrow0\) proves
\[
 \|K^{1/2}e^{itF}\psi\|
 \le\|K^{1/2}\psi\|+|t|\sqrt{C/2}\,\|\psi\|.              \tag{5}
\]
It also holds for unnormalized vectors. Equations (1)--(4) therefore give
energy square-root drift at most \(45/2\) for either sign. Source time
translates \(\alpha_s(U_X)=e^{isH}U_Xe^{-isH}\) and adjoints obey the
same bound for arbitrary source times \(s\).

For a word \(W\) of \(k\) such unitaries, triangle inequalities give
\[
 \|K^{1/2}W\psi\|\le\|K^{1/2}\psi\|+{45k\over2}\|\psi\|,
 \qquad
 \operatorname{Tr}\rho W^\dagger K W
       \le(\sqrt\epsilon+45k/2)^2.                         \tag{6}
\]
The trace inequality follows by applying the vector inequality to the
columns of \(\rho^{1/2}\), with the Hilbert--Schmidt triangle inequality.
Thus it uses the selected state directly. No variational vector is inserted
in the expectation.

There is a sharper averaged bound for complete random-sign histories.
For one pulse,
\[
 \tfrac12\sum_{\sigma=\pm1}U_{X,\sigma}^\dagger K U_{X,\sigma}
       =K+A_X,\qquad
 A_X=\sum_b(\cos\delta_{Xb}-1)T_b,\qquad
 \|A_X\|\le2025/4.                                        \tag{7}
\]
For the member-zero corner, \(\|A_X\|\le180\). Averaging every sign
and retaining every sign record in a \(k\)-pulse history gives final
mean excess at most \(\epsilon+k(2025/4)\), even with arbitrary source
time evolutions between pulses. Each original sign sequence has
probability exactly \(2^{-k}\). No branch is normalized or discarded.
A static classical origin and arbitrary idle quantum reference can be
included throughout, since (5) and (7) hold on every origin block and
after tensoring an identity.

## 3. Joining the existing filtered probes

The [fixed triangular probes](../particle_primary_field_algebra/PROOF.md)
are contractions \(B\) with \(\|[H,B]\|\le d=2/T\). Their adjoints
have the same bounds. A finite functional-calculus estimate gives
\[
 \|[K^{1/2},B]\|\le2\sqrt d.                               \tag{8}
\]
For completeness, use
\(K^{1/2}=\pi^{-1}\int_0^\infty t^{-1/2}K(K+t)^{-1}\,dt\).
Below a splitting point \(s\), bound the commutator with
\(K(K+t)^{-1}\) by two. Above it use
\(\|[K(K+t)^{-1},B]\|\le d/t\). The two integral bounds are
\(4\sqrt s/\pi\) and \(2d/(\pi\sqrt s)\). Taking \(s=d/2\)
gives \((4\sqrt2/\pi)\sqrt d<2\sqrt d\). When \(d=0\), the
commutator vanishes by functional calculus. No lower spectral gap is used.

A mixed word of \(k\) local phase unitaries and filtered contractions with
fixed widths \(T_j\) consequently obeys
\[
 \|W\|\le1,\qquad
 \operatorname{Tr}\rho W^\dagger K W
 \le(\sqrt\epsilon+D_W)^2,\qquad
 D_W={45k\over2}+\sum_j2\sqrt{2/T_j}.                       \tag{9}
\]
Time translates and adjoints are allowed at every position. Finite sums
with absolute coefficient sum \(c\) and common drift bound \(D\) have
norm at most \(c\) and final-energy bound
\(c^2(\sqrt\epsilon+D)^2\). These estimates allow one common state on
the algebra containing both menus.

## 4. Complete original response and common time dynamics

For stationary \(\rho\), the original transfer measure of a contraction
word is
\[
 \mu_{n,W}(d\omega)=\sum_{i,j}p_i|W_{ji}|^2
                        \delta_{E_j-E_i}(d\omega).
\]
It has total mass at most one, exactly one for a unitary word. Since
\(|E_j-E_i|\le(E_j-E_0)+(E_i-E_0)\), (9) implies
\[
 \int|\omega|\,d\mu_{n,W}
       \le(\sqrt\epsilon+D_W)^2+\epsilon.                 \tag{10}
\]
This is uniform for each fixed word. The finite downward-record theorem
`Lean/QFT/GroundRecordTail.lean`, with the actual same-source ground
floor and selected mean excess, also gives
\[
 \mu_{n,W}(({-\infty},-\eta])\le\epsilon/\eta
       \quad(\eta>0).                                    \tag{11}
\]
Thus the source-marginal weighted time modulus is
\[
 \|\alpha_t(W)-W\|_{\rho_n}^2
 \le\min\{4,\,2|t|[(\sqrt\epsilon+D_W)^2+\epsilon]\}.      \tag{12}
\]
This is not operator-norm continuity. Initial complete energy measurement
preserves the stationary source marginal and generally changes its
correlations with the retained idle reference. No stationarity of the
joint purification is assumed.

Use the bounded sequence C* algebra and norm-null quotient construction
from the owning filtered-algebra proof. Adjoin the countable local phase
menu and its rational source-time translates. Choose one subsequence for
all mixed rational-time star-word moments. Positivity and normalization
pass to the limit, with the null ideal and GNS left ideal distinguished
as in that proof. The rational time action is isometric and leaves the
selected limiting state invariant.

Equation (12), followed by finite polynomial approximation, makes its
canonical rational-time GNS implementers strongly continuous. They
therefore extend uniquely to a strongly continuous real time group fixing
the cyclic vector. For each dense word the spectral measures converge to
its continuous limiting measure. Equation (11) removes their negative
spectral support, so the canonical generator is nonnegative on the
whole GNS space. This uses the general GNS interface in
`Lean/QFT/ColimitGNSRepresentation.lean`; it does not assert that that
module formalizes this entire source-sequence construction. The coherent
finite-colimit hypotheses of `ColimitSelectedStateGNS.lean` are not assumed.

Real multitime correlations are controlled as well. For a contraction
factor \(A\) with energy-root drift \(r_A\), a vector of energy-root norm
\(e\), and a time change \(h\),
\[
 \|(\alpha_h(A)-A)\psi\|
 \le\sqrt{2|h|}\,(2e+r_A\|\psi\|).                        \tag{13}
\]
This follows by inserting \(e^{ihK}Ae^{-ihK}\), using
\(\|(e^{ihK}-I)\xi\|\le\sqrt{2|h|}\|K^{1/2}\xi\|\),
and applying (5) or (8). Apply (13) to each suffix of a mixed word and
telescope; every suffix has the bound (9). The resulting uniform
Hölder modulus extends the common rational-time word limits to arbitrary
real times and makes convergence uniform on compact time menus.
Fixed finite coherent Kraus words, retained finite pointer densities and
fixed complete instruments formed by uniform functional calculus inherit
these moment limits. The growing full reference and origin register have
no asserted trace-norm limit. There is no effective subsequence rate.
The real time group acts on the represented von Neumann algebra; no
point-norm continuous action on the unfiltered sequence C* algebra is
claimed. Spatial continuity and canonical momentum of these unfiltered
phase fields require their own attachment.

## 5. An actual ultraviolet operator-norm obstruction

Let \(\Omega_0\) be the filled negative-beta determinant and constant
quantum links on the full primary physical carrier. It is the physical
[expansion anchor](../particle_primary_spectral_scale/PROOF.md), not the
selected preparation or an eigenvector of \(H\). Each centered cell
number annihilates it, so \(U_X\Omega_0=\Omega_0\). The literal Wilson
hopping on one bond creates orthogonal particle/hole packets of squared
norm one per internal row, summed over both directions and both beta
spin components. The five native representation blocks contain fifteen
internal rows in total, giving total norm squared fifteen.

Different bonds have orthogonal occupation/character packets. A pair
acquires opposite phases on the two orientations, with common squared
modulus \(|e^{i\delta}-1|^2=2(1-\cos\delta)\). Non-hopping terms
commute with the unitary. Every resulting native link character is present
at every selected cutoff. Therefore, with \(B=3L^2(L-1)\),
\[
 {1\over V}\sum_X\|[H,U_X]\Omega_0\|^2
       ={30B k_n\over a^2V},\qquad
 k_n=\sum_X(1-\cos\delta_{Xb}).                            \tag{14}
\]
The exact member-zero value is 180. For every member the tent lower bound
\(k_n\ge2a^2D_R\) gives
\[
 {30B k_n\over a^2V}\ge{60B D_R\over V}
 =40\left(4R+{2\over R}-{1\over R^3}-{1\over2R^5}\right)
 \ge160R={160\over a}.                                    \tag{15}
\]
Consequently no capacity-independent absolute bound on every
\(\|[H,U_X]\|^2\) exists for this apparatus. This does not lower-bound
the selected-state weighted second moment: the anchor is absent from
that expectation. It does not prove collapse or survival of the source.
It explains why (5), which controls energy forms rather than absolute
commutator norms, supplies a useful common-history estimate.

## 6. The physical boundary

Uniform form control does not imply an inelastic limiting response.
Take the finite two-level control
\(H_a=\operatorname{diag}(0,1/a)\),
\(F_a=\pi\sqrt a\,\sigma_x\), its exact ground state and
\(a=4^{-n}\). Then
\(\|[F_a,[H_a,F_a]]\|=2\pi^2\) is uniform, while the unitary response is
\[
 (1-\sin^2(\pi\sqrt a))\delta_0+
       \sin^2(\pi\sqrt a)\delta_{1/a}.
\]
Its mean tends to \(\pi^2\), its second moment diverges, and its weak
limit is entirely elastic. Its phase gates have roots-of-unity eigenvalues.
This tests the implication and is not the primary source.

For the primary unitaries a selected energy-weighted tail bound, populated
fixed-energy contrast, or a suitable weighted current covariance estimate
is required for nonelastic survival. The present construction supplies
common time histories and positive GNS energy on the enlarged source
algebra. It supplies no physical momentum for the unfiltered fields,
Lorentz covariance, surviving interaction, species or GeV rest mass.

## 7. A uniform finite source-marginal modular response

The [actual polynomial modular bound](../particle_primary_gravity_premises/PROOF.md#2-exact-modular-slopes-of-the-selected-polynomial)
applies to the same selected density, without replacing it by an exact
Gibbs state. Its preparation has Taylor order \(m\ge1\), scaled argument
\(0<y_*\le1/2\), cooling parameter \(\beta\), and
\[
 \delta={y_*^m\over m!T_m(y_*)}
 \le{y_*\over m+y_*}\le{1\over3}.
 \tag{16}
\]
The first bound retains just the top two positive Taylor terms. The
selected cooling plan cannot choose order zero: its spectral width and
\(\beta\) have product greater than one, so the order-zero state-error
upper bound is one, exceeding the permitted target at most \(1/4\).
For every source unitary \(U\), the existing finite passivity proof gives
\[
 \beta(1-\delta)W(U)\le D(U\rho U^\dagger\|\rho)
       \le\beta W(U),\qquad
 W(U)=\operatorname{Tr}(U\rho U^\dagger-\rho)H.             \tag{17}
\]

Retain every classical origin/sign pair \(r=(X,\sigma)\) with its
original probability \(p_r=1/(2V)\), and compare the two source-marginal
classical--quantum states
\[
 \Omega_{CQ}=\bigoplus_r p_r U_r\rho U_r^\dagger,
 \qquad \tau_{CQ}=\bigoplus_r p_r\rho.
\]
Their block relative entropy is exactly
\(D(\Omega_{CQ}\|\tau_{CQ})=\sum_rp_rD(U_r\rho U_r^\dagger\|\rho)\).
The already source-derived local-phase mean satisfies
\(0<m_*\le\sum_rp_rW(U_r)\le M_*\), where
\(m_*=(16/9)h_*>0.06149\) and \(M_*=14175/32\). Equations (16)--(17)
therefore supply the capacity-uniform finite response
\[
 {2m_*\over3}\le {D(\Omega_{CQ}\|\tau_{CQ})\over\beta}
       \le M_* .                                         \tag{18}
\]
The receipt also gives the sharper actual-member lower bound
\((1-\delta)m_n\), with the exact polynomial \(\delta\).

This relative entropy concerns the source marginal and its complete
classical record. The arbitrary idle quantum reference is retained by the
apparatus, but is traced out for this particular entropy comparison. A
relative entropy between the full purified input and output can be
infinite and is not equated with (18). The parameter \(\beta\) is the
declared cooling normalization, with no geometric temperature or diamond
interpretation. No regional area law, local gravitational kernel or
Einstein dynamics follows. Nor does (18) pass the unbounded energy mean
to the GNS limit: the ultraviolet countercontrol in section 6 explains
why that requires an additional energy-weighted tail bound.
