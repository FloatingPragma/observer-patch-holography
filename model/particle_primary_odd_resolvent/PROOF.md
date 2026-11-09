# Full complementary resolvent of a physical odd composite

On the primary eight-cell Wilson source, the physical weak-pair/singlet
composite has a controlled full-complement self-energy at interaction
coefficient one. Exact fourth parent moments bound its actual resolvent;
an independent cubic polynomial certifies prepared spectral weight above
\(1/26\) in the source-energy interval \([-722.22108,-714]\). One complete
clock retains more than \(1/27\) in the widened reported interval. The
authenticated finite density floor transfers this complete final-clock
event to a strictly positive original selected-state event. Neither result
isolates a pole or supplies a uniform physical particle residue.

## 1. Same source and a specified physical analysis vector

Use member zero of the [primary contract](../particle_common_contract/PROOF.md):
all native modes and quantum links, the physical Gauss space, degree-one
product cutoff, original selected polynomial density, idle reference,
energy zero and source clock. Let \(\mathcal H_\Delta\) be the complete
reducing physical sector with
\(\Delta=(0,2,0,0,1)\) in native charge order \((-4,-3,1,2,6)\).
Its [full-source sector bound](../particle_primary_sector_threshold/PROOF.md)
admits the outward rational floor

\[
 H_\Delta\ge L I,\qquad L=-18055527/25000=-722.22108.     \tag{1}
\]

The authenticated floor is approximately \(-722.2210758583\), strictly
above (1). No projection to the minimal-defect space is made in (1).

The normalized vector \(\phi\) has the two weak positive-beta modes
in one spin row and the charge-six mode in one spin row occupied at a
single corner. All other matter is the filled negative-beta anchor and
all links are constant. This is the physical local weak epsilon
contraction times the charge-six singlet. Its negative increment is
the corresponding three-hole vector. Every native block is retained.
The [exact composite calculation](../particle_primary_composite_transport/PROOF.md)
gives, with \(b=-711\) and \(K=H_\Delta-b\),

\[
 \langle K\rangle_\phi=0,\quad
 v=\langle K^2\rangle_\phi=132301/735,\quad
 d_3=\langle K^3\rangle_\phi=369869271/42875.             \tag{2}
\]

This vector is an analysis vector unless the declared paired-memory
preparation is applied. That preparation keeps the displaced selected
source state and every reference correlation in one physical memory,
with the compensating opposite increment in another. It changes the
source marginal explicitly. Prepared probabilities are assigned to this
vector. Section 4 separately proves a finite final-energy probability for
the trilinear success output on the original selected density, with its
own density-floor factor.

## 2. Exact fourth parent moment

Separate the complete parent operator as \(H_{\rm par}=H_0+T+B\),
where T is all native Wilson hopping and B is the signed magnetic
multiplication. The vector obeys \(H_0\phi=b\phi\). Its first complete
parent action is in the actual product cutoff. If \(P_{\ell}\) denotes
that cutoff, then

\[
 K^2\phi=P_{\ell}(H_{\rm par}-b)^2\phi,
 \quad\langle K^4\rangle_\phi
       \le\|(H_{\rm par}-b)^2\phi\|^2.                  \tag{3}
\]

The inequality retains every omitted link and matter direction. The
right side is computed exactly, rather than estimated with a global
perturbation norm.

### Conditional Slater and Haar calculation

Fix all link coordinates g. The matter vector \(\phi\) is a Slater
vector. In each native block its occupied projector is \(P_-\otimes I_d\)
at every cell, with one entire positive-beta spin row added at the
chosen cell in the weak and charge-six blocks. For the hole vector one
negative-beta spin row is removed instead. These spin projectors commute
with every native R(g). Let k(g) be the complete one-body hopping.
The finite Slater identity is

\[
 \langle e^{tT(g)}\rangle_\phi
             =\det(I-P+P e^{t k(g)}).                   \tag{4}
\]

The conditional first moment is zero. Unitarity of each native link and
the literal Wilson spin blocks give conditional variance 12 per native
dimension, hence 180 for the full representation. This variance is
independent of g. Bipartite parity gives conditional third cumulant zero.
Expanding the logarithm of (4), the fourth conditional cumulant is

\[
 \begin{split}
 \kappa_4(g)={}&\operatorname{Tr}(Pk^4)
 -4\operatorname{Tr}(PkPk^3)-3\operatorname{Tr}(Pk^2Pk^2)\\
 &+12\operatorname{Tr}(PkPkPk^2)-6\operatorname{Tr}((Pk)^4).
 \end{split}                                               \tag{5}
\]

The one-body matrix is block diagonal in the actual native irreducible
blocks. These traces therefore add across blocks before Haar averaging.
No statistical independence of their representations on a common link
is assumed. Since the conditional variance is the constant 180, the
cross term in the fourth moment is exactly \(3(180)^2\).

Every trace in (5) is a sum of closed four-step walks. Of the 168
oriented cube walks, 120 reduce by adjacent reversals. Their native
link products cancel to identity; their traces contribute the dimension
d. The other 48 traverse a plaquette. Each has a link appearing once,
whose nontrivial native matrix entry has zero Haar mean. They vanish.
The exact spin trace sums per dimension, in the order displayed in (5),
are

\[
\begin{array}{c|ccccc|c}
 &1&2&3&4&5&\kappa_4\\\hline
\text{anchor}&72&36&48&24&15&-18\\
\text{one added spin row}&153/2&81/2&213/4&225/8&147/8&-18\\
\text{one removed spin row}&135/2&63/2&177/4&165/8&99/8&-18.
\end{array}                                               \tag{6}
\]

The code checks both spin choices at all eight cells. The full result is

\[
                   \langle T^4\rangle=3(180)^2-18(15)=96930.
                                                               \tag{7}
\]

An independent literal 32-mode character calculation retains all
occupations and all twelve Haar exponents after two hopping actions.
For the anchor, one added particle and one hole it gives variance 12
and fourth moment 414. This checks the scalar native block and the CAR
signs independently of (4). It does not replace the other native blocks;
their exact contribution is supplied by the matrix-walk proof.

### Magnetic and electric terms

T and B commute on the complete parent link-function space. Odd total
hopping expectations vanish by bipartite parity. The existing exact
native-character integrals give

\[
 \langle B^2\rangle=1/735,\qquad
 \langle B^4\rangle=41/6002500,
 \quad\langle(T+B)^4\rangle=581831145041/6002500.          \tag{8}
\]

In (8) the mixed term is \(6(180)\langle B^2\rangle\), because the
conditional hopping variance is constant. The admitted first-action
packets give

\[
 \|(H_0-b)(T+B)\phi\|^2=2219349/5.
\]

In the cross term with \((T+B)^2\phi\), hopping/loop contributions
have an unused plaquette edge and integrate to zero. The surviving
magnetic cubic cross term, with its actual electric denominators, is
\(-188/18375\). Consequently (3) gives the full-space bound

\[
 \langle K^4\rangle_\phi\le M_4
       =9738478674383/18007500
       \simeq540801.259163.                              \tag{9}
\]

This is the exact parent fourth moment and an upper bound after the
degree-one compression. It is not asserted to equal the compressed
fourth moment. The bound is over 56 times smaller than the inherited
global-norm second-action estimate. For the specified Slater composite
and its translated counterpart the transport remainder therefore gives

\[
 0<|t|\le3/(2M_4)\simeq2.77366\,10^{-6}
 \ \Longrightarrow\
 |\langle\phi_y,e^{-itH}\phi_x\rangle|^2\ge |t|^6/256.  \tag{10}
\]

No late-time transport rate or invariant particle band is inferred.

## 3. Full complementary resolvent at coefficient one

Define the rank-one P=|phi><phi| in the complete physical sector,
Q=I-P, \(D=QH_\Delta Q|_{Q\mathcal H_\Delta}\), and
\(u=QH_\Delta\phi\). Q includes all other minimal states and every
other physical matter/link state; it is not a small packet complement.
Equations (1), (2) and (9) imply

\[
 D\ge L,
 \quad\|u\|^2=v,
 \quad\mu=\langle u,Du\rangle/v=b+d_3/v,
 \quad\operatorname{Var}_{u/\sqrt v}(D)\le
 \tau=M_4/v-v-(d_3/v)^2.                              \tag{11}
\]

Numerically \(\mu\simeq-663.0743148254\) and
\(\tau\simeq527.5560705844\). For real z<L set

\[
 \Sigma(z)=\langle u,(D-z)^{-1}u\rangle,
 \quad G(z)=\langle\phi,(H_\Delta-z)^{-1}\phi\rangle
                      ={1\over b-z-\Sigma(z)}.        \tag{12}
\]

The second identity is the exact finite Feshbach identity. No interaction
expansion or omitted-space invariance is used. Here are explicit
two-sided bounds on its full self-energy. Put
\(e=L-z>0\), \(d=\mu-L>0\), \(d_0=b-L>0\). Then

\[
 {v\over\mu-z}\le\Sigma(z)\le
 U(z)=\min\left\{
 {v(ed+\tau)\over e[d(e+d)+\tau]},
 {v d_0\over v+e d_0},\ {v\over e}\right\}.             \tag{13}
\]

The lower bound is convexity of the scalar inverse applied to the
spectral measure of D in u/sqrt(v). For the first upper bound, take
c=mu+tau/d. The quadratic polynomial interpolating 1/(q-z) at L and
tangentially at c majorizes it on q>=L, since its difference is

\[
 { (q-L)(q-c)^2\over (L-z)(c-z)^2(q-z)}\ge0.           \tag{14}
\]

The quadratic coefficient is positive, so replacing the actual variance
by its cap tau can only increase its expectation. The two-point measure
on L and c with mean mu and variance tau evaluates that expectation,
giving the first term in (13). This proves the Gauss-Radau bound from
the actual source moments, including the variance-upper-bound direction.

For the second term, positivity of H_Delta-L gives, by the positive
Schur complement and a limit from below,
\(\langle u,(D-L)^{-1}u\rangle\le d_0\), with zero spectral weight
on any singular kernel. The function r/(1+er) is increasing and
concave for r>=0. Applying it to r=1/(q-L) in the same normalized
measure gives the second term. This also guarantees U(z)<b-z.
Thus

\[
 {1\over b-z-v/(\mu-z)}\le G(z)
                       \le{1\over b-z-U(z)}.           \tag{15}
\]

All constants are source-derived and the action has coefficient one.
For example, at z=L-8 the enclosure is approximately
\(0.06045815<G(z)<0.07034152\); at z=L-16 it is
\(0.04028075<G(z)<0.04186508\). The receipt retains exact rationals
and more frequencies. These are resolvent expectation bounds, not
positions of poles.

For a complete finite readout, extend
\((L-z)(H_\Delta-z)^{-1}\) by zero on other occupation sectors. It is
a positive physical contraction. Its square-root success operator and
the full complementary square root give success probability
\((L-z)G(z)\) in the explicitly prepared vector. This uses functional
calculus of the original H, not evolution under QHQ or a replacement
clock. Exact finite apparatus services are admitted; no efficient
compiler or executed full spatial histogram is claimed.

## 4. A low-energy spectral window with original prepared weight

Let \(a=3\), \(r=v/a\), \(d_0=b-L\), and
\(p(k)=(-a-k)(k-r)^2\). Equations (2) give the positive expectation

\[
 \langle p(K)\rangle=v^2/a-av-d_3
                 =13236625211/8103375>0.                \tag{16}
\]

On k>=-a the polynomial is nonpositive. On the permitted low interval
[-d0,-a] it is decreasing and at most
\((d_0-a)(d_0+r)^2\). Integration against the complete spectral
measure proves

\[
 \Pr_\phi\{H_\Delta\in[L,-714]\}
 \ge{v^2/a-av-d_3\over(d_0-a)(d_0+r)^2}>1/26.           \tag{17}
\]

One complete original-H clock, with resolution 1/64 and total tail plus
implementation allowance 1/512, therefore reports
\([L-1/64,-714+1/64]\) with original prepared weight greater than
1/27. Every clock output, compensating physical memory and displaced
source/reference correlation is retained. The lower reported endpoint
is above the authenticated global ground upper. The interval can contain
many lines; this argument proves neither an isolated branch nor a
threshold onset. Its input is the specified prepared bare composite,
not the original polynomial density's odd-probe output.

### A finite final-energy event on the original selected density

The [authenticated physical density floor](../particle_spatial_modular_response/PROOF.md#3-physical-gauss-support-and-the-boundary-algebra)
on this same primary source is

\[
 \rho\ge cP_G,\qquad c=2^{-N},\qquad
 N=\lceil2\beta(U_\rho-L_\rho)\rceil+B_{\rm dim}.       \tag{18}
\]

Here \([L_\rho,U_\rho]\) is the declared source preparation spectrum
bound, and \(2^{B_{\rm dim}}\) bounds its ambient dimension. The code
authenticates all three quantities through the common-contract regional
consumer and retains c in factored form. The positive polynomial density
satisfies (18) directly; no approximation error is subtracted from c.

Let \(\Omega\) be the normalized filled-beta, constant-link physical
anchor and A the fixed local weak-pair/singlet creation monomial from the
[complete selected-state apparatus](../particle_primary_selected_odd_response/PROOF.md).
The literal full-CAR calculation gives \(A\Omega=\phi\), up to the
recorded unit phase. The physical anchor is only a matrix-element witness;
it does not replace the original input density. Positivity gives

\[
 A\rho A^\dagger\ge cA P_GA^\dagger
                    \ge c|\phi\rangle\langle\phi|.     \tag{19}
\]

Let \(\widehat E_J\) be the positive event effect of the complete
implemented final source-energy clock used in (17), for its widened
reported interval J. Its prepared expectation exceeds 1/27, including
the finite clock's actual tail and implementation allowance. Applying
(19) to this same effect proves

\[
 \Pr_\rho\{\text{local creation success, final clock in }J\}
 =\operatorname{Tr}(\widehat E_J A\rho A^\dagger)>c/27.  \tag{20}
\]

For the retained uniform menu of four spin choices at each of V cells,
one specified success row has Kraus operator \(A/\sqrt{4V}\). Its
original event weight is therefore greater than \(c/(108V)=c/864\)
at member zero. The complete source-clock error is inside the prepared
expectation before multiplication by c. Subtracting a fixed clock error
after that multiplication would not prove (20).

This experiment records only the final energy after the original-state
creation intervention. No initial-energy recording, energy-transfer law
or success normalization enters (20). The loaded compensating memory,
failure outcomes and untouched quantum reference are all retained. The
operator inequality (19) concerns the source marginal. It does not
assert domination of a pure source-reference density by the prepared
source-reference output. The factored finite-volume c supplies no usable
uniform weight under refinement, isolated line or physical mass.

## 5. Why the retained-penalty complement test does not isolate this band

This limitation concerns a specific full-space inequality family. Let
P_min contain all constant-link, minimal-three-defect vectors in either
odd sector, and Q_min=I-P_min. For the spatially uniform native onsite
retentions r_j and common electric retention s, the inherited local
Schur method has

\[
 H\ge-720-\mathcal L+K_r,\quad
 K_r=3\sum r_jD_j+(105s/2)\sum C_\ell,
 \quad\mathcal L=12\sum c_j+6/7.
\]

The beta vacuum forces each scalar shift
\(c_j\ge d_j/(105(1-s)C_j/2-\delta)\). Add a physical neutral
particle-hole singlet in an unused cell to a minimal odd composite.
This constant-link vector lies in Q_min, with penalty at most
\(9+6r_*\), where r_*=min r_j. It supplies an upper bound on the
minimum penalty available to this comparison, even if every local Fock
loss c_j is solved exactly.

The [exact two-case scalar-loss bound](../particle_primary_branch_resolution/PROOF.md#7-an-exact-obstruction-to-retained-penalty-isolation)
gives

\[
 \begin{array}{c|c|c}
 r_*\text{ range}&\mathcal L\text{ strict lower}&
 \text{comparison-floor strict upper}\\\hline
 [0,1/2]&184218/24871&-720+12-184218/24871<-715.4069\\
 [1/2,1)&2867416/274505&-720+15-2867416/274505<-715.4457.
 \end{array}                                               \tag{21}
\]

These ceilings lie below the authenticated interacting weak and color
trial uppers, approximately -714.7851 and -714.7739. Thus this retained
scalar-penalty comparison cannot put the whole bare complement above
those trial values. Equation (21) is not an upper bound on the actual
complement's energy and does not disprove an isolated physical band.
A substantially better trial or a connected/dressed comparison would
change that test. The resolvent bounds (13) are actual coefficient-one
consequences despite this limitation; their domain z<L supplies no pole
exclusion inside the low window (17). No uniform dimensionless source
gap, relativistic mass shell or GeV mass is proved.
