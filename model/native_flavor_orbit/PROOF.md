# Flavor orbits of the native-current memory

## One current and the complete embedding group

The source is the twelve-current spin/exterior realization in
[`native_matter_attachment.py`](../common_source/native_matter_attachment.py),
with its selected fifteen-state internal space, three-copy family factor,
canonical kinetic form `I45`, and both components of each Higgs invariant line.
The auxiliary spin lift and the selected family factor are declared architecture.
This theorem does not attach that architecture to the distinct radial-defect
source or identify its singular values with interacting physical poles.

The source-computed primitive charge has five distinct eigenvalues
`-4,-3,1,2,6`, of multiplicities `3,2,6,3,1`. It is an exact linear combination
of the twelve currents. A commuting operator therefore has no off-diagonal
charge blocks. On each block, the complete current-commutator equation has
rank `8,3,35,8,0` in `9,4,36,9,1` complex unknowns. The identity lies in each
kernel, so each kernel consists exactly of the scalar matrices.
`algebra.commutant_certificate` recomputes these ranks over the exact source
number field, using all twelve currents.

It follows that the complete commutant on the fifteen-state space is `C^5`.
Tensoring the three-dimensional gauge-trivial family factor gives

\[
 \{G_a\otimes I_3\}'=\bigoplus_{q\in\{-4,-3,1,2,6\}}M_3(\mathbb C).
 \tag{1}
\]

Its unitary group is exactly `U(3)^5`. These are all embeddings that preserve
the specified currents and canonical kinetic inner product on this same
45-dimensional space. An isometric embedding into a larger space, a changed
current identification, or a changed kinetic form requires its own map.

Let `S_q` denote the five family unitaries. In the all-left CAR convention the
symmetric pairing kernel transforms by `K -> S^T K S`. The directed source
Higgs block is `conjugate(Y)`, and its opposite block is its transpose. Thus
the physical coefficient transformation is

\[
 Y_u'=S_1^\dagger Y_u\overline{S_{-4}},\qquad
 Y_d'=S_1^\dagger Y_d\overline{S_2},\qquad
 Y_e'=S_{-3}^\dagger Y_e\overline{S_6}.                 \tag{2}
\]

The same `S1` acts on both quark components because they belong to one weak
doublet. Tests check (2) on all 45 rows and both Higgs components with complex,
nondiagonal family unitaries. A control rotating only one weak component
fails the actual current equations.

Consequently every channel singular value is invariant. The left squared
matrices `Hu=Yu Yu†` and `Hd=Yd Yd†` undergo one common similarity by `S1`.
Their spectral-projector overlaps, hence squared mixing moduli when their
eigenvalues are distinct, are invariant. The polynomial

\[
 I_{\rm CP}=\operatorname{tr}[H_u,H_d]^3                 \tag{3}
\]

is invariant without a nondegeneracy assumption. In particular, the fixed
feedback law that attaches the same real positive `Y` to all three channels
has identical singular spectra and `Hu=Hd` in every admissible embedding.
Its hierarchy restriction cannot be removed by a change of family frame.
For degenerate eigenvalues, individual basis vectors inside their eigenspaces
carry no separately invariant mixing angle.

The complete Higgs-line calculation leaves three arbitrary complex `3x3`
coefficients, of total complex dimension 27. Equation (1) describes changes
of family basis; it does not select these coefficients. The exact complex
control in `coefficient_freedom_certificate` has nonzero (3) while satisfying
the same invariant Higgs tensor factorization. Thus the CP constraint below
is a consequence of the specified memory coefficients, not of the native
gauge representation.

## The complete current-gated memory law

[`native_gauge_family_records.py`](../common_source/native_gauge_family_records.py)
specifies a different readable source law. Its trace-orthonormal current frame
has central positive Casimir `C`; the complete current instrument has Kraus
operators `M0=cos(theta sqrt(C))` and
`Ma=i F_a sin(theta sqrt(C))/sqrt(C)`. All thirteen outcomes are retained.
Its binary coarse graining records jump or no jump, with jump effect
`sin^2(theta sqrt(C))`. Both Higgs components supply invariant left and right
domains, on which the effect is scalar. The two-leg binary probabilities
therefore factorize for every state on those domains, including entangled
states.

The unchanged native C/E raw record has 128 labels and event factors
`p_native=1/12`, `q_native=1/6`. The fine two-leg current instrument has 169
outcome pairs, giving 21,632 fine native/current labels per old memory input.
The source records their binary jump/no-jump grouping. This combined binary
record has 512 labels per old memory input, including zero-weight labels.
Both groups contain every fine outcome assigned to them; none is selected
away. Its effective flags have

\[
 p={1\over12}\sin^2(\theta\sqrt{c_L}),\qquad
 q={1\over6}\sin^2(\theta\sqrt{c_R}),                  \tag{4}
\]

where `(cL,cR)` is `(25/8,4)` for up, `(25/8,3)` for down, and `(9/8,3)`
for charged leptons. The angle, coupling, independent legs, reset port 2,
memory controller and coherent feedback pulse are declared implementation
data. Reading them determines their consequences without selecting which
implementation describes nature. Their addition is a change from passive C/E
repair alone. The retained default is `theta=pi/5`; the existing geometric
catalogue also permits `pi/3` and `pi/2`.

The twelve-state memory applies the source half-turn `gC` on flag C, then
`gE` on a joint flag. Other records and core destinations reset to port 2.
Every old input is retained. Put

\[
 a=p(1-q),\quad b=pq,\quad D=1-a^2+ab.
\]

For `0<=p<1, 0<=q<=1`, the reset part has probability `1-p>0`, so it is a
strict contraction on zero-sum signed measures. The stationary probability
is unique. Direct substitution into the complete transition gives

\[
 \pi_2={1-a\over D},\quad \pi_3=a\pi_2,\quad
 \pi_{10}={ab\pi_2\over1-a^2},\quad \pi_9=a\pi_{10},
 \qquad\pi_i=0\text{ otherwise}.                     \tag{5}
\]

Here `D>0` and `1-a^2>0`; all entries are nonnegative and sum to one.
Zero stationary entries do not remove transient memory inputs or records.
The source geometric rank-one projectors satisfy `P9=P10`. The stationary
unconditional pulse coefficient is therefore exactly

\[
 Y(p,q)={(1-a)P_2+a(1-a)P_3+abP_{10}\over D}.          \tag{6}
\]

The four source duty policies retain every branch. After tracing the current
auxiliary legs and using the declared binary grouping, the reduced
memory/feedback Kraus factors are
`sqrt(w_r)|F(m,r)><m| tensor exp(-i h chi_D(r) X tensor P_F(m,r))`,
with identity on disabled pulses. Summing all 6,144 reduced-channel labels
gives identity on the entire old-memory and feedback-target input space.
The scalar effect identity on each current leg justifies this reduction,
including states entangled with a reference; it does not identify the
fine current-output states with scalars. The underlying fine instrument has
259,584 labels when all twelve old-memory inputs are included. The first
derivative at zero pulse step, averaged in (5), is

\[
 Z_D=\sum_{m,r}\pi_mw_r\chi_D(r)P_{F(m,r)}.            \tag{7}
\]

These are first-order stationary coefficients. The finite-step channel need
not be unitary, and identifying an accumulated finite-memory protocol with
the averaged generator requires its separate averaging error bound.
The actual source theorem in `native_family_record_memory.py` supplies such
a bound on its stated product-initial-state protocol.

For duties `always,C,E,CE`, the exact traces are `1,p,q,pq`; the rank on
`0<p,q<1` is respectively `3,3,2,2`. In particular
`Z_C=Y-(1-p)P2`, which generally differs from `pY`. Duties E and CE have
projector support precisely `{2,10}`. Positivity follows from (7); the rank
is the dimension spanned by its positively weighted projector ranges.
There is no division by the duty probability and no conditioning on success.

## Source-generated mixing and the CP constraint

For the actual real geometric frame, direct multiplication of (6) gives

\[
 [Y(p,q),Y(p,r)]_{02}=
 { (\sqrt5-1)p^2(q-r)
   [p^2(1-q)(1-r)+q+r-1]\over10 D(p,q)D(p,r)}.        \tag{8}
\]

For `0<p<=1/12` and `q,r` in `[0,1/6]`, the square-bracket factor is at most
`1/144+1/3-1=-95/144`. Thus unequal right gates give noncommuting matrices.
At `theta=pi/5`, the up and down left gates agree, whereas their right
arguments obey `0<sqrt(3)pi/5<2pi/5<pi/2`. Strict monotonicity of sine on that
interval proves `q_up>q_down`. Both coefficients are positive definite.
Their squares cannot commute: commuting positive matrices have commuting
positive square roots. The quark mixing is therefore nontrivial.

Every coefficient (7) is real symmetric, for every real angle and every
retained real record weight in this projector class. This statement also
holds for arbitrary real linear combinations of the real source projectors.
Consequently `Hu,Hd` are real symmetric, their commutator is real
antisymmetric, and its cube has zero trace. Equation (3) vanishes exactly.
For the ordinary three-family charged current this is the CP-conserving
real-coefficient class; no inference about a separate strong theta angle
or an added Majorana sector is made. Equation (2) proves that complex family
basis changes cannot evade this obstruction. A different coherent control
that changes the coefficient class requires a different source law.

For the unconditional policy, positive coefficients have singular-value sum
one in each channel. Policy C has equal up/down sums because their p agrees.
The event-only policies retain an exact zero singular direction per channel.
These are joint action constraints, not fits to observed particle pole masses.

## Certified spectrum and mixing calculation

The calculation retains all four duty policies at the source default angle.
Each spectrum has three charged coefficients, all three family values in
each coefficient, and every row of the 45-dimensional neutral-Higgs kernel.
The joint action assembly uses the declared static mean-feedback service in
`native_record_yukawa_field.py`. The source `native_family_protocol.py`
distinguishes that service from its finite correlated-memory experiment.
No equality of these two finite processes is assumed or proved here.
The source Higgs partial isometries have unit nonzero modulus. Hence each up
or down singular value occurs six times, each lepton value twice, and three
minimal neutral directions have zero singular value. The ranks are 42 for
always/C and 28 for E/CE. The resulting nullities are 3 and 17.

Every reported interval uses rational arithmetic. Machin's identity
`pi=16 atan(1/5)-4 atan(1/239)` and alternating arctangent sums enclose pi.
Integer square-root brackets enclose the rational Casimirs. On `[0,2]`, the
sine Taylor sum through degree `2n-1` has absolute remainder bounded by
`x^(2n+1)/(2n+1)!`; all source arguments lie in this interval. Outward decimal
rounding enlarges these enclosures. Rational interval evaluation then bounds
every coefficient entry. Exact support ranks supply its zero roots.

For each positive root, deterministic rational bisection uses interval LDL
congruence of `Z-xI`. A pivot interval excluding zero supplies its sign;
the number of negative pivots is the inertia. The certified counts at each
root bracket are `i` and `i+1`, and brackets are disjoint. No numerical
eigensolver admits or selects an interval.

For distinct eigenvalues, the exact spectral projector is

\[
 P_i(Z)=\prod_{j\ne i}{Z-\lambda_jI\over\lambda_i-\lambda_j}.
\]

Interval evaluation of `tr(P_i(Zu)P_j(Zd))` encloses every squared mixing
modulus. Intersecting with `[0,1]` uses the independent projector-probability
identity. This includes the isolated zero directions in the rank-two
policies. Independent high-precision diagonalization of the literal parent
coefficients checks the enclosures but does not define their admission.

## Complete finite histories and the generator grading

Write `X6=X tensor I3`, `Z6=Z tensor I3`, where X and Z are the Pauli flip
and sign matrices of the declared two-component feedback register. Complex
conjugation is in these source bases, and `Theta=Z6 K`. Every real source
projector, every real pulse step h, and every disabled pulse obey

\[
 U_j(h)=e^{-ihX\otimes P_j},\qquad
 [U_j,X_6]=0,\qquad Z_6\overline{U_j}Z_6=U_j.                 \tag{9}
\]

Products preserve these identities in their original order. They hold for
every retained history of all four duties, without a stationary-input or
commuting-projector assumption. The containing unitary group is
`U=I tensor B+iX tensor C` with real B,C; in the X eigenbasis it is
`diag(V,conjugate(V))`, with `V in U(3)`. This group description does not
assert reachability of every V by the specified stochastic controller.

If -1 is absent from the spectrum, the principal generator
`H=i Log(U)/T`, `T>0`, is Hermitian. Conjugation respects the principal
logarithm on this domain and reverses the factor i. Therefore

\[
 [H,X_6]=0,\quad Z_6\overline H Z_6=-H,\qquad
 H=I\otimes iA+X\otimes S,\quad A^t=-A,\quad S^t=S,
 \quad A,S\text{ real}.                                   \tag{10}
\]

To derive the last statement, write a commuting Hermitian H as
`I tensor C+X tensor D`. The antiunitary identity makes C purely imaginary
and D real. Hermiticity gives the stated transpose symmetries. The sufficient
condition `sum_j abs(h_j)<2` follows from
`norm(U-I)<=sum_j norm(U_j-I)<=sum_j abs(h_j)<2`, which excludes -1.

The reduced memory/feedback Kraus operators are real memory matrix units
times the histories (9), with nonnegative real square-root weights. Thus
the joint channel intertwines Theta conjugation and left and right
multiplication by X6. The identities hold on correlated memory/feedback
inputs and with an untouched reference, on which the antiunitary uses
ordinary conjugation. Retaining record labels in a real orthogonal basis
preserves them. The actual fine current Kraus operators contain complex
current matrices; covariance on those output legs requires a separate
current-leg antiunitary. The scalar-effect reduction proved above justifies
tracing them, without postselection. A correlated input need not define a
state-independent reduced channel from its memory marginal alone.

If these reduced channels converge to a continuous unitary group near time
zero, their channel identities pass to the limit. Differentiation gives
`[H,X6]=0` and `H+Theta H Theta^-1=c I`; subtracting `c I/2` gives (10).
This statement concerns a continuous group, not an arbitrary logarithm of
one endpoint. Indeed, with
`A0=[[0,1,0],[-1,0,0],[0,0,0]]`, the Hermitian matrix
`2pi X tensor iA0` exponentiates to I but has an imaginary off-diagonal
block. Its eigenvalues in units 2pi are `-1,-1,0,0,1,1`. The principal
logarithm of that endpoint is zero.

## Retained one-step and finite-history noise

For an independent stationary memory preparation, group identical pulses
into weights `(w0,w2,w3,w10)`, with w0 the disabled-pulse identity. The exact
complete record sums give, with `pi2=(1-a)/D`,

\[
\begin{array}{c|cccc}
 &w_0&w_2&w_3&w_{10}\\ \hline
\mathrm{always}&0&\pi_2&a\pi_2&ab/D\\
C&1-p&\pi_2-(1-p)&a\pi_2&ab/D\\
E&1-q&q-ab\pi_2&0&ab\pi_2\\
CE&1-b&b-ab\pi_2&0&ab\pi_2
\end{array}                                                \tag{11}
\]

These are sums of nonnegative source weights and sum to one. Direct source
evaluation gives `tr(P_i P_j)=1/5` for distinct i,j in `{2,3,10}`. For pulse
step epsilon and `c=cos(epsilon)`, consequently

\[
 t_0=\operatorname{tr}(U_j)/6=(2+c)/3,\qquad
 t_1=\operatorname{tr}(U_i^\dagger U_j)/6=(7+8c)/15.          \tag{12}
\]

The normalized Choi state is `J=sum_h w_h |U_h>><<U_h|/6`. With r=1-w0,
its exact purity defect is

\[
 \delta=1-\operatorname{tr}J^2
 =2w_0r(1-t_0^2)
  +(r^2-w_2^2-w_3^2-w_{10}^2)(1-t_1^2).                    \tag{13}
\]

For `0<p,q<1` and `epsilon not in 2pi Z`, delta is strictly positive for
all four duties: always has distinct positive active weights, and each
other duty has positive idle and active weights. For any unitary channel V,
using the unhalved diamond norm,

\[
 \|\Phi-\mathcal V\|_\diamond
 \ge 2(1-\sqrt{1-\delta})\ge\delta.                       \tag{14}
\]

Indeed the normalized maximally entangled input gives a Choi output whose
overlap with any pure state is at most `lambda_max(J)<=sqrt(tr J^2)`.
The two-outcome measurement onto that pure state gives trace distance at
least twice one minus that overlap. This is an explicit reference-assisted
witness, not a claim about every pure input. The receipt encloses (13) by
rational intervals at epsilon=1/4 for all twelve source-channel/duty pairs.

For N>=2, start with any memory state independent of the feedback/reference.
The first instrument dephases old-memory coherences. Keep every subsequent
record and its original order. Pairwise contributions to `1-tr(J_N^2)` are
nonnegative, allowing a lower bound from two grouped histories. For always
duty choose all C=false, and C=false until the last step followed by
C=true,E=false. Their weights are `(1-p)^N` and `(1-p)^(N-1)a`. The first
N-1 steps reset to port2 and give the same prefix; the final pulses are U2
and U3. Therefore

\[
 1-\operatorname{tr}J_N^2
 \ge 2a(1-p)^{2N-1}(1-t_1^2)>0.                           \tag{15}
\]

For duties C,E,CE use the common inactive prefix C=E=false, of per-step
weight `c0=(1-p)(1-q)`. Compare the all-inactive history with the same
prefix followed by an active event of respective weight
`d_C=a,d_E=(1-p)q,d_CE=b`. From port2 their last pulses are U3,U2,U2,
including the actual core reset after the joint event. Thus

\[
 1-\operatorname{tr}J_N^2
 \ge 2d_Dc_0^{2N-1}(1-t_0^2)>0.                           \tag{16}
\]

Grouping here bounds the complete sum; it does not select the physical
output. Equations (14)--(16) exclude equality to any unitary target channel
for each finite N in the stated domain. They do not give an asymptotic
noise floor. At epsilon=2pi k all pulses are I; boundary gate choices can
also make the channel unitary. The always-duty weak-pulse averaging bound
in the parent is consistent with (15). Neither a unitary dilation retaining
the apparatus nor reconstruction of a coefficient by tomography is excluded.

## Channel logarithms and the separate coherent correction

The complete twelve-memory classical/quantum operator space has dimension
`12*6^2=432`. Output memory ports0 and1 are always empty, so its one-step
map has rank at most360. A finite matrix exponential is invertible; this
complete map has no matrix logarithm. Removing these outputs is a change
of domain. Even on the four recurrent ports, identity feedback blocks
realize the transition matrix with characteristic polynomial

\[
 (\lambda-1)(\lambda+a)(\lambda^2-a^2+ab).                 \tag{17}
\]

For a>0 the eigenvalue -a excludes the standard principal logarithm whose
domain avoids the negative real axis. Other branches and continuous-time
embeddability on a changed domain are not classified here.

A useful coherent diagnostic is the Hilbert--Schmidt projection of a
random-unitary channel onto traceless Hamiltonian commutators. In dimension
d=6 it is `-i[H_coh,.]`, with

\[
 H_{\rm coh}={i\over2d}\sum_h w_h
 [\overline{\operatorname{tr}U_h}U_h
       -\operatorname{tr}(U_h)U_h^\dagger].                \tag{18}
\]

This follows by taking the superoperator Hilbert--Schmidt pairing with
`-i[B,.]` for traceless Hermitian B; two such commutators have pairing
`2d tr(BH)`. It is basis independent and depends only on the channel, but
is a tangent projection, not its logarithm or exact evolution generator.

There is a nonzero diagonal correction even after summing all records.
Start at port2, use always duty, `p=1/12,q=1/6`, and
`cos(epsilon)=99/101,sin(epsilon)=20/101`. The four two-step projector
histories `(2,2),(2,3),(3,2),(3,10)` have weights
`(1-a)^2,a(1-a),a(1-b),ab`. In (10) for (18), putting
`t_d=(10+16c+4c^2)/5`, direct multiplication gives

\[
 A={t_dc(c-1)\over6}
 \{a(b-a)[P_3,P_2]+ab[P_{10},P_3]\},\qquad
 A_{02}={1656589(\sqrt5-1)\over899081864640}>0.              \tag{19}
\]

The separate iA cannot be discarded while identifying the finite process
with a static mass service. The off-diagonal S remains real; if it is
assigned the same canonical Dirac coefficient role, its invariant (3)
vanishes. This conditional statement is not a physical CP theorem for the
entire history algebra. In the X basis the X=+1 block `K=S+iA` is an
arbitrary Hermitian matrix in the containing symmetry class. For example,

\[
 K_u=\operatorname{diag}(1,2,3),\qquad
 K_d=\begin{pmatrix}1&1&i\\1&2&1\\-i&1&3\end{pmatrix}
 \quad\Longrightarrow\quad
 \operatorname{tr}[K_uK_u^\dagger,K_dK_d^\dagger]^3=-34560i.
                                                               \tag{20}
\]

Both lift to generators (10). Assigning K as a Yukawa matrix would also
replace the other X block `iA-S` by `-iA-S`, and hence changes the source
operator. Equation (20) refutes an inference of universal zero CP from
the symmetry alone; it is neither a source-reachability claim nor a
physical Yukawa assignment. The antiunitary in (9) is a spectral grading
in the auxiliary feedback space, not an identified physical CP operation.

## Scope

The result classifies the full same-current canonical-kinetic family orbit
and the stated complete-record memory coefficients. It supplies distinct
charged spectra and nontrivial real mixing from an admitted current-gated
source, and excludes charged-current CP violation in that real coefficient
class. The finite-history theorems retain a separate coherent correction,
give exact nonunitarity witnesses, and distinguish principal history
generators from channel logarithms. They supply no physical clock,
interacting continuum, repair-protected
particle identification, matching to measured poles, or selection of the
implementation describing nature.
