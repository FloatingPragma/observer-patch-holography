# Source work, modular response and the geometric premises

The selected primary polynomial state admits a quantitative comparison
between physical source work and modular response. For every physical
source unitary, its relative-entropy cost lies between two explicit
multiples of the work supplied on the unchanged source clock. Regional
restriction gives an entropy-work inequality without identifying the
regional modular generator with a local Hamiltonian. A neutral beta-row
flip on member zero supplies strictly positive work, at least 18 source
energy units, and a nonzero modular cost. Proper geometry and an Einstein
area law do not follow from these finite statements.

## 1. The source, support and reused finite results

Use [the primary contract](../particle_common_contract/PROOF.md): the open
Wilson action at coupling one, product link cutoff, full native matter
carrier, physical Gauss projection \(P_G\), selected polynomial density
\(\rho\), canonical conjugate reference with zero Hamiltonian, and
\(\alpha_t(A)=e^{itH}Ae^{-itH}\). Operators and logarithms below act on
the physical support. No Gibbs density replaces the preparation.

[Primary stress](../particle_primary_stress/PROOF.md) supplies the exact
partition \(H=\sum_xh_x\), physical local currents \(i[h_x,h_y]\),
complete regional flux and source-derived norm bounds. Its continuity
law instantiates [CoarseEnergyWard](../../Lean/QFT/CoarseEnergyWard.lean).
[Spatial modular response](../particle_spatial_modular_response/PROOF.md)
supplies the actual regional restriction, boundary-center support,
faithfulness floor and finite entropy first law. Those identities are
reused, rather than inferred from spacetime geometry.

## 2. Exact modular slopes of the selected polynomial

Write the actual preparation on \(P_G\) as

\[
 Q(H)=T_m(y(H))^s,\quad T_m(y)=\sum_{k=0}^m y^k/k!,\quad
 y(E)=\frac{\beta(U-E)}{2s},\quad
 0\le y\le y_*\le\frac12.
 \tag{1}
\]

Here \([L,U]\), \(m\ge1\), \(s\), and \(\beta\) are the authenticated
member's preparation parameters. The source polynomial proof establishes
(1). On the physical carrier \(\rho=Q^2/Z\), so

\[
 K(E)=\log Z-2s\log T_m(y(E)),\qquad
 K'(E)=\beta\left(1-\frac{y(E)^m}{m!T_m(y(E))}\right).
 \tag{2}
\]

The ratio \(y^m/T_m(y)\) is increasing for \(y\ge0\), because
\(mT_m-yT'_m=\sum_{k=0}^{m-1}(m-k)y^k/k!\ge0\).
The exact rational number

\[
 \delta=\frac{y_*^m}{m!T_m(y_*)}<1
 \quad\hbox{therefore gives}\quad
 \beta(1-\delta)\le K'(E)\le\beta.              \tag{3}
\]

At member zero, \(\beta=22272\), \(m=11\),
\(s=536870912\), and \(y_*=1990647/7340032\).
The code computes (3) by rational arithmetic, without evaluating the
large polynomial matrix or approximating its eigenvalues. The exact
member-zero result satisfies \(\delta<1.116\times10^{-14}\).

Integrating (3) from the source spectral upper bound gives the exact
global comparison \(K=\beta H+cI+R(H)\), where
\(0\le R\le rI\), \(r=\beta\delta(U-L)\). After removing the scalar
phase at the midpoint of that range, the source unitaries
\(\rho^{it}=e^{-itK}\) and \(e^{-it\beta H}\) differ in norm by at
most \(|t|r/2\). Their complete channels, including any untouched
reference, therefore differ in output half-trace distance by at most
\(\min\{1,|t|r/2\}\). With the convention \(\alpha_s(A)=e^{isH}Ae^{-isH}\),
the modular observable flow \(\rho^{it}A\rho^{-it}\) corresponds to
source parameter \(s=-\beta t\). The channels dual to these observable
flows use the adjoint unitaries and have the same bound. At member zero
\(r<3.25\times10^{-6}\).
This is a global clock comparison on the actual preparation, without a
regional modular-flow or proper-time identification.

Let \(V\) be any unitary on the physical carrier, applied only to the
source. Define its cyclic source work and modular cost by

\[
 W(V)=\operatorname{Tr}(V\rho V^\dagger-\rho)H,
 \quad D_V=D(V\rho V^\dagger\|\rho)
         =\operatorname{Tr}(V\rho V^\dagger-\rho)K.
 \tag{4}
\]

The last equality uses unitary invariance of the source entropy. The
source-passivity argument in
[field state response](../particle_field_state_response/PROOF.md) applies
because its eigenpopulations decrease with \(E\). For clarity, order
\(E_1\le\cdots\le E_d\) and \(p_1\ge\cdots\ge p_d\).
For any increasing function \(f\), write its diagonal values as a scalar
plus positive increments of upper spectral projectors. Every lower
rank-\(k\) spectral projector has
\(\operatorname{Tr}(P_kV\rho V^\dagger)\le\sum_{i\le k}p_i\).
Summing the increments gives
\(\operatorname{Tr}(V\rho V^\dagger-\rho)f(H)\ge0\).
This also follows from the doubly stochastic matrix \(|V_{ij}|^2\).

Apply this finite rearrangement to the increasing functions
\(K(E)-\beta(1-\delta)E\) and \(\beta E-K(E)\). Then

\[
 0\le\beta(1-\delta)W(V)\le D_V\le\beta W(V).  \tag{5}
\]

Equation (5) concerns the actual polynomial state. It neither assumes
exact Gibbs equilibrium nor identifies modular time with the source
clock. The constant \(\beta\) is its cooling parameter, not a geometric
diamond normalization or a Newton coefficient.

## 3. Regional entropy bounded by local source work

Suppose \(V=V_R\otimes I\) in the ambient owner-cell factorization and
\([V,P_G]=0\). The supported regional density
\(\sigma=\operatorname{Tr}_{\bar R}\rho\) and its physical boundary
sectors are those of the spatial modular construction. The support is
preserved by \(V_R\). Relative-entropy monotonicity under this actual
partial trace and (5) give

\[
 0\le D(V_R\sigma V_R^\dagger\|\sigma)
      \le D_V\le\beta W(V).                    \tag{6}
\]

For the complete two-branch intervention with Kraus operators
\(\sqrt{1-t}I,\sqrt tV\), retain its coin and reference. The regional
source marginal is \(\sigma_t=(1-t)\sigma+tV_R\sigma V_R^\dagger\).
The existing finite first law becomes

\[
 S(\sigma_t)-S(\sigma)+D(\sigma_t\|\sigma)
       =tD(V_R\sigma V_R^\dagger\|\sigma)
       \le t\beta W(V).                        \tag{7}
\]

Both left terms are nonnegative: entropy concavity and equal branch
entropies give the first, and relative-entropy positivity gives the
second. Thus (7) is a floor-independent upper bound on their sum. It can
be combined with the inherited quadratic remainder bound on the actual
regional support. No equality \(K_R=\beta H_R+\text{constant}\) is used.

There is a local quantitative instance on every selected member. Set
\(A=h_0\), with support the actual owner support of that cell energy,
\(\|A\|\le k\), and \(\|[H,A]\|\le c\), where \(c\) is the complete
flux cap for the energy-index region \(\{0\}\). The supporting marginal
region can contain several owner cells; it is distinct from that energy
index set. For \(V=e^{i\theta A}\), stationarity gives \(W'(0)=0\),
and \(|W''(\theta)|\le2kc\). Integration and passivity yield

\[
 0\le W(e^{i\theta A})\le\theta^2kc,
 \quad D(V_R\sigma V_R^\dagger\|\sigma)
       \le\beta\theta^2kc.                     \tag{8}
\]

The receipt evaluates the inherited choice \(\theta=1/1000\). It also
specifies the rational choice
\(\theta_*=1/[4(k+\beta c)]\). Writing \(x=\beta c/k\), the right
side in (8) is \(x/[16(1+x)^2]\le1/64\). This choice has the same
inverse-energy units as the source duration. The small bound comes from a decreasing intervention
angle, not a uniform nonzero continuum response. All coefficients and
supports are regenerated from the same member.

## 4. A nonzero work and modular-response witness on member zero

On the native \(Q_6=6\) block at every cell apply the Dirac matrix
\(\gamma_1\), with identity on the other four native blocks and on all
links. Its second quantization \(V_{\rm flip}\) is a physical even
unitary, preserves all native occupation counts, and commutes with every
Gauss constraint. It anticommutes with the block's beta matrix. Let
\(D_6\) count its positive-beta particles and negative-beta holes.
The exact onsite work is

\[
 \Delta H_{\rm onsite}=\frac6a(2V-D_6).          \tag{9}
\]

Electric and magnetic terms are unchanged. On one native dimension-one
bond, the Wilson hopping one-body matrix has singular values \(1,1,0,0\)
between the endpoints. Its full-Fock Hermitian bond norm is \(2/a\).
The difference of its unitary conjugate and itself has norm at most
\(4/a\), with all occupations retained. With \(B\) internal bonds,

\[
 \frac{6(2V-\langle D\rangle)-4B}{a}
       \le W(V_{\rm flip})\le\frac{12V+4B}{a}.  \tag{10}
\]

The actual [occupation bound](../particle_primary_occupation_bound/PROOF.md)
gives \(\langle D\rangle<5\) at \(a=1,V=8,B=12\). Thus
\(18<W(V_{\rm flip})\le144\). The exact rational receipt retains the
stronger lower bound, approximately 18.05293682. Equation (5) gives a
strictly positive source modular cost with no faithfulness-weight
subtraction and no change of preparation.

The Hermitian work observable \(V_{\rm flip}^\dagger H V_{\rm flip}-H\)
has norm at most 144. The existing complete binary-observable service
with effects \((I\pm\Delta H/144)/2\) therefore has original positive
record probability at least \(9/16\). Its labels \(\pm144\) estimate
mean cyclic work. They are not an asserted trajectory-work spectrum.
A uniform complete-output half-trace error \(e\) subtracts at most \(e\)
from this event bound. Root services and the work-supplying intervention
are declared finite apparatus operations; no large primary run is claimed.

## 5. Retained reference and coin are different entropy domains

All operations above can be tensored with the identity on the canonical
reference and implemented with complete retained records. The entropies
in (4), (6) and (7) are nevertheless source marginals. The canonical
source-reference input is pure. For the witness in (10), positive work
implies that the source marginal changes, so its transformed canonical
joint vector is not collinear with its input. Consequently

\[
 D((V_{\rm flip}\otimes I)|\Psi\rangle\langle\Psi|
       (V_{\rm flip}^\dagger\otimes I)
       \,\|\,|\Psi\rangle\langle\Psi|)=+\infty. \tag{11}
\]

Its source cost in (5) is finite. Equation (11) is an actual-source
counterexample to extending the marginal work bound to the full retained
purifier. No reference is discarded from the program; a partial trace
specifies the observable functional whose bound is proved.

The retained classical coin plus region has entropy
\(H_2(t)+S(\sigma)\); its derivative at zero is singular. The coin plus
full source-reference has entropy \(H_2(t)\). Their entropies are not
the differentiable marginal functional in (7). No expression here is
identified with generalized gravitational entropy or proper area.

## 6. An exact boundary obstruction to the periodic coframe transfer

[Geometric response](../particle_geometric_response/PROOF.md) declares a
coframe family using periodic unitary translations \(T_i\), and obtains
the flat identity \((I-T_i)^\dagger(I-T_i)=2I-T_i-T_i^\dagger\).
The actual primary matter shifts are open partial isometries. Their link
entries are unitary in the parent coordinate algebra, but they have no
matter image across an exterior face. Put \(\Pi_i=I-T_i^\dagger T_i\).
Then direct substitution of those open shifts into the periodic formula
instead gives

\[
 (I-T_i)^\dagger(I-T_i)=2I-T_i-T_i^\dagger-\Pi_i,
 \quad h_{\rm primary}-h_{\rm naive\ coframe}
       =\frac{\beta_D}{2a}\sum_i\Pi_i.         \tag{12}
\]

There are \(L^2\) missing-face sites per direction. In the full native
carrier beta has thirty positive and thirty negative modes per cell.
Thus the second-quantized difference has norm exactly
\(45L^2/a\), equal to 180 at member zero. This value also holds on the
physical Gauss carrier: filling all negative-beta rows at every cell,
with constant links, attains the negative extremum. Each filled native
15-row block is its top exterior power; the native connected gauge
representation has determinant one because its generators are traceless.
The corresponding filled positive rows attain the other extremum.
This is the [physical singlet construction](../particle_primary_capacity/PROOF.md)
in the beta basis. Both vectors are physical, and no invariance of their
span under hopping is asserted. Product link compression leaves the
boundary term unchanged.

Equation (12) does not contradict the periodic theorem's domain. It
rules out its unmodified use on this open primary. One can declare a
boundary-completed coframe family, but that declaration and its response
must be specified and tested. Neither choice supplies proper geometry.

## 7. Exact scope of the Einstein theorem interfaces

The accompanying [premise map](premises.json) records the owning theorems
and their actual domains. The finite source energy Ward identity, supported
modular first law and (5)-(8) are established on this primary state and
clock. Their mathematical outputs do not supply the following inputs of
[the matched Einstein theorem](../matched_einstein_limit/PROOF.md):
observer-derived proper metric and clocks; action-matched symmetric
Lorentz stress; physical diamond modular kernel; fixed-proper-volume area
response with a positive common coefficient; generalized-entropy balance;
and a common shrinking-ball refinement with errors \(o(\ell^4)\), matched
curvature, Ward and Bianchi limits.

`EinsteinPremiseLink.lean` proves the entropy pairing on a finite simplex
while taking the split \(K=2\pi B+Z\) as an argument. Defining \(B=K/(2\pi)\)
would supply an algebraic name, without the physical diamond kernel.
`SmallBallPremises` explicitly takes stationarity, continuum kernel and
fixed-volume area as fields. `QuantitativeSmallBall.lean` proves tensor
decoding and conditional error transport. `SourceOrderEinsteinComposition`
uses its source Gram quotient and supplied tensors; it does not identify
the primary Wilson state or prove the physical geometric premises.
The boundary and purifier counterexamples above concern concrete invalid
attachments, rather than defects in those conditional theorems.

These gravity premises are separate from a first finite source-energy or
GeV comparison. No gravity input is required by the source's finite mass
readout theorem.
