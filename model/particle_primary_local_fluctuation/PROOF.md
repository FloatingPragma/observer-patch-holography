# Local Weyl records with a uniform energy moment

The source is the selected open Wilson family with sixty CAR modes per
cell, all quantum link coordinates, product word cutoff, zero-external-charge
Gauss constraints, polynomial state \(\rho_n\), source Hamiltonian \(H_n\)
and idle conjugate reference. No physical parameter or state is changed.
Write \(a=2^{-n}\), \(R=1/a\), \(L=2R^2\), \(V=L^3\), and
\(\epsilon=2^{-n-4}\). A fixed energy threshold below is in source-clock
units, or physical units with one fixed declared clock conversion.

The [selected hopping estimate](../particle_primary_energy_tightness/PROOF.md)
supplies
\[
 -a\langle H_{\rm hop}\rangle_{\rho_n}/V\ge h_n\ge h_*>0.03459 .
 \tag{1}
\]
This uses a variational energy upper bound to constrain the actual
selected state. It does not replace that state by the trial vector.
The source also satisfies
\(\operatorname{Tr}\rho_n(H_n-E_{0,n})\le\epsilon\).

## 1. A bounded local phase apparatus

Let \(N_y\) count all sixty native occupations at cell \(y\), and set
\(n_y=N_y-30I\). These gauge-invariant operators commute and have norm
at most thirty. For an integer displacement \(k\), define the periodic
tent \(t_R(k)=\max(1-d_L(k,0)/R,0)\), where \(d_L\) is circular distance
on the list of \(L\) coordinate addresses. Put
\[
 f_R(y-X)=\prod_{j=1}^3t_R(y_j-X_j),\qquad
 U_{X,\sigma}=\exp\!\left(i\sigma\pi a\sum_y f_R(y-X)n_y\right),
 \quad \sigma=\pm1.                                      \tag{2}
\]
The classical address \(X\) is uniform on the \(V\) source cells; the
sign is uniform and independent. Both are retained and used throughout
the same history. The complete instrument has Kraus operators
\[
                   K_{X,\sigma}=U_{X,\sigma}/\sqrt{2V}.    \tag{3}
\]
There is no success selection. Every source occupation, reference and
classical label is retained.

Only the apparatus coefficients in (2) wrap around the coordinate list.
The Hamiltonian has no wrapping bond, its Gauss constraints are unchanged,
and no periodic momentum is introduced. For an address separated from
every face by at least \(R\), the profile is an ordinary local tent of
physical radius one. The boundary profiles can have separated pieces
on opposite faces; Section 4 removes these profiles with a quantified
change of observation domain.

Every phase coefficient in (2), divided by \(\pi\), is an integer
multiple of \(R^{-4}\). On the occupation basis the unitary consists
of roots of unity of order \(2R^4\). It is a finite product of local
number-phase gates and commutes with the link cutoff and all Gauss
projections. Its norm is one, even though the generator's norm grows.
No unbounded field has been declared a contraction. An encoder error
of at most \(e\) radians in each nonzero cell coefficient produces
\[
 \|U-\widetilde U\|\le30(2R-1)^3e.                       \tag{4}
\]
The same bound, capped at one, controls the complete instrument's
half trace distance for arbitrary retained quantum references.
It contains no inverse probability of a selected outcome.

## 2. Exact original work and the discrete gradient

Onsite, electric and magnetic terms commute with (2). A hopping
bilinear \(c_y^\dagger S R(U)c_z\) acquires the scalar phase
\(\exp(i\sigma(\phi_z-\phi_y))\), with
\(\phi_y=\pi a f_R(y-X)\). Averaging the retained signs cancels the
sine term. If \(T_{yz}\) is the complete Hermitian hopping bond,
including its \(1/a\), its contribution to the mean energy change is
\[
       -(1-\cos(\phi_z-\phi_y))\langle T_{yz}\rangle .
 \tag{5}
\]
The full-Fock bond norm is at most \(30/a\). The parent identity
survives the source compression because (2) commutes with its projector.

For any actual interior bond in a chosen direction, summing the
cosine factor over all addresses gives the same number \(k_n\).
This is translation of a finite list of apparatus coefficients; it
does not use state translation invariance. The symmetric tent gives
the same \(k_n\) in all three directions. Therefore the exact original
signed first moment of the complete response is
\[
 M_{1,n}=\frac1{2V}\sum_{X,\sigma}
       \operatorname{Tr}\rho_n(U_{X,\sigma}^\dagger H_nU_{X,\sigma}-H_n)
       =-{k_n\over V}\langle H_{\rm hop}\rangle_{\rho_n}.  \tag{6}
\]

The one-dimensional finite sums are
\[
 \sum_k t_R(k)=R,\quad
 S_R:=\sum_k t_R(k)^2={2R^2+1\over3R},\quad
 \sum_k(t_R(k+1)-t_R(k))^2={2\over R}.
 \tag{7}
\]
They also hold at \(R=1,L=2\). Hence the squared-gradient sum in
one direction for the three-dimensional tent is
\[
                   D_R={2(2R^2+1)^2\over9R^3}.           \tag{8}
\]
Since \(|a\Delta f_R|\le1\),
\(2(a\Delta f_R)^2\le1-\cos(\pi a\Delta f_R)
\le5(a\Delta f_R)^2\). At level zero \(k_0=4\) exactly;
at other levels the rational enclosure is
\[
                       2a^2D_R\le k_n\le5a^2D_R.         \tag{9}
\]
Combining (1), (6), (9) and
\(\|H_{\rm hop}\|\le30B/a\), \(B=3L^2(L-1)\), proves
\[
 0<{16\over9}h_*
 \le {2aD_R}h_n\le M_{1,n}
 \le {150B\over V}aD_R\le {14175\over32}.                 \tag{10}
\]
For \(n\ge1\), the penultimate upper bound equals
\(100(4+2/R^2-1/R^4-1/(2R^6))\), whose maximum on
\(R=2,4,\ldots\) is \(14175/32\). Level zero uses
the exact \(k_0\) and has upper bound \(180\).
The strict positive bound in (10) refers to the signed energy-weighted
original response, not a probability normalized within a band.

## 3. Unfiltered probability tightness and its exact boundary

In a physical energy basis the measure is
\[
 \mu_n(d\omega)=\frac1{2V}\sum_{X,\sigma,i,j}
 p_i |(U_{X,\sigma})_{ji}|^2\,\delta_{E_j-E_i}(d\omega).
 \tag{11}
\]
It is a probability measure, with total original mass one.
The source preparation is stationary. The reference purification need
not be stationary under its idle-reference joint Hamiltonian.
Equation (11) can be read with two complete source energy instruments,
retaining the initial record. That initial measurement preserves the
source marginal and generally changes the joint purification.

Completeness and the selected ground excess give
\[
 \int_{\omega<0}|\omega|\,d\mu_n
 \le\frac1{2V}\sum_{X,\sigma,i}p_i(E_i-E_0)
                    \langle i|U_{X,\sigma}^\dagger U_{X,\sigma}|i\rangle
 \le\epsilon.                                           \tag{12}
\]
This reuses the [ground-record bound](../../Lean/QFT/GroundRecordTail.lean)
in its first-moment form. Thus (10) supplies
\[
 \int|\omega|\,d\mu_n\le{14175\over32}+2\epsilon
 \le {14179\over32},\qquad
 \mu_n(|\omega|>E)\le\min\!\left(1,{14179\over32E}\right).
 \tag{13}
\]
These are uniform original-probability bounds for the actual local
unitary experiment. They use neither a temporal filter nor a division
by a small successful weight. The correlations satisfy
\[
 |F_n(t)-F_n(s)|\le{14179\over32}|t-s|,\qquad
 F_n(t)=\int e^{it\omega}\,d\mu_n(\omega).                  \tag{14}
\]
Consequently subsequences have weak probability limits with continuous
characteristic functions and total mass one. By (12) every such limit
is supported on nonnegative energies. This is the one-response measure
conclusion, not a common field representation for all apparatuses.

The lower first moment in (10) does not exclude a purely elastic limit.
The finite two-level control
\[
 \mu_a=(1-a)\delta_0+a\delta_{1/a},\qquad a\downarrow0,
 \tag{15}
\]
is generated by a unitary rotation of the ground vector. It has total
mass one, first moment one, uniform probability tightness and
\(\mu_a\Rightarrow\delta_0\). Thus passing the positive first-moment
floor to the limit requires an additional uniform-integrability bound,
or a directly populated fixed positive-energy window. This control
tests the inference; it is not a replacement for (11).

## 4. The same retained bulk-origin protocol

Use the authenticated
[retained-origin cube](../particle_primary_origin_momentum/PROOF.md)
\(F_n=[b_n,L_n-b_n)^3\). For \(n\ge1\), \(b_n\ge R\), so
every profile with \(X\in F_n\) is unwrapped and supported in the
ordinary source box. Restricting the classical choice to this
predeclared menu is a different complete preparation of the address;
it does not condition the quantum experiment on success.

Let \(f_n=|F_n|/V\), \(d_n=1-f_n\). On a common classical label
register, the full and restricted complete channels have half trace
distance exactly \(d_n\), for every input source-reference state.
The same upper bound holds after any common complete energy
readout. It tends to zero along the source sequence.

For each fixed address, the sign-averaged work is bounded above by
\(90k_n/a\), from the bond norm and the three directional cosine sums.
For \(n\ge1\) this is at most \(2025/4\); level zero has bound \(180\).
If \(m_n\) is the lower value computed in (10), the retained-origin
measure therefore obeys
\[
 M^{F}_{1,n}\ge
       \max\!\left(0,{m_n-(2025/4)d_n\over f_n}\right),
 \qquad
 \int|\omega|\,d\mu^F_n\le{2025\over4}+2\epsilon.          \tag{16}
\]
The nonnegative lower bound also uses passivity: the selected
populations are nonincreasing functions of energy, so a unitary cannot
lower their mean energy. This is the finite rearrangement inequality,
with the exact polynomial populations from the
[pair-channel proof](../particle_primary_pair_channels/PROOF.md).
The code uses the member's sharper address cap.
The explicit level-32 receipt has a positive lower value in (16).
No matrix of that large carrier is allocated.

The total-variation comparison gives the same subsequential energy
limits for the two origin menus. It does not identify these Weyl
families with the regular spatially smeared probe menu in a momentum
representation, or transfer survival through a further spatial smear.
It makes no trace-norm convergence claim for the growing reference
register itself.

## 5. A finite populated interval and complete clock errors

A finite second moment follows without the electric cutoff norm.
The exact tent gradient satisfies
\(\sum_{y,j}|\Delta_j f_R(y)|=6R^2\).
Using \(|e^{it}-1|\le|t|\) gives
\[
 \|U^\dagger HU-H\|\le180\pi R^2,\qquad
 \int\omega^2\,d\mu_n\le324000R^4.                       \tag{17}
\]
At level zero the twist is local number parity. For the actual
corner star the neighbour Gram matrix has eigenvalues
\((3\pm\sqrt3)/2\), each twice per native row. The sum of its
singular values gives the complete Fock norm, including every
occupation sector. With fifteen native rows,
\[
 \|U^\dagger HU-H\|^2
 \le16(15)^2(3+\sqrt6)<19620.                             \tag{18}
\]
Parent-coordinate estimates descend by compression; the unitary
commutes with that projection.

For any probability response with signed first moment at least \(m>0\)
and second moment at most \(S\), split its positive first moment below
\(m/4\), inside \([m/4,4S/m]\), and above \(4S/m\).
The first and last pieces are at most \(m/4\) each. Cauchy-Schwarz
then proves
\[
              \mu([m/4,4S/m])\ge m^2/(4S).               \tag{19}
\]
The receipt evaluates (19) on the actual source. Its interval is broad
and capacity dependent; no line or rest frame is inferred.

The [complete two-energy clock](../particle_primary_joint_spectroscopy/PROOF.md)
is used with resolution \(m/16\), single-clock tail \(w/16\) and
implementation error \(w/16\), where \(w=m^2/(4S)\).
Two tails and both complete clock errors cost at most \(w/4\).
The reported interval enlarged by \(m/8\) at each endpoint therefore
retains at least \(3w/4\) original probability. All \(2VJ^2\) fine
records, the sign, address, initial energy and reference are retained.
The exact phase service in (2) is declared; an approximate phase
encoder additionally costs (4). None of these large source or clock
register plans is described as an efficient execution.

This supplies actual bounded spatial fluctuations with uniform original
energy tightness and certified finite response. A nonzero limiting
excitation residue, a common momentum attachment, an isolated branch,
and its rest-energy identification remain separate physical premises.
The declared SI conversion can convert an admitted energy interval;
the construction supplies no source-derived particle mass in GeV.
