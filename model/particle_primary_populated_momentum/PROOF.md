# Populated energy records in a scanned spatial analyzer

Use member zero of the [primary contract](../particle_common_contract/PROOF.md):
the open eight-cell Wilson/Gauss carrier, all sixty CAR modes per cell,
its selected polynomial density, idle conjugate reference and source
Hamiltonian clock. The construction changes a neutral apparatus phase.
It changes no source action, preparation, clock or energy zero.

The [pair-channel proof](../particle_primary_pair_channels/PROOF.md)
supplies complete spatial/internal instruments and positive discrete
Walsh records. The [checkerboard calculation](../particle_primary_staggered_response/PROOF.md)
supplies an actual Wilson cancellation in one of those records. The
construction here calculates the response to a continuously scanned
finite spatial phase, including a neighborhood of that cancellation.
It bounds original joint spatial/energy outcomes, without identifying
the spatial label with conserved momentum.

## 1. One complete phase-scanned experiment

Write \(x\in\{0,1\}^3\), and use the physical local pair components
\(C_{xa}\), \(a=0,1,2,3\), from the pair-channel proof. Their scalar
component is
\[
 C_{x0}=2^{-1/2}\sum_{\alpha=1}^2 B_{x\alpha},\qquad
 B_{x\alpha}=c^\dagger(p_{x\alpha})c(m_{x\alpha})
\]
in the dimension-one native block of integer charge six. Both endpoints
of each local bilinear are at the same source site, so its gauge charge
is zero. It is CAR even and preserves all global native occupations.

Put \(C_{0,xa}=C_{xa}\),
\(C_{1,xa}=(I-C_{xa}^\dagger C_{xa})^{1/2}\). A uniform coherent
32-dimensional address, the controlled complete columns and analyzer
\[
 W_{(k,a),(x,b)}=8^{-1/2}(-1)^{k\cdot x}
             e^{-i\kappa\cdot x}\delta_{ab},
 \qquad k\in\{0,1\}^3                                      \tag{1}
\]
give the 64 success/failure Kraus operators
\[
 K_{b,k,a}=32^{-1/2}\sum_{x,c}W_{(k,a),(x,c)}C_{b,xc}.
 \qquad \sum_{b,k,a}K_{b,k,a}^\dagger K_{b,k,a}=I.            \tag{2}
\]
Equation (1) is the unchanged Walsh analyzer multiplied by a unitary
diagonal phase. For rational \(\kappa_i/\pi\), its phases are algebraic
roots of unity in the declared finite matrix service. General real
phases describe an analytic family; their finite approximations have the
explicit error below. No irrational input oracle is assumed.

For a successful scalar record, let \(q=\kappa+\pi k\). Then
\[
 K(q)=\frac1{16\sqrt2}\sum_{x,\alpha}e^{-iq\cdot x}B_{x\alpha}.
                                                               \tag{3}
\]
The finite phases depend on \(q\) modulo \(2\pi\) per coordinate.
The receipt retains both the analyzer convention and its signed
representative in \([-\pi,\pi)\). Wavevector units are radians per
source-cell spacing, with \(a=1\) in this member.
All other internal records and all failures are retained with the source
and arbitrary quantum reference. A successful scalar event is not
renormalized by its probability or by its analyzer multiplicity. An
incoherent random choice of site removes the spatial cross terms and
does not implement (1).
Parseval also gives the scalar subinstrument bound
\[
 \sum_kK_{0,k,0}^\dagger K_{0,k,0}
     =\frac1{32}\sum_xC_{x0}^\dagger C_{x0}\le\frac14I.     \tag{3a}
\]
It controls the union of all eight successful scalar rows with no
individual-row normalization.

## 2. The actual Wilson dependence on the spatial phase

Let \(D_6\) count beta-positive particles and beta-negative holes in
this native block. The exact CAR identities and full-Fock quasispin
bound in the [coherent response proof](../particle_primary_odd_residual/PROOF.md)
extend to unit-modulus coefficients:
\[
 \|K(q)\|^2\le9/64,\qquad
 [K(q)^\dagger,K(q)]=(16-D_6)/512.                         \tag{4}
\]
Local phase conjugation preserves the quasispin norm. Bilinears at
different sites commute, giving the second identity independently of
the phase.

In the actual beta chart write the hopping one-body operator as
\(h=\left(\begin{smallmatrix}A&B\\B^\dagger&D\end{smallmatrix}\right)\).
If \(W_q\) is the diagonal cell-phase matrix, the raising block is
\(R_q=\left(\begin{smallmatrix}0&W_q\\0&0\end{smallmatrix}\right)\), and
\[
 [h,R_q]=\begin{pmatrix}-W_qB^\dagger&AW_q-W_qD\\0&B^\dagger W_q\end{pmatrix}.
                                                               \tag{5}
\]
The literal Wilson diagonal blocks on an oriented edge are \(-I_2/2\)
and \(I_2/2\). Hence the creation coefficient on \(x\to y\) is
\(-(e^{-iq\cdot x}+e^{-iq\cdot y})I_2/2\), multiplied by the actual
unit-modulus native link character. There are eight oriented edges and
two spin rows per axis. Pointwise on the full Haar carrier,
\[
 \|AW_q-W_qD\|_{\rm HS}^2
       =16 c(q),\qquad c(q)=\sum_{i=1}^3\cos^2(q_i/2).     \tag{6}
\]
No link is replaced by its identity value. The number-preserving
particle/hole remainder has trace-zero one-body matrix of norm at most
\(3/2\), since multiplication by \(W_q\) preserves the inherited
off-diagonal block norm bound.

The full-Fock creation and number-sector inequalities from the coherent
response proof therefore give
\[
 L_q^\dagger L_q\le16c(q)(1+D_6/2),\qquad
 N_q^\dagger N_q\le36D_6,
 \quad [H,K(q)]-6K(q)=(L_q+N_q)/(16\sqrt2).               \tag{7}
\]
The first bound follows by singular-value decomposition and the
hard-core-pair creation estimate on every occupation sector. For the
second, on \(D_6=j\) use
\(\|N_q\|\le(3/2)\min(j,32-j)\) and
\(\min(j,32-j)^2\le16j\). These are form bounds, including blocked
and fully occupied states. The product link projection commutes with
the probes and \(D_6\), so compression and Gauss restriction preserve
them.

For any \(t>0\), (7) yields the source-derived centered moment bound
\[
 \int(\omega-6)^2d\mu_q(\omega)
 \le M_t(c,d):=\frac{(1+t)16c(1+d/2)+(1+1/t)36d}{512},    \tag{8}
\]
where \(d\) is the authenticated selected bound on
\(\operatorname{Tr}\rho D_6\), and
\(\mu_q=\sum_{i,j}p_i|K(q)_{ji}|^2\delta_{E_j-E_i}\) is its original
stationary transfer measure. The actual
[native allocation](../particle_primary_branch_resolution/PROOF.md)
proves \(d<46/25\). Equation (4) gives total original success weight
at least \((16-d)/512\). Consequently
\[
 \mu_q([6-r,6+r])\ge\frac{16-d}{512}-\frac{M_t(c,d)}{r^2}.
                                                               \tag{9}
\]
This is an inequality on the selected source density, with no auxiliary
trial vector substituted for it.

## 3. Populated bands throughout the scan

For every spatial phase, \(c(q)\le3\). With \(t=6/7\), \(r=21/4\)
and \(d<46/25\), (9) gives
\[
 \mu_q([3/4,45/4])>7057/1317120>5/1024.                 \tag{10}
\]
This controls every scalar record in every declared rational-phase
analyzer, including zero spatial phase. It does not establish a rest
frame for that record.

If \(|q_i-\pi|\le1/4\) in each coordinate, then
\(c(q)\le3/64\). Taking \(t=7\), \(r=3\) in (9) gives the
finite open neighborhood bound
\[
 \mu_q([3,9])>391/44800>1/128.                          \tag{11}
\]
Thus the narrower populated band is stable under an actual variation of
the detector phase; it is not confined to one exceptional exact mask.

The executable configuration is \(\kappa=(\pi/16,\pi/16,\pi/16)\),
using order-32 phase gates. Let \(s=k_1+k_2+k_3\). Since
\(\sin^2(\pi/32)<(11/112)^2\),
\(c(q)<3-s+s(121/12544)\). Exact rational use of (9) gives:

| Shifted axes \(s\) | Scalar records | Transfer band | Original weight per record exceeds |
|---|---:|---|---|
| 0 | 1 | \([3/4,45/4]\) | \(7057/1317120\) |
| 1 | 3 | \([1,11]\) | \(30027/3920000\) |
| 2 | 3 | \([3/2,21/2]\) | \(80081/8467200\) |
| 3 | 1 | \([3,9]\) | \(6883/705600\) |

The Young parameters are respectively \(6/7,1,3/2,9\). The bounds
sum because these eight successful records are disjoint outcomes in
one complete apparatus. Their corresponding energy-band union has
original weight greater than \(1367869/20580000>0.06646\).
Different interval widths are bounds, not measured dispersion curves
or separate species assignments.

## 4. Actual spatial transformations and complete energy clocks

Using centered positions \(r_x=x-(1/2,1/2,1/2)\) changes (3) only
by a known scalar phase. The exact physical cube action from the
pair-channel proof gives
\[
 U_R\widetilde K(q)U_R^\dagger=\widetilde K(Rq),\qquad
 U_I\widetilde K(q)U_I^\dagger=-\widetilde K(-q).          \tag{12}
\]
Here \(R\) is a proper cube rotation and \(I\) is physical inversion.
Both commute with the source Hamiltonian and its polynomial density.
Their original ideal energy-transfer measures therefore obey
\(\mu_q=\mu_{Rq}=\mu_{-q}\), with the correspondingly transformed
apparatus. The retained canonical purification has the matching
source/conjugate-reference covariance. No source-only invariance of
that purification is used. Spatial covariance classifies these finite
operator records; it supplies no continuum spin assignment.

Use the [complete source clocks](../particle_primary_joint_spectroscopy/PROOF.md)
with resolution \(\eta=1/64\), per-clock tail allowance \(1/32768\)
and per-clock implementation allowance \(1/32768\). Their event loss
is at most \(1/8192\). Each interval in the table expands by \(1/32\)
at both endpoints. All \(64J^2\) energy/probe histories, initial pointer
and quantum reference survive. The initial energy apparatus preserves
the stationary source marginal and can change its reference correlations.

If the analyzer's diagonal-phase implementation differs in operator
norm by \(e_\phi\), its complete output half-trace error for arbitrary
retained references is at most \(e_\phi\). This follows from the
controlled-column isometry and the single diagonal unitary in (1).
There is no site-count, outcome-count or successful-normalization factor.
Coordinate phase errors \(\delta_i\) give
\(e_\phi\le\sum_i|\delta_i|\). Subtract the clock and phase errors
from each original event bound. For the union event, subtract them once,
not once per row. Exact declared phases therefore give reported union
weight greater than \(1/16\). The receipt retains stronger values.

The finite theorem `Lean/QFT/GroundRecordTail.lean` transports each
ideal band \([l,u]\) to final excitation \([l,u+2]\), at original
weight loss at most \((9/64)(1/16)/2=9/2048\). Every table row retains
positive weight. This bounds energy above the same source ground;
it is not a rest-energy measurement.
For their union, (3a) improves the discarded initial-excitation tail to
\((1/4)(1/16)/2=1/128\). The same original event therefore has final
excitation in its row-dependent interval \([l,u+2]\) with weight
greater than \(1367869/20580000-1/128>15/256\). This deduction uses
one initial energy cut shared by the complete instrument.

## 5. The momentum boundary and the independent control

The [relative-observation construction](../particle_primary_origin_momentum/PROOF.md)
has canonical commuting generators \(H\ge0,P\) on its regular limiting
domain. The [coherent-offset construction](../particle_primary_momentum_transfer/PROOF.md)
attaches sinc response kernels to that joint spectral measure, with
actual finite aperture and clock errors. Those results do not turn
the present unsmeared member-zero operator into a populated component
of that limit. Its seed, spatial regularization and refinement domain
must be joined with surviving original weight. A checkerboard label
at fixed lattice spacing moves to wavevector of order \(1/a\) on
refinement, so the table proves no fixed physical-energy or momentum
survival.

The finite open source has cube rotations and inversion, but no proved
commuting continuous translations. `Lean/QFT/PeriodicMomentum.lean`
uses a supplied commuting translation and cannot create it. Equation
(3)'s label is an operational Fourier-aperture setting. Its effects are
generally neither projections nor conserved operators. Cooling alone
does not bound initial momentum, as the explicit control in the owning
momentum proof demonstrates. No finite spatial label is promoted to a
rest frame, relativistic dispersion or a GeV particle mass.

An independent control executes a complete physical six-dimensional
two-endpoint CAR/link sector. It retains four matter modes, link flux
\(-1,0,1\), both Gauss constraints and every admitted occupation. Its
Wilson spin block is \(\left(\begin{smallmatrix}-1&i\\i&1\end{smallmatrix}\right)/2\),
with onsite coefficient three and electric coefficient \(63/4\).
The stationary full-rank control density is proportional to \((40I-H)^2\),
with its actual idle purification. Four complete analyzer records and
144 energy histories check normalization, inversion covariance and
retained-reference completeness. Its effects fail idempotence and
commutation with \(H\); its coherent and incoherent laws differ. This
control tests the readout logic and is not substituted for the primary
source or selected density. The primary proof uses (4)--(9) and its
authenticated source occupation bound; it executes no large primary
matrix, spectrum or histogram.
