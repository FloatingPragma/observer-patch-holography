# Weighted source energy, stress and finite boost response

## 1. Reused theorems and actual source attachment

Use exactly the selected cube source, its ground Omega_B, internal charge Q
and Hamiltonian H_B from `native_ordered_ground` and `source_energy_effort`.
The densities e_x are the same half-bond plus staggered-field densities in
`native_spin_clock_join/energy.py`. That package proves, on any finite
graph,

\[
 H_B=\sum_xe_x,\qquad j_{xy}=i[e_x,e_y]=-j_{yx},\qquad
 \dot e_x+\sum_yj_{xy}=0.
\]

The integer source matrices here instantiate those exact premises for all
eight cells. No substitute field action, thermal state or zero-field vacuum
is used. `spin_stress_response` supplies the corresponding driven-source
work identities. The general exponential flow and its commutator
derivative are the finite results of `Dynamics/PrivateInner.lean`;
`StoneConverse.lean` fixes a generator only up to scalar until the common
clock and reference phase are supplied. Those inputs are supplied here by
the same H_B and E0 used in the energy experiment.

For cube axis i put r_i(x)=(2 bit_i(x)-1)/2 and define

\[
 B_i=\sum_x r_i(x)e_x,\qquad
 P_i=i[H_B,B_i]=\sum_{x<y}(r_i(y)-r_i(x))j_{xy}.
\]

The last equality is the inherited weighted Ward identity. These source
operators commute with Q. Each B_i has norm at most 5/2 and each P_i norm
at most 2. Here P_i denotes the source energy-flux candidate, with the
declared lattice and clock units. No identification with spatial
translations or physical particle momentum is included in that notation.
The existing bipartite field-current cancellation gives an additional
attachment: the H_B and H_J energy fluxes P_i are exactly the same operator.
The producer checks this identity on all three cube axes. Their force
difference is the explicit i[H_B-H_J,P_i] torque. The finite-box H_J force
is nonzero as well. These operator identities retain the prepared state;
using a different clock on its histories requires the common-state bounds
of `native_spin_clock_join`, rather than a zero-field vacuum substitution.

There are actual conserved finite translation readouts on this source.
For each even coordinate shift a in {3,5,6}, let T_a permute cube sites
by bitwise xor with a. Exact source matrices show T_a^2=I, [H_B,T_a]=0,
[Q,T_a]=0, T_a e_x T_a^dagger=e_(x xor a), and the same transport for
the local coherent creation A_(q,x). The unique selected ground has
T_a eigenvalue+1. The principal Fourier logarithm is
K_a=pi(I-T_a)/2, with exp(iK_a)=T_a and [H_B,K_a]=0. This directly
instantiates the general polynomial-conservation theorem in
`QFT/PeriodicMomentum.lean`, without importing its separate CAR model.
The energy flux is different: T_a reverses P_i on either flipped axis,
whereas K_a commutes with T_a and H_B. These even shifts form the finite
group Z2 times Z2; their logarithms supply no three independent continuum
momenta or physical momentum normalization.

## 2. Local stress and transport of energy dipoles

The local densities

\[
 p_x^i=\tfrac12\sum_y(r_i(y)-r_i(x))j_{xy},\quad
 \tau_{xy}^i=i([e_x,p_y^i]-[e_y,p_x^i]),\quad
 f_x^i=i[e_x,P_i]
\]

satisfy sum_x p_x^i=P_i, tau_xy=-tau_yx and

\[
 \dot p_x^i+\sum_y\tau_{xy}^i=f_x^i,\qquad
 \sum_xf_x^i=i[H_B,P_i].
\]

These are Hermitian operators with finite source support. For a nearest
neighbor graph, e_x is supported in its radius-one ball, j_xy vanishes
when dist(x,y)>2, and p_x has support in the radius-three ball. Consequently
the stress label vanishes when dist(x,y)>4. The finite cube calculation
checks all three columns and every cell without a state restriction.
The total force is nonzero, so it cannot be renamed a closed antisymmetric
internal flux. `QFT/LatticeMomentumBalance.lean` proves that last finite
sum implication. Its Weyl/contact matrix attachment is separate; only the
general finite sum theorem is reused with these spin-source premises.

`QFT/CoarseEnergyWard.lean` proves push-forward continuity and composition
for any additive-group-valued quantities. It applies directly to these
matrix-valued currents. The forced counterpart applies to each real
matrix pairing, preserving f rather than setting it to zero. The executable
checks the eight coordinate block maps for all three source columns.
For a finite time interval, the energy increment is e_x(t1)-e_x(t0) and
the transported current is its time integral. Integrating the operator
identity gives the same block theorem without replacing that integral by
an endpoint current.

For any block map c and chosen representatives R_a, let

\[
 E_a=\sum_{c(x)=a}e_x,\quad
 D_i=\sum_x(r_i(x)-R_{c(x),i})e_x.
\]

Then B_i=sum_a R_ai E_a+D_i and

\[
 P_i=\sum_{a<b}(R_{b,i}-R_{a,i})J_{ab}+i[H_B,D_i],
 \quad J_{ab}=\sum_{c(x)=a,c(y)=b}j_{xy}.
\]

The subcell dipole D_i is required. It is nonzero for twelve of the
twenty-four tested axis/block pairs. For successive block maps, its exact
composition is D_(fine,coarse)=D_(fine,middle)+sum_a(R_a-S_d(a))E_a.
Thus the source provides a refinement rule for energy polarization beyond
coarse totals. No decay of that dipole in a continuum scaling is assumed.

The same local creations A_q of `source_energy_effort` give a concrete
response of these columns. Their unconditioned initial accelerations
<A_q^dagger i[H_B,P_i] A_q> are positive. The certified values are about
0.1383566404 for q=+1 and 0.4410041545 for q=-1, equal on the three cube
axes. These are energy-transport responses of the finite created channels,
not a physical dispersion relation or a claim of conserved momentum.

## 3. Exact boost defect and a complete source observation

Set L=H_B-E0 I and

\[
 D_i^{\rm boost}=i[B_i,P_i]+L.
\]

This is the defect in the second commutator required by a normalized
Lorentz energy/boost pair. The first relation i[B_i,L]=-P_i holds by
construction. An exact coordinate reflection followed by the source's
global pi spin rotation preserves H_B and reverses B_i. Uniqueness of the
selected ground therefore gives <B_i>=0. Direct source certification gives

\[
 \langle B_i^2\rangle>0.02499935348,\qquad
 \langle B_i L B_i\rangle>0.06557717030.
\]

In the actual ground,
<D_i^boost>=-2<B_i L B_i>, with certified upper bound below -0.13115434068.
For an arbitrary positive candidate speed c, use B_i/c and P_i/c. The
corresponding ground defect is -2<B_i L B_i>/c^2 and is nonzero. A speed
choice cannot fix this selected finite ground. Independently, the trace
of i[B_i,P_i]+L equals Tr L>0, excluding an exact full finite positive-energy
Lorentz algebra for any nontrivial source. These are finite statements;
neither excludes an appropriately controlled continuum limit.

The weighted generator is operational. Its 32 nonzero Pauli terms follow
from the actual e_x, with coefficient absolute sum 7/2. A native Hadamard
pointer controls exp(-i eta B_0), approximated by a symmetric product with
n subdivisions. Its product error is at most (7/2)^3 |eta|^3/(3n^2).
The source ground-preparation prefix, all literal controls and every
first failure are included. On the true ground the ideal minus probability
obeys

\[
 p_-\ge \eta^2\langle B_0^2\rangle/4
          -\eta^4\|B_0\|^4/48.
\]

This follows from the scalar cosine remainder and the spectral theorem.
The source reflection also makes the characteristic real; with the small
declared eta it is positive, so nonzero p_- implies failure of ground-ray
invariance, rather than just an arbitrary constant phase. At eta=1/1024
and n=4 the complete 11,821-pulse experiment starts from blank, retains
all 131,072 work rows and has unconditioned minus probability certified
above 5.1235e-9. No success normalization is applied. The measurement shows
that this finite weighted-energy flow is a nontrivial source operation;
it is not identified as an operational spacetime boost.

## 4. Quantitative same-source boost-history reduction

For one axis let X_eta=e^(i eta B) X e^(-i eta B). The existing finite
inner-flow derivative gives L_eta'=-P_eta and
P_eta'=-L_eta+D_eta. Variation of constants yields exactly

\[
 L_\eta-L\cosh\eta+P\sinh\eta
   =-\int_0^\eta\sinh(\eta-s)D_s\,ds,
\]
\[
 P_\eta-P\cosh\eta+L\sinh\eta
   =\int_0^\eta\cosh(\eta-s)D_s\,ds.
\]

Hence the respective operator errors are bounded by
(cosh|eta|-1)||D|| and sinh|eta| ||D||. Duhamel and telescoping give a
complete-history halftrace bound

\[
 \min\{1,\|D\|[(\cosh|\eta|-1)\sum_j|t_j|
                      +\sinh|\eta|\sum_k|s_k|]\}.
\]

Here the histories use the two stated operator clocks with durations t_j
and s_k, and arbitrary identical intervening complete channels, records
and untouched references. No local-observable or success-only distance
replaces that bound. The P clock is a controlled source operator; its
attachment to spatial translations is a separate requirement.
If a duration depends on a retained record, its absolute value in the sum
is the largest controlled-block duration at that stage, as in the existing
common-state clock theorem. The bound therefore includes finite feedback
without replacing it by an average duration.

The first component has an eta-squared error, but the second has an
eta error. Sending a boost step to zero cannot discard the accumulated
second-commutator defect and infer finite rapidity covariance. A continuum
application must control these defect matrix elements on the actual
propagated observer states, retain energy-dipole refinement, identify the
spatial action and preserve the common clock. `WardLimitManifest.lean`
requires these source, refinement, chart, weak-limit and residual
attachments before its limit theorem yields a physical Ward identity.
The finite-rotation/one-boost cone theorem consumes operational covariance
of the whole history/readout family. The algebraic celestial map alone
does not discharge it.
