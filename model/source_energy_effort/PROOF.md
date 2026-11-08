# Work and action in the protected source

## 1. Source and complete charge instrument

Use the selected eight-spin H_B source with J=1 and B=1/4. Its full integer
eigensystem certificate and unique ground Omega_B are those of
`native_ordered_ground`. The catalogue proves that the conserved charge Q
is total internal Sx, with the Marshall signs included in its matrix
representation, and that Q Omega_B=0. The physical spin frame is obtained
by the exact same Marshall diagonal used in `spin_modular_geometry`.
No Hubbard clock or infinite ordered state is substituted.

At site 0 let Z=Z_0 and A_q=Pi_(q X_0)Z, q=+-1. Then

\[
 [Q,A_q]=q A_q,\qquad A_+^\dagger A_++A_-^\dagger A_-=I,
 \qquad A_++A_-=Z.
\]

The source instrument is Z followed by one coherent X_0 pointer coupling.
It preserves coherence within each created charge sector. Reading all
individual X spins and grouping their outcomes would be a different,
dephasing instrument. This executable prepares the ground from blank and
executes the local instrument, retaining both pointer outcomes, reserved
pointer labels and every first failure.

The finite local weights below depend on this selected finite-B source.
They are separate from collective operators in an ordered continuum limit.

## 2. When records suffice to measure energy

Let F_r be a complete finite POVM on a fixed source space. A real score
s_r gives an unbiased energy estimator for every input state exactly when

\[
 H=\sum_r s_r F_r.
\]

Necessity follows because equality of expectations for every density
matrix implies equality of Hermitian operators; sufficiency is linearity.
More generally, the full record distribution determines mean energy for
every state exactly when H lies in the real span of the effects. If it
does not, take a nonzero Hermitian D orthogonal to that span but with
Tr(HD) nonzero. Since sum F_r=I, Tr D=0. For sufficiently small delta,
rho_+-=I/d+-delta D are densities with identical record laws and different
energy expectations. Repeated samples cannot distinguish these laws.

Every fine bare C/E branch on the protected one-particle-per-site code is
I/sqrt(2), with raw repair record 0 and the retained coin label. The producer
checks all 12,288 combinations of eight-spin basis input, cube seam, law and
coin. Linearity proves the assertion on coherent states and entangled
references. Induction gives the same result for any schedule whose choices
depend on these records. Every history effect is a scalar multiple of I.
The source Hamiltonian is noncentral: the all-up Z product state has
expectation 0, and the alternating Z product state has expectation -6.
Therefore any functional of only these repair histories is energy-blind.
This result constrains the choice of an energy observable; it does not
prevent measurement with the source's coherent controls.

## 3. Local creation work and a complete energy estimator

For a prepared ground and any bounded creation A define its unconditioned
excitation moment

\[
 W_A=\langle A^\dagger(H_B-E_0)A\rangle
     =\langle A^\dagger[H_B,A]\rangle.
\]

The Hermitian source observable
B_A=A^dagger H_B A-{A^dagger A,H_B}/2 has this same expectation in the
ground. Its finite support follows directly from the local source action.
For A=A_q, the executable constructs its exact integer numerator and
bounds its norm by a row norm. If the normalized rational ground proposal
has vector error epsilon, its W_q error is at most 2||B_q||epsilon; its
channel probability error is at most 2epsilon. Positive lower probability
bounds permit interval division to obtain conditional mean energies.
There is no discarded channel or normalization of the complete experiment.

Since [H_B,Q]=0 and Q Omega_B=0, the vectors A_+Omega_B and A_-Omega_B are
energy-orthogonal. Thus

\[
 W_++W_-=\langle ZH_BZ-H_B\rangle.
\]

The operator on the right is explicitly
-sum_(y~0)(X_0X_y+Y_0Y_y)/2+X_0/4 in the physical frame. This is the
expected system-energy work of the local unitary creation. It is positive
even though a subsequent bare protected repair consumes no system energy.

For a direct measurement, write the exact source action as
H_B=h0 I+sum_(j=1)^44 a_j P_j, with h0=-3, bond coefficients 1/4 and
selection-field coefficients+-1/8. Six source fair coins choose 64 settings;
the additional 20 settings have coefficient 0 and retain both pointer labels.
Setting j measures P_j and scores h0+64 a_j r, r=+-1. The effect identity in
Section 2 holds exactly. The producer checks this full matrix identity
against the actual source H_B, including its constant.

The source protocol first creates A_q, then chooses a setting and measures
the Pauli string. It retains all 262,144 work rows of each two-pointer
instrument. Prefix success amplitudes are never divided by their norm.
Prefix failures, selector failures, all conditional first failures and
all 256 setting/charge/energy pointer labels constitute the complete output.
Every failure has declared score 0. Every successful energy score subtracts
the same certified ground reference E_ref, with error delta_E.

Let eta be the sum of the complete creation, selector and largest
conditional measurement errors, including the selected-ground error. Then
the charge weights have error at most eta. With M the largest absolute
excitation score, each unconditioned energy moment has error at most
2M eta+delta_E. Here eta<6.65e-13 and M<33. The measured moments agree with
the independent exact commutator intervals. This is an unbiased first-moment
measurement of H_B in each created sector, not a spectral projection of H_B.
The conditional means supply no particle pole or relativistic rest mass.

## 4. Action from source controls and its clock readout

The full native source supplies twelve anti-Hermitian six-level generators
G_p. A literal integrated pulse has G=sum_p u_p G_p, including all twelve
binary-rational controls and every complement level. For the full input
rho of that pulse, its signed action is

\[
 a_G(\rho)=\int_0^1\operatorname{Tr}(e^{sG}\rho e^{-sG}\,iG)ds
          =\operatorname{Tr}(\rho\,iG).
\]

For a measured tape, sum this quantity over every surviving prefix, using
its subnormalized input density. Halted failure branches incur no later
pulse action. This gives an operational source functional, computable from
the specified matrices and observable via their controlled phase. It is
signed and depends on the state. The producer evaluates exact Q(sqrt5)
intervals for literal Z, inverse Z, X and I pulses on the full six-level
blank, each with the same pulse count. Their actions differ. Constant
phases cannot be removed from a controlled clock.

For the specified common clock L=H_B-E0 I, the action of an autonomous
interval is

\[
 \mathcal A_\rho(T)=\int_0^T\operatorname{Tr}(\rho(t)L)dt
                  =T\operatorname{Tr}(\rho L),\qquad
 \dot E= i\operatorname{Tr}(\rho[H_B,H_B])=0.
\]

Hence energy is the relative generator-action rate for this clock, while
ongoing system power is zero. Inserting bare C/E repairs does not change
either statement on the protected code. Realizing the controlled clock
with a native pulse tape has a separate implementation cost. No equality
between unread bare C/E dynamics and H_B time evolution is asserted.

The action rate has an independent complete source readout. A native
Hadamard pointer controls a symmetric product for exp(-it(H_B-E_ref)),
including the scalar ground-phase correction, then reads its sine
quadrature. On charge q the signed pointer moment equals
<A_q^dagger sin(tL)A_q> in the ideal instrument. With ||L||<=28,

\[
 \left|t^{-1}\langle A_q^\dagger\sin(tL)A_q\rangle-W_q\right|
 \le t^2 28^3/6.
\]

If S is the sum of absolute Pauli coefficients, the symmetric-product
operator error is S^3 t^3/3. Add all literal-angle, full-generator,
preparation and phase-reference errors. A complete output error nu changes
the signed moment by at most 2nu, giving total energy error
2nu/t+t^2 28^3/6. At t=1/4096 the actual 4,907-pulse source history gives
error below 3.59e-4 and independently recovers both unconditioned work
moments. No continuum derivative of a rounded binary control is assumed;
this is a finite-time measured estimate with an explicit remainder.
Specifically, nu includes t delta_E for the controlled ground-reference
phase, so the displayed derivative error includes 2 delta_E. The scalar
reference is retained in both the pulse word and the target characteristic.

A separate complete native Z^2 control returns every ideal input to itself
and has zero ideal system work, with four primitive pulses and absolute
angle sum 2pi. Pulse count and absolute control effort therefore do not
equal this state energy or its work. They measure a different resource.

## 5. Finite thermodynamic meaning and physical boundary

For any complete instrument, the system energy balance is
Delta E=sum_r Tr(H K_r rho K_r^dagger)-Tr(H rho). It is an exact finite
identity and retains every record. On the pure selected ground the global
logarithmic reference is not faithful. The actual two-site restriction is
faithful and supplies K_A=-log rho_A through `spin_modular_geometry`.
For local Z creation its regional entropy is unchanged by unitary
conjugation, but

\[
 \Delta\langle K_A\rangle
   =D(Z\rho_A Z\|\rho_A)>0.
\]

The source certificate bounds this relative entropy independently of the
positive H_B work. The modular logarithm error contributes 2kappa and the
reference-vector error contributes 16epsilon. This uses the actual source
reference and no Gibbs state or temperature. It discharges the finite
modular identity on this marginal, not an identification K_A=beta H_B.

Physical apparatus energy, memory-reset heat, a reservoir temperature and
the autonomous origin of the controlled field clock require their own
source constructions. A fair repair record can require memory even when
the protected system energy is unchanged. The finite created-channel
energies here are separate from collective continuum charge excitations,
mass shells, laboratory masses and an interacting relativistic limit.
