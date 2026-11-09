# Koide balance and the capacity-dependent response

The [complete coherent attenuator](BARE_EFFECTIVE_CONTROLLER.md) implements
a balanced bare coefficient and a deformed effective action, with its
attenuation, phase and gain included in the configured apparatus. It retains
all loss records and has a full-CAR weak-feedback attachment. This makes
the two action laws compatible at different stages. They cannot both give
the exact same canonically normalized mass triple when the response
coordinate u is positive. Deriving the attenuation law from primitive
repair, or interpreting it as a physical kinetic or self-energy correction,
requires a separate source calculation.

## The two specified coefficients

Use the same regular three-cycle R, family frame and phase delta. Put

```math
T_delta=(e^{i delta}R+e^{-i delta}R^2)/sqrt(2),
B=(I+T_delta)^2,
D=(I+r T_delta)^2,\qquad r=e^{-chi}.
```

R is unitary, R^3=I, and T_delta is Hermitian. The
[balanced native analyzer](../standard_model_source_map/KOIDE_SOURCE.md)
gives Y=(5/108)B. The
[capacity law](../charged_lepton_response/P_DEPENDENCE.md#the-source-law-and-its-capacity-variable)
gives a positive common multiple of D, where

```math
chi=(u^2/512)/(1-u/77),\qquad
delta=2/9+[(P-phi)^2/(21*pi)]/[1-u/27-u^2/135].
```

The [response fixed point](../../code/particles/leptons/derive_charged_mcpr_completion_conditional.py#L157-L198)
computes chi from its stipulated register sources and feedback. The
[mass readout](../../code/particles/leptons/derive_charged_mcpr_completion_conditional.py#L201-L239)
then inserts exp(-chi) into the root amplitude. These operations are not a
calculation of the fermion kinetic form or its self-energy.
The [conditional incidence calculation](../../code/particles/leptons/derive_charged_face_incidence_conditional_theorem.py#L279-L287)
also distinguishes bare balance from endpoint attenuation and requires a
physical map between those stages.

Assume I+T_delta and I+r T_delta are strictly positive. Their eigenvalues
are 1+sqrt(2) cos(delta+2*pi*j/3) and
1+sqrt(2)r cos(delta+2*pi*j/3). The trigonometric sums give

```math
Q(B)=2/3,\qquad Q(D)=(1+r^2)/3<2/3\quad(chi>0),
Q(M)=sum_j sigma_j(M)/(sum_j sqrt(sigma_j(M)))^2.
```

A common scale and left/right unitary family frames preserve Q. They
therefore cannot identify these two triples. The same-phase comparison
isolates the amplitude change. The separate
[calibrated phase](../charged_lepton_response/CALIBRATED_SANITY.md) is
chosen using the muon/electron ratio and differs from delta(P,u); it cannot
be substituted for the capacity law's phase. Comparing that calibrated
coefficient with D requires accounting for both changes.

## An exact possible correction

Since B and D commute and are positive, the kinetic matrices

```math
Z_L=Z_R=B D^{-1}
```

are positive and obey

```math
Z_L^{-1/2} B Z_R^{-1/2}=D.                         (1)
```

This follows by simultaneous diagonalization: each eigenvalue b is divided
by b/d. A second exact possibility keeps canonical kinetics and changes
the mass coefficient by

```math
E=D-B=2(r-1)T_delta+(r^2-1)T_delta^2.              (2)
```

Common positive normalizations can be restored separately. Equations
(1)-(2) show that a family-dependent kinetic or mass correction can connect
the two shapes. They obtain the correction by solving for the desired
output; they do not establish that a source interaction generates it.
The inverse construction requires strict positivity. At a zero root it
does not define an invertible kinetic form on the full family space.

A passive field-coordinate change supplies no such effect. Starting from
canonical kinetics and mass B, invertible field changes S_L,S_R produce
kinetic forms S_L^* S_L, S_R^* S_R and mass S_L^* B S_R. Write their polar
decompositions S_L=U_L P_L and S_R=U_R P_R. Canonical normalization then
returns U_L^* B U_R, with the original singular values. Holding B fixed
while imposing the kinetic matrices in (1) changes the action. Physically,
it changes how the three families propagate relative to their mass term.

## What the finite source theorems supply

The [complete C-duty inverse](../mcpr_native_memory_bridge/PROOF.md#exact-inverse-of-the-complete-memory-coefficient)
realizes the whole P-dependent coefficient by two admissible rate choices
and a gain. It retains disabled records and preserves canonical kinetics
under its gauge-trivial family frame. Its
[correlated weak-pulse theorem](../mcpr_native_memory_bridge/PROOF.md#correlated-weak-pulse-dynamics)
gives a finite-error action realization. These are positive constructions
of the specified response; they supply no balanced-to-deformed matching
theorem for one fixed controller.

The [free field theorem](../charged_lepton_response/RELATIVISTIC_MASS.md#exact-mass-shell-of-the-native-quadratic-action)
gives m_j=|v_H| sigma_j(Y) with canonical kinetics. It contains no hidden
chi correction. The [ordered-memory result](../charged_lepton_response/SOURCE_PHASE_SELECTION.md)
has an identity-block orientation potential and a simultaneous dissipator.
Its orientation generator commutes with every coherent coefficient, so
conjugation by that generator does not change its phase or balance.

### The fixed native diffusion

The dissipator requires a separate calculation. The exact
[source replay](memory_diffusion.py) uses all 128 native records for each
of twelve old memory inputs, at p=1/12, q=1/6 and unit gain. For a duty's
enabled transition E, full transition F and stationary law pi, put

```math
G=(I-F+pi 1^T)^{-1}-pi 1^T,\quad w=E pi,
K=diag(w)/2-ww^T/2+E G diag(w),\quad V=K+K^T,
L_V(B)=-(1/2) sum_ij V_ij [P_i,[P_j,B]].
```

This is the proved [complete all-lag coefficient](../particle_memory_response/PROOF.md#summing-every-lag),
with V positive semidefinite. It acts identically on the family factor of
I tensor B and X tensor B, since each pulse Hamiltonian is X tensor P_i.
Set H=R+R^2 and J=i(R-R^2). The following exact entries of
[L_V(H),R], in the native family frame with zero-based indices, are nonzero:

| Duty | Entry | Value |
|---|---|---|
| always | (0,1) | -33292836(sqrt(5)-1)/2151685171 |
| C | (0,0) | -8729149(sqrt(5)-1)/4800041280 |
| E | (0,0) | 159298835(sqrt(5)-1)/23040198144 |
| CE | (0,0) | 443064679(sqrt(5)-1)/691205944320 |

Thus the diffusion leaves the circulant algebra for every specified duty.
The code also computes a nonzero [L_V(J),R] for each duty. Since P_i, V
and R are real, the H commutator is real and the J commutator is purely
imaginary. For C_theta=aI+b(cos(theta)H+sin(theta)J), their real and
imaginary parts cannot cancel when b cos(theta) is nonzero. This includes
the strictly positive chamber |theta|<pi/12 of the balanced source.
There is therefore no literal scalar harmonic attenuation semigroup on
that coefficient in this frame. The same conclusion holds for the actual
Y_theta=C_theta^2/3: its H coefficient is
(2ab cos(theta)+b^2 cos(2theta))/3>0 throughout that chamber.

This calculation concerns the stationary unread h G1 correction in the
declared pulse clock. It does not exclude other rates or controllers,
an isolated-time coincidence, an isospectral match after a family rotation,
or a physical self-energy. At fixed gain and target time the correction
vanishes as h tends to zero. It derives no relation between its rate and
chi=(u^2/512)/(1-u/77).

Even a different channel that damped the root would require an action-level
check. For a unital completely positive map Phi and Hermitian C,
Phi(C^2)-Phi(C)^2 is positive semidefinite. Consequently, evolving the
Yukawa effect C^2/3 and squaring an evolved root are different operations.
The tests verify the nonzero infinitesimal variance
L_V(H^2)-L_V(H)H-H L_V(H)=sum_ij V_ij[P_i,H]^*[P_j,H]
for all four duties, as well as an independent three-ray reconstruction
and the full six-dimensional Nambu double commutator. No particle pole
statement follows from these channel moments.

The [complete C-duty response](NATIVE_RESPONSE.md) gives a basis-independent
finite-step restriction: a specified pure input has strictly positive
purity loss at every pulse angle outside 2*pi*Z, for all 0<p,q<1.
Its unread working-register channel therefore differs from every unitary
mass-only replacement. This includes both admissible inverse-rate branches
without fitting a kinetic matrix to their target spectrum. The source's
flagged dilation is retained. Its weak-pulse limit has the proved mean
coefficient and unchanged common kinetic H0; the finite correction has
both the ordering potential and diffusion, with the stated clock and state.

### Independent response matching

The [mass-response theorem](../charged_lepton_response/MASS_RESPONSE_BRIDGE.md#finite-stability-theorem)
bounds the singular values of Z_L^{-1/2}(M_0+E)Z_R^{-1/2} when actual
kinetic and correction bounds are supplied. It computes neither Z_L,R
nor E. Constant kinetic forms must multiply the full Lorentz-covariant
kinetic terms; normalizing only the time derivative does not establish a
common propagation speed.

For a common-source identification, independently extract those kinetic
forms and the mass response on the declared state and clock, retaining
the complete observation records. The concrete test is the equality

```math
Z_L^{-1/2}(B+E)Z_R^{-1/2}=A(P,u)D                 (3)
```

up to specified current-preserving family frames, or a certified error
bound for its singular values. A(P,u)>0 includes the declared common
normalization. Deriving Z or E from the target right-hand side alone is
insufficient. Observable electron thresholds and unstable-lepton
resonances also require the same theory's infrared and response matching.
Neither exact Koide balance of a bare action nor the finite P-dependent
implementation supplies those physical corrections by itself.

## Constructed maps, excluded replacements and physical matching

The [attenuation controller](BARE_EFFECTIVE_CONTROLLER.md) supplies a
complete sixteen-flag analyzer with C_r=aI+r b(zR+conjugate(z)R^2).
Its coherent weak feedback realizes Y_r=C_r^2/3, and its full-CAR lift
preserves the source currents with the stated Higgs covariance and archive
error. Supplying r=exp(-chi), the P-dependent phase and the common gain
therefore constructs the full specified action without postselection.
The [general charged-field attachment](CHARGED_FIELD_ADMISSION.md) accepts
every positive Hermitian Y, including that coefficient and g Z_C. With
the specified kinetic action and Gaussian reference it gives three free
Dirac families per retained valley, with finite controller and compact-
momentum lattice errors. These are constructions of configured fields.

The [three-input calibration](CALIBRATED_TRIPLET.md) makes this existence
statement explicit for the supplied electron, muon and PDG2026 tau central
values. The inverse Fourier transform of their positive square roots gives
`r^2=3Q-1`, a phase and a common action scale. Here `r<1`, so the complete
attenuator realizes the triplet with all sixteen flags retained. Its
certified `10^-36 GeV` rational realization error is not a physical theory
uncertainty. All three masses are inputs. In the distinct balanced Koide
comparison only the muon/electron ratio fixes the phase and tau is held out.
Neither construction derives the controller from primitive repair dynamics
or substitutes for the independent kinetic/self-energy test (3).

The exact restrictions have narrower domains. A common unit or passive
field change cannot change the Koide quotient. The fixed native memory
diffusion does not literally attenuate the C3 harmonic. Its finite C-duty
working-register channel has the positive purity-loss witnesses in
[NATIVE_RESPONSE.md](NATIVE_RESPONSE.md); this excludes an exact unitary
mass-only replacement of that channel. The same proof's fourth-trace
identity excludes a three-mass relativistic free spectrum for the retained
finite-step ordering Hamiltonian alone. It does not identify that
Hamiltonian with the complete noisy source. The charged-field attachment
also identifies the precise default-spectrum, mirror and massless-spectator
obstructions, without excluding other configured actions or interactions.

The attenuator's r, phase, gain and coherent services are readable inputs.
The primitive C/E weights admit several complete attenuators and do not
by themselves select the stipulated chi fixed point. A physical comparison
also needs the common interacting state, clock, charged observable and
kinetic/self-energy matching, with thresholds and resonances distinguished.
The configured action construction supplies no measured mass, natural-
particle identity or defect-to-species assignment by itself.
