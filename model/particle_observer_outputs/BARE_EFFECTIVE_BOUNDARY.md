# Finite interaction response, coefficient reconstruction and particle mass

The literal finite quantum-Higgs source determines a complete nonlinear
response from its chosen family coefficient. Both addition lines of a
resolved channel, together with the neutral transition, recover that
coefficient without using the charged-sector energy offset. The result is
an observable coefficient ratio on the same finite clock. It is not a
relativistic particle-mass ratio.

The owning executable is [`bare_effective.py`](bare_effective.py). It
consumes the actual [charged response](CHARGED_HIGGS_RECORDS.md), including
its original probabilities and all channel multiplicities. It rejects an
invisible line, unresolved ordering or negative reconstructed squared
coefficient. No measured mass enters this calculation.

## Which energies are outputs of these experiments

Set `G=4 Tr(Y†Y)>0`, `d=2 epsilon_F+epsilon_H` and
`Delta=sqrt(d²+4G)`. The [neutral source](HIGGS_RECORDS.md) has sector
ground `E0=(d-Delta)/2`. Its scalar-number correlator has transition energy
Delta and original residue `G/Delta²`. Complete two-time projective
records, with their actual backaction, recover Delta from
`J(2t)/(4J(t))=cos²(Delta t/2)` on the stated unaliased window. This last
ratio is a probability/time relation, not a ratio of particle masses.

The [thermal construction](THERMAL_CLOCK.md) gives the same gap for any
reproducible input with nonzero empty-state weight. Ground cooling is
unnecessary for this observable. For `Y=diag(1,2,3),d=1`, its complete
energy menu is `-7,8,1`, with ranks `1,1,305`; ground-relative values are
`0,15,8`. The connected ground scalar-number correlator has only the
transition line 15. Its full ordered correlator also has an elastic
zero-frequency term. The faithful
Gibbs state populates all three spectral blocks, so a complete energy
detector can compare gaps 7, 8, 15 and their dimensionless ratios. Degenerate
zero frequencies are not massless particles.

The modular identity `-log(rho_beta)=beta H+log(Z)I` fixes
`sigma_s=alpha_(-beta s)` in the declared conventions, on all 307-state
operators. It supplies no new spectrum or laboratory scale. The chosen
temperature and faithful reference do not turn this equilibrium experiment
into a Lorentz-invariant vacuum.

For a squared singular value `s_j=sigma_j(Y)²`, let `kappa=1` for the
left-doublet channel and `kappa=2` for the charged-singlet channel. The
same finite ground and compensated probe have addition energies

\[
 \omega_{j,\pm}=\epsilon_F+
  \frac{\Delta\pm\sqrt{\Delta^2-4\kappa s_j}}2,
 \qquad
 \omega_{\rm rem}=\frac{d+\Delta}2-\epsilon_F.        \tag{1}
\]

Their complete residues follow from the literal source Gram matrix. All
1,614 dark addition directions are retained. The source's local quantum
scalar couples a fermion pair; it has not been replaced by a fixed Higgs
expectation. Therefore (1) is a different observable from the quadratic
action coefficient `v_H sigma_j(Y)`.

Changing epsilon_F at fixed d adds `epsilon_F Q*` to the source. It shifts
addition lines together and removal oppositely, while every neutral output
is unchanged. Some permitted choices produce negative energies relative
to this constrained neutral ground. This is a cross-sector reference issue,
not a prediction of a negative particle mass. The charged probe's reference
energy is part of the stated apparatus convention.

## A coefficient inverse that cancels the charged offset

Suppose the original weights of both addition lines in a channel are
strictly positive and those two energies are resolved. Their separation is

\[
 \delta_j=\omega_{j,+}-\omega_{j,-}
         =\sqrt{\Delta^2-4\kappa s_j},\qquad
 s_j=\frac{\Delta^2-\delta_j^2}{4\kappa}.              \tag{2}
\]

The additive epsilon_F cancels exactly. For two channels with the same
kappa and positive denominator,

\[
 \frac{\sigma_i(Y)}{\sigma_j(Y)}=
 \sqrt{\frac{\Delta^2-\delta_i^2}{\Delta^2-\delta_j^2}}. \tag{3}
\]

One positive clock/energy rescaling multiplies all gaps by the same factor;
it cancels in (3). It cannot repair a wrong ratio of the corresponding
source coefficients. The theorem recovers the coefficient of this finite
vertex; a separate physical matching is required before calling (3) a
lepton ratio.

For unresolved channel degeneracies, combine their weights and retain the
original multiplicity. Let `P_lambda` be a spectral projector of the source
Gram operator on its one-fermion field index. A channel-resolved addition
readout uses `conjugate(P_lambda)=P_lambda^T` on the retained dual charged
reference. This transpose matters for a complex family coefficient. These
projectors commute with the corresponding gauge representation. Such a
projector is an additional specified apparatus operation,
not a postselection that renormalizes a favored line. If `s_j=0`, its
upper addition residue vanishes on this ground. The ground experiment
therefore does not observe both lines needed by (2); the executable rejects
that channel instead of inferring a hidden peak. A different admitted
preparation or complete coefficient apparatus can test it separately.

The actual diagonal fixture returns coefficient ratios `1:2:3` in both
internal channel types. The coherent analyzer at `z=1` (phase angle zero),
pass amplitude `r=1` returns
`1:1:(34+24sqrt(2))`, including its repeated pair. These are different
configured sources. Neither inverse derives a target-fitted phase or
attenuation. The uncalibrated fixture contains no such inputs; it does not
erase the ancestry of a different calibrated construction.

For `epsilon_F=0`, rationalizing (1) gives

\[
 \omega_{j,-}=
 \frac{2\kappa s_j}{\Delta+\sqrt{\Delta^2-4\kappa s_j}}. \tag{4}
\]

At small `kappa s_j/Delta²` its leading dependence is
`kappa sigma_j(Y)²/Delta`. Thus even the leading local charged energy is
quadratic in the family singular coefficient. Treating it as the free
fixed-Higgs mass would change the observable. The exact response is
available; a guessed radiative correction is unnecessary for this finite
calculation.

## Propagating finite measurement errors

Let the certified intervals be
`Delta in [D_-,D_+]`, `omega_- in [l_-,l_+]`,
`omega_+ in [u_-,u_+]`, with `D_->0` and `u_--l_+>0`.
Put `delta_-=u_--l_+`, `delta_+=u_+-l_-`. Then

\[
 \frac{D_-^2-\delta_+^2}{4\kappa}
 \le s_j\le
 \frac{D_+^2-\delta_-^2}{4\kappa}.                    \tag{5}
\]

This rectangular enclosure needs no independence assumption. The code
requires a nonnegative lower endpoint and strictly positive original
probability lower bounds for both lines. Integer-square enclosures give
the square-root interval; division by a strictly positive reference
interval gives the coefficient ratio. Complete detector, preparation,
clock, sampling and spectral-resolution errors belong in the supplied
intervals. Equation (5) does not include an uncalculated continuum or
particle-matching error.

## The prepared reference is a sector ground

The distinction has an exact ambient-carrier control. In the same diagonal
fixture, take `u=C_1† Omega` with the scalar register set to its vacuum,
where `C_1†` is one weak component of the literal pair vertex. This vector
has `Q*=2`. Replaying the actual CAR signs gives
`||u||²=28` and `||W u||²=1764`. The normalized trial span of u and Wu has
Hamiltonian compression

\[
 \begin{pmatrix}0&\sqrt{63}\\\sqrt{63}&1\end{pmatrix},
 \qquad
 E_{\rm ambient,min}\le\frac{1-\sqrt{253}}2<-7.        \tag{6}
\]

The span is not invariant; (6) is a Rayleigh upper bound, not an extra
exact sector spectrum. It proves that the prepared 307-sector ground is
not the full 786432-state carrier ground. It does not identify a lower
physical Gauss-singlet vacuum: the trial is gauge charged, and a
compensating reference and its energy define another observation domain.
The constrained preparation and all finite readouts above remain valid.
An interacting spatial extension must specify and justify its own common
state on the sectors that its hopping actually connects.

## Reuse and the shortest physical continuation

The finite Born, instrument and Gibbs theorems apply to these complete
matrices. Repeating those derivations supplies no missing spatial momentum.
The [finite spectral readout](../spectral_readout/README.md) provides
positive Fourier-clock effects and alias/error bounds when controlled
evolution is supplied. The [matrix response construction](../matrix_spectral_response/README.md)
provides the full-residue and Schur-response method; its cube spectrum and
state are not substituted for this Higgs source. The
[free Dirac theorem](../charged_lepton_response/RELATIVISTIC_MASS.md)
gives the configured quadratic mass shell, while its
[mass-response bound](../charged_lepton_response/MASS_RESPONSE_BRIDGE.md)
propagates actual kinetic, self-energy and matching estimates when those
are supplied. The [Lorentz mass-shell theorem](../../Lean/Geometry/MassShellKinematics.lean)
is kinematics with a declared mass parameter. None of these statements
identifies (1) as the physical spectrum.

The nonredundant continuation for one observable is:

1. Specify one spatial Hamiltonian containing the chosen local vertex,
   charge-preserving transport and the claimed gauge dynamics. Choose its
   common state, regulator and clock. The [gauge-coordinate construction](FIELD_ROUTE.md)
   supplies a bounded finite route and cutoff estimates; its one-cell
   rotor is a different experiment. Preserve every sector connected by
   the actual hopping, rather than fixing local Q*=0 by deletion.
2. Calculate the same gauge-invariant or compensated probe's spatial and
   temporal response, with original residues and finite-error control.
   The neutral scalar-number channel avoids a charged dressing apparatus
   and is the smallest finite observable here. A neutral reference or
   charge assignment alone does not establish its dispersion.
3. Establish a stable relativistic branch, or the appropriate controlled
   threshold/resonance, on that same state and a stated volume/cutoff family.
   Reuse the [original-weight threshold theorem](../spin_shell_readout/THRESHOLD_READOUT.md)
   only after supplying its near-edge spectral-weight premise. Constant
   kinetic whitening alone does not fix the spatial speed. Chiral matter
   assignments require the corresponding regulator/mirror control.
4. Attach the observed spin, currents and interactions to one physical
   species, then compare its mass or energy ratio with an independent
   measurement in the matching convention. A scalar conversion line is not
   automatically the Higgs; a zero frequency is not a transverse photon;
   internal lepton-type weights are not three measured charged leptons.

A comparison of one well-defined observable does not require a proof of
architectural uniqueness, a complete Standard Model catalogue or Einstein
dynamics. It does require these same-source propagation and identification
statements. The finite construction provides no independent laboratory
lepton, photon or Higgs mass at its stated domain.

```sh
python -m model.particle_observer_outputs.bare_effective
python -m pytest -q model/particle_observer_outputs/test_bare_effective.py
```
