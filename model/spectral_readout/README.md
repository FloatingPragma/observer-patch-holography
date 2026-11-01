# Complete finite-resolution energy readout on the original ground response

The positive spectral measure of the original interacting ground admits a
complete, positive, finite-resolution energy measurement. Finite Fourier
samples determine all its outcome probabilities. The source calculation
below evaluates this measurement on the **actual previously certified
78-line cube measure**, retaining every outcome. It does not replace the
measure by fitted poles, label filter centers as masses, or infer the
infinite-ground numerical spectrum from the cube.

For the selected infinite ground, the same construction has weak spectral
convergence and strong convergence of its **compressed resolvents** at
fixed spacing. Exact sum rules give numerical energy-tail bounds. A
finite-to-infinite numerical comparison still requires the ground-state
discrepancy that the existing compactness construction does not compute.
The code leaves that term unknown unless explicitly supplied as a
conditional premise.

## Positive spectral measure and exact sum rules

Use the already proved isometry `V`, signed doubled generator
`L=L_omega directsum -conjugate(L_omega)`, and contact map `J` from
[the spatial response theorem](../retarded_spatial_locality/README.md).
For a unit vector `f` in the one-particle lattice space define

```
mu_f(B)=<Vf,1_B(L)Vf>,       phi_f(t)=integral exp(-itE) dmu_f(E).
```

This is a probability measure on real **signed addition/removal energy**.
The positive generator in each summand describes excitation energy above
the ground; the signed convention is the retarded response convention.
For the cube down-addition probe there is no removal component, and all
78 energies are positive.

There is a complete charge-tagged interpretation. For a normalized CAR
mode, the Kraus pair `K_add=c(f)^dagger`, `K_remove=c(f)` obeys
`K_add^*K_add+K_remove^*K_remove=I`. In each output branch measure the
**nonnegative excitation energy** of the positive ground generator, and
retain the addition/removal label. Assign its sign only when forming the
response-frequency variable. Negative removal frequency therefore never
means negative physical excitation energy. An exact two-mode control has
addition cost 3 and removal cost 2, each with probability one half; its
signed response has frequencies 3 and -2. This mathematical instrument is
not claimed as a new executed native long-time measurement.

In the selected hidden-SU(2) invariant infinite ground,
`V*J=J*J=rI`, with `1/8<=r<=1/2`, and `LV=Vh+UJ`. Consequently

```
integral dM(E)=I,
integral E dM(E)=h+UrI,
integral E^2 dM(E)=h^2+2Urh+U^2rI,
```

where `M(B)=V*1_B(L)V`. The second-moment identity is a quadratic-form identity
obtained from `(LV)*(LV)`; it does not assume that `Vf` lies in `Dom(L^2)`.
For every unit `f`, the first two scalar moments are finite and

```
mu_f(|E|>=R) <= <f,(h^2+2Urh+U^2rI)f>/R^2.
```

For a one-mode local probe, `<f,hf>=0` and `||hf||^2=3c^2/2`, because
there are six incident half-strength hopping terms. At `c=1/8,U=1`,

| Quantity on the selected infinite ground | Exact interval or bound |
| --- | --- |
| Mean signed energy | `[1/8,1/2]` |
| Second moment | `[19/128,67/128]` |
| Variance | `[17/128,35/128]` |
| Probability of `|E|>=2` | at most `67/512` |
| Probability of `|E|>=4` | at most `67/2048` |
| Probability of `|E|>=8` | at most `67/8192` |

These are state-uniform conclusions for the specified selected-ground
class, not computed individual poles. For a smeared probe with
`||hf||<=K`, the sharper uniform second-moment bound is
`K^2+UK+U^2/2`; the simpler `(K+U)^2` also holds. Thus any continuum
smearing construction with uniformly bounded kinetic norm has uniformly
tight energy measures. This statement does not retain a nonzero
connected four-point function in a continuum limit or prove a mass shell.

## A complete positive finite-resolution POVM

Choose an integer `M>=2` and a time step `tau>0`. For
`k=0,...,M-1`, put

```
K_k(E)=| M^(-1) sum_(j=0)^(M-1) exp[-ij(tau E-2pi k/M)] |^2,
q_k=integral K_k(E) dmu_f(E).
```

Every `K_k` is nonnegative. Finite root-of-unity orthogonality gives
`sum_k K_k(E)=1` pointwise: in the squared sum, `j-l` lies strictly
between `-M` and `M`, so only `j=l` survives the sum over outcomes.
Consequently the effects `K_k(L)` form a complete POVM and the `q_k`
are a normalized probability vector. No outcome is dropped or
postselected. Its finite Fourier formula is

```
q_k=1/M + (2/M) sum_(n=1)^(M-1) (1-n/M)
                     Re[exp(2pi i kn/M) phi_f(n tau)].
```

A finite system with a supplied controlled evolution can implement this
POVM by preparing a uniform clock register, applying the controlled
powers of `exp(-i tau L)`, applying the Fourier analyzer and measuring
all clock outcomes. Alternatively its probabilities are determined by
both quadratures at the `M-1` displayed times. This is a measurement
construction, not a claim that the infinite GNS generator is an available
native gate or that those additional native executions have occurred.

Assign principal energy centers
`E_k=2pi k/(M tau)` for `k<M/2`, and
`E_k=2pi(k-M)/(M tau)` otherwise. These are **readout centers**, not
asserted eigenvalues. The energy period is `2pi/tau`: energies differing
by that amount have identical sampled data. Increasing `M` alone cannot
remove this aliasing.

## Resolution and alias bounds

Let `0<tau Delta<=3<pi`, and use circular distance modulo `2pi/tau`.
For angular distance `d>=tau Delta`,

```
K_k(E)<=1/[M^2 sin^2(d/2)]<=pi^2/(M^2 d^2).
```

Summing the worst case over all `M` outcomes gives
`10/(M tau^2 Delta^2)`, since `pi^2<10`. A sharper bound sums the
actual grid tails. Write `D=M tau Delta/(2pi)`. Along either directed
arc the scaled distances form a unit-spaced sequence whose first entry
is at least `D`. For `D>1`, each effect is at most `1/(4q^2)`, and

```
sum_(j>=0) 1/(D+j)^2 <= 1/(D-1).
```

There are at most two arcs. Hence the circular error probability is at
most `1/[2(D-1)]=pi/(M tau Delta-2pi)`. Using `pi<22/7`, the code
returns the minimum of one, the first bound, and
`22/(7M tau Delta-44)` when the denominator is positive.

To control ordinary real energy, exclude true energies within `Delta`
of a Nyquist edge. The rational interior `|E|<=3/tau-Delta` lies inside
that safe region. If the second moment is at most `C`, its complement
has probability at most `C/(3/tau-Delta)^2`. This alias term is zero
when an independently certified entire support lies in that interior.
Adding it to the circular bound controls
`P(|E_readout-E_true|>=Delta)`. This probability refers to the joint
measure `K_k(E)dmu_f(E)`, so it also bounds every interval-bin error:

```
mu_f(I eroded by Delta)-error <= P(E_readout in I)
                            <= mu_f(I enlarged by Delta)+error.
```

For a bounded Lipschitz function `g`, the expectation error is at most
`Lip(g) Delta+2||g||_infinity error`. Thus taking `tau` small enough to
control the tails and then `M` large enough constructs any desired
finite weak energy resolution. A pointwise spectral density or an
isolated pole is not needed, and neither is concluded.

## Executed source calculation on all 78 finite spectral lines

`finite.py` authenticates the complete prior source and its exact moments
`1,1/2,131/256`. It evaluates the Fourier expression using only rational
arithmetic, outward dyadic intervals, and bounded Taylor/trigonometric
identities. Argument halving and double-angle identities handle long
times. An independent direct sum of squared Fourier amplitudes computes
the same probabilities and the actual joint resolution error.

For `tau=1`, the complete 16-outcome readout is approximately:

| Readout center | Probability |
| ---: | ---: |
| 0 | 0.483266606026 |
| 0.392699081699 | 0.083595687954 |
| 0.785398163397 | 0.025139421111 |
| 1.178097245096 | 0.358720295560 |
| 1.570796326795 | 0.011795465706 |
| 1.963495408494 | 0.003949923856 |
| 2.356194490192 | 0.002344166059 |
| 2.748893571891 | 0.001729576088 |
| -3.141592653590 | 0.001459060757 |
| -2.748893571891 | 0.001364393614 |
| -2.356194490192 | 0.001399012671 |
| -1.963495408494 | 0.001574751882 |
| -1.570796326795 | 0.001971304202 |
| -1.178097245096 | 0.002833022115 |
| -0.785398163397 | 0.005037196832 |
| -0.392699081699 | 0.013820115573 |

The receipt retains outward rational intervals for every probability and
center, for `M=8,16,64`. At `Delta=1/2`, the actual joint error lies in
`[0.100284517427,0.100284517428]`,
`[0.055193158231,0.055193158232]`, and
`[0.020209031433,0.020209031434]`, respectively. The finite source has no
alias term here because all its energy is below `1.85<3-Delta`.
Negative reported centers are finite-filter leakage from this positive
source; the bin centered at zero does not prove a zero-energy particle.

The already executed complete native Ramsey experiment at time one
agrees with the first harmonic used here, within its retained error.
Higher harmonics through times `7,15,63` are source calculations in this
packet. There are no new native long-time executions or finite-shot
samples, and no physical energy-unit calibration.

## Finite-volume, native and statistical errors

Let `epsilon` bound the complex error in every sampled characteristic
function. The Fourier formula bounds one bin's error by
`(M-1)epsilon/M`. For an event containing `b` bins, use
`min(b,M-b)(M-1)epsilon/M`, capped by one; the empty and complete events
have exactly zero error.

Independent noisy Fourier coefficients need not reconstruct nonnegative
probabilities. The bound capped by one applies to two actual probability
laws, or to a scalar event estimate after clipping it to `[0,1]` using
`clip_event_estimate`. Clipping cannot increase its error from the true
event probability. The code does not assert that a raw noisy Fourier sum
is already a probability law, or normalize away any native failure outcome.
The complement bound uses the exact Fourier identity that all raw bins
sum to one before this scalar clipping.

For normalized local CAR probes, the original averaged-dynamics bound
contributes at most twice its operator error to the anticommutator
characteristic function. Add the separately supplied native complete
instrument error, the statistical error, and

```
epsilon_state = max_(sampled t)
 |nu_L({tau_t(a),a^dagger})-omega({tau_t(a),a^dagger})|.
```

The state discrepancy tends to zero along the selected subsequence but
has no supplied numerical rate. `comparison_contract` returns an unknown
total when this input is absent; a supplied number is marked explicitly
as a conditional premise. It never substitutes zero for that missing
information. It also does not identify the finite cube ground with a
large-box ground.

For binary Ramsey quadratures, probability error at most `e` in each
quadrature gives complex characteristic error at most `2sqrt(2)e<3e`.
With `N` independent shots per quadrature, Hoeffding and a union bound
give failure probability at most `4(M-1)exp(-2Ne^2)` over all sampled
times. The calculator bounds this rationally without generating samples.
Native control errors refer to complete instruments, including failure
outcomes, rather than success-conditioned distributions.

## Spectral convergence on the actual infinite ground

Construct each finite comparator as the direct sum over its common
translated full-ground labels, including the conjugate removal summand.
It has an isometry `V_L` from the same lattice one-particle space and a
self-adjoint signed generator `L_L`. This uses one label throughout each
correlation, not a product of separately averaged states.

The existing local-dynamics convergence and weak-star state convergence
give convergence of every local characteristic function, uniformly on
compact time intervals, along the same selected ground subsequence.
Normalization and the uniform second-moment bound
`||L_L V_L||<=sqrt(3)c+U` give tightness. Therefore the positive scalar
spectral measures converge weakly, for each fixed local probe; density
extends this to all square-summable probes. All continuous finite-resolution
readouts above converge. A sharp interval probability converges when its
endpoints carry no mass in the limiting probe measure. Weak convergence
alone cannot certify an atom or passage through such an endpoint.

One can strengthen the fixed-spacing resolvent statement without
identifying the varying GNS Hilbert spaces. Put
`G_L(z)=V_L^*(z-L_L)^(-1)V_L`, acting on the common one-particle lattice
space. Its matrix entries converge to `G(z)` in either half-plane. The
full-CAR finite-range bound applies uniformly to every translated finite
Hamiltonian, including its zero exterior. Hence the exponentially decaying
kernel bound from the spatial response proof holds uniformly for these
finite comparators. For finitely supported input vectors, truncate the
output to a large finite region: the tail is uniformly small by that
bound and the remaining finitely many entries converge. Uniform operator
bounds then give `G_L(z)f->G(z)f` for every `f`. The adjoints converge by
the same argument in the other half-plane.

This is **strong convergence of the compressed resolvents on the common
one-particle space**. It is not an assertion of ordinary strong-resolvent
convergence of `L_L` on a single unconstructed common GNS Hilbert space.
It supplies no convergence rate beyond the explicit dynamics terms and
the retained unknown state discrepancy. Together with the exact positive
kernels and calculated finite distributions, it closes a controlled
finite-resolution energy-readout interface; the continuum particle and
gravity identifications remain separate requirements.

## Replay

```
.venv/bin/python -W error -m pytest -q model/spectral_readout
.venv/bin/python -W error -m model.spectral_readout --write
.venv/bin/python -W error -m model.spectral_readout --check
```

The cold receipt binds this package, its proof fragment, and the transitive
source manifests of the prior exact spectrum, native first harmonic, and
infinite-ground response theorem. Tests include aliasing, a discontinuous
sharp-projection counterexample, charge-sign interpretation, negative raw
statistical estimates, and independent positive-kernel evaluation.
