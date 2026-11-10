# Phase spectroscopy on the interacting defect-wave source

## The measured quantity and the unchanged source

Let U be the [joint defect-wave tick](../PROOF.md), with the same link words,
wave registers, five controlled repairs and clock. This experiment reads the
spectral phases of U. It does not define a Hamiltonian by relabelling a
relaxation operator, replace U by an effective generator, or assign a rest
mass to a clock label.

The source begins in `|g=1> tensor A|port0,internal3>`. The complete interior
Gauss instrument has positive and negative outcomes with probabilities
`1/2,1/2`. The negative outcome is retained as a preparation failure. The
positive outcome proceeds through the experiment below. Write psi for its
normalized state and mu for its spectral measure under U:

`c_t = <psi,U^t psi> = integral exp(i t theta) d mu(theta)`.

The initial preparation is an owned finite control. It is not the faithful
stationary reference or a proposed vacuum. All statements about its spectral
weight retain that preparation domain.

## Complete clocks and a filter that removes exact flat phases

For an integer M define the complete Fourier clock

`K_l = (1/M) sum_(t=0)^(M-1) exp(-2 pi i l t/M) U^t`, `0<=l<M`.

Orthogonality of the finite characters gives

`sum_l K_l* K_l = (1/M) sum_t U^(-t) U^t = I`.

Its implementation prepares a uniform clock register, coherently applies
U^t conditional on the clock value, then reads every inverse-Fourier pointer
outcome. Each controlled U is compiled from the same owned gate product.
Controlled source gates and the clock register are declared optional apparatus;
they do not change the idle law or its tick duration. The executed clock has
M=4, so its only Fourier phases are `1,i,-1,-i`.

Before this clock apply the complete instrument

`F_even=(I+U^2)/2`, `F_nonflat=(I-U^2)/2`.

Unitarity gives `F_even* F_even+F_nonflat* F_nonflat=I`. If `U v=+v` or
`U v=-v`, then `F_nonflat v=0`. Thus the active outcome excludes both exact
flat eigenphases, including any auxiliary register modes at those phases.
It does not assume that all other modes represent particles.

All eight fine success-branch Kraus maps are `K_l F_b P_+`, where P_+ is
the Gauss projection. The ninth terminal map is P_-, the retained Gauss
failure. Their effects sum to

`P_- + P_+ sum_b F_b* (sum_l K_l* K_l) F_b P_+ = I`.

Every successful branch commutes with the interior gauge transformation,
since U does. The same proof tensors with the identity on every untouched
reference. The finite instrument theorem holds on the complete source space.
The implementation executes rational complex amplitudes, using separate exact
real and imaginary parts to avoid floating-point arithmetic.

The controlled clock requires four clock states, one filter bit and one Gauss
pointer bit; there are nine terminal labels because the Gauss failure stops
the experiment. These apparatus costs are separate from the source register
and its exact-arithmetic storage. The maximum source power actually executed
is U^5, arising from the quadratic filter and four clock values.

## Exact populated phase response

Literal execution on the carrier gives

| t | c_t |
|---|---|
| 0 | `1` |
| 1 | `0` |
| 2 | `55/72` |
| 3 | `-9571/54000` |
| 4 | `279031/540000` |
| 5 | `-19868809/48600000` |

The active filter probability, conditional on the positive Gauss preparation,
is

`p_nonflat = (1-Re c_2)/2 = 17/144`.

Its original-trial probability is `17/288`; the Gauss failure, complementary
filter and every Fourier bin remain in the
[complete receipt](results.json). The histogram is computed from actual
output vectors and independently checked against the moment formula

`p_l = 1/M + (2/M^2) Re sum_(t=1)^(M-1) (M-t) exp(-2 pi i l t/M) c_t`.

This displays a populated, non-flat frequency response of the interacting
program. It is stronger than assigning spectral labels to unpopulated modes.
It supplies no isolated particle pole, conserved momentum or physical energy
normalization.

An exact spectral support bound follows without diagonalization. Let
`x=sin^2(theta)`, so `E_mu x=17/144`. If `w=mu(x>=1/16)`, then

`17/144 <= (1-w)/16+w`, hence `w>=8/135`.

The unnormalized active measure is `d nu=x d mu`, and its weight on this
same event satisfies

`nu(x>=1/16) >=17/144-1/16=1/18`.

Consequently at least `1/36` of all original trials lies in the active
spectral event `|sin(theta)|>=1/4`, with the Gauss success weight included.
At least `8/17` of the active conditional spectral weight lies there. This
event is separated from both 0 and pi by `arcsin(1/4)` in circular phase.
It is a support bound for this preparation, not a uniform mass gap.

Replacing the tick by identity makes the non-flat filter identically zero.
That separately labelled control tests whether the filter is merely imposing
a positive answer. The actual positive weight follows from U's executed
second moment.

## A disturbance of an explicit stationary reference

The [background preparation](../background/PROOF.md) supplies a normalized
positive-Gauss state Omega_0 satisfying `U Omega_0=Omega_0`. Its entries are
uniform over the ten link words, twelve ports and sixty addresses, with
internal component 0 and parity 0. This is an explicit pure stationary
reference; stationarity does not identify a physical vacuum or ground energy.

Apply the source's complete five-outcome face measurement to Omega_0. Each
output `P_f Omega_0` has squared norm `1/5`. Apply both filters and the complete
four-bin clock to every face outcome. The forty retained records sum to one
per prepared-background trial. Conditional second moments and active weights
are:

| Face index | Conditional c_2 | Active weight per background trial |
|---|---|---|
| 0 | `107/120` | `13/1200` |
| 2 | `29/30` | `1/300` |
| 5 | `107/120` | `13/1200` |
| 16 | `29/30` | `1/300` |
| 17 | `641/720` | `79/7200` |

The active weights sum to `283/7200`. Every active output is orthogonal to
Omega_0: `F_nonflat* Omega_0=0`, and the polynomial clock commutes with the
filter. Thus a complete gauge-invariant apparatus prepares and reads
nonstationary disturbances relative to a stationary state of the same source.
This reference-relative statement uses no energy-log branch or species name.
The background preparation has its own retained failure records; the forty
histories here are conditional on that preparation's success and do not
replace its original-trial probabilities.

## Reused theorem interfaces and resolution

The construction reuses the complete Fourier instrument from
[joint spectroscopy](../../particle_primary_joint_spectroscopy/PROOF.md)
and the source-independent summed Fejer bound in
[spin-shell readout](../../spin_shell_readout/PROOF.md). Their Hamiltonians,
momentum operators, states and species assignments are not transferred.
For a circular turn margin d and M bins, the scalar tail bound is

`ell(M,d)=min(1,1/[2(Md)^2]+1/[2Md])`.

The existing exact planner selects M=1024 for margin `1/64` turn and
requested tail error `1/16`, giving the sharper bound `17/512`. This is a
resource/error plan for an extended apparatus; its clock alone uses powers
through 1023, and including the quadratic filter requires powers through 1025.
The retained experiment
executes four bins; it does not claim the 1024-bin histogram was computed.
The tail bound is circular. An unwrapped energy requires an independently
selected generator and alias-free energy interval.

The finite theorem
[`finite_unitary_powers_have_no_limit`](../../../Lean/QFT/FiniteUnitaryScatteringNoGo.lean)
also applies: this nontrivial finite U cannot have an ordinary operator limit
as its powers tend to infinity. This does not exclude scattering defined
through relative dynamics or a controlled growing family. It excludes using
settling of this fixed finite U as the scattering proof.

## Why phase readout alone does not select a mass

For `U=sum_a exp(-i theta_a) P_a`, all operators

`H_n=(hbar/tau) sum_a (theta_a+2 pi n_a) P_a`

with integer n_a give the same U. Every experiment built from integer U
ticks and fixed source instruments has the same record law for those
generators. Even positivity and a common ground vector need not select a
branch: the two positive matrices

`H_0=diag(0,pi/2)`, `H_1=diag(0,5pi/2)`

both exponentiate to `diag(1,-i)` at unit time, yet their positive gaps
differ by a factor of five. Their polynomial-in-U stationary references also
agree. This two-state countercontrol is a logical boundary; it is not an
alternative claimed spectrum of the full defect-wave program.

The explicit stationary references give a stronger countercontrol on the
actual joint source. Let H_0 be its spectral logarithm with eigenvalues in
`[0,2 pi)`, so `exp(-i H_0)=U`, and let `P_0=|Omega_0><Omega_0|`. For positive
integers k define

`H_k=H_0+2 pi k (I-P_0)`.

The two terms commute, and `exp(-2 pi i k(I-P_0))=I`, so every H_k produces
the same actual source tick. Every H_k is positive and has the same unique
ground vector Omega_0. All orthogonal vectors have energy at least `2 pi k`;
the independently checked orthogonal stationary state Omega_1 has exactly
that energy. Thus k=1 and k=2 give exact gaps `2 pi` and `4 pi`, respectively,
on the same source, state, stationary reference and tick. These are nonlocal
spectral-log countercontrols, not proposed physical generators. A local
continuous-clock or refinement condition may exclude them; positivity,
stationarity and a unique ground vector alone do not.

A physical energy attachment needs a selected continuous clock/generator or
a controlled small-step family, an appropriate reference and energy zero,
and populated momentum-dependent spectral responses. Relativistic dispersion,
spin/statistics, interactions and a species dictionary have separate content.
Unit conversion does not choose an energy branch or identify a particle.
