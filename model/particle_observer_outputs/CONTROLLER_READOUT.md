# Read the coherent controller from complete finite spectra

The configured coherent source has an observable attenuation and phase class.
Its complete interacting spectral records recover the three family
coefficients. Their moments recover the pass amplitude `r` and `cos(3 theta)`.
The derivation uses the declared source normalization, or the positive-root
trace identity to remove one common scale. It contains no measured particle
mass. Reading a configured controller does not select that controller through
source dynamics, and does not erase the ancestry of a target-calibrated input.

## Reused source and record construction

The [native coherent analyzer](../standard_model_source_map/KOIDE_SOURCE.md)
gives the exact branch weights `(10,5,5,2,50)/72` and a regular three-cycle
`R`, with `R^3=I`, `Tr R=Tr R^2=0`. The
[complete attenuator](../particle_candidate_atlas/BARE_EFFECTIVE_CONTROLLER.md)
specifies a pass amplitude and unit phase while retaining all sixteen flags.
Use its literal coefficients

\[
a=\sqrt5/6,\qquad b=a/\sqrt2,\qquad
C=aI+rb(zR+\bar zR^2),\qquad Y=C^2/3,\quad z=e^{i\theta}.
\tag{1}
\]

The [charged source construction](CHARGED_HIGGS_RECORDS.md) and
[coefficient inverse](BARE_EFFECTIVE_BOUNDARY.md#a-coefficient-inverse-that-cancels-the-charged-offset)
use the same finite vertex, state, neutral reference and clock. For every
resolved charged channel with both original addition weights positive,

\[
\sigma_j(Y)^2=\frac{\Delta^2-(\omega_{j,+}-\omega_{j,-})^2}{4\kappa},
\quad \kappa=1\text{ or }2.\tag{2}
\]

The charged additive offset cancels. Multiplicities are retained: a family
eigenvalue of rank `m` appears `4m` times in the left-doublet channel and `2m`
times in the charged-singlet channel. The two channels independently recover
the same three coefficients. The complete probe retains every addition and
removal event and the conjugate charged reference. Formula (2) selects which
record coordinates enter an estimator; it does not condition the original
experiment on an addition event or discard its other outcomes.

[`controller_readout.py`](controller_readout.py) composes the owning
`response_lines` and `channel_reconstruction` functions. It consumes their
full-CAR Gram theorem without repeating the 1650-column source calculation.
The executable control at `r=1,z=1` retains the repeated family eigenvalue;
its full original probe probabilities sum to one.

## Exact moment inverse

The source's [moment proof](../particle_calibration_selection/PROOF_PHASE_ATTENUATION.md#exact-coherent-moments)
gives, with `T1=Tr Y` and `T2=Tr Y^2`,

\[
T_1=\frac5{36}(1+r^2),\qquad
T_2=\frac{25}{3888}+\frac{25r^2}{648}+\frac{25r^4}{2592}
    +\frac{25\sqrt2\,r^3}{1944}\cos(3\theta).
\tag{3}
\]

Consequently,

\[
r^2=\frac{36T_1}{5}-1,\qquad
\cos(3\theta)=
\frac{1944\left(T_2-25/3888-25r^2/648-25r^4/2592\right)}
     {25\sqrt2\,r^3}\quad(r>0).
\tag{4}
\]

The code rejects `r^2` outside `[0,1]` or a cosine outside `[-1,1]`. At `r=0`,
`Y=(5/108)I` and `T2=25/3888`. Every phase gives the same `Y`; these spectral
records contain no phase information. This statement concerns `Y` and its
interacting spectral readout. It makes no claim that every possible coherent
observable of the analyzer's retained loss arms is phase independent.

For `r>0`, the unordered spectrum determines

\[
\phi=\tfrac13\arccos(\cos 3\theta)\in[0,\pi/3],\qquad
\theta=\pm\phi+2\pi k/3\pmod{2\pi}.
\tag{5}
\]

Phase reversal exchanges two family coordinates; the `2pi/3` shifts permute
the three coordinates. An unordered spectrum therefore cannot distinguish
these settings. Calling `phi` the absolute phase requires an additional
choice of principal chamber. An oriented, phase-sensitive apparatus would
be needed to remove this remaining equivalence.

## Positive-root admission and the common-scale nuisance

The full-spectrum inverse admits the positive-`C` branch. In addition to the
moment constraints it checks all three positive coefficients and

\[
\sum_j\sqrt{3\sigma_j(Y)}=3a.\tag{6}
\]

Checking only the first two moments is insufficient to authenticate an
arbitrary triple as this source family. Equation (6) is the missing branch
constraint. For completeness, set `c_j=sqrt(3 sigma_j)`. Equation (6) places
these roots in the plane `sum c_j=3a`; equation (3) places them on the circle
`sum(c_j-a)^2=6r^2 b^2`. Every point on that circle has the form
`c_j=a+2rb cos(theta+2pi j/3)`, up to a family permutation. Thus the plane,
radius and positivity conditions give membership in this analytic source
family, rather than only agreement with two moments. At `r=0` the circle
reduces to its centre. Controllers with a negative `C` eigenvalue exist;
they are outside this inverse's declared branch.

The source parser admits exact rational `r` and exact rational real and
imaginary parts of `z`. The general inverse describes the analytic family;
its matrix controls use parameters admitted by that parser. The nonzero
control is `r=1/2,z=(63+16i)/65`, giving
`cos(3theta)=201663/274625`. No parser extension is inferred from the analytic
identities. The optional analytic example `theta=pi/12` is not an admitted
literal parser input.

An absolute application of (4) needs the normalization `a=sqrt(5)/6` in the
same source clock as (2). If the reconstructed coefficients are
`y_j=s sigma_j(Y)` with a known positive common scale `s`, divide by `s`
first. Physical seconds or GeV are unnecessary.

There is also a ratio-only route, conditional on the same positive-root
source family. Equation (6) fixes the otherwise unknown common scale:

\[
s=\frac{(\sum_j\sqrt{y_j})^2}{3a^2},\qquad
Q=\frac{\sum_j y_j}{(\sum_j\sqrt{y_j})^2},\qquad
r^2=3Q-1.\tag{7}
\]

The normalized coefficients `y_j/s` then enter (4). This is source
normalization from a stated structural premise, not an identification of
arbitrary measured masses with `Y`. The
[finite Koide theorem](../../Lean/ObserverPatchHolography/KoideCirculant.lean)
supplies `koide_quotient_phase_independent`; here its quotient is
`Q=(1+r^2)/3`, so `Q=2/3` corresponds to `r=1`. That phase-independent identity
cannot choose a phase. The second moment supplies phase-class information.

## What this closes and what it requires

The configured controller is reconstructible from its admitted complete
finite source spectrum, with the stated discrete phase ambiguity. No
particle-data input is needed for this reconstruction. This is compatible
with [classical-record nonidentifiability](../particle_source_selection/PROOF.md):
the interacting coherent apparatus supplies additional operational data
beyond the classical committed record law.

The source architecture may choose its controller. If a controller was
chosen using target particle data, observing that choice through its
resulting spectra preserves that calibration ancestry. Neither (2) nor (4)
derives the target choice from repair dynamics. The finite coefficients also
need their spatial, relativistic and physical-observable matching before
they can be called measured lepton masses.

These exact tests do not turn finite-error spectral estimates into an exact
phase estimate. The owning coefficient inverse provides conservative
intervals for source coefficients. A phase interval must propagate those
intervals through (4) and retain visibility and common-clock errors. The
`r^-3` factor makes phase recovery ill-conditioned near `r=0`; `acos` is
also ill-conditioned at its endpoints. The module exposes exact source
controls and rejects inadmissible inputs; it does not claim a finite-error
phase certificate from interval midpoints.

Run `python3 -m pytest -q model/particle_observer_outputs/test_controller_readout.py`.
Run `python3 -m model.particle_observer_outputs.controller_readout` for the
machine-readable exact certificate.
