# What two material cameras see when the source forms freely

One mobile excess meets two marked heads on the 80-cell Interval Echo board. Four units of free evolution leave correlations in the heads after the mobile is removed. At the gateway between those heads, a vacancy and an excess both see the correlations, and each camera gives the same image when its emitter and receiver trade places. Yet they disagree even on the **sign** of the image-area response calibrated with ordinary exchange energy. Near the fully formed source, the vacancy reports **+29.66** ordinary units per unit of correlation intervention and the excess reports **−1.97**. Transverse axis conventions, fixed linear distortions and probe-specific scalar gains cannot explain the full matrix discrepancy. These are model measurements of material beam images. Reading them as one proper curvature requires an additional instrument or a derivation of how the two probes couple.

## How the same source and its comparison are made

Use the existing 80-cell federation with strength-eight marked contacts at
$(0,1)$ and $(38,39)$. Prepare the usual positive mobile at port zero and
the registered starting head coins. Evolve the entire board for four units,
then remove the mobile without selecting its position. Call the resulting
two-head density $\rho$. Keep the left and right heads from **different**
independent copies of the same encounter to produce
$\tau=\rho_L\otimes\rho_R$. Both stocks have exactly the same one-head
states, even reversal and parked energy 16. A source-side coin mixes them as
$C_s=(1-s)\tau+s\rho$. The receiver receives no preparation labels.
The incoming mobile, two material marks, initial coins, exposure, removal
and later camera appointments are controlled. The interaction during formation
is free. Translation of the stored heads and a measurement of their inertia
are outside this experiment. [Formation and partial trace](../../../full_algebra/native_precision.py) · [Original mobile experiment](../../../../../../model/experiments/interval_echo_einstein_station/mobile_collar_source/README.md)

The head entropies are $S(\rho)\simeq0.45146624$ and
$S(\tau)\simeq0.49319119$ nats. The encounter stores
$I(L:R)=S(\tau)-S(\rho)\simeq0.04172496$ nats of mutual information.
For $K_X=-\log X$, the entropy slope in the direction $\rho-\tau$
at the product is exactly zero because $K_\tau$ is the sum of its one-head
generators. At the formed source the slope is
$-I(L:R)-D(\tau\Vert\rho)\simeq-0.29384843$ nats. The area responds at both
endpoints. The projected camera area cannot supply a universal measurement of
this modular slope on this preparation family. [Source logarithm and enclosures](certificate.py) · [Independent numerical source readout](measure.py)

## Two independent material probes and an ordinary-energy calibration

Each camera sends one chosen probe charge from a pixel on one five-port
ring. After four units, it resolves the opposite ring and compares the
two arrived coins with an **independent** source preparation, using the two
registered local Bell detectors. All starts enter the ledger, including
missed pixels and failed coincidences. Both species, excitation and vacancy,
use the same physical comparison rails and source state. Their charge is
set by the sender. A third calculation mixes their raw pixel probabilities
with equal weights before computing the image.
[Ports, storage and read gates](../reciprocal_screen/live_reference/README.md)

For the independent ordinary calibration, unread left-sign dephasing raises
the odd-reversal fraction $f$ from zero, preserving the conditional two-coin
state. The parked meter sees $16+f$ in native energy units. Enroll every
parity outcome in the coin test. The original parity-agreement veto would
weight even/odd hits by $(1-f)^2,f^2$, making the first-order normalized
image response zero at $f=0$. Here reversal controls pulse timing only.
The stored target coin is in its even storage register after transfer.
The compiled Bell effect projects onto the same coin pair for both
reference parities. With $K_{ji,q,\pi}$ the native flight block, the actual
matrix of hit probabilities is

\[
P_{ji}(s,f)=\tfrac14\sum_{\pi=e,o}w_\pi
 \operatorname{tr}(C_s^{\mathsf T}K_{ji,q,\pi}C_s K_{ji,q,\pi}^\dagger),
\quad (w_e,w_o)=(1-f,f).
\]

The native generator is real symmetric, so the reverse block is its
transpose. Transposition and cyclicity of the trace give the same hit
probability in reverse for each populated sector. This proof uses the
registered stocks with the same conditional coin state in both sectors.
Source and reference are independent draws, including their reversal
lotteries. [Native parity tests](../reciprocal_screen/live_reference/apparatus.py) · [Parked meter](../../../null_operator/measure.py)

## A comparison that survives changes of station frame

Use the existing two-dimensional pentagonal bearings $b$ on each local
ring. For uniform incident pixels the centered image matrix is
$J=2\operatorname{Cov}(b_{\rm det},b_{\rm emit}\mid\mathrm{hit})$.
Keep the two source derivatives $F=\partial_s J$, $O=\partial_f J$
on the **same** source and camera. Since the input and reference are
independent, $P(s,0)$ is quadratic in $s$ and $P(s,f)$ is linear in
$f$. Exact endpoint slopes follow from the three stocks at $s=0,1/2,1$.
Writing these probability matrices as $P_0,P_h,P_1$, the slopes are
$4P_h-P_1-3P_0$ at the product and $3P_1-4P_h+P_0$ at the formed
stock. Both point in the increasing-$s$ direction, with one-sided limits
at the physical endpoints. The independent ordinary slope is
$2[P(s,1/2)-P(s,0)]$. Differentiating the image after constructing these
probability slopes retains the throughput normalization.
[Image and native bearings](../../../beam_response/measure.py) · [Interval propagation](bounds.py)

The area-equivalent ordinary units are
$(\partial_s\det J)/(\partial_f\det J)$. Both numerator and denominator
acquire the same factor under independent invertible linear changes of
emitter and receiver axes. This ratio also applies to unsigned image area
at a nonsingular baseline. For the full matrix, $J\mapsto LJR$ gives
$O^{-1}F\mapsto R^{-1}(O^{-1}F)R$. The same pencil is obtained after
adding any fixed matrix offset to $J$. Its determinant keeps its sign
even if each instrument has a different nonzero scalar gain on either
derivative, since a scalar enters its determinant squared. A reversal of
the viewing direction gives the same determinant and trace. The flight
duration and enrolled port apertures are held fixed. Changing either is
a different measurement.

| Gateway | Equivalent ordinary units near the product, vacancy / excess | Equivalent ordinary units near the formed source, vacancy / excess |
| --- | --- | --- |
| $(1,38)$ | $+1.5283$ / $-14.2920$ | $+29.6578$ / $-1.9670$ |
| $(3,192)$ | $[0.005,0.031]$ / $[0.133,0.187]$ | $[-0.002,0.024]$ / $[0.088,0.139]$ |
| $(10,56)$ | $[-0.368,-0.367]$ / $-0.02063$ | $[-0.289,-0.287]$ / $-0.01467$ |
| $(22,44)$ | $[0.038,0.039]$ / $[-8772,-3236]$ | $[0.158,0.160]$ / $[1373,1733]$ |
| $(27,36)$ | $[-2.56,0.05]$ / unresolved | $[-5.01,-1.14]$ / unresolved |

These are the five distinct gateways incident to the two source-owner cells.
Ranges are outward enclosures from the native propagator. Single numbers are
rounded. Weak ordinary calibration slopes can leave large or unresolved ratios.
Comparing different gateways alone allows spatial variation. The decisive
comparison uses the **same** gateway, two identical source interventions and
separately calibrated probes. At the formed source the pencil determinants
are $-1.27585$ and $+3.88597$ at $(1,38)$. The same strict sign disagreement
occurs at $(22,44)$. Pooling charge labels before forming an image produces
a well-defined protocol-dependent response. At $(1,38)$ it reads $-3.27140$
ordinary units near the formed source. These observations constrain this
camera as a candidate common metric. They do not determine continuum
curvature, an ordinary-matter coupling to it, or a galactic mass density.
[All five native tables, invariants and pooled images](certificate.json)

## Numerical and sampling checks

The native certificate uses the integer formation column and four-tick
Chebyshev propagator. The existing formation-vector error is below
$10^{-20}$ before the last entry rounding. Normalization and partial trace
give trace-norm error below $10^{-18}$ for $\rho$ and $2\times10^{-18}$ for
$\tau$. Each four-by-four kernel block has operator error below $10^{-19}$.
Rounding the source entries to $2^{-80}$ adds below $10^{-23}$ in trace norm.
With $\|K\|\le1$ and the Bell factor $1/4$, replacing both source copies and
both kernel factors changes a probability by less than $2\times10^{-18}$.
The archived decimal centers include a further rounding allowance, giving
radius $4\times10^{-18}$. Interval arithmetic carries the resulting pixel
boxes through both derivatives. A denominator interval containing zero is
reported as unresolved. [Integer wave calculation](certificate.py)

Two independent, seeded **synthetic** acquisitions retain the complete
five-pixel plus miss multinomial counts in each direction. The first
uses $2^{48}$ starts per emitted pixel and stock, four stocks and Bonferroni
binomial intervals. The held run uses $2^{60}$ starts per pixel and stock
in 4,096 independent batches, adding the ordinary calibration at the formed
source. For an observed hit count $X$ the held mean enclosure has half-width
$\sqrt{2XL}+2L$, with $L=\log(2\cdot2500/0.01)$. The standard two binomial
Bernstein tails and a union bound cover all 2,500 held probabilities with
at least 99% confidence. The code also verifies that every deterministic
native interval lies within its archived count interval. At the central
gateway the held forward and reverse area-equivalent intervals both lie
inside $[29.54,29.78]$ for vacancy and $[-1.99,-1.94]$ for excess. The
forward pencil determinant intervals lie within $[-1.277,-1.275]$ and
$[3.815,3.959]$. Some outer mounts resolve only in the integer calculation.
[First acquisition](results.json) · [First confidence intervals](bounds.json) · [Independent held acquisition and intervals](held.json)

From the repository root, with one numerical thread:

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=code:.
pkg=matched_einstein_limit.source_implication.full_geometry.station_balls.formed_optics
python -m "$pkg.measure"
python -m "$pkg.bounds"
python -m "$pkg.certificate"
python -m "$pkg.held"
python -m pytest -q --pyargs "$pkg"
```
