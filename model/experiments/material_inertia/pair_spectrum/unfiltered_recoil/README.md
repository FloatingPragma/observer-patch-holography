# Long-flight recoil with every final port retained

A common impulse gives a positive late centre displacement of the freely moving excess–vacancy pair on the 20-cell Interval Echo board. On launch owners 0 and 5 the signed displacement lies in [0.01666, 0.01720] in the normalized owner coordinate. On owner 10 it lies in [0.01959, 0.02016]. Each interval includes the largest contribution the unresolved final states and their interference can make. These are numerical enclosures of the **mean of the complete terminal address distribution**. The older short-flight experiment also supplies actual simulated count banks. [Protocol](PROTOCOL.md) · [Full numerical receipt](results.json) · [Prior counted measurements](../impulse_response/README.md)

An observer prepares the sharp pair and selects the original lowest doublet. In the three runs, about 0.0122%, 0.0262%, and 0.307% of preparations pass. Each accepted launch receives one of two opposite address phases, flies freely, and ends with a read of both port addresses. The observer assigns each address its original graph-harmonic weight and averages the two. Subtract the two impulse-arm means to obtain the signed displacement. All terminal addresses contribute, including pairs on separate owners. The per-original-launch difference with a zero score for initial rejection is the quoted conditional displacement times the reported initial admission probability. The clock unit, harmonic, initial spectral filter, and impulse are supplied instruments. [Source and graph coordinate](../impulse_response/measure.py)

At impulse magnitude 1/4 and dwell 32,768 the complete-address mobility is [1.0170, 1.0498] × 10⁻⁶ on owners 0 and 5 and [1.1962, 1.2301] × 10⁻⁶ on owner 10, in harmonic-coordinate units per impulse and native exchange time. An independent calibration setting uses half the impulse and half the dwell. Four times its signed displacement predicts the held displacement to within the chosen 10% tolerance on every launch. The ratios lie in [0.9276, 1.0643] on owners 0 and 5 and [0.9367, 1.0558] on owner 10. A separate changed law halves all free propagation rates after launch while the preparation, impulse, clock and detector stay fixed. Every changed-law mobility interval lies strictly below its native partner. Complete propagation of both 57,360-component vectors at the original short dwell of 1/16 gives a negative signed displacement on each launch. The sign depends on the observation interval.

The signed-load impulse at the long dwell gives mobility intervals containing zero, with absolute upper bounds below 2.49 × 10⁻⁸. This limits the unfiltered response to a supplied material phase gradient. The charged-probe contact experiment measures net signed-load cancellation plus common scattering. The electrical interpretation requires an independently sourced and calibrated field. Each launch has a reproducible finite response in its declared spatial coordinate. The distinct mobilities across launch settings, the earlier failure of a scalar dispersion forecast, and the need for a calibrated distance delimit the physical mass inference. OW-19's complete local assessment remains partial. [Charged probe and prior dispersion controls](../README.md) · [Official row](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L121)

## Bounding the terminal reader

Let P be the exact projector onto the 24-line moving band and Q its complement. Let y be either signed post-impulse launch and U the full evolution. Write v=UPy and w=UQy. For the address reading X, with norm at most one,

$$
|\langle Uy,XUy\rangle-\langle v,Xv\rangle|
\leq 2\|QXv\|\,\|Qy\|+\|Qy\|^2 .
$$

The first term accounts for interference. The bound uses QXv, which is small on these localized launches, and requires no assumption on phases or propagation within Q. All signs and both measurement settings get separate bounds. The half-rate calculation uses the same spectral projectors.

The calculation reruns the parent Schur counts and the 24-vector Ritz calculation. The projected evolution error in norm is at most the projector error plus dwell times the Ritz residual by Duhamel's formula. This also bounds the error in QXv once the independent projector error is included. Operator norms of the two-column matrices bound both launch vectors and every superposition in the rest doublet. Replacing a unit exact rest vector by its normalized projection into the approximate rest doublet costs at most twice the rest projector error in norm, or eight times this error in the signed difference of the two means. Its coefficient direction differs from the numerical sharp-port preparation by at most twice the rest projector error divided by the square root of the numerical admission. On the two-dimensional numerical subspace the constant part of the signed-response matrix cancels, so the final error is bounded by its spectral width times that coefficient difference. All spectral bounds use floating-point residuals with an added 10⁻¹² allowance and the parent spectral separation estimates. The result is a numerical certificate under these estimates. These intervals have no sampling component. The comparison at short dwell evaluates each final address probability explicitly as a check. [Bound and complete short-flight calculation](measure.py) · [Independent exact small-system test for interference](test_measure.py)

From the repository root, replay the certificate and run its interference check with

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python -B -W error -m model.experiments.material_inertia.pair_spectrum.unfiltered_recoil.measure
.venv/bin/python -B -W error -m pytest -q model/experiments/material_inertia/pair_spectrum/unfiltered_recoil
```
