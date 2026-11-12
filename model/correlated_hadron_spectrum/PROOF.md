# Correlated retained hadron fits

The input is the retained 32-stream, 128-configuration quenched Wilson
table on the 12³×24 source at beta 5.7. No gauge configurations are
regenerated. All four kappas, two sink choices and the original primary
and alternate time windows are retained: 48 fits, each with 32 outer
delete-one-stream refits. Pion and rho fits use one periodic cosh;
nucleon fits use one exponential. These are finite fit models, not
complete spectral reconstructions.

## Paired covariance and fixed rules

Let j_i be the retained delete-one-stream correlator means and jbar their
average. Reconstruct compatible stream summaries x_i=n jbar-(n-1)j_i.
Then their mean is jbar and deleting stream i gives j_i exactly in real
arithmetic. The separately rounded recorded central mean is used in the
central fit; its discrepancy from jbar is retained. The largest-coordinate
relative discrepancy must be at most 1e-11. A uniform shift of all x_i
does not change any covariance.

For m retained stream summaries the covariance of their mean is

    C = sum_i (x_i-xbar)(x_i-xbar)^T / (m(m-1)).

For the full n-stream sample this equals the original jackknife expression
(n-1)/n sum_i (j_i-jbar)(j_i-jbar)^T. For every outer deletion, compute
this matrix on the remaining 31 stream summaries and refit that replica.
All paired labels survive through hadron masses and Wilson-potential scales.
This reconstructs the information present in the reduced stream means;
it supplies no missing within-stream autocorrelation or raw gauge data.

The fixed rule standardizes C by its diagonal into a correlation matrix R.
There must be at least p+2 streams for p time coordinates, strictly positive
diagonal entries, full numerical rank at relative eigenvalue tolerance
1e-12, and condition number at most 1e8. Failed windows are rejected without
a pseudoinverse, shrinkage, diagonal fallback or omitted replica. The rule
is stated in code before the retained fits are evaluated and contains no
experimental hadron mass. These thresholds define numerical admission,
not a theorem controlling finite-sample covariance estimation error.

For model vector f(m), the conditional best amplitude is

    A(m) = f(m)^T C^-1 y / (f(m)^T C^-1 f(m)).
    chi²(m) = (y-A(m)f(m))^T C^-1 (y-A(m)f(m)).

The positive metric proves this amplitude is the unique least-squares
minimum at fixed m. The mass search keeps the original interval [0.02,3],
uses 150 grid points and refines every grid-bracketed local minimum, with
both endpoints retained. This is a deterministic bounded numerical
optimization, not a proof of a global continuous minimum. Boundary and
nonpositive-amplitude replicas are explicitly listed. No chi-square
confidence or goodness-of-fit distribution is asserted for the estimated
covariance. Original alternate windows and sinks give fit sensitivity;
they are not selected by agreement with experimental masses.

## Common mass and scale chain

The original pion-positive, relative-error-at-most-25% admission determines
the chiral fit points. Every original kappa fit remains in the receipt,
including the inherited unresolved lightest nucleon. The original chiral
forms and their diagonal mass-error weighting are retained as a distinct
approximation. They do not become fully correlated errors-in-variables
regressions merely because the first-stage time fits are correlated.

Every outer replica refits the mass chain with its matching Wilson-potential
scale. The light-point condition uses the imported ratio c=0.1780. With
the supplied rho anchor R=0.77526 GeV,

    m_pi=c R,  m_rho=R,  m_N=R (a m_N)/(a m_rho).

With the alternative conventional string anchor S=0.44 GeV,

    m_rho=S (a m_rho)/sqrt(a² sigma),
    m_pi=c m_rho,  m_N=S (a m_N)/sqrt(a² sigma).

The receipt retains the joint six-output replica covariance for both
calibrations. Zero conditional pion/rho variance under the rho anchor is
algebraic shared-scale cancellation, not zero physical uncertainty.
Crossing every admitted hadron variant with the original potential variants
gives a finite fit-sensitivity envelope. The unchanged diagonal-time chain
is independently recomputed and compared with its frozen receipt.

No result bounds sea-quark, volume, continuum, calibration or physical
source/clock-transfer error. Covariance-aware inference changes the
finite-ensemble fit; it is not an independent mass prediction. Synthetic
correlated cosh/exponential data recover known coefficients, while singular
and excessively conditioned matrices are rejected by the same rule.
