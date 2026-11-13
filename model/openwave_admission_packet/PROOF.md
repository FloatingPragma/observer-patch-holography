# Finite native field replays and their physical domain

All tests use the same exported 180-coordinate Nambu symbol

`H(k)=H0+sum_i Gamma_i sin(k_i)`.

The original native-carrier frame identifies each charged reducing block.
Its exact Clifford identities give
`H(k)^2=(m_j^2+sum_i sin(k_i)^2) I` on that block. Both energy signs and every
other coordinate remain in the complete action. A positive eigenvector
defines a heralded preparation of probability one half. The complementary
failure event is retained with probability one half.

The executable is the full original Taichi stencil on a periodic cube of
side six. A separate NumPy eigendecomposition of the full symbol restricted
by the proved frame supplies an eigenvector. The code verifies its residual
against the full 180-coordinate symbol. It then compares every coordinate
of the propagated spatial field to the analytic plane-wave phase. The
stability test uses event parameter 2; the rest and nonzero-momentum tests
use 3/4. All use 80 fourth-order Runge-Kutta steps. The original bound
`T^5 M^5/(120 n^4)` applies to exact-arithmetic integration on the stated
symbol-norm domain. Numerical roundoff is checked diagnostically and is not
covered by that truncation bound. An energy shift 0.03 fails the same
full-vector comparison by more than 100 times its integration bound.

The three commands support bounded parts of criteria 4, 11 and 30. Free
finite norm persistence is not nonlinear localized-particle stability.
Three charged block coefficients do not prove exactly three interacting
physical minima. The sine symbol is a finite lattice dispersion; it is not
an exactly Lorentz-invariant physical continuum. The source contains its
other Nambu, valley and massless coefficients. No incompatible construction
is used to fill a missing condition of these three tests.

The source export uses an analyzer phase selected with charged-lepton
target exposure and the electron calibration. The mass spectrum is therefore
a calibrated comparison for the admission audit. The unitless dispersion
and stability results are consequences of the declared finite action.
They are not compared with an unexposed observational target. The supplied
tests alone do not meet the genuine-prediction admission gate. An official
admission, independent blind certification and cell award are separate
maintainer actions; no such action is performed by this package.
