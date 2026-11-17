# Orthonormal continuum cells and an actual moving CAR algebra

For a cube of coordinate side h centred at x_a, let f_a(x)=h^(-3/2) 1_cube(x). Distinct tiled cubes are orthonormal in the continuum one-particle Hilbert space; tensor each with all four spinor components. Their span is an exact finite CAR subalgebra. The source is the continuous proper-cutoff Dirac projector P(p), with a specified periodic or antiperiodic spin structure on a continuum torus. The numerical grid enumerates observer cells and Fourier quadrature; it is not a postulate that the metric changes in discrete steps.

The Fourier form factor is w_h(p)=product_j sinc(p_j h/2). With the local beta product P_beta=(I-beta)/2 separated exactly, the covariance is

C_ab = P_beta delta_ab + (h^3/V_box) sum_p w_h(p)^2 [P(p)-P_beta] exp(i p dot (x_a-x_b)).

The sum is over every continuum torus momentum, not only a first Brillouin zone. Regrouping those momenta into a base grid plus aliases is an identity. The implementation sums successive alias cubes for the decaying difference; the exact P_beta contribution is never truncated. Alias-shell convergence is measured. This is still an approximation to the exact continuum covariance, and the finite cell span is not the full regional CAR algebra. Subcell observables remain in the source but are outside this declared observer algebra.

For the moving observer use actual normalized scaled functions

f_(a,epsilon)(x) = exp(3 epsilon/2) f_a(exp(epsilon) x).

Both cube centres and widths shrink by exp(-epsilon). This is induced by a unitary coordinate pullback on the one-particle space. A vector field agreeing with the dilation on all observed cubes can be smoothly tapered outside a larger neighbourhood, preserving the same observed first-order CAR map. No random selector or entropy-dependent radial adjustment is involved.

Writing d=x_a-x_b and C_reg=C-P_beta delta, the derivative at epsilon=0 is

delta C_shape(d) = integral h^3 exp(i p dot d) [
    (-3 w_h^2 - p dot grad_p w_h^2)(P-P_beta)
    - i(p dot d) w_h^2(P-P_beta)] d^3p/(2pi)^3.

The torus uses the corresponding complete Fourier sum. This is computed independently from the coframe response. The literal finite scaled-cell covariance is also evaluated, and its finite entropy derivative agrees with the analytic derivative.

For a coframe mode sigma(x), proper volume changes by 3 integral_union sigma dx. An inward dilation with coefficient average_union sigma changes the coordinate volume by the opposite amount. For a Fourier mode the average includes the exact cell sinc factor, as well as the sum over centres. Thus proper volume is fixed to first order for the actual union of cubes. This union approximates a geometric ball only with its explicitly stated shape discretization error; no equality to a perfect sphere is assumed.

In infinite flat space, simultaneous uniform metric g=e^(2epsilon)delta and the scaled observer cells are related by an exact coordinate isometry. Changing integration variable p to exp(-epsilon)p makes the joint covariance equal to the original covariance, at every cell size. The cancellation is a coordinate identity, not the absence of nonuniform curvature response. On a fixed-period torus, uniform physical rescaling changes the entire box and is not an isometry, so a finite-box response remains. That response is measured rather than subtracted.

The massless antiperiodic box avoids the isolated zero-mode filling ambiguity while preserving cubic symmetry. A separate periodic small-mass control is retained. These are explicit continuum spin structures, not species mass calibrations.
