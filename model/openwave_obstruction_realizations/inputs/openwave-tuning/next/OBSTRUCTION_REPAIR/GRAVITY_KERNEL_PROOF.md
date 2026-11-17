# The continuum comparison and the missing clock response

This is an independent mathematical comparison for the regulator worker.
It neither supplies the microscopic entropy kernel nor assumes that its
equilibrium calculation has passed. `gravity_kernel_review.py` checks the
identities symbolically.

In three spatial dimensions let g_ij=exp(2σ)δ_ij, linearized around flat
space. For a radius-r ball, the metric changes its coordinate-region area
and volume by 2∮σ dA and 3∫σ dV. At the round reference sphere the boundary
shape contribution satisfies δA_shape=(2/r)δV_shape. Therefore the
first-order fixed-proper-volume area response, independent of the detailed
first-order boundary displacement, is

    δA_V = 2∮σ dA − (6/r)∫σ dV.

For w_r(x)=(r²−|x|²)/(2r) inside the ball and zero outside, distributional
differentiation gives Δw_r=δ_boundary−(3/r)1_ball. Thus δA_V=2∫w_r Δσ.
For a Fourier component q, the vacuum-area entropy kernel is

    K_r(q) = −2η |q|² W_r(q),
    W_r(q) = 4π r⁴ [(3−z²)sin(z)−3z cos(z)]/z⁵,  z=|q|r.

W_r(0)=4πr⁴/15, whereas K_r(0)=0. Constant coframe changes cannot supply
the required curvature response. Finite-ball windows have zeros away from
q=0; dividing by a window near a zero needs its own error control.

If, on the same actual source, the matter modular kernel is independently
found to be 2πW_r(q)T00(q), stationary entropy would give
Δσ=−πT00/η. Only then does a localized energy E, with the stated decaying
boundary condition, give σ=E/(4ηr). The common window cancels across ball
radii. This explains why agreement between independently measured kernels
at different radii is a meaningful test. Substituting these continuum
formulas for the kernels to be measured would instead assume the result.

The curved-coframe check σ=a|x|² gives δA_V=16πar⁴/5. Its scalar curvature
is R=−24a at first order, so δA_V/A=−Rr²/30, as required. A fixed-coordinate
area calculation would give 8πar⁴ and would test a different variation.

The spatial result alone does not establish attractive gravity or slow
proper clocks. For a static metric

    ds²=−exp(2φ)dt²+exp(2σ)δ_ij dx^i dx^j,

Γ^i_00=exp(2φ−2σ)∂iφ. If the lapse remains one (φ=0), an initially resting
massive test body has zero initial geodesic acceleration, independently of
the spatial σ profile. Its proper clock is not redshifted. Deriving the
weak-field relation φ=−σ requires additional temporal/spatial field
equations or the relevant boosted entropy conditions on the actual
source. Installing that relation after a successful spatial solve would
leave the requested dynamical mechanism unproved.
