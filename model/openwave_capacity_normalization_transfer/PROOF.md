# Cosmological normalization and a validated long-wavelength metric response

The domain is the declared flat Einstein model with conserved ideal dust,
ideal radiation without anisotropic stress and the reconstructed canonical
capacity scalar. Its present matter and radiation fractions, Hubble value,
event-horizon choice and scalar attachment are supplied. The finite observer
capacity has no proved transfer to this horizon, scalar stress or physical
clock. Architectural choices are permitted; uniqueness of the construction
is not required for the following calculation.

## Two normalization freedoms

Write H(a)=H0 E(a) and use the reduced Planck mass M_Pl. For every positive
lambda the transformation

    H0 -> lambda H0, rho_i0 -> lambda squared rho_i0,
    L_horizon -> L_horizon/lambda, k -> lambda k

preserves the dimensionless background, the canonical scalar divided by
M_Pl, the scalar equations normalized by H squared, k/H0 and the complete
dimensionless perturbation equations. The future-radius integral
L=a integral_a^infinity da'/(a' squared H) transforms by 1/lambda as well.
The event-area law rho_dark=3 M_Pl squared/L squared is compatible with this
entire family. It therefore does not determine a member's absolute H0 or
present density 3 M_Pl squared H0 squared Omega_dark0.

Independently, the linear perturbation system is homogeneous: multiplying
every perturbation by A, including A=0, gives another solution. A transfer
coefficient normalized to a primordial potential is independent of A.
Neither that coefficient nor the dimensionless background selects a
primordial fluctuation amplitude or its statistics. A native area/clock
relation and a native state-selection/preparation law could constrain these
freedoms; such relations are not supplied by the declared background.
Choosing H0 or a density from observation is an empirical calibration.

## Exact reduction of the regular adiabatic mode

Let N=log(a), h=H_N/H, s=phi_N/M_Pl, b=s_N,
A=V_phi/(H squared M_Pl) and M=V_phiphi/H squared. The background identities
are b=-(3+h)s-A and h=-(3 Omega_m+4 Omega_r+s squared)/2. For constant R set

    T=R(H/a) integral_0^a da'/H(a'),  Phi=R-T,
    U_m=U_r=T, delta_m=-3T, delta_r=-4T,
    y=sT, z=bT+s(Phi+hT).

Then T_N=Phi+hT and Phi_N=-T_N. Substitution in all seven dust/radiation/
scalar evolution equations and the Hamiltonian constraint gives exactly
zero residual. The executable checks every identity symbolically, including
the momentum constraint. In the scalar second derivative, use
b_N=-h_N s-(3+h)b-Ms+2hA; the h_N terms cancel.
This also reuses the exact constraint-propagation theorem of the parent
gravitational attachment.

R=3/2 fixes the asymptotic radiation-era potential to one. For an exact
power law H=C a^(-p), T=R/(p+1) and Phi=Rp/(p+1). Radiation gives T=1/2,
Phi=1; dust gives T=3/5, Phi=9/10 under this same curvature normalization.
The variables U_i denote the regular limit of aH theta_i/k squared; no
division by zero is executed. At strict k=0 the comoving matter/radiation
densities and scalar displacement vanish. The solution is the residual
adiabatic coordinate mode. Its regular long-wavelength metric coefficient
is not a measured density-growth mode or a CMB observable.

## A rational quadrature certificate

Set t=-log(a), x=sqrt(Omega_dark) and r=Omega_r/(Omega_m+Omega_r).
The existing validated background supplies outward intervals for

    x_N=x(1-x squared)(1+r+2x)/2, r_N=-r(1-r),
    h=-(3+r)/2+(1+r)x squared/2+x cubed.

On the unit square, -2<=h<=0. Both terms in

    h_N=[(1+r)x+3x squared]x_N+(1-x squared)r(1-r)/2

are nonnegative. The exact inequality x(1-x squared)<=2/5 follows from
its maximum 2/(3 sqrt(3)) and 4/27<4/25. Hence x_N<=4/5,
h_x<=5 and 0<=h_N<=5(4/5)+(1/2)(1/4)=33/8.

For F(t)=exp(-t)/E(t), differentiation gives

    F_t=(h-1)F, F_tt=[(h-1) squared-h_N]F.

Thus F is positive and decreasing, and F_tt/F lies in [-25/8,9]; in
particular |F_tt|<=9F. For the same-density constant-Lambda control,
h_N is one half the variance of dilution exponents 3,4,0 with weights
Omega_m, Omega_r, Omega_Lambda. These weights are nonnegative and sum to
one, so 0<=h_N<=2 by the range-four variance bound. The same |F_tt|<=9F
certificate applies.

On a step of length d, the trapezoidal error is therefore at most
3 d cubed F_left/4. The code uses an outward upper bound for F_left,
retains every endpoint interval and sums all step errors. The background
logit uses the original monotone second-order Taylor map and its global
third-derivative remainder; no floating background solution supplies these
bounds. Exponentials, square roots and logs use the existing rational
enclosures. Every interval operation rounds outward to thirty decimal
places, and the background endpoint also retains its binary rounding.

The integration stops at t=12. Since E(a)>=sqrt(Omega_r0)/a squared,
the omitted improper integral satisfies

    0 <= integral_0^(a_min) da/E(a)
      <= a_min cubed/[3 sqrt(Omega_r0)].

This bound is below 8.1 times 10^-15 for the supplied inputs. At a requested
redshift, the rational lower enclosure of log(1+z) is inserted into the
grid. The log-interval width times the left F upper bound encloses the
remaining integral endpoint shift. The logit interval is enlarged by four
times that width, using its global derivative bound. H at the requested
redshift uses the exact rational a=1/(1+z) and this enlarged density interval.
These operations certify the final metric coefficient, rather than only a
quadrature on an approximate background.

At step 1/128 the present normalized capacity potential is enclosed in
[0.734179279095,0.734257299186]; the same-density Lambda control is enclosed
in [0.709108311479,0.709181069590]. Their ratio lies in
[1.035249403259,1.035465650732]. The corresponding finite-k matter-transfer
ratios are different quantities. Both their physical roles and their
normalizations are kept separate.

An independent step-1/64 calculation overlaps the refined intervals, whose
widths are smaller. Four inherited full seven-equation integrations at
k/H0=10^-4 and 2 times 10^-4 give present potentials inside the validated
zero-k intervals. Their comoving matter amplitudes scale approximately as
k squared and their constraints are checked without projection. Those ODE
runs are floating diagnostics; they supply no rigorous nonzero-k error and
do not replace the interval proof.

The result is conditional on classical Einstein gravity, the supplied
background fractions and canonical scalar. Native horizon/capacity and
stress acquisition, a physical amplitude, primordial statistics, neutrino
free streaming, recombination and the Boltzmann transfer are outside this
calculation. No independent cosmological normalization or match with Nature
is asserted.
