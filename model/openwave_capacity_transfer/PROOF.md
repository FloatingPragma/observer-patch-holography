# Linear gravitational attachment of the canonical capacity scalar

## Declared model and conventions

Take the retained flat background, its supplied present matter/radiation/dark
fractions, and its canonical scalar reconstruction. Couple conserved ideal
dust and ideal radiation (zero anisotropic stress) to classical Einstein
gravity. In conformal Newtonian gauge the two metric potentials agree and are
denoted Phi. This is a declared effective theory. It is not a derivation of
quantum spacetime, the scalar from native capacity records, or a primordial
state. The standard perturbation equations use the conventions of
[Ma and Bertschinger](https://arxiv.org/abs/astro-ph/9506072); the separation of
regular growing perturbations from freely supplied initial modes is also
discussed by [Doran et al.](https://arxiv.org/abs/astro-ph/0304212).

Write N=log(a), m=Omega_m, r=Omega_r, s=phi_N/M_Pl, b=s_N,
h=H_N/H=-(3m+4r+s^2)/2, A=V_phi/(H^2 M_Pl), M=V_phiphi/H^2,
Q=k^2/(aH)^2, U_i=aH theta_i/k^2, y=delta_phi/M_Pl, z=y_N.
The background satisfies b+(3+h)s+A=0. The inherited exact scalar theorem
gives 2 <= M <= 6 and rest-frame propagation speed squared one.

## Full scalar Einstein/fluid system

With P=Phi_N, the momentum constraint is

    P = -Phi + (3m U_m + 4r U_r + s y)/2.

Energy-momentum conservation and the perturbed Klein–Gordon equation give

    delta_m,N = -Q U_m + 3P,       U_m,N = h U_m + Phi,
    delta_r,N = -(4/3)Q U_r + 4P,  U_r,N = (1+h)U_r + delta_r/4 + Phi,
    y_N = z,
    z_N = -(3+h)z - (Q+M)y + 4sP - 2A Phi.

Define D=3m delta_m+3r delta_r+s z-s^2 Phi+A y, and the 00 constraint

    C = Q Phi + 3(P+Phi) + D/2.

Use m_N=(-3-2h)m, r_N=(-4-2h)r, Q_N=-2(1+h)Q,
A_N=M s-2h A, and the background scalar equation. Direct differentiation and
polynomial cancellation, checked exactly with SymPy, prove

    C_N = -(3+2h) C,
    Phi_NN+(4+h)Phi_N+(3+2h)Phi
      = (r delta_r+s z-s^2 Phi-A y)/2.

Thus the spatial trace follows, and C=0 initial data satisfy the remaining
Einstein equation for their interval of existence. The exact propagation is
`C(a)=C(a0)*(a0/a)^3*H(a0)^2/H(a)^2`. No numerical constraint projection is
performed. A deliberately nonzero initial C tests this propagation law.

## Stable variables and regular initialization

The numerical state uses Delta_m=delta_m+3U_m, Delta_r=delta_r+4U_r and
dU=U_r-U_m, avoiding cancellation of order-one superhorizon terms. Substitution
in the equations gives

    Delta_m,N = -Q U_m + 6r dU + (3/2)s(y-s U_m),
    Delta_r,N = Delta_r -(4/3)Q U_r - 6m dU + 2s(y-s U_r),
    dU_N = h dU + Delta_r/4.

The regular k=0 growing adiabatic mode, normalized to unit radiation-era
potential, is determined by

    T = (3/2)*(H/a)*integral_0^a da'/H(a'),
    Phi = 3/2-T, U_m=U_r=T, y=sT,
    z=bT+s(Phi+hT), delta_m=-3T, delta_r=-4T.

Indeed T_N=Phi+hT and Phi_N=-T_N satisfy the momentum constraint; the
background identities give C=0. A supplied finite-k initial correction
`delta_r -> delta_r-2Q Phi/(3r)` enforces the 00 constraint exactly. It is a
controlled early-time preparation, not a derived primordial spectrum. Moving
a0 from 1e-7 to 1e-8 tests the dependence on this preparation. The background
integration stops at a=1e-14; H(a)>=sqrt(Omega_r0)/a^2 bounds the omitted
dimensionless initial integral by less than 1e-14 for the runs supplied.

## Numerical evidence and limits

DOP853 evolves k/H0 in {0.1,1,10,100}, with independent tolerance and start-time
controls and 601 sampled constraint evaluations per trajectory. The retained
results report all tolerances, solver counts, residuals and transfer values.
An independent pure-dust Einstein–de Sitter control recovers Phi=1,
U_m=2/3 and Delta_m=-2k^2 a/3, including growth rate one.

The present matter transfer ratios to the same-density Lambda control are
approximately 0.97585, 0.97001, 0.91494 and 0.90822. Transfer normalization
means a linear response coefficient; it does not posit a physical perturbation
of unit amplitude. The numerical evidence is floating point, with no rigorous
interval ODE certificate. It supplies neither a CMB/large-scale-structure
likelihood nor a new empirical match. Free-streaming neutrinos, recombination,
anisotropic stress and the full Boltzmann hierarchy remain required for such
a comparison. The exact constraint identities close the linear gravitational
attachment of this declared scalar/fluid model only.
