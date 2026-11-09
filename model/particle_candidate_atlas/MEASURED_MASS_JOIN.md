# Necessary matching budgets for the configured charged actions

The complete C-duty inverse and the generic positive-family full-CAR theorem
realize the specified capacity-dependent coefficient as a finite feedback
mean and a specified free Dirac action. This calculation propagates its four
certified source-coordinate boxes into dimensionless masses and ratios. It
also gives necessary sizes of corrections if that action is to reach the
quoted charged-lepton ratio windows. Here rho_match denotes matching error, distinct from the MCPR root amplitude
rho=sqrt(2) exp(-chi). No correction is inferred to occur in
the source, and no measured value is used to alter its parameters.

## Original action normalization

For each admitted source pair P,u, retain the owning response law

    chi=(u^2/512)/(1-u/77),
    delta=2/9+[(P-phi)^2/(21 pi)]/[1-u/27-u^2/135],
    d_j=[1+sqrt(2)exp(-chi)cos(delta+2pi j/3)]^2,
    kappa=(u/50)/(1+u/31+u^2/310),
    A^3=exp(3kappa)/[2*6^14*d_0*d_1*d_2],
    sigma_j(Y)=A*d_j=m_j^(0)/|v_H|.

The roots are in the owning positive chamber and are sorted without changing
their original amplitudes. A is the stipulated determinant normalization,
not an electron-mass calibration. Rational Taylor enclosures and integer
cube-root bounds preserve the full input box. The two complete-memory inverse
branches implement the same three eigenvalues with different admitted rates
and gain. No probability is divided out in their physical pairing controller.

The forward and measured-endpoint coordinates are supplied decimal inputs.
The two closure variants are certified roots of their declared census and
scale equations. Their conditional source construction is distinct from a
proof that primitive repair selects those choices. The law's target-informed
architectural ancestry is retained.

## Ratio distortion required by a reference window

Let r=(1,r_mu,r_tau) be the source ratios, and let [L_i,U_i] be specified
positive reference intervals, including L_0=U_0=1. Suppose componentwise
corrections multiply the three masses by positive factors f_i. A common
unit cancels, and the output ratios belong to the reference intervals only
if f_i/f_0 belongs to [L_i/r_i,U_i/r_i]. The smallest possible spread of
these three factors is exactly

    D_min(r)=max(1, max_i(L_i/r_i)/min_j(U_j/r_j)).             (1)

The lower bound follows from the largest lower endpoint and smallest upper
endpoint. It is attained by choosing points in each interval within that
minimal spread. It retains the shared electron denominator; separately
adjusting two ratios without that shared factor is unnecessary.

For source boxes r_i in [l_i,u_i], a rigorous necessary lower bound is

    D_required=max(1, max_i(L_i/u_i)/min_j(U_j/l_j)).           (2)

Replacing u_i by l_i in the numerator and l_j by u_j in the denominator
gives an upper enclosure for D_min throughout the box. Rectangular interval
independence can weaken the lower bound; it cannot create a false exclusion.

The finite mass-response theorem supplies the conditional distortion budget

    D_available = (1+rho_match)/(1-rho_match) * (1+eta)/(1-eta)
      * sqrt[(1+eps_L)(1+eps_R)/((1-eps_L)(1-eps_R))].         (3)

Thus D_available<D_required excludes that particular matching hypothesis.
Here eta bounds ||E||/sigma_min(M0), eps_L/R bound complete Lorentz-covariant
kinetic forms, and rho_match is a separate observable-matching bound. None is
computed physically by this calculation. In a matching-only model,
rho_match >= (D_required-1)/(D_required+1). For equal left/right kinetic bounds
and no other correction the same expression is necessary for eps_L=eps_R.
For a correction on only one kinetic form the necessary bound is
(D_required^2-1)/(D_required^2+1).

## A sharper additive self-energy condition

With canonical kinetic forms and no separate matching error, assume only
||E|| <= eta*m_e^(0). Weyl's singular-value inequality then permits
|m_i'-m_i^(0)| <= eta*m_e^(0), not eta*m_i^(0). Writing
t=m_e'/m_e^(0)>0, matching requires simultaneous interval intersections

    (r_i-eta)/U_i <= t <= (r_i+eta)/L_i.

The exact minimal normalized additive budget in this unrestricted constant
matrix class is

    eta_min(r)=max(0, max_(i,j)
                    [L_j*r_i-U_i*r_j]/[L_j+U_i]).             (4)

Necessity follows by comparing every lower and upper endpoint. Conversely,
their nonempty intersection supplies t, and diagonal corrections in the
source singular basis attain the allowed shifted masses. This sharpness
statement concerns an allowed matrix perturbation, not a calculated source
self-energy or a proposed phenomenological fit. For source intervals, use
l_i and u_j in the numerator for a necessary lower bound, with the opposite
endpoints for an upper enclosure. The regime eta<1 preserves the lightest
positive coefficient.

## Meaning of comparison and symmetry

CODATA2022 ratio windows are quoted standard-uncertainty intervals, not
theory errors or joint confidence regions. Exact agreement with estimated
central values is not required. Zero lower correction budget means only that
the action ratio boxes can meet those windows. It establishes no interacting
particle identification. A nonzero lower budget states the minimal size of
a hypothesized correction within the stated matrix class; it gives no
evidence that the source produces that correction.

A common mass rescaling preserves all ratios and Koide Q, so it cannot
repair a disjoint ratio. Family-dependent kinetic and Higgs-mass coefficients
are compatible with the declared gauge representation: the family factor
is gauge trivial, a left kinetic form acts on the entire lepton doublet,
and the charged mass tensor transforms with the Higgs.
The owning mass-response theorem gives exact examples in which such allowed
corrections change Q. Symmetry permission does not calculate their magnitude
or select a correction preserving Q=2/3. For a deformed source with
Q<2/3, an additional demand for exact balanced Koide requires at least the
corresponding Q distortion, even if its measured-ratio windows overlap.

Electron matching concerns a charged infrared threshold under a specified
state and dressing. Muon and tau matching concerns a resonance or inclusive
line shape. No three isolated real mass atoms, or common unrestricted
electric-sector gap, are inferred. The physical source clock, interacting
state, kinetic normalization, self-energy, mirror treatment and observable
definition remain explicit matching inputs.

## Owning evidence

- `model/charged_lepton_response/p_dependence.py` and its authenticated four
  source branches; `model/mcpr_native_memory_bridge/PROOF.md` and the two
  complete inverse branches per capacity input.
- `model/pixel_mass_closure/closure.py:mass_coefficients`, the complete
  Higgs-unit normalization, with conditional root receipts.
- `model/particle_candidate_atlas/CHARGED_FIELD_ADMISSION.md`, the general
  positive-family finite full-CAR and specified free-field attachment.
- `model/charged_lepton_response/MASS_RESPONSE_BRIDGE.md`, finite whitening,
  additive perturbation, charged threshold and resonance distinctions.
- CODATA2022 primary constants: https://physics.nist.gov/cuu/pdf/all.pdf.
