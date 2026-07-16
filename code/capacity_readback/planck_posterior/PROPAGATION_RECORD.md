# Planck posterior propagated through the Λ→N map

Executes the pending item named in GAP-A5 (CL-3/CL-4): "the joint cosmological
posterior propagated through the same map." Artifact:
[planck_lambda_to_N_propagation.json](planck_lambda_to_N_propagation.json),
generator [propagate.py](propagate.py).

## Map and reduction

Λ = 3Ω_Λ H₀²/c², N_Λ = 3π/(Λ ℓ_P²) = πc⁵/(Ω_Λ H₀² ħG). The measurement side
reduces to one sufficient statistic, Ω_Λh²: ln N_Λ = const − ln(Ω_Λh²) − ln G.
The bridge side N_CRC = 3.5321315434×10¹²² is certified at relative width
1.6×10⁻²⁵ and enters as a point. c and ħ are exact by SI; G contributes
2.2×10⁻⁵ relative, negligible but carried in quadrature.

## Uncertainty model

Gaussian in ln(Ω_Λh²), built from published Planck 2018 VI base-ΛCDM marginals
by three routes with independent degeneracy assumptions: (1) (Ω_Λ, H₀) pair
with correlation ρ ∈ [0.90, 1.00]; (2) (H₀, Ω_mh²) pair with ρ ∈ [−1.00, −0.90];
(3) the Ω_mh³ degeneracy eliminating Ω_mh². The full PLA chains are not in the
repo; the spread across routes and sweeps bounds the model error at ±0.05 in z.

## Results

| combination | N_Λ | gap | z |
|---|---|---|---|
| TT,TE,EE+lowE+lensing (the combination the corpus consumes: N_Λ = 3.313×10¹²² back-solves to its centrals) | 3.3126×10¹²² | 6.63% | **2.37–2.48** |
| +BAO (sensitivity case) | 3.2633×10¹²² | 8.24% | 3.82–3.95 |

1. The circulating "about 2.5 one-dimensional Planck standard deviations"
   survives propagation: the correct sufficient-statistic treatment gives
   z ≈ 2.4–2.5 for every plausible correlation. The prose figure's implicit
   fully-correlated linear sum happened to approximate the right answer
   because Ω_Λ and H₀ are in fact near-degenerate in the Planck posterior;
   that assumption is now stated rather than silent.
2. The likelihood-combination choice is material: +BAO moves the comparison to
   ≈3.9σ. Per the program's own anti-fitting rules the combination must be
   frozen before CL-3 becomes evaluable; the defensible default is the one the
   ledger already consumed (TT,TE,EE+lowE+lensing).
3. Everything here is conditional on CP-1–CP-3 and the readback map F (CL-7).
   This artifact upgrades the *measurement-side* uncertainty of the comparison
   from prose to computation; it does not promote the row.

## Remaining refinement

A definitive version samples the PLA `base_plikHM_TTTEEE_lowl_lowE_lensing`
chains directly through the map (removes the Gaussian and degeneracy
approximations; expected shift ≲0.05 in z).
