# CP-1 corrected-balance candidates

Companion to [G2_GAP_1_COUPLING_THEOREM.md](G2_GAP_1_COUPLING_THEOREM.md). Status:
diagnostic scan with declared menu; no ledger row moves. Generator and artifact:
`proof/nclosure/` at repo root (scan.py, nclosure_scan.json); measurement-side
sigma from `proof/planck_posterior/planck_lambda_to_N_propagation.json`.

## Setting

The certified arithmetic of the coupling theorem is exact under two independent
recomputations, so the 6.6 percent offset between the conditional bridge
capacity and the Λ-located capacity, if physical, sits in a premise. Working
hypothesis of this scan: the Λ-located capacity is correct and the balance
condition CP-1 carries a missing term. The required correction is exact:

```
Delta = 6*pi/(P_fwd*alpha_U) - ln(N_Lambda/pi) = 0.06417
```

(2.28×10⁻⁴ relative in the exponent), at the Planck 2018 TT,TE,EE+lowE+lensing
centrals with propagated sigma_lnN = 0.0266. Discrete-count alternatives are
excluded by magnitude: m_rep = 22 shifts the exponent by 8 percent, an
11-port read by 9 percent, the alpha_U mode choice by 6×10⁻⁴ relative in N.

## Scan

Menu of 40 one-term structural corrections; 14 survive at |z| < 1. The Λ readout
at Planck precision selects nothing inside this set; proximity carries no
weight at this menu density. The three candidates with structural readings:

| candidate | z | corrected closure |
|---|---|---|
| seed 15π/16 | 0.014 | N = (15π/16)·exp(6π/(P·α_U)); seed π(1 − 1/16), with 16 readable as β_EW² or as d² at spacetime dimension d = 4 (d² − 1 = 15 is also the small-ball denominator of the Einstein branch) |
| balance term (π/24)·α_U per port | 0.016 | X/12 = π/(2Pα_U) − (π/24)α_U; equivalently α_U·X = 6π/P − (π/2)α_U², equivalently a relative depth balance Γ_EW = t_tr·(1 − P·α_U²/12) |
| half-period per tick pair | 0.048 | X = 6π/(Pα_U) − π/48 |

Two-loop-shaped terms are excluded by the scan: α_U·ln(1/α_U) lands at z = 2.5
and α_U·ln(1/α_U)/(2π) at z = −1.6, so an RG-truncation reading of the offset
fails, and the surviving corrections are structural rather than
higher-loop.

Under the +BAO likelihood combination the target moves to Delta = 0.0795 and
the survivor set changes; the likelihood-combination freeze recorded in the
propagation artifact is prerequisite to any scoring.

## Discharge obligations

Each candidate faces the CP-1 gate in corrected form: a counting or geometric
derivation from declared screen structure, with the same three named splices
the tick-projection certificate scope note leaves open (the port/adjoint
identification, the pixel-area weight per port read, the port/channel
identification). A derivation landing one of these terms closes CL-7 in
corrected form and moves CL-3 from a 2.4σ conditional offset to a sub-0.1σ
conditional landing at declared menu price M = 40, still conditional on CP-2
and the F construction. A derivation landing the uncorrected balance keeps the
2.4σ offset and returns the burden to the Λ side. Either derivation decides;
nothing else in this scan does.
