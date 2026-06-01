# OPH/BD Open Computation Audit Digest

Date: 2026-06-01

Source status: transcribed from a user-provided Pro handoff. The referenced
`sandbox:/mnt/data/...` workbook, report, result certificate, and README are
not mounted in this local workspace, so this digest preserves the reported
result and adds a local reproduction script for the low-energy target equations.

Local reproduction script:

```bash
python3 code/particles/uv/oph_bd_open_computation.py \
  --out code/particles/runs/uv/oph_bd_open_computation_results.json
```

## Selected Representative

The handoff fixes the selected representative as:

```text
BD_{n=1}^{OPH}
= Bouchard-Donagi E8 x E8 heterotic SU(5) Standard Model,
  one-Higgs stratum.
```

Reported selection rationale:

- Calabi-Yau threefold with Z2 fundamental group.
- Invariant SU(5) bundle.
- Observable SU(3) x SU(2) x U(1).
- Three generations.
- No exotic matter.
- A region with precisely the MSSM plus one Higgs pair.

Claim boundary: this is an encoded candidate-sieve result and low-energy
threshold target reduction. It is not yet a raw cohomCalg/Sage/Macaulay2
recomputation of the full Bouchard-Donagi geometry and bundle data.

## OPH Electroweak Inputs

```text
alpha_2(mZ) = 0.03377843630219015
alpha_Y(mZ) = 0.010131601067241624
v           = 246.76711732749683 GeV
m_H         = 125.1995304097179 GeV
m_t         = 172.35235532883115 GeV
alpha_3(ref)= 0.1184
M_SUSY      = 1 TeV in the unification proxy
mZ(ref)     = 91.187978085123 GeV, inferred from the reported proxy values
```

## Reproduced Low-Energy Targets

Using the local reproduction script:

```text
y_t_OPH                 = 0.987745211164
lambda_H_OPH            = 0.128706603202
lambda_MSSM_tree_max    = 0.0689737254093
Delta_lambda_min        = 0.0597328777922
m_h_tree_max            = 91.6524602856 GeV
```

Interpretation:

```text
BD_{n=1}^{OPH} needs a moderate positive stop/Higgs threshold,
not an exotic threshold.
```

## Stop-Threshold Proxy

The script uses:

```text
Delta m_h^2 ~= 3 m_t^4 / (2 pi^2 v^2)
  [ log(M_S^2/m_t^2)
    + (X_t^2/M_S^2)(1 - X_t^2/(12 M_S^2)) ].
```

Representative reproduced values:

| tan beta | X_t / M_S | Required M_S, one-loop proxy |
| ---: | ---: | ---: |
| 2 | 0 | 3046.266 GeV |
| 2 | sqrt(6) | 679.714 GeV |
| 5 | 0 | 1191.830 GeV |
| 5 | sqrt(6) | 265.933 GeV |
| 10 | 0 | 968.658 GeV |
| 10 | sqrt(6) | 216.137 GeV |
| 50 | 0 | 901.608 GeV |
| 50 | sqrt(6) | 201.176 GeV |

These are not spectrum-generator outputs. They are target bands for the
eventual BD soft-term and SUSY-breaking calculation.

## Gauge-Unification Proxy

With SU(5)-normalized hypercharge,

```text
alpha_1 = (5/3) alpha_Y
```

and `M_SUSY = 1 TeV`, the reproduced one-loop SM/MSSM proxy gives:

```text
log10(M_U/GeV)                 = 16.0815812332
M_U                            = 1.20665e16 GeV
alpha_U^{-1}                  = 26.0176807530
alpha_3_pred(mZ)              = 0.111511310319
Delta(alpha_3^{-1})_needed    = -0.521754255407
```

This is the string/GUT threshold target in the 1 TeV proxy:

```text
Delta_string/GUT(alpha_3^{-1}) ~= -0.52.
```

## Operator Audit

The local script checks hypercharge neutrality for:

```text
Required Yukawa terms:
  Q Hu uc
  Q Hd dc
  L Hd ec

Dangerous MSSM-charge-allowed operators:
  L L ec
  L Q dc
  uc dc dc
  Q Q Q L
  uc uc dc ec
```

The latter are gauge-allowed by the MSSM charge algebra and therefore require
a string selection rule, discrete symmetry, or geometric vanishing. The handoff
reports the intended paper-level statement as:

```text
OPH selects the branch; the BD/BCD geometry supplies vanishing RPV trilinears.
```

This statement still needs source-level citation verification before promotion
into the compact paper.

## Encoded Candidate Sieve

Reported encoded audit set:

| Candidate | Score | Verdict |
| --- | ---: | --- |
| BD_{n=0}^{SU(5),Z2} | 7/8 | rejected: no Higgs pair |
| BD_{n=1}^{SU(5),Z2} | 8/8 | selected |
| BD_{n=2}^{SU(5),Z2} | 7/8 | rejected/nonminimal: extra Higgs pair |
| Braun-He-Ovrut-Pantev SU(4), Z3 x Z3 | 4/8 | backup witness only |
| generic Spin(32)/Z2 heterotic | 0/8 | rejected as minimal OPH class |

Digest statement:

```text
BD_{n=1}^{SU(5),Z2} is the unique selected candidate inside the encoded audit set.
```

## Remaining Gates

The handoff reduces the open problem to exact target equations, not to a solved
BD moduli point. Remaining gates:

| Gate | Needed input | Status |
| --- | --- | --- |
| Cohomology audit | explicit X, V, equivariant action | structural witness known; raw compute not run |
| Yukawa/superpotential audit | sheaf cup products, harmonic representatives, instanton data | qualitative BCD witness reported; numeric textures not run |
| Gauge thresholds | heavy spectrum, Kahler/complex moduli, string scale | target deltas computed; string thresholds not run |
| SUSY breaking/decoupling | hidden sector, mediation, soft terms | low-energy target computed; UV soft terms not run |
| Moduli isolation | full map F_{BD,n=1}(m) | equations defined; full map unavailable here |

Final target system:

```text
F_{BD,n=1}(m_star) = O_OPH
```

with low-energy target components:

```text
y_t(m_star)                              = 0.987745211164
lambda_H(m_star)                         = 0.128706603202
Delta lambda_SUSY(m_star)                = 0.0597328777922
Delta(alpha_3^{-1})_string/GUT(m_star)   ~= -0.521754255407
```

## Immediate Next Steps

1. Import the raw workbook, JSON certificate, and original script once they are
   available outside `sandbox:/mnt/data`.
2. Verify the BD and BCD source claims against the original papers before any
   public wording change.
3. Add explicit raw BD geometry/bundle data in a computational format accepted
   by cohomCalg, Sage, or Macaulay2.
4. Run cohomology and equivariance checks as receipts, not prose claims.
5. Treat the low-energy target script as a threshold target generator until the
   full string-threshold and soft-term calculations exist.
