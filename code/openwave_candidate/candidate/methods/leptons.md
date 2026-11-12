# Declared lepton register and fixed empirical comparison

Let `m=(tail-head)/sqrt(2)`, `p=(tail+head)/sqrt(2)` on one seam,
with a two-state record and two address orientations. Define

```text
K = |p,1><m,0| + |m,0><p,1|
Q_o = I - K^2 + (-1)^o K
W = Q R_orientation
V|load,record,o> = (-1)^(record*o)|load,record,o>
U0 = T (VW tensor VW tensor I3)
U(theta) = exp(i theta J_C) U0.
```

`T` shifts the neutral three-state register forward or backward when the
first seam's record is one. The preparation is the tensor product of
`(tail+sqrt(2) idle)/sqrt(3)`, blank records and uniform orientations.
The coherent analyzer is the equal superposition of three orthogonal flags:
first seam idle/second seam agrees, and first seam repairs at either
orientation/second seam idle. Every vector lives in the full `12*12*3=432`
space. The tick is represented exactly as the integer matrix `4 U0`.

The compressed matrix is

```text
M(theta) = (sqrt(2)/3) C(theta)
C(theta) = I/sqrt(6) + (exp(i theta) R + exp(-i theta) R^-1)/sqrt(12)
lambda_k proportional to 1 + sqrt(2) cos(theta+2 pi k/3).
```

The declared mass-read convention uses `lambda_k^2`, with the electron,
muon and tau labels assigned to characters `k=1,2,0`. The phase is the
exclusive read count `P(C repairs+ xor E agrees)`. Its interpretation as
phase per unit time consumes the declared weak-round and clock construction;
this portable calculation evaluates the tick and read count, not that full
continued apparatus. There is no spatial boundary-value minimization here.

| Equation or observable | Candidate implementation | Pinned public source |
|---|---|---|
| Seam toggle, event and orientation record | [`object_tick`, leptons.py:31](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/leptons.py#L31) | [native toggle identity](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/lepton_register_round/kick.py#L53) |
| Complete tick, preparation and analyzer | [`register_tick` and `source_and_analyzer`, leptons.py:49](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/leptons.py#L49) | [owning register code](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/lepton_register_round/register.py#L98) |
| Exact compression and complete read count | [`compression_coefficients` and `exact_read_counts`, leptons.py:86](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/leptons.py#L86) | [register proof](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/lepton_register_round/PROOF.md) |
| Squared-character ratios and reference residuals | [`spectrum` and `comparison`, leptons.py:133](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/oph_openwave/leptons.py#L133) | [fixed comparison inputs](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/lepton_cycle/experimental_inputs.json) |

Candidate links identify the OPH package in its own repository. No
OpenWave model directory or maintainer-awarded cell is assigned.

## Computed results and gates

The exact integer checks require `4U0.T @ 4U0 = 16 I432` and the stated
coherent compression coefficients. Exhausting all 432 output coordinates
gives flag probabilities `1/12,1/6,1/72`, hence phase `2/9`; the extracted
off-diagonal/diagonal squared ratio is `1/2`. Both values feed the
high-precision spectrum evaluation. Direct matrix diagonalization supplies
a separate numerical check of the three-character formula.

The declared readout yields `mu/e=206.7703159727293886` and
`tau/e=3477.4728371045985323`. The measured electron-energy anchor gives
`0.00051099895069`, `0.10565941449590446` and `1.7769849708134271` GeV.
This is one calibrated scale, not three independent GeV predictions.

The frozen CODATA2022 direct muon/electron comparison is
`206.7682827 +/- 0.0000046`. The residual is `+9.833581354` ppm and
`+442.0158107` reference standard uncertainties. Therefore a claim of exact
experimental agreement fails. No physical-theory uncertainty or mechanism
for the residual is supplied. No channel-error correction is fitted to it.

The frozen PDG2024 tau comparison is `1776.93 +/- 0.09 MeV`. The predicted
ratio has residual `+30.93583508` ppm. Keeping the unspecified cross-edition
covariance within its positive-semidefinite range gives a reference-only
standardized discrepancy between `0.61078304` and `0.61079060`. These are
marginal comparison statements, not a joint fit, likelihood or probability
that the physical theory is correct. The code uses 70 decimal-digit working
arithmetic; displayed decimal values are numerical evaluations of exact
algebraic/trigonometric expressions, not interval certificates.

The missing-orientation control keeps the full tick unitary but changes the
compression to a quadrature form and is rejected by the compression gate.
Deleting one output row is rejected by the full unitarity gate. Thus a
unitarity check alone does not validate the mass-read construction.

## Architectural choices and omissions

The two independent seam objects, their prepared amplitudes, coherent
analyzer, orientation record, three-state representation, weak-round clock,
square-of-amplitude readout, particle-label assignment and electron anchor
retain their declared or calibrated status. The public parent additionally
names the double-load, toggle, two-cell, flux and return-test premises.
These are not nine independently measured parameters, and no counting
argument here proves that the empirical agreement is parameter-free.
Particle masses and Koide relations were known during architecture design.

Three Fourier characters do not establish exactly three physical mass
minima. The code does not compute spatial localization, interacting lepton
poles, independent absolute units, radiative corrections or a finite-round
mass error. OpenWave's full electron-rest-energy and lepton-spectrum
criteria are therefore not marked satisfied.

**Scope**: exact reconstruction of the declared 432-state seam register and
retrospective comparison of its squared-character readout, with one empirical
GeV anchor; no independent physical mass spectrum or OpenWave cell is inferred.

## Reproduction and review

```sh
python -m oph_openwave.leptons
python -m unittest discover -s tests -p 'test_leptons.py'
```

Minimal inspection set: `oph_openwave/leptons.py`, `tests/test_leptons.py`,
this note, and the common run receipt. The producing implementation was
also compared against both oriented and unoriented literal `432x432` ticks
in the pinned RER source; every integer entry agrees. The [independent adversarial audit](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/evidence/source_parity_audit.json) checks literal source parity and a separate Hermitian eigensolver, detects changed tick/orientation controls, and shows that changing the phase preserves Koide while changing the ratios. The [audit script](https://github.com/FloatingPragma/observer-patch-holography/blob/main/code/openwave_candidate/candidate/audits/source_parity.py) supplies the computation. Conditional coefficients are confirmed; physical mass identification is not established.
