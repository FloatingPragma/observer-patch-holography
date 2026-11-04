# Particle families, masses and defects

The conditional derivation connects the finite Standard Model representation
and three-family theorems to the same-current Higgs tensors and a canonical
field action. It gives three charged-lepton families, six quark flavors,
three minimal massless neutrino fields, the gauge fields and one Higgs scalar.
The [proof](PROOF.md) derives the permitted action and mass formulas;
the [formal bridge](FORMAL_BRIDGE.md) identifies the finite theorem inputs.
Source selection, the interacting continuum and a microscopic defect
identification are separate premises.

## The conditional particle spectrum

With `H=(0,v/sqrt2)`, the charged masses are the singular values of
`v Yu/sqrt2`, `v Yd/sqrt2` and `v Ye/sqrt2`. The calculation constructs the
full 45-component mass kernel from the retained Higgs tensors and checks
its electric neutrality and fermionic pairing. Arbitrary exact complex
Yukawa matrices, including zero and repeated singular values, are accepted.
Spin, color and antiparticle multiplicities are kept separate from flavor.

The default uses the repository's empirical common-`M_Z` calibration. These
are **running or auxiliary tree-action coefficients in GeV**, not measured
pole masses or independent predictions:

| Fields | Calculated mass coefficients in GeV |
|---|---|
| electron, muon, tau | 0.000487789, 0.102760, 1.74551 |
| up, charm, top | 0.00123594, 0.624250, 169.890 |
| down, strange, bottom | 0.00270177, 0.0537480, 2.86240 |
| three minimal neutrino fields | 0, 0, 0 |
| photon, eight gluon components | 0 |
| W+, W- | 80.8480 each |
| Z | 92.2264 |
| Higgs radial mode | 131.280 |

The bosonic formulas are `mW=gv/2`, `mZ=v sqrt(g^2+gY^2)/2` and
`mH=sqrt(2 lambda)v`. The stated tree potential is a comparison prescription
at the calibrated running scale. Loop matching is required to compare its
coefficients with observed poles. The minimal neutrino action does not fit
oscillations; the proof states the separate dimension-five mass extension.
The [complete calculation](results.json) retains all parameter provenance,
exact mass expressions, family ratios and multiplicities.

## Protected defects and interactions

The [defect catalogue](DEFECTS.md) distinguishes integer hedgehog winding,
Hopf textures and the symmetry-protected lattice index. It states their
precise protection conditions and the source-to-field maps required by the
massive-family hypothesis. Massless fields have no defect assignment.

An additive conserved defect charge must also survive every admitted
interaction. With topologically trivial minimal neutrinos and the full
charged-current couplings, all three charged leptons share a charge, as do
the three up flavors and the three down flavors. Neutral conjugate-pair
vertices force the Z and Higgs charges to zero. Thus distinct protected
windings cannot simply label the three generations under these premises.
Excited levels within a common sector and neutral composites require their
own dynamical construction. No observed particle-to-defect assignment is
claimed by the catalogue.

## Reproduce

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -W error -m model.conditional_standard_model --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -W error -m pytest -q model/conditional_standard_model
```

The receipt binds the finite theorem sources, calibration and topology
proofs. Its passing status concerns the scoped algebra and calculation;
it accepts no physical completion gate.
