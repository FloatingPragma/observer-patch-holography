# Particle-simulation receipts

This directory publishes the two machine-readable receipts used for the
finite simulator statements in *Deriving the Particle Zoo from Observer
Consistency*.

- `k1_64k_completion_audit.json` records the two completed 65,536-patch
  calibration-null runs and binds their configurations and principal arrays
  by SHA-256. The simulator contains no quark, Higgs, or Yukawa coupling; the
  reported values are not quark-mass evidence.
- `poft_transport_emission_targeted_20260712.json` records the direct
  three-label permutation-transport assay on two 4,096-patch and two
  65,536-patch states. All direct POFT emission receipts are false. The assay
  excludes the natural direct carrier tested there, while leaving a future
  source-derived complex lift open.

`manifest.json` binds the archived receipt bytes. The receipts also contain
SHA-256 hashes of their upstream configurations and arrays, but those
upstream run arrays are not mirrored in this compact package. Consequently
this is public custody for the printed report values and verdict boundaries,
not an independent raw-array replay archive. Neither receipt records the
simulator commit that produced it, and that exact producer revision was not
recoverable from the retained files. Recomputing the assays therefore requires
both the matching upstream arrays and identification of the producing revision
in the public
[OPH-FPE simulator](https://github.com/muellerberndt/oph-physics-sim).
