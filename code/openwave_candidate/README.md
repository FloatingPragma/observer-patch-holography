# Finite OPH comparison computations

This bundle supplies portable, headless computations and independent
checks of explicitly declared finite OPH constructions. Its scientific
source basis is public RER commit
`3bcacf1276ddb3ead35be2b9adcd63b6643e2d35`. OpenWave dependency compatibility
is tested against commit `ef55175d04737e73c92f33f9f9cf3d17c593da54`.

Start with the [package and commands](candidate/README.md). Each method
note starts with equations, links their implementation and lists the
physical quantities not computed. The
[definition-only dispersion comparison](certification/README.md) retains
two independent implementations, full polynomial comparison and
altered-input controls.

The [run receipt](evidence/portable_run.json) records the five calculations.
The [public-clone check](evidence/clean_public_reproduction.json) records
installation from the pinned public revision, all 19 package tests and
exact agreement of all five scientific outputs in an isolated environment.
Supporting audits retain the exact charge lattice, full finite spin
kernel, original repair branches, calibration ancestry, source parity and
negative quark-splitting result. The export manifest authenticates the
public bundle's files.

**Scope**: finite mathematical computations and conditional readouts;
no OpenWave admission or official cell, complete physical spectrum,
physical Yukawa derivation, annihilation calculation or cosmological
abundance prediction is inferred. The private planning/scoreboard
documents are not part of this bundle.
