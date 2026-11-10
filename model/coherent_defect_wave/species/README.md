# Shared-link occupation extension

The coherent defect-wave clock has an exact fermionic occupation extension.
All wave occupations use one common ten-state link register. The one-wave
restriction equals the parent clock exactly; the full exterior algebra supplies
CAR, and the wave-number-controlled repair preserves the same interior Gauss
constraint. A normalized constrained two-wave preparation has an exact complete
five-face readout. This is a declared finite fermionic architecture, with no
physical spin, relativistic dispersion or measured particle identification.

The spatial occupation circuit uses the exterior lift of a dressed one-body
hop between neighboring cells. Its links stay in their original cells, both
Gauss constraints hold, and a two-wave detector records cell occupancies
0, 1, 2 with probabilities 1/4, 1/2, 1/4. The gate differs from partial SWAP
of two permanently occupied tensor registers. That distinction prevents an
unsupported transfer of the tensor-wave dispersion calculation.

A complete two-cell readout combines the actual whole-cell translation with
the same spatial CAR clock. A localized constrained excitation gives four
momentum/phase records with probabilities `1/120, 59/120, 1/4, 1/4`.
The clock also moves the retained defect flux with positive probability.
The measured clock bins are finite phase responses, without a chosen energy
logarithm or a physical mass assignment.

```sh
python -m model.coherent_defect_wave.species --check
python -m pytest -q model/coherent_defect_wave/species
```

The [proof and premise crosswalk](PROOF.md) explain the shared-link lift, its
architecture choices and the conditions needed to attach finite field theorems.
The [candidate table](CANDIDATES.md) separates protected flux, wave channels and
the fifteen-state auxiliary Standard Model representation. The
[exact receipt](results.json) includes all detector outcomes and code ancestry.

The single-wave clock commutes with the full internal `U(3) x U(2)` action.
Its three identical triplet copies cannot acquire different spectral masses
from that clock alone. The internal charge pattern of a declared exterior
representation is compatible with one Standard Model family. Physical fermion
statistics belong to the separately declared wave Fock algebra; even auxiliary
occupation masks do not establish them.

The parity-controlled multi-wave repair and the presence-controlled repair
coincide on one wave but produce different two-wave records. The receipt
retains that comparison as evidence that the one-wave law does not uniquely
fix its many-wave extension. Neither extension uses observed particle masses.
