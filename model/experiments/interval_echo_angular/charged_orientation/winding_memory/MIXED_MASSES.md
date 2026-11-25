# Positive formation lines on a mixed winding guide

The registered law permits a mobile insertion below 0.05694 native energy units. Its energy depends on the prepared pattern of blank and unit-contrast contacts. Subtracting the unit-filled material with the very same heads gives the same spectrum for a vacancy and its charge conjugate. All private references on the route are erased.

| Physical guide | Lowest positive formation line | Phase-zero edge for the repeated motif | Width of the lowest band |
|---|---:|---:|---:|
| 15 links, seven even marks separated by blanks, one double blank | 0.034436490 | 0.033961403 | 0.00093454 |
| 144 links, nine periods of 16, a triple blank per period | 0.034423786 | 0.034423112 | 0.00002756 |
| 144 links, three periods of 48, a triple blank per period | 0.003821664 | 0.003821511 | 0.00000104 |
| 144 links, strict even alternation | 0.000411524 | 0 exactly | lower band touches zero |

The 15-port census covers all 32,767 nonempty subsets of unit marks at a fixed launch on the physical cycle. Identifying identical lifted paths and their reversals gives 8,383 cases per homogeneous contrast parity. Five even cases and zero odd cases have a positive eigenvalue below the cutoff. The smallest such level lies in $[0.034436489,0.034436491]$. The omitted all-blank ring has energy at least 16. Among the 128 individual parity assignments to the seven marks of the displayed 15-link motif, three give a positive eigenvalue below the cutoff. The other nonuniform parity layouts are outside this census.

The cutoff follows from the one-defect ceiling of 198 and the laboratory mass ratio $1776.93/0.51099895069=3477.3652619$. An electron in a one-defect family would need native formation energy at most 0.056939661234 if the tau respects the ceiling. The muon/electron comparison uses 206.7682827. All laboratory values enter after constructing the native spectra. [One-defect bound](../../../../candidate_laws/interval_echo/openwave/row_11_spectrum/marked_bands/material.py) · [PDG tau listing, p. 1](https://pdg.lbl.gov/2025/listings/rpp2025-list-tau.pdf) · [CODATA electron and muon values](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf)

## Which gap survives increasing the guide length?

For an even unit mark, after host subtraction the enabled pair contributes $-4|x+y|^2$. A blank contributes $8(|x|^2+|y|^2)+|x-y|^2/2$. Substituting $-|x+y|^2=|x-y|^2-2(|x|^2+|y|^2)$ writes the period form as a positive weighted difference sum plus an on-site coefficient $8(d_{\rm blank}-d_{\rm mark})$. For strict alternation each coefficient vanishes. A constant wave then has zero energy. A recurring double or triple blank produces positive coefficients at the extra blank sites. A null vector would have to be constant and vanish there, so the repeated motif has a positive lower band edge. Taking entrywise magnitudes reduces the difference form, placing the minimum at zero Bloch phase. This is a phase of translation through a prepared guide. Conversion to momentum in observer spacetime requires a separate spatial calibration.

The physical head can reverse once in either direction. With the direction initially fixed, the reachable configurations unwrap to a finite path. A blank lets the charge travel without writing a head. At each end of the path the next mark blocks another turn. The code constructs this lift from the exact native blocks with diagonal and hopping coefficients $(8.5,-0.5)$ on a blank, $(-4,-4)$ on an enabled even mark and $(-4.5,-4.5)$ on an enabled odd mark. A blocked marked hop has zero relative energy. Signed endpoint currents and the same-head subtraction fix these entries for both charges. Integer determinant recurrences count eigenvalues strictly below the laboratory cutoff and handle zero roots separately. The physical bottom levels and the positive phase-zero edges in the table have rational enclosures of width $2\times10^{-9}$.

## Local source, receiver and mass comparison

The long ring traverses 48 cell owners. Removing the remote compensator and all marked guide edges leaves a biconnected graph for each displayed motif. Each write can borrow an ordinary reference while holding the other endpoint fixed, as in the [winding preparation](README.md). Addressed full transfers carry a sharp vacancy to the middle of a period. At release the entire masked guide evolves freely. The detector one period away reads its own vacancy. On the period-16 guide the finite bottom line has source weight 0.002400 and source-to-receiver pole amplitude 0.004566. At one inverse bandwidth, the free release gives 861 arrivals in 32,768 simulated launches, with exact law probability 0.026524. The broadband count includes all lines. The pole weight is calculated separately from the native generator. A frequency-resolved observer measurement on the same prepared line would test that pole separately.

The tested periods on the long ring are 12, 16, 18, 24, 36, 48 and 72. Anchoring the lowest line in each spectrum and testing every other line yields no simultaneous laboratory ratio match. For period 48 the closest ratios are 205.18547 and 4323.61135, respectively 0.7655% below the muon and 24.3358% above the tau. Each level describes an excitation of the specified loaded guide. The calculation supplies neither an intrinsic rest mass in MeV nor a physical spin on those levels.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.interval_echo_angular.charged_orientation.winding_memory.mixed_masses --write
```

[Certified levels, apertures and counts](mixed_masses.json) · [Executable sweep and exact enclosures](mixed_masses.py)
