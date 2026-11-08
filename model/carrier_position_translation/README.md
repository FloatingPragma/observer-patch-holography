# Position changes carried by native records

The actual rank-three carrier readout obeys an exact operator covariance
law on every recorded repair branch. It holds for coherent inputs and an
untouched reference, and retains every outcome. The repair maps are partial
translations with state-dependent acceptance. Their complete uniform mean
contracts the readout. The protected spin field's total-occupation readout
is constant on its spin sector, so this detector cannot locate those field
excitations. [PROOF.md](PROOF.md) gives the operator identities, uncertainty
bounds and precise boundary of the proposed position–momentum identification.

The certificate checks all 180 native branches and 4,680 source columns,
the complete adaptive two-event observer, and the actual protected cube's
256 spin states and 12,288 repair columns. It uses the existing source
geometry and repair definitions without introducing a spatial lattice or
translation generator.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.carrier_position_translation --check model/carrier_position_translation/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q --strict-markers --strict-config model/carrier_position_translation
```
