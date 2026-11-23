# Interval Echo on the spherical screen

The row 4 experiment asks whether a prepared adjacent pair persists when every seam of the 20-cell sphere acts freely. The observer occupies one cell and records one of its 30 occupied seams or a quiet outcome. Two-load states with separated ports and states outside the cell both occur in the model. The exchange-only arm uses the same initial pair, sphere and detector. The experiment identifies a bound excitation of the specified law. The electric charge, spin and particle species require separate measurements. [Interval Echo specification](../../model/candidate_laws/interval_echo/SPEC.md) · [Public criterion](https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/MODELS.md#L104)

OpenWave awarded **⚠️ Partial** on the adopted Interval Echo pin. **Full-credit resubmission is on HOLD.** The portable contact instrument and the all-time adjacency bound strengthen the evidence for the pair. The adopted specification supplies the coupling tested in the release. [Resubmission assessment and clearing condition](row_04/CHECKLIST.md) · [Maintainer decision](https://github.com/openwave-labs/openwave/pull/745#pullrequestreview-6005365373)

The written-field extension supports a strong candidate for the electron’s charge and field sector. Its independent flux and remote-probe tests use a different law from the original two-load release. Spin, exchange and formation-energy tests determine the scope of that interpretation. For the pinned matter-wave study, both routes connect the same gateway ports through the same 194 cells. The shortest internal routes put 62 levels in the interval that is empty on the repeated four-port route. [Electron candidate and related rows](../../model/candidate_laws/interval_echo/openwave/electron_field/README.md) · [Matter-wave route comparison](../../model/candidate_laws/interval_echo/openwave/row_30_route/README.md)

## Reproduce the author measurement

Use Python 3.12 with NumPy 2.3.5 and SciPy 1.18.1. From this directory run

```sh
cd row_04
python -B row_04_check.py --board FEDERATION_LEVEL0.npz
python -B row_04_window.py --board FEDERATION_LEVEL0.npz
python -B native_check.py
```

At native exposures 1, 7 and 37, 4,048, 4,052 and 4,028 of 4,096 simulated terminal records per setting show a pair in the launch cell. The control gives 1,780, 1,182 and 824. The saved receipt includes every simulated trial in draw order as compressed address bytes. The appointments are supplied exposures. Each trial has one terminal read. A separate calculation bounds the ensemble probability of a local pair above 94% throughout the interval from 0 to 37. [Frozen local design](row_04/ROW_04_DESIGN.md) · [Raw records](row_04/row_04_receipt.json) · [Window calculation](row_04/row_04_window.json)

The third command uses the containing RER checkout to compare both generator matrices entry by entry with the pinned law, and both terminal branches with the 99-state native echo. To exercise the scientific rejection, append `--remove-contrast` to the first command. The local persistence gate rejects that alternative on the same board. The registered platform reproduction uses the recorded author settings and discloses prior knowledge of binding. An optional subsequent held release uses a different source and exposures. Its selector is the first commit on the first-parent history of official OpenWave main containing the M10.3 task. [Held release program](row_04/registered_release.py)

OpenWave's M10.3 merge selected seam (87, 92) in cell 7 for the three registered exposures, 3, 11 and 29. On 4,096 simulated releases at each appointment the original registered program counted 4,051, 4,056 and 3,939 local pair records. The exchange-only arm counted 1,797, 1,078 and 940. All four registered bounds passed. The selector fixes the launch independently of the author. Pair binding at the original strength was known before registration. Reproduce from the row directory with `python -B registered_release.py bb697927ff37446130d4c6e8c68126b6a3eac00d --output /tmp/oph-registered-row4.json`, using a fresh output path. [Public registration merge](https://github.com/openwave-labs/openwave/commit/bb697927ff37446130d4c6e8c68126b6a3eac00d) · [Raw held release](row_04/registered_release.json) · [Comparison with the energy derived from repair records](../../model/experiments/interval_echo_consensus/modular_standard.json)

The public repair simulator illustrates one patch comparing loads and retaining a repair record in a related mean-repair model. The sparse two-load simulation and its terminal seam record are reproduced by the commands above. [Public repair view](https://simulation.floatingpragma.io/patches?scene=repair)
