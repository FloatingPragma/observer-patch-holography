# Protected winding and quantized native current

The compact texture has integer shell degrees \(-1,+1\), preserved by
every record of its specified local repair law. A unit degree cannot be
split into two equal integer degrees. Independently, the actual native
current has spectrum in one-third integers on the entire fermion carrier
and integral spectrum on colour singlets. Every admitted field record
preserves those sectors, so a half-unit charge eigenvalue is excluded by
the actual dynamics. This proof requires no identification of winding
with current.

The same texture permits neutral and charged fermion occupations. Thus
winding alone cannot determine the native current on the unchanged
carrier. Physically identifying that current with observed electric
charge remains a separate dictionary problem; the quantization and
history-preservation statements hold for the declared source.

The [proof](PROOF.md) distinguishes this winding from the simulator's
integer face register. The [receipt](results.json) retains every stochastic
proposal, an exact guard sweep, a fresh 200-sweep source run and the
current-attachment counterexample. Arithmetic statements are also checked
by [Lean](Protection.lean).

```sh
python -m model.openwave_charge_protection --write
python -m model.openwave_charge_protection --check
python -m pytest -q model/openwave_charge_protection
# From the Lean directory:
lake env lean ../model/openwave_charge_protection/Protection.lean
```
