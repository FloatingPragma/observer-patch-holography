# A finite thermometer made from native records

Three stopped native coin prefixes, `0`, `10` and `11`, have probabilities
`1/2`, `1/4` and `1/4` on every one of the 78 source configurations. Their
lengths define a nonconstant energy `E=(1,2,2) log(2)` in units of native
prefix cost. Comparing `b` freshly generated prefixes selects the reference
`pi_b(h)=exp(-b E(h))/Z_b` when they agree. The energy is fixed before that
comparison, and `pi_b(0)/pi_b(1)=2^b` identifies the integer parameter.

The [proof](PROOF.md) identifies the exact finite record-register kernel,
its stationary reference, dimensionless entropy balance and resource cost.
At `b=1`, at most two native attempts implement the exact heat bath. For
`b>1`, a bounded retry controller keeps its incoming register on timeout;
it preserves the same reference but is not an idempotent heat bath. Exact
resampling uses almost-sure retry with finite expected cost and no finite
worst-case bound. Every failed comparison and its fine source records are
retained. Fresh prefix records can be generated on one continuing carrier;
no separately prepared replicas are required.

The [certificate](certificate.py) checks the prefix law on all 78 inputs,
replays finite comparison histories, and verifies exact register stationarity,
detailed balance, supported fluctuation identities and the entropy split.
Copying one prefix instead of generating independent records is a retained
computed comparison. The construction supplies the mathematical finite-core
premises for this declared observer register; it does not identify physical
energy, temperature, erasure work or the native occupation Hamiltonian.

```sh
python3 -m model.native_record_thermodynamics.certificate --check
python3 -m pytest -q -W error model/native_record_thermodynamics
```
