# Complete lepton register cycles

The native sixty-address source executes both the phase round and the
mass-read round with all four return outcomes retained. A declared fresh
preparation and SWAP reset preserve every rejected source state and its
quantum correlations. The repeated complete channels have explicit error
bounds and resource counts. Their controlled limits generate the declared
current and squared-register-eigenvalue rates; physical pole masses and
an independent energy unit require a separate attachment.

[The proof](PROOF.md) gives the full instruments, the rank-three nonreturn
memory, stationary reference families and parameter sensitivities.
[The receipt](results.json) records exact identities and native numerical
checks. The muon/electron residual is reported without changing a parameter.

```sh
python3 -m model.lepton_cycle --check
python3 -m pytest -q model/lepton_cycle
```

The source is the existing
[lepton register](../lepton_register_round/PROOF.md). Its frozen receipt and
tests are unchanged. Fresh double loads, preparation gates, controlled kicks,
inverse ticks, return apparatus and growing archive storage are declared
services, with their error allowances stated explicitly.
