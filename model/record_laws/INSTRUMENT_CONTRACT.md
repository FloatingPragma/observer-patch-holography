# The coherence table and the recorded instrument

The four requirements in the [classification](README.md#theorem-1-every-admissible-law-is-a-coherence-table)
constrain the **unread channel** `Phi`. They imply that it admits a
definite-current Kraus representation. They do not force an arbitrary
implementation of that channel to commit those current labels.

For native C, choose its two branches `A=K_(coin=0,current=0)` and
`B=K_(coin=0,current=1)` and replace them with

```text
L_plus  = (A+B)/sqrt(2),
L_minus = (A-B)/sqrt(2).
```

All other branches stay fixed. Direct expansion gives
`L_plus rho L_plus* + L_minus rho L_minus* = A rho A* + B rho B*`.
Thus completeness, the classical transition law and every property of the
unread channel are unchanged. Each new branch, however, combines different
currents. This is a change of record measurement basis before commitment;
it is not classical relabelling of already committed currents.

The [exact control](instrument_contract.py) checks equality on all 6,084
matrix units of the actual 78-state carrier. Put `x=(loads at 0,2)` and
`y=(loads at 1,2)` on seam `(0,1)`. Then `A|y>=|y>/sqrt(2)` and
`B|x>=|y>/sqrt(2)`. For the input `(x+y)/sqrt(2)` the two new outcome
probabilities are `1/2,0`; for its incoherent mixture they are `1/4,1/4`.
Even for the basis input `x`, the original moving-record probability is
`1/2`, while the new plus-record probability is `1/4`. The plus label
cannot be decoded as one current. The source channel is the same in every
case.

Consequently the classification is a classification of channel coherence
tables, with a compatible current instrument available for each table.
Statements about current histories additionally use that instrument: its
selective maps have definite current, with the prescribed diagonal branch
probabilities. Extra coin or endpoint labels need their own refinement
contract. A general record basis is not covered by the current-history
indistinguishability statement.

The repository's executed C/E instruments already satisfy these additional
conditions; their existing predictions and receipts are unchanged. The
counterexample is a finite mathematical instrument, not a new physical
implementation of coherent record access. It specializes the existing
[complete-record contract](../OBSERVER_RELATIVE_RECOVERY.md#complete-records)
to the channel classification and requires no federation, scheduler or
physical clock.

```sh
python -m model.record_laws.instrument_contract --check model/record_laws/instrument_results.json
python -m pytest -q -W error model/record_laws/test_instrument_contract.py
```
