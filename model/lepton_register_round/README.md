# The charged-lepton generation register

Two seam objects on one face of the carrier and a neutral register carrying
the face third-turn form one 432-state tick. Read on the Fourier analyzer of
the flags "E agrees", "C repairs +" and "C repairs -", the tick compresses
to `(sqrt2/3) C(theta)`, the balanced circulant with `|b|/a = 1/sqrt2`. The
balance and the phase are read counts: the hold and the two oriented repairs
weigh `1/6`, `1/12`, `1/12`, and two independent objects give the exclusive
read count `1/12 + 1/6 - 2/72 = 2/9`.

The package checks this on the public seam program in exact arithmetic,
together with the weak round that turns the count into the phase of the
oriented current (`2/9` per unit clock time in the limit of many rounds), the
unit of the declared clock, the kick as the commutator of an oriented toggle
with the tick of the pair, the counts on the sixty-address carrier, and the
selection of the pair on one port by a class-3 flux. A source of the double
load and the lepton scale are open.

From the repository root:

```sh
python -m model.lepton_register_round --write
python -m model.lepton_register_round --check
python -m pytest -q -W error model/lepton_register_round
```

[PROOF.md](PROOF.md) states each check and separates exact computation from
declared premises. [results.json](results.json) holds the tables, the 80 named
checks and the input hashes.
