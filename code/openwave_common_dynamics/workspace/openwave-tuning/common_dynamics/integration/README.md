# Joint field, particle and clock audit

The additive program reads the position, full finite field and a native clock
after the same massive matter-field evolution. It separates the clock wire
from the position flags and accounts for all reset and export failures.
The audit independently reconstructs the finite source laws and checks the
field-cutoff and unit conventions.

The scientific domain and error composition are in [PROOF.md](PROOF.md).
[INTEGRATION.json](INTEGRATION.json) contains the executed audit and ordered
block-program digest. [CRITERIA_AUDIT.md](CRITERIA_AUDIT.md) states the connection
required for each of the 31 OpenWave criteria.

`verified.complete_program` and `verified.error_contract` reject changed
producer bytes. Run `python audit.py`, `python units_review.py`,
`python verified.py`, `python current_review.py`, then `python freeze.py` from this directory to
reproduce the package. The full primitive expansion is generated lazily;
the certificate preserves exact ordered blocks and repetition counts.
Its finite checks do not establish a quantum Maxwell vacuum, observed
particle spectrum or official platform validation.
