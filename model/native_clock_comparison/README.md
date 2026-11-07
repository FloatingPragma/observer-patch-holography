# Native clock comparison

The complete native Y measurement rejects identifying the scalar-record
transfer generator with the working carrier's clock while preserving its
stationary reference. Exact logarithm bounds give a Y mean greater than
1/2000 at transfer parameter 1/10. Every modular-clock parameter gives zero
on that reference. Source preparation and complete detector errors leave a
half trace separation greater than 1/5000 under the stated candidate-service
accuracy allowance.

The transfer continuation is a tested hypothesis, not an executed native
service. Its record-Hilbert ground direction differs from the working mixed
state. Reference-preserving unitary encodings cannot remove the obstruction;
general encodings or different preparations require their own comparison.
Bare repair records alone are insensitive to interleaved diagonal phases
after tracing the terminal working carrier, including arbitrary untouched
references. Reading declared clock specifications is allowed.

[The proof](PROOF.md) gives the finite-time bound, complete-output error
accounting and exact domains. The receipt imports the executed analyzer's
source-bound certificate and recomputes the error coefficients used here.

```sh
python -m model.native_clock_comparison.packet --write
python -m model.native_clock_comparison.packet --check
python -m pytest -q -W error model/native_clock_comparison
```
