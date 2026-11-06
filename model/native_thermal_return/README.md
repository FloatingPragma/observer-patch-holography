# Native current readout for thermal return records

The [current-only controller](../native_clock/configuration.py) implements
the full-configuration comparison used by the native closed-history thermal
construction. An initial census restores and reports the binary source
configuration. Its signed-current records then supply the local field word.
A terminal census reports the body endpoint for comparison. Every census
failure, waiting timeout and different-state result remains an outcome.

The [proof](PROOF.md) and [certificate](certificate.py) provide two separate
bodies: regional incident-clock histories, and fixed global-attempt histories
with full configuration words. The latter has no waiting cap. Both charge
2648 initial-restoration attempts and 2520 terminal-readout attempts at the
stated precision, with complete public marginal error at most `17/2^29`.
A thirty-seam covering sweep and a bounded native warm-up replace a supplied
uniform starting law; their preparation error and attempts are explicit.
Conditioned return probabilities include their normalization error.

For the incident-clock construction, the source transfer determines
`H=-log T` before the selected loop length supplies `beta=L`. The complete
attempt variant attaches the same readout to the separate transfer `P` and
its reconstruction. These are finite Gibbs descriptions of selected native
history ensembles. The construction supplies no physical energy scale,
Kelvin temperature, laboratory clock or coherent native Gibbs preparation.
The comparison is the stated public marginal, with private tapes physically
retained and ignored by subsequent control.

```sh
python3 -m model.native_thermal_return.certificate --check
python3 -m pytest -q -W error model/native_thermal_return
```

The producer freshly recomputes both 66-column census success tables, checks
the native current inverses, compares capped incident histories with an
independent endpoint oracle, and enumerates the complete-attempt branch law.
It rechecks the inherited source pins for the readout, transfer and mixing
certificates. It does not rerun every upstream numerical campaign.
