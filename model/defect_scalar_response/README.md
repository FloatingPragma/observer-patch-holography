# Scalar response of native charge sectors

Changing the cost of a scalar excitation changes the ground polarization
and the complete energy response of native charge-removal probes. This
package computes that susceptibility exactly for the quark and lepton
sectors of the same interacting local source. It sums every energy
outcome in an exact second-moment rule and keeps the
one-half pump failure probability. No measured particle data enter.

The [proof](PROOF.md) specifies the scalar intervention and the charge-class
interpretation; the [receipt](results.json) gives the exact results.
Relating this source-clock susceptibility to physical Yukawa measurements
requires a field-amplitude intervention and a physical observer dictionary.

```sh
python -m model.defect_scalar_response --write
python -m model.defect_scalar_response --check
```
