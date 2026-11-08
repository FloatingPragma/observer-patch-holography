# Local charged preparations and particle comparison

The ordered source has bounded preparations in compact spatial regions
with a certified positive energy-momentum response. Its configuration
probabilities give nondegenerate finite creation Gram matrices. Filtering
with the same source clock yields separated, nonzero equal-sign charge
configurations whose entire energy is positive. These results provide
inputs for a particle comparison; a source dispersion and controlled
incoming and outgoing limits are required for scattering.

[The configuration proof](PROOF.md), [local preparation bounds](LOCAL_CREATORS.md),
and [positive multipacket construction](POSITIVE_MULTIPACKETS.md) state the
analytic domains and complete-record probabilities. The executable also
checks a complete local instrument on the actual selected eight-site
source. Its nilpotent local creator provides a control against transferring
the continuum Gram theorem to a finite matrix.

Run from the repository root:

```sh
python -m model.spin_multipacket_comparison --check model/spin_multipacket_comparison/results.json
python -m pytest -q model/spin_multipacket_comparison
```

The receipt replays the exact source operators and analytic bound producers.
It binds the source proof inputs and the retained ground certificate. The
continuum conclusions are analytic theorems; the finite matrix experiment
is not a continuum numerical simulation. No physical species, mass,
spin-statistics relation or scattering cross section is assigned.
