# Complete spectral records of the native quantum-flux ground family

The full208-state native source has an exactly verified charge-conjugation
symmetry and a certified N=2 ground sector. Its two adjacent-charge sector
thresholds coincide in the interval
`[1.1077131829806603,1.1077131833590101]`, in source energy units.
These are finite source thresholds. Their interpretation as observed
particle masses requires a relativistic spatial field construction.

The eight-channel addition/removal apparatus retains all four CAR modes,
both signs, every initial occupation and the quantum flux register.
References compensate charge and fermion parity. Combining this instrument
with the polynomial cooling family and the same source clock supplies
complete16-outcome, two-bin energy-transfer records. Their probabilities
are exact rationals enclosed by outward dyadic intervals in the receipt.

The finite theorem `Lean/QFT/GroundRecordTail.lean` proves that mean
ground-energy excess epsilon bounds the original weight of energy
decreases of at least eta by `epsilon/eta`. The executable supplies its
stationary state and complete instrument premises. No spectral gap or
successful-branch normalization is required. A separate finite-resolution
plan retains Fourier-filter leakage, aliasing and source-clock error.

The [proof](PROOF.md) states the reference apparatus, source domain,
GNS limit premises and exact arithmetic. Floating Lehmann weights are
available as diagnostics; they are separate from the exact certificate.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -W error -m pytest -q model/particle_ground_spectral_records
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m model.particle_ground_spectral_records --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m model.particle_ground_spectral_records --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m model.particle_ground_spectral_records --diagnostics
```
