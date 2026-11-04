# Relativistic dispersion of the supplied lattice masses

The periodic mass source has exact moving one-body bands. All eight lattice
valleys are retained, giving two Dirac copies at each rest energy
`|m0+|m||` and `|m0-|m||`. Their continuum dispersion is
`E²=v²|p|²+M²`, with a proved finite-spacing error. The bounded radial
source supports the pure-spin branch in wave packets arbitrarily far from
its pinned core. The [proof](PROOF.md) bounds the packet boundary and the
varying radial texture on that same operator.

For the radial bulk parameters `m=64,v=1,a=1/64`, the positive energies are:

| Model momentum | Lattice bulk energy | Continuum energy |
| --- | ---: | ---: |
| `(0,0,0)` | 64 | 64 |
| `(1,0,0)` | 64.00781138753645 | 64.00781202322104 |
| `(1,1,1)` | 64.02343130344798 | 64.02343321003646 |

These decimals display rationally enclosed values. The lattice shell error
on `|p_j|<=1` is at most `1/4096`. A specified compact exterior packet
has continuum eigenvector residual below `0.000606`; its analytic bound
requires `16*10^18` spin-site rows and does not claim their execution.
The numerical control executes 11,664 rows including its boundary.
The receipt also evaluates the supplied `a=1` radial source: at `(1,0,0)`
its exterior bulk energy is `64.00553158453005`. The smaller spacings are
declared refinements of that one-body action; the finite interacting core
theorem is not transported along them.

Matching model light speed one to `299792458 m/s` fixes
`length_unit=c_SI*T`, where `T` is seconds per Hamiltonian time unit.
Its energy unit is `hbar/T`. The [units calculation](CALIBRATION_PROOF.md)
computes the conversion coefficients and accepts an explicitly identified
same-generator physical clock. It supplies no absolute `T`; absolute GeV
mass entries therefore remain null.

The localized interacting charged gap is distinct from these exterior
one-body bands. Motion of the defect, interacting relativistic poles,
physical particle identification and a physical time calibration require
additional source results.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -W error -m pytest -q model/source_mass_shell
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -W error -m model.source_mass_shell --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -W error -m model.source_mass_shell --check
```

The [receipt](results.json) keeps checks of exact source matrices, all signed
branches and zero-mass multiplicities, rational intervals, source hashes,
and numerical admission decisions. Its strict replay checks the complete
calculation; no raw platform-dependent floating residual defines an interval.
