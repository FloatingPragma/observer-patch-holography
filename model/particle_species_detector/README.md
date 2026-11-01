# Identifying the model's propagating Weyl particles

This detector identifies **eight operational Weyl flavors**, labeled by
the eight lattice valleys `b in {0,1}^3`, with both signed energy branches.
The labels are attached to a concrete native-source Fourier/occupation
instrument, not assigned from an energy fit. On the stated periodic
`8×8×8` lattice the measured family has

\[
N=1,\qquad (-1)^N=-1,\qquad
E_\pm=\pm\sqrt2/16\simeq\pm0.0883883476483.
\]

These are fermionic one-particle modes of the supplied CAR field. The
positive branch has four flavors of each cone helicity `±1/2`. All eight
have continuum mass zero. They are model particles with operational
momentum, energy and helicity labels; no electron, photon or other
empirical species name is inferred.

```sh
.venv/bin/python -W error -m pytest -q model/particle_species_detector
.venv/bin/python -W error -m model.particle_species_detector --write
.venv/bin/python -W error -m model.particle_species_detector --check model/particle_species_detector/results.json
```

## The identified family and the original field

Keep the original periodic action
`H=dGamma[-c sigma.sin(k)]+U sum_x n_x0 n_x1`, with `c=1/8` and `U>=0`.
The one-particle sector is exactly invariant and contact is zero throughout
it. Every one of the 1,024 physical CAR modes remains in the computation.
The detector does not replace them by ten positional qubits.

Use lexicographic site order and component order,
`mode=2*((x*8+y)*8+z)+spin`. For each valley `b`, choose momentum labels
`ell=(4b_x,4b_y,4b_z+1)`, with `k=pi ell/4`, and either component
`r=+1` (spin label zero) or `r=-1` (spin label one). The normalized mode is

\[
\psi_{b,r}(x,s)=512^{-1/2}e^{ik\cdot x}\delta_{s,(1-r)/2},
\qquad E_{b,r}=-r(-1)^{b_z}\sqrt2/16.
\]

`algebraic_modes()` represents every amplitude exactly as
`(A+sqrt(2)B)/32` with Gaussian-integer `A,B`. Exact matrix multiplication
proves normalization, mutual orthogonality and `H psi=E psi` for all
sixteen columns. The generating rule is independently checked against
the existing periodic source matrix on its `4^3` domain.

Consequently the rank-one effects obey `H Pi_label=E_label Pi_label`,
and `P_out=I-sum_label Pi_label` commutes with `H` on the full one-particle
space. These are actual spectral energy labels of the source Hamiltonian.
The other 1,008 Fourier/component records remain present. They do not all
diagonalize `H`, and the detector does not assign the family energy formula
to them.

The exact total directional particle current is
`J_z(k)=-c cos(k_z) sigma_z`. On this family it has the same eigenvalue as
`H`: `±sqrt(2)/16`. Every forward z bond carries signed current
`±sqrt(2)/8192`. This is a nonzero propagation-current witness even though
an individual plane wave has uniform density. Current is stationary on
these modes; `[H,J_z]` is nonzero at generic momenta, retained as a negative
control. No global conservation claim is made for that component.

Define `s_j=(-1)^b_j`, `chi=-s_x s_y s_z`, and the proper valley frame
`R_b=-chi diag(s)`. The cone-spin helicity on the selected positive-z
offset is `s_x s_y r/2=chi sign(E)/2`. This is the calibrated spin frame of
the [massless continuum construction](../relativistic_sector_audit/README.md).
The raw detector measures momentum and **CAR component occupancy**; its
helicity label uses this explicit frame map. Bare `sigma_z/2` differs from
that helicity for some valleys. Physical Lorentz spin and electromagnetic
charge are not silently attached.

| Flavor `b` | Momentum labels `ell` | Positive-energy component | Positive cone helicity |
| --- | --- | ---: | ---: |
| 000 | (0,0,1) | 1 | −1/2 |
| 001 | (0,0,5) | 0 | +1/2 |
| 010 | (0,4,1) | 1 | +1/2 |
| 011 | (0,4,5) | 0 | −1/2 |
| 100 | (4,0,1) | 1 | +1/2 |
| 101 | (4,0,5) | 0 | −1/2 |
| 110 | (4,4,1) | 1 | −1/2 |
| 111 | (4,4,5) | 0 | +1/2 |

The other component at each momentum is its negative-energy branch with
opposite cone helicity. The eight flavors are distinguished by measured
valley momentum; their equal energy magnitude alone would not distinguish
them. They are not eight independently postulated microscopic fields.

## A constructive complete detector

Let `F[k,x]=exp(-2pi i k.x/8)/sqrt(512)`, separately on each component.
The detector applies the fermionic Gaussian extension `Gamma(F)` and
copies every output occupation to a fresh pointer. This yields all binary
occupation strings, not only the one-hot outcomes. Vacuum, all other
one-particle records, multiple occupations and every source/export abort
remain in the output dictionary. On `N=1` its 1,024 effects are
`F* |ell,spin><ell,spin| F`. On other number sectors the same measurement
reports free-mode occupations, not total interacting energy.

The `native_fft` module gives an explicit three-axis FFT8 network:
4,608 two-mode butterflies and 768 bit-reversal transpositions. Every
nonadjacent operation has an adjacent FSWAP route; the determinant phase
on double occupation is retained. A butterfly is a second-mode phase
followed by the existing hopping word with coefficient `-i` and duration
`pi/4`. An exact calculation in `Z[zeta_8]/(zeta_8^4+1)` checks all 64
one-axis entries and all 64 Gram sums. Tensoring the three transforms
establishes the complete mode map. Full four-column native butterfly
checks include the vacuum and two-occupation scalar.

One Fourier transform uses 1,260,800 adjacent FSWAPs and 32,909,824
formal native pulses. Prepare a selected mode by exciting one local
occupation from the full vacuum with the native `Y(-pi/2)` word, followed
by the inverse network. Continue with the original `H`, then apply the
forward network and the 1,024 complete occupation pointers. The complete
preparation-and-detection word has 65,822,721 formal pulses. The original
Hamiltonian continuation is a separate specified model-law service; those
pulse counts do not include a synthesized field evolution.

Exact literal-angle and source-control certificates give a full-output
error below `4.153e-9` for each Fourier transform. Including the excitation
and all pointers gives

\[
\epsilon_{\rm source}<8.31\times10^{-9}.
\]

All first-export failure branches contribute to this norm. The bound
holds on the full Fock space and arbitrary untouched references, except
that the stated excitation preparation starts from the supplied vacuum.
Interposed original-H evolution is unitary on every number sector, so
it does not amplify the preparation error. No claim that the literal
native map is exactly Gaussian or number preserving is used.

For an ideal prepared member of the catalog, detector output is its exact
label with probability one at every model time. Half trace distance of
the complete native and ideal outputs is at most the vector error.
Therefore the **unconditional** probability that the specified native
preparation/evolution/detector reports its correct label is at least
`1-8.31e-9`. The wrong-label, vacuum, outside, multiple-occupation and
abort probabilities together are at most this error. Nothing is divided
by a survival probability.

## Execution, counts and the comparison boundary

The replay executes the ideal Gaussian network on all 1,024 one-particle
rows, prepares all sixteen columns, evolves them under the original
source matrix at time eight, and executes the detector. It retains every
output probability. The complete-vector discrepancy from the exact
spectral prediction is about `1.94e-15`; the Gram discrepancy is about
`3.90e-15`. These floating diagnostics are separate from the rigorous
native-source bound. They are **not** a direct execution of a
`2^1024`-row literal-native state or its complete abort array.

A PCG64 illustration with seed 20261101 samples 65,536 records from the
exact uniform mixture of all sixteen ideal labels. Every address count,
including zero outside counts, is retained. This is model simulation,
not laboratory data or a literal-native word run. Hoeffding and a union
bound give simultaneous probability error at most `1/64`, except with
probability below `32(3/8)^32`; the source error is separate. The rigorous
correct-label theorem above does not depend on this sampled illustration.

For the finite periodic refinement `L=8*2^m`, `a=2^-m`, `c=(1/8)/a`,
the physical model period stays eight and the offset momentum stays
`p_z=pi/4`. Its exact one-particle energy is

\[
E_\pm(a)=\pm\frac{1}{8a}\sin(a\pi/4)
\longrightarrow\pm\pi/32\simeq\pm0.098174770425,
\qquad m_{\rm continuum}=0.
\]

The receipt encloses these energies rationally for five spacings. This
finite-torus statement follows directly from the same eigenvalue formula;
it does not assume an infinite-lattice wave-packet limit. The finite-spacing
mass-shell defect is a regulator correction, not a tachyonic mass.

The empty-sector reference is not the ground state; the negative branch
has been retained and is not relabeled as a physical antiparticle. No
stable interacting vacuum, physical charge, positive massive rest energy
or identification with the protected cube's neutral response is claimed.
The separate retrospective physical comparison records what the measured
model quantum numbers do and do not match; it fits no energy scale.

## Comparison with measured particles

The fixed comparison baseline is PDG 2025, not a claim to use a later edition.
Its target data are read only by `comparison.py`; they do not enter the
field, detector, energy calculation or source parameters.

| Measured particle | PDG mass in MeV/c² | These continuum modes | Comparison status |
| --- | ---: | ---: | --- |
| Electron | 0.51099895000 ± 0.00000000015 | 0 | Rejected direct identification |
| Muon | 105.6583755 ± 0.0000023 | 0 | Rejected direct identification |
| Tau | 1776.93 ± 0.09 | 0 | Rejected direct identification |

The primary sources are the PDG [electron](https://pdg.lbl.gov/2025/listings/rpp2025-list-electron.pdf),
[muon](https://pdg.lbl.gov/2025/listings/rpp2025-list-muon.pdf), and
[tau](https://pdg.lbl.gov/2025/listings/rpp2025-list-tau.pdf) listings.
A common positive conversion of units cannot map any of these zero
continuum masses to a positive charged-lepton mass. Using the equal finite
energies as masses would instead give all ratios one, and would also
confuse kinetic energy with rest mass.

Masslessness alone is compatible with the photon, but the present modes
have cone-helicity magnitude `1/2`, whereas the [photon](https://pdg.lbl.gov/2025/listings/rpp2025-list-photon.pdf)
has spin one. A spin-preserving direct identification therefore fails.
For three such massless modes all mass-squared differences are zero.
This disagrees with the nonzero splitting ranges in the selected
three-neutrino fit in [PDG Table 14.7](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-neutrino-mixing.pdf):
`Delta m21² in [6.92,8.05]×10^-5 eV²/c⁴` and
`Delta m32² in [2.388,2.531]×10^-3 eV²/c⁴` for the stated normal-ordering
analysis. These are reported fit ranges, not exact mathematical bounds.
A single massless neutrino is not excluded by this test, nor is a photon
constructed as a composite or in a different sector.

The result is therefore an operational identification of the model's
Weyl flavors, with **zero direct identifications with measured physical
species**. This is a scoped comparison of this constructed sector, not a
no-go theorem for other OPH sectors or mass-generation mechanisms.
