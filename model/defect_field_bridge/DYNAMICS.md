# Source-executed dynamics of the colored two-particle sector

This experiment combines the actual native C/E repair, its ordered-cover
source implementation, a colored fermionic exchange code, and the static
quasiparticle Hamiltonian in `model/defect_field_benchmark/quasiparticles.py`.
The calibrated matter inputs are explicitly supplied running parameters at
the stated scale. The calculation is a forward realization and observation
test, not a derivation of those inputs from repair counts.

For internal space V, the encoded colored state is

```
E_col psi = (|0> psi - |1> SWAP psi)/sqrt(2).
X_order E_col = E_col (-SWAP).
```

Thus the existing source-implemented ordering flip induces the actual
two-fermion `-SWAP`, including its internal-label permutation. Both the
encoder and its inverse have literal source words: seven controlled swaps
between two seven-bit registers, an ordering-bit rotation, and its scalar-
corrected Z operation. Each local controlled swap is checked on all eight
inputs. The experiment executes the complete 15-qubit preparation vector,
including every padded address. It then contracts the native source
instrument with this state using exact finite sums of tensor factors. No
Schmidt truncation, successful-branch normalization, or probability clipping
is performed. The inverse encoder retains both ordering outcomes. This last
map is an ideal algebraic contraction, checked against its nine local source
factors, rather than a direct execution on the entire joint register.

The positive-mode Bogoliubov basis identifies a static quadratic Hamiltonian
with `E_vac + sum_a E_a gamma_a^dagger gamma_a`. Its two-quasiparticle sector
is invariant. Attaching that 90-mode internal space to each of the twelve
native ports is an explicit spatial replication: the underlying CAR space
has 1080 modes. This packet uses its hard-core two-quasiparticle subspace,
with distinct occupied ports, not the original 90-mode Fock space. This is
distinct from the bare two-fermion sector, which the
pairing term need not preserve. Electric charge also preserves the chosen
positive-mode space. Broken gauge generators and a changing Higgs background
can have anomalous quasiparticle blocks and are outside this execution.

The designed scheduler uses the actual native current `j` to program

```
t_j = t*(1+j/4),
U_j = exp[-i*t_j*(E+E_vac/2) + i*kappa*j*Q].
```

The joint branch is the colored native partial injection followed by
`U_j tensor U_j`. The spatial replication would sum twelve local vacuum
constants. We choose a common energy-reference counterterm `-11*E_vac`, so
this packet retains the displayed single `E_vac`. This declared reference
does not affect fixed-two-particle Ramsey gaps, and any extension to other
particle-number sectors must retain its separately specified vacuum law.
The joint branch is not obtained by
taking an exterior power of a nonunitary Kraus operator. Completeness follows
from the complete native effects and the full unitary phase gates; the actual
output tree checks it numerically. This scheduler and clock-unit dictionary
are explicit designs, not an energy inferred from the number of repairs.

`dynamics_controls.json` freezes the programmed values rounded once to
`1e-9 GeV`. Replay never recomputes a new tape by rounding platform-dependent
eigenvalues. Instead it verifies the frozen values against the fresh complete
quasiparticle spectrum with numerical admission `1e-8 GeV`. All 90 modes and
their charge multiplicities remain; eigensolver vectors are not persisted as
program inputs. The source executes all 128 columns of each seven-qubit phase
gate, with identity on the padding. The quantum output factors retain all
complex phases and all padded addresses.

The neutral test pair contains one opposite-charge partner and a coherent
superposition of two equal-charge, distinct-energy quasiparticle modes.
After native repair and evolution, a record-selected Ramsey read follows the
transported internal label. It retains all C eight and E sixteen raw native
slots, both decoded ordering outcomes, and all `128*128` internal pointer
outcomes. The independent predicted selected probability is
`p_native*cos(t_j*(E_b-E_a)/2)^2`. This is an observed **energy difference**;
it is not a measurement of an absolute rest mass or the vacuum offset.
The latter remain specified action data in this experiment.

A fresh neutral pointer then measures the same particle's specified electric
charge without a superposition of different total charges. Its source word
is a Ramsey rotation, a controlled phase `exp(i*(Q+pi/2))`, and the inverse
rotation. All 256 input columns are executed and compared with that operator.
For a charge eigenvalue `q` its zero-output probability is `(1-sin(q))/2`;
this is distinct for each of the seven retained values from -1 to +1.
Both pointer outcomes and their quantum outputs are retained in every raw
native branch. The word is QND because each pointer Kraus operator is a
function of Q. This is a specified charge readout, not an electromagnetic
field-response experiment: the common EM phase cancels on the neutral pair.

The saved arrays store source factors, internal factors and complete pointer
probabilities. They determine every unnormalized post-measurement source
vector. Complete probabilities are evaluated as sums of squared amplitudes
in bounded blocks; this combines executed source operators with the ideal
decoder and ideal computational-pointer contractions, rather than directly
simulating one giant register. Every decoder factor and the four-column
pointer-copy gate is independently source-executed. The nine decoder factors
and fourteen pointer copies each satisfy numerical operator admission
`epsilon=2e-7`. Telescoping bounds their combined operator discrepancy by
`delta=(1+epsilon)^23-1 < 4.601e-6`; for a normalized input the complete
probability total-variation discrepancy is at most `delta+delta^2/2`.
These are analytic composition bounds conditional on the numerical local
checks, not outward-rounded numerical certificates. The source
instrument's arbitrary-reference guarantee supplies its action on the
internal registers. Native preparation, internal preparation, phase gates,
decoder, Ramsey operation and all fourteen internal pointer copies are
counted, including the charge probe and its fresh neutral pointer. The
complete saved probability shape per raw branch is
`(2 decoded ordering, 2 charge pointer, 128 first internal, 128 second internal)`.
Proper-code complement outcomes are retained with exact ideal
probability zero by the invariant proper-code source law; their slots are
included in the resource accounting. Floating execution is not an interval
or hardware proof.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.defect_field_bridge.dynamics_packet --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/defect_field_bridge/test_dynamics.py
```

The native twelve-port positions are not identified here with physical
momentum or spacetime, and this fixed-background sector does not implement
pair-creating interactions or a dynamical vacuum. It supplies a complete
same-process colored exchange and measured spectral-response experiment;
those further physical identifications must use this same carrier and law.
