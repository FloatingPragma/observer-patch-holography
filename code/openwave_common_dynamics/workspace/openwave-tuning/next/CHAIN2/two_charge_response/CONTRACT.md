# Execution contract

Run from the repository root with Python, NumPy, SciPy, and SymPy:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 openwave-tuning/next/CHAIN2/two_charge_response/native_contract.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 openwave-tuning/next/CHAIN2/two_charge_response/run.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 openwave-tuning/next/CHAIN2/two_charge_response/mobile_pair.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 openwave-tuning/next/CHAIN2/two_charge_response/edge_echo.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 openwave-tuning/next/CHAIN2/two_charge_response/full_v_echo.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 openwave-tuning/next/CHAIN2/two_charge_response/inverse_transport.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 openwave-tuning/next/CHAIN2/two_charge_response/finalize.py
```

Inputs are local immutable copies of the three producer freezes plus an
explicit import annex. Replay verifies their hashes before using them.
It never refreshes pins against producer directories. The source ID is
`native-local12-mass-fermion-plaquette-v3-proposed`. The command uses no
reserved acceptance configurations, changes no producer files, and
awards no official completion row.

The field is one shared set of 288 Z16 registers. Both occupied matter
modes are in the same fermionic Fock carrier. Local coins use the
exterior square, including two distinct modes at one port, and transport
retains every canonical anticommutation sign. The 896 nonphysical matter
modes and all work/complement sectors belong to the native source;
their ideal probability is zero on the declared clean physical input.

The preparation starts with 3,203 blank native source carriers, creates
two orthogonal orbitals, initializes the definite field basis
`F=P0+Py`, and retains every reset timeout. Its shared Gauss background
is `G=2delta0`. The field-string and orbital choices are declared
inputs. The detector uses 3,200 fresh native ready pointers to read all
matter occupations and electric bits, preserving every failure and
padding outcome. Source work bits are 3200, 3201, 3202. The one-particle
address-parity detector is excluded.

`NATIVE_CONTRACT.json` supplies emitted preparation words, exact rational
operator and reset bounds, full-Fock evolution counts, readout pulse
counts, and complete terminal-law TV bounds. All bounds retain an
arbitrary original spectator reference. They do not infer independent
reset probabilities or erase native global phases. Clock operations
and clock-pointer relocation are outside this consumer; no physical
clock error is borrowed by relabeling a macrostep.

The complete compressed laws store canonical occupied mode pairs,
288 field labels, complex amplitudes, and probabilities. Omitted
entries have exact zero in the declared representation. The last
current is a deterministic function of the terminal modes. Recorded
norms are unconditional; no outcome is postselected or renormalized.
Separate-field controls explicitly have two field registers and are
source substitutions.

The development limit is three macrosteps and 120,000 blocks of
16-by-16 complex128 coefficients. This gives at most 491,520,000 bytes
of final block-array payload, excluding keys, orbit work arrays,
intermediate states, Python overhead, native replay columns, and
compressed record buffers. Resource logs report actual support sizes.
Exceeding a limit raises `ResourceRefusal` before a law is accepted.
Binary64 residual checks are numerical checks, not interval arithmetic.

The dark-orbital comparison with different positions is classified as
a charged field-preparation susceptibility. At fixed initial field,
the exact source gives zero added-dark-charge effect on the probe.
The two-mobile experiment instead changes couplings at one identical
complete preparation and reports unordered-pair and current response.
Neither result establishes a Coulomb law, a physical mass scale,
reciprocal labelled acceleration, a continuum limit, or stable matter.

`FULL_V_ECHO.json` supplies an additional preparation and apparatus for
the identical full source at both supplied masses. Inverse local coins
prepare two held/moving event qubits on an edge outside face0; one full
source step is followed by the literal native inverse transport,
one-body phase calibration and event analysis. The full sixteen-label
plaquette field factor is retained in all four event branches. This
complete64-outcome interference law isolates a connected two-charge
phase. Its inverse words, resets and full readout are charged in the
native contract. `ECHO_PROOF.md` states the factorization and its scope.
`inverse_transport.py` streams every literal inverse pulse and binds its
chronological JSONL digest and full938,112-pulse count; materializing a
large native-word array is unnecessary.
