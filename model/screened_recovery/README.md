# Screened `q=1` code: compiled continuing recovery

This packet closes the source compilation gap for the sharp `q=1` qubit from
the [particle screening packet](../particle_screening/README.md). It emits a
nine-qubit encoder, the complete record-selected recovery menu, and a two-step
C/E continuation using the existing [source compiler](../common_source/ENCODED_MEMORY.md).
The result is **actively serviced quantum memory** on the declared carrier.
It is not an autonomous localized particle, rest energy, or a physical lifetime.

## Finite theorem and actual model map

Fix the classical core `c=(0,1)` and the native total-two carrier. The
screening packet supplies the two columns

\[
V|0\rangle=(|0,2\rangle+|1,3\rangle+|0,4\rangle+|1,6\rangle)/2,
\quad
V|1\rangle=(|0,3\rangle+|1,2\rangle+|0,6\rangle+|1,4\rangle)/2.
\]

Each source pair is a nonedge and has `q=1`. The prior
[all-seam flagged-code theorem](../common_source/ENCODED_MEMORY.md#optimal-native-carrier-code)
applies because every retained C/E raw branch is injective on these supports
and has the same support count `m_h` for both logical columns. This packet
checks the actual endpoint maps for all 720 raw slots: 90 C and 118 E slots
have positive code weight. Including the fair coin, `p_h=m_h/8`, and the
weights sum to one for every selected seam. The exact branch identity is

\[
V^\dagger K_h^\dagger K_hV=p_h I_2.
\]

The source encoder starts with an unknown qubit in addresses zero and one,
and eight blank code qubits. Two `Y(-pi/4)` source pulses make four equal
syndrome terms on wires 3 and 4. An emitted full-register Gray permutation
maps those terms to the eight declared physical code addresses. For a
positive branch, another emitted permutation sends its two disjoint native
image sets to the first `m_h` syndrome positions. Inverse `m_h`-term mixing
and the encoder give a complete 512-row recovery unitary `R_h`, with

\[
R_h K_h V=\sqrt{p_h}V.
\]

The support sizes are 1, 2, 3 and 4. The size-three mixer uses the existing
real-state source compiler; the others use zero, one or two `Y` pulses.
Every emitted encoder and inverse recovery Gray schedule is checked as an
integer action on all 512 addresses against independently enumerated native
images. The inherited route certificate checks all 30 physical seams on
all 78 native states and the source padding. Actual floating source words
then execute all 208 positive recovery branches on both logical columns,
retaining every source row. A two-attempt C/E path additionally executes the
native source instrument with all raw slots before selecting and recovering
the branch. Fresh source replay operator discrepancies must pass the fixed
`1e-8` numerical admission ceiling; backend-specific values are diagnostics.

For any arbitrary entangled reference and any finite number `t` of attempts,
induction on the branch identity preserves the joint logical/reference state
after each selected recovery. A classical policy may select later C/E laws
and seams from earlier complete flags. The flag probability at each step is
independent of the unknown logical amplitudes. This exact composition assumes
an immediate recovery after every native attempt and supplied, accurate
source services. It does not give a physical error rate or a calibrated time.

## Core transport and service cost

The native record determines a co-moving core label. For uniform seam choice,
it changes with probability `1/15` per attempted repair under either C or E;
it becomes a nonedge with probability `1/30`. The source recovery returns the
*fixed* code `V`. If the transported core label is retained, that code has
`q=1` probability `29/30` averaged over a uniform attempt. The protocol
therefore archives the old and new core labels and resets the classical label
to `(0,1)` after correction. It uses a fixed 15-bit status and label record
per attempt. This is an explicit metadata reanchoring service, not physical
motion or a free autonomous repair. The resulting code is sharp `q=1` with
respect to the reset label. No native quantum scratch register is added.

The complete resource ledger counts source pulses per invocation, even when
host calculation caches a matrix. Initial encoding costs **65,484** native
source pulses. With a uniformly selected native seam, the exact expected
native repair, routing and selected recovery cost is `37,499,681/60` pulses for C
and `36,870,437/60` for E per attempted repair. The respective maxima for
one positive branch are 799,631 and 783,290 pulses. The fixed archive uses
23 C or 24 E classical bits per attempt: five seam-selector bits, one coin,
two or three raw-record bits, and the 15-bit core-reset entry. The initial
core label uses seven bits. Source proper-code load/export and classical
control are apparatus services declared by the compiler; the core
reset has no assigned source pulse or physical work cost here.

The retained source composition admission ceiling is `1e-6` through 100
corrected attempts for either law, including the initial encoder. Every
replay recomputes the floating product diagnostic and checks this ceiling;
the receipt does not retain a historical backend norm as a portable fact. Its product
bound counts every emitted full-register macro and source primitive and takes
the maximum over orthogonal retained classical control blocks. This is not a
directed-rounding or hardware fault certificate. It concerns successful
processor blocks and the complete native raw control menu; it does not add
processor-failure density or classical-reset fault certification. Exact finite composition
holds for all finite attempt counts only under ideal supplied services.

## Particle gate

The active code improves one-port screening compared with the separate q=1 control, but occupies four exterior ports and has full-exterior trace
distance at least `4/9` from the transported reference. It needs recovery
and core reanchoring after every attempted native event. No autonomous
localized state/observable, invariant bond, statistics phase, measured rest
energy, gravitational stress, or observer-clock lifetime is identified.
CP-03 is open. The service does give a concrete component for testing
future candidate particle dynamics and cost/energy hypotheses.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/screened_recovery/test_protocol.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.screened_recovery.experiment --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.screened_recovery.experiment --diagnostics
```

The [campaign](campaign.json) pins its predecessor commit, input hashes and
exact science result. The [receipt](screened_recovery_results.json) binds
the imported source closure, this packet's code, tests and documentation,
and retains fixed rational numerical admission intervals. The verifier owns
fresh execution and compares the canonical receipt exactly. These intervals
summarize freshly passed numerical gates, not rigorous interval arithmetic.
The optional diagnostics command prints current backend observations without
claiming that historical residuals were reproduced. No physical data comparison
is added.
