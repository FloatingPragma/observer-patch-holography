# Reading a neutral rest phase

OpenWave records OW-03 as untested at [head `0df73dc6`](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L437). The full local criterion is partial. On the unchanged Interval Echo seam a coherent superposition of filled endpoints and a symmetric excess–vacancy pair has a relative frequency of 32 exchange units. The primitive echo, free phase contacts and sharp memory/load pointers cannot distinguish it from its phase-erased mixture. A supplied block-mixing Ramsey analyzer resolves the phase. This limits that neutral-creation route to a rest clock. [Protocol](PROTOCOL.md) · [Proof and scope](PROOF.md) · [Complete records](results.json)

## Which additional operation is needed?

The three filled, zero-contrast seam words form a reducing subspace of every primitive operation. Every finite adaptive program built from them has zero probability to create a neutral pair from that subspace. It also misses coherence between that subspace and its complement. This extends to uniform filled federations with blank or zero-contrast seam memories. Nonuniform reference matter and additional coherent source/analyzer operations require their own tests. The deposited 17-appointment experiment retains all 99 joint memory/load outcomes and outcome-dependent second contacts. Pure and dephased native probabilities agree within 10⁻¹⁵. The same pair census detects probability 1/2 under a changed, charge/swap-preserving local law with a filled-to-pair mixing term, against zero natively. The ideal extra analyzer resolves `(1 ± cos(32t))/2`; it supplies the operation absent from this grammar.

## What an energy reference can fix

For an entire apparatus in one total-load sector, adding λ times its total charge to the Hamiltonian changes no complete transcript probability. The proof covers balancing references, arbitrary charge-conserving instruments, untouched quantum references and classical feedback with outcome-dependent waits. Every coherent arm must share its elapsed physical interval, and the same shift applies to the reference and controller. Coherently applying the shifted sample Hamiltonian against literal idle is a conditional drive, which the deposited sensitivity control distinguishes by 0.707106781 in its plus-port probability.

The filled-host assignments `(vacancy, excess) = (35,67)` and `(51,51)` give the same allowed records and the same neutral formation energy of 102. Any estimator from these records has worst expected absolute error at least eight exchange units on those two conventions. A physical inertial-mass or relativistic reference could select an additional interpretation. The neutral sum alone cannot select either charged formation energy. These are identifiability boundaries, with no exclusion of a clock on a different declared preparation.

## Replay

From the RER root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_03_clock.neutral_reference.measure
python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_03_clock/neutral_reference
```

The replay reads `results.json` and fails on a difference; `--output PATH` explicitly creates a separate first record. Nine tests check the exact local echo block, complete adaptive readout, supplied analyzer, uniform-versus-conditional shift, balancing reference, minimax bound invalid-probability rejection, the public replay command and a charge/swap-preserving changed law. Every appointment has 262,144 simulated enrollments. The raw numerical probabilities are retained. A checked normalization residual below 2×10⁻¹² is assigned to an explicit complementary bin for sampling and recorded as an L1 correction bound; no negative probability is clipped and no probability vector is rescaled. No physical clock calibration or laboratory energy scale is established by this equality-at-every-reading result.
