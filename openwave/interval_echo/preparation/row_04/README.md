# Row 4 · Free persistence of a spectroscopy excitation

**Official mark ⚠️ Partial, 25 November 2026.** The written test is “the particle state persists in free evolution.” M10.3 selected contrast strength 8 with pair persistence in view. M10.4 proved persistence at every strength on its eleven-write material. At strength ½ and exposure 37 the same pair and detector gave 0.995945 on the written board and 0.440929 on blank memories. At strength 8 both boards clear the historical 0.94 gate. Review accepted the computation and retained partial credit because the eleven writes were selected for persistence. [M10.3 decision](https://github.com/openwave-labs/openwave/pull/746) · [M10.4 decision](https://github.com/openwave-labs/openwave/pull/754) · [Official cell](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#L438)

**Scope**. One load dressed by native arrow records on the 20-cell Interval Echo board. The detector tests localization of that excitation at the first independently appointed read after release. Each appointment uses a fresh preparation. The claim is persistence, and the smallest object is one load with its coherent arrow history. Its ports remain connected to other readings through the native law.

## Finite implication from the launch record

Let $D=(|a\rangle-|b\rangle)/\sqrt2$ for either of the two enrolled signed histories, $P_D=|D\rangle\langle D|$, and $P_T$ the terminal detector effect for the two physical tips. The exact signed local reader at launch has hit effect $P_D/8$. Let $p$ be its hit probability per original start, including the commissioning refusals. The flight generator satisfies $HD=(\kappa/2)D$ after removal of the common seam energy. The proof checks the two columns over the complete reachable configuration space, including exits. The detector inequality $P_T\succeq P_D$ and preservation of $P_D$ give

$$
\Pr(\text{two-tip terminal reading at }t)\ge
\operatorname{tr}(P_D\rho_{\rm release}) = 8p .
$$

The source observer retains the prefix bit, all reverse-query outcomes and the final center read. Another observer compares its local arrival records on the two tips. The trace of the released ensemble is the commissioning acceptance probability. Both observers use every original start as their denominator. A lower bound on $p$ closes the launch premise on this construction. Subtract the read uncertainty before multiplying by eight. A direct persistence read tests the inference on independent preparations. [Native selective effects](../../../../model/candidate_laws/interval_echo/openwave/row_04_spectral/phase_reader.py) · [Full-column and flight checks](../../../../model/candidate_laws/interval_echo/openwave/row_04_spectral/experiment.py)

## Choice ledger

| Choice and how it was fixed | Dependence and disposition |
| --- | --- |
| Oriented material and port-zero launch | The public donor `44f7680df` specified ascending, reverse-local and descending orders to compare reachable transport histories. Ascending permitted complete 32-configuration spectroscopy. The alternatives do not even reach the two enrolled histories. The transport question and tractable instrument fixed this deciding material. [Original arrow task](https://github.com/FloatingPragma/observer-patch-holography/blob/44f7680df55b4e68b89c301896df3e10237f046d/model/candidate_laws/interval_echo/openwave/row_30_arrows/PROTOCOL.md) |
| Both signed analyzer columns and their source | The donor applies a 90% joint wavelength and frequency criterion to all 32 spatial columns and 31 distinct energy levels. Both matching columns, 14 and 15, enter the release. The replay uses the donor's original pulse arrays. The whole clock-transfer study qualifies all 32 words; it uses a contact-transfer meter independent of the persistence read. The separate staggered reading and the failed donor rest and speed tests are disclosed. [Original spectroscopy task](https://github.com/FloatingPragma/observer-patch-holography/blob/44f7680df55b4e68b89c301896df3e10237f046d/model/candidate_laws/interval_echo/openwave/row_30_selective/PROTOCOL.md) · [Clock transfer](../../../../model/candidate_laws/interval_echo/openwave/row_04_spectral/CLOCK_TEST.md) |
| Coupling and flight | The exact flight implication applies at every real strength. The carried preparation works on the sixteen tested nonzero strengths with its clock calibration. Zero has no transfer clock and is inconclusive for this launch. Fixed historical waits fail at detuned strength. No coupling is selected by a persistence reading. [Complete release terms and outcomes](../../../../model/candidate_laws/interval_echo/openwave/row_04_spectral/DONOR_RELEASE.md) |
| Commissioning and terminal detector | Addressed isolation, matched local fuel, pulse service and retry budget are supplied. The 631 fair commissioning operations finish with probability $(1-2^{-N})^{631}$ at cap $N$. At cap 8 the exhaustion probability is 0.91539. At every cap of 16 or more the commissioning loss is below 0.00959, so the computed conditional bound clears the donor's 0.98 benchmark with failures counted. The donor's signed histories fix the detector tips. The changed-law return can miss them. The separate five-query source was developed from those known rays and is excluded from the tested launches. [Commissioning operations](../../../../model/candidate_laws/interval_echo/openwave/row_04_spectral/source.py) |

The comparison addresses the written-material route in the maintainer decision. The donor predates the M10.4 layout. The public stability criterion and historical pair work were known when the donor was published. The donor's spectroscopy result was known when it was selected for this transfer. The submitted claim would depend on the donor's own selection record. The supplied commissioning procedure changes throughput and its accepted branch reinstates the donor's material before the original pulse word starts. [Maintainer's independent-source route](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/openwave/xperiments/m10_oph/research/findings/m10_4_maintainer_review.md#5-the-two-decisions)

## Result and replay

The minimum computed all-start bound on 32 resolved launches exceeds **0.99999950**. At exposure 37 the largest dephased return is 0.877116, rounded up, and the largest rate-matched random-swap return is 0.427295. The exact inequality covers the read at 150 and beyond it. The simulated witness draws give a minimum simultaneous 99% lower bound of 0.9941931. This bound refers to sampling the computed probabilities, and the exact implication carries the time-domain conclusion. [Complete launch and terminal receipt](../../../../model/candidate_laws/interval_echo/openwave/row_04_spectral/donor_release.json)

From the OPH root, run:

```sh
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_spectral.donor_release
```

**HOLD.** The experiment replaces the adjudicated object, memory and reader. A candidate nomination, review of the independent provenance of this choice on the exact outgoing version, registered in-platform replay and release clearance are required. [Version and checks](../REVIEW.md)
