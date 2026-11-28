# A separate probe on the same fresh arrival

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

Use the clock-bypassed local write on the four-fan 320-cell route. Keep native bridge rate 4 and wait 30.8125. The receiver at 1045 records the signed load and the face-0 word pair 334/365 in the positive chart, without distinguishing the two words. For a vacancy the word pair is 177/146 in the complementary chart. This is the same gate used by the [fresh reference witness](../../../fan_path/local_launch/WITNESS.md). The independent probe is enabled after that record, at the unchanged collar on faces 1044–1051–1054 and 1044–1054–1055. Read the probe at exposure 2 on 1054–1055. Enroll all sixteen added-head backgrounds and all four source and probe sign combinations. [Original collar protocol](../PROTOCOL.md)

## The arrival instrument and unchanged collar

The actual free propagation produces amplitudes indexed by the upstream word, the local receiving word and the Bell-half input. After selecting the terminal port and the two words, let $A_e$ be the resulting two-by-two matrix for upstream word $e$. An unread reference gives the unnormalized local density $\sum_e A_e A_e^\dagger/2$. The two reference Y readouts give $\sum_e A_e|\bar y_\pm\rangle\langle\bar y_\pm| A_e^\dagger/2$. Their sum is the unread-reference density. All these traces count every original launch. The replay checks the full two-qubit Pauli records of the reconstructed gate against the witness archive and verifies charge conjugacy with complemented inside and outside words. The opposite source chart uses the gauge that gives the same calibrated reference axes.

Restricted to the two archived code columns, the face-pair projector is a scaled isometry. It retains half the probability for every code input and preserves their relative amplitudes after normalization. The actual terminal detector projects directly on the face words. A preceding projection onto the entire archived code would remove additional arriving amplitudes. The executable records the effect of that extra operation separately and uses the direct face-gated density in every collar comparison.

The collar embeds both accepted words as existing source columns at hub 1045, with the same external word on each side of the local density. A face circulation reverses three heads while preserving the incident count at every vertex. The two columns therefore share the same charge and head-balance sector on every external background. In the removed-source control, an isolated full incoming bridge pulse transfers the source backward while changing only the bridge head. The pulse is a swap on every accepted microscopic pair. The source stays parked on the preceding fan, the bridge closes and the local density passed into the collar stays the same. The parked-source arm instead keeps the source at 1045 and retains its original contact diagonals. The hub degree cost is zero and the count of enabled source hops is fixed by incident head balance. This arm equals the removed-source read at every time because the remaining difference is a scalar.

When the connector ring is cut, the original fan and the outside probe terms act on separate factors. Each remaining move and rate depends only on its own factor's word and load. The external factor starts identically in the mobile-source and removed-source arms, including on a local fan correlated with the reference or upstream words. Its probe probabilities therefore agree exactly. The executable independently evolves all five arms and checks the predicted cut and parked equalities. It uses the archived Taylor enclosures for free flight and collar evolution; the reported allowances are numerical bounds for computed probabilities.

The collar probe and the reference analyzer share a causal gate on fresh copies. The entanglement witness includes a terminal face analyzer on its own cohort. The exterior probe changes the local evolution after the gate, so measuring its response establishes no claim that the probe preserves entanglement or that a single terminal measurement reveals both observables.

## Records at probe exposure 2

The ideal clock-bypassed stock has exact zero clock refusals and zero source-write refusals. The charge-and-face gate accepts **1.584370%** of original launches. Reference Y+ and Y− give **0.702405%** and **0.881965%** joint acceptance per original launch, respectively. Every other original launch is a charge-or-face refusal. Within the accepted cohort, the probe has distinct crossing and no-crossing outcomes.

At each fixed background let $\Delta_{s,p}$ be the mobile-minus-removed probe probability for source sign $s$ and probe sign $p$. The part odd in both signs is $(\Delta_{++}-\Delta_{+-}-\Delta_{-+}+\Delta_{--})/4$. The following means use all sixteen head backgrounds with equal weight. The ranges keep each head background fixed across the four signs before comparing the resulting contrast.

| Source analyzer report | Mean contrast per original launch | Mean contrast given gate | Range given gate |
| --- | ---: | ---: | ---: |
| Ignored | $-1.9680573\times10^{-5}$ | −0.124217% | −0.447254% to +0.002237% |
| Y+ | $-8.8862018\times10^{-6}$ | −0.126511% | −0.456182% to +0.002362% |
| Y− | $-1.0794371\times10^{-5}$ | −0.122390% | −0.440143% to +0.002138% |

The per-original contrast allowance is at most **$1.774\times10^{-10}$ in probability**. The conditional contrast allowance is at most **$2.541\times10^{-8}$ in probability given gate**, or $2.541\times10^{-6}$ percentage points. For every reference treatment, added-head word 3 gives the negative extremum and word 11 gives the positive one. In the matching-sign example at background 9, permitting source motion reduces crossing by about 1.0910 percentage points conditional on the gate with the reference ignored. A force law would require measurements of the response across distance and comparison with an independent coupling standard. This experiment measures an independent probe response to the captured material's local motion on the supplied head backgrounds. [Complete probability laws](results.json)

The upstream two-word arrival in the witness experiment has a nonseparable receiving-word/reference record. The subsequent probe comparisons resolve diagonal probe outcomes and, on the reference-resolved cohorts, a single reference axis. Those data leave post-probe reference coherence unmeasured.

Replay from the repository root.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.field_charge_records.electric_collar.fresh_arrival.measure
~~~
