# Ordinary energy of the released excitation during its drive

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

The same lowest mode that contributes to the distal absorption after the parent falls to ground has a light **ordinary contact-energy average across its physical pulse period**. At the original calendar this average is 0.066486 on 23 dimers and 0.033796 on 45. The independently calibrated light cutoff lies between 0.097012 and 0.097017. The common half-duty control gives average energies above 0.49. [Protocol](PROTOCOL.md) · [Complete records](results.json) · [Original release and absorption](../README.md)

### Which excitation is calibrated

The original separate parent starts excited. The daughter mode is measured on sister cohorts after the same local parent-ground read, using native holds and their inverse to read and restore its lowest one-mode branch. The original receiver tag combines the distal local absorption with a terminal vacuum. The orientation-phase operator is linear in normal modes up to conserved parity, so a successful absorption to the vacuum has a one-mode contribution. All but the lowest one-mode frequency sit above the light cutoff. On 23 and 45 dimers there were 28,871 and 8,608 successful joint tags per $2^{26}$ simulated originals, respectively. The computed bright-only leakage upper bounds are $1.415\times10^{-5}$ and $5.082\times10^{-8}$. If $p$ is the joint tag probability and $epsilon$ that leakage bound, the lowest-mode mass at the absorber is at least $(\sqrt p-\sqrt\epsilon)^2$ whenever $p>\epsilon$. This follows by subtracting the largest allowed bright amplitude from the accepted amplitude and using the detector's contraction bound. The count and leakage records therefore identify a contribution of the same low mode to the parent-loss event, even with bright interference. [Primary detector and leakage theorem](../PROOF.md)

### Two contact layers, one independent phase test

Keep the row 20 route and schedule. The local reference is the same coupled vacuum. From fresh admitted specimens, observers read each active contact using its checked native spectrum at the start of its layer. Its energy expectation remains constant throughout that layer. The first link has half the other links' dwell, so it contributes half their time weight. Independently scale both existing dwells by $\lambda$ with their ratio fixed. If $A$ and $B$ are the integrated generators of the two layers and $U(\lambda)=e^{-i\lambda B}e^{-i\lambda A}$, the derivative of a continuous eigenphase of $U$ equals the expectation of $A$ at the start and $B$ at the opening of its layer, using the convention $U|\psi\rangle=e^{-i\theta}|\psi\rangle$. Taking the mode-minus-vacuum derivative and dividing by the original period gives the same cycle-averaged energy as the independent contact readings. The program verifies this identity numerically on the exact signed aperture and checks its gap to neighboring modes.

| Guide and drive | Ordinary average, simultaneous 99% | Lowest stroboscopic frequency |
| --- | ---: | ---: |
| 23, original | $[0.06640310,0.06658988]$ | 0.06680803 |
| 45, original | $[0.03360267,0.03395899]$ | 0.03396663 |
| 23, half common link duty | $[0.50029911,0.50111820]$ | 0.50064667 |
| 45, half common link duty | $[0.49331187,0.49605720]$ | 0.49469179 |

The intervals include the propagated numerical allowances and simultaneous uncertainty from seeded simulated contact draws. Each contact bank represents $2^{42}$ original launches and includes both failed heralds. Lowest-mode enrollment succeeds on about 0.487% and 0.299% under the original calendar. Shuttering the parent link makes the parent-ground outcome impossible. In the inherited full detector, quarter, half and full first-link openings all exceed the bright-only bound. The medium duty and the parent contact play different roles.

The low mean includes the time-weighted sum of much larger layer contrasts. The first switch changes the ordinary daughter-minus-vacuum energy by about −20.293 and −20.790. An eigenmode that returns every cycle has zero net switch work relative to the returning vacuum over a full cycle. Energy records made during actual source depletion and absorption need their own complete work comparison. The measured mean gives a positive light-energy calibration for the powered material mode. It supplies no ordinary energy line on the unpowered host. [Ordinary release from the same enrolled branch](../README.md#removing-the-drive-and-reading-ordinary-energy)

### What fixes the charge

The original flight retains the amount-two guard on each 02/20 dimer. A full receiver read gives its baseline charge on all admitted native records. The changed test opens the last link's exchange in just one link layer on the same daughter, keeping every other address, dwell and guard. The other port of the receiving dimer has no link in that layer. If the two active ports have unlike extreme loads, their actual blank native pulse exchanges them with probability $\sin^2(\delta_B/2)$. The terminal dimer then reads either −2 or +2 relative to its filled standard. The original guarded probability is zero, while the changed-law probabilities of a nonzero terminal reading are $8.8950\times10^{-7}$ and $9.0533\times10^{-7}$ per admitted daughter. The neighboring dimer receives the compensating charge. This tests local neutrality in the held domain, with the total charge fixed in both arms.

With an additive conserved charge, the whole daughter's net charge follows by subtracting the unchanged parent charge and the charge in the rest of the material from the original total. Charge transfer by the detector must be included in this bookkeeping. A simultaneous neutral reference read on the complement before absorption would allow local charge fluctuations within the daughter. The original terminal vacuum check follows a long absorption aperture and resolves the guarded guide. The charged replacements test the pointer's calibration on separate stocks. A release with admissible charged alternatives requires the corresponding joined complement record or a detector-equivalence proof.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_09_neutral.critical_emission.cycle_energy.measure
~~~
