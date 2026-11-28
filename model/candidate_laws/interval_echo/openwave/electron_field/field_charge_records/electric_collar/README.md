# Separate material probe on the asymmetric terminal fan

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

The terminal charge-and-code record from the original asymmetric three- or four-fan source triggers a local probe. The probe is a second excess or vacancy on a written seam with both endpoints outside the source fan. It reports where its own material unit goes. It never rereads the arriving unit's gateway. This is an addressed experiment on the written-ring extension. [Protocol](PROTOCOL.md) · [Original two-energy preparation](../../fan_path/README.md)

| Terminal hub | Connecting and outer faces | Probe seam |
| ---: | --- | --- |
| 1011 | 1009–1016–1017 and 1009–1013–1017 | 1009–1013 |
| 1045 | 1044–1051–1054 and 1044–1054–1055 | 1054–1055 |

Both endpoints of each probe are two intact-board hops from its hub. The executable retains all sixteen added-head words, both source signs, both probe signs, the two original inputs, and the two inherited source appointments. The probes move under an autonomous local generator after the apparatus opens. Both material units have separate routes and the probe has a supplied one-seam enclosure. The probe's bare contact is independently calibrated at written rate 4, including the eight-unit degree detuning on 1009–1013. [Complete records](results.json)

## The fixed-source null

On every sharp history admitted by the original terminal read, the inverse incoming bridge contact admits the returning source transfer. Its two-rail term has equal diagonals and rate 4. An isolated full contact of duration $\pi/8$ moves the source to the preceding fan and changes only the incoming bridge head. The captured internal head word stays the same. The removed unit is then parked away from the probe with the gateway shuttered. This realizes the neutral sham without reconstructing an unread field from simulator data.

Keep the source at the terminal hub in the other arm and retain the diagonal of its material contacts. The diagonal counts the allowed hops, four exchange units apiece. The number is the incident head count for an excess or the incident tail count for a vacancy. All ring flips keep that count fixed at the hub. A probe hop uses a distinct seam and also keeps it fixed. The two arms therefore differ only by a scalar from the held source. The two histories predict identical probe probabilities at every later local wait. The scalar also cancels in a neutral-subtracted probe energy comparison in this clamped experiment. This null holds for either sign and for an arbitrary retained coherence in the admitted local field words.

Allowing the source to move inside the fan changes the question. The first added face is the only term acting on both the old fan and the outer probe's written support. If it is removed, the remaining generator factors over those two supports even with both material contacts open. The same-field removed-source and mobile-source probe probabilities must then agree at every wait. Each trial starts with the exterior background prepared independently, so this implication also covers a captured field correlated with an unobserved upstream fan.

## What changes when the source moves

The most resolved banks change probe crossings by about two percentage points among captured trials. For example, with an excess source and excess probe at hub 1011, the probe crosses on **29.8654%** of captures with the source mobile and **31.8713%** after source removal. When the connecting ring is cut, it crosses on **26.2199%** of captures in both cases. The full launch denominator gives a mobile-minus-removed difference of about $-1.40031\times10^{-5}$ in this bank.

The sign test compares the four signed-source and signed-probe interventions on the same added-head word. For each sign setting let $\Delta_{q,p}$ denote the recorded probe probability per original launch with the source mobile minus its matched same-field removed-source probability. The component odd in both signs is $(\Delta_{+,+}-\Delta_{+,-}-\Delta_{-,+}+\Delta_{-,-})/4$. Using the same source appointment 27.0625, input 2 and probe exposure 2, it changes sign when the outer head preparation changes on either physical route.

| Fans | Added-head word giving negative response | Negative component per launch | Added-head word giving positive response | Positive component per launch | Error per component |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 5 | $-3.728419\times10^{-6}$ | 9 | $+1.443237\times10^{-6}$ | $<5.96\times10^{-10}$ |
| 4 | 13 | $-1.396017\times10^{-6}$ | 6 | $+3.130851\times10^{-6}$ | $<9.37\times10^{-10}$ |

Conjugating both material signs and every oriented written head preserves the complete generator and detector law. The comparison above holds the added word fixed across the four sign banks before changing that background. It resolves an independent material response to the captured source's motion through the collar. The static source null, head-sensitive response and equal two-hop endpoints leave an electric force law unresolved. The original signed material unit and the extra probe were both supplied, so an electric unit selection remains a separate test. OpenWave rows 1 and 20 stay on HOLD.

## Replay

The code encloses the two admitted energies and the corresponding finite eigenframe, reconstructs the inherited capture amplitudes for the ideal binary band pass, and propagates the collar in dyadic steps with Taylor remainder and binary64 roundoff bounds. The documented finite yes-and-reset operator differs from that pass in norm by less than $1.469\times10^{-21}$. Every later transfer or measurement is contractive, so this difference is within the carried allowance on original-launch events. The record gives the probability of terminal acceptance, the probability of subsequent probe transfer, and both refusal rates. Subtracting transfer from terminal acceptance gives terminal acceptance with the independent probe still at its starting port. [Finite source instrument](../../fan_rotation/COHERENT_FILTER.md)

Run from the repository root:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.field_charge_records.electric_collar.measure
~~~
