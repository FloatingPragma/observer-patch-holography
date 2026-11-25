# One thermal reference for a measured dimer, its repair and its refining family

## Protocol

Use the one-double protected dimer sector and the existing energy decoder. Supply equal address cohorts on the first 20-owner board. At each address prepare both stationary modes in separate cohorts. Decode the incoming mode into rail A or B, run a half exchange on AB and read B against vacant C with the registered balanced-pair pointer, then undo the pointer. The proposal is now sharp. Attempt the same native BC filter that independently produced the one-appointment reference. Read C against vacant D and undo. A dark read accepts; undo the energy decoder and publish the retained mode and the unchanged address. A bright read places the mobile on C. A whole BC exchange restores it to B, ready to try another AB proposal on the same owner. The finite acquisition stops after four unsuccessful attempts, publishing timeout. All completed and failed starts count.

Freeze 2,048 starts per old mode at each of the 20 owners. Keep the accepting depth and sign in every bank, as well as timeouts. Test the complete geometric prediction, the source's independently filtered frequency, input independence, and its failure when the filter or proposal shutter is parked. Verify that the corrected bright branch really restores B and that the accepted state is the indicated stationary mode. Retain the random-time rule and the finite cap as different measurement contracts. Use only the registered finite material matrices with one numerical-library thread. The units of the clock, gap and multiplier are internal.

## What each address returns

Decoding sends the energy-32 ground mode to A and the energy-33 excited mode to B, up to phases. A half AB contact puts half the probability on each; the B/C binary pointer resolves the proposal. The BC filter transmits a proposed A with probability one and a proposed B with probability one half. The three disjoint outcomes of one attempt therefore have weights 1/2 (ground), 1/4 (excited) and 1/4 (leak), for **both** old signs. After a leak the mobile is sharp on C; the full BC contact puts it on B with certainty. The next attempt has the same three probabilities. After a dark read the inverse decoder returns an energy eigenstate. A fresh local decoding and balanced pointer read identifies the returned sign. The full embedded source matrix certifies each pointer, including undo. Measurements and correction pulses stay on the original owner. [Native transaction](surface.py)

For accepting attempt t≥1, the probabilities of terminal ground and excited reports are respectively (1/4)^(t−1)/2 and (1/4)^(t−1)/4. The genuine timeout probability after L attempts is 4^(−L). The indefinitely callable service terminates almost surely with mean 4/3 attempts. The terminal probabilities are 2/3 and 1/3, independent of both the old sign and the number of attempts to acceptance. The finite-cap receipt has an explicit third terminal category. Reporting timeouts as either sign would change the measured kernel.

| Source or repair arm | Ground | Excited | Unresolved or refused | Excited among admitted |
| --- | ---: | ---: | ---: | ---: |
| Independent one-appointment source | 20,462 | 10,181 | 10,317 filter refusals | 0.33225 |
| Repair starting from ground | 27,164 | 13,650 | 146 | 0.33444 |
| Repair starting from excited | 27,209 | 13,596 | 155 | 0.33319 |

All 20 owners supplied both old signs. Across 81,920 repair starts the 301 four-attempt timeouts compare with 320 expected. Parking BC instead gives 40,766 excited outputs in 81,920 starts with no timeout. Parking the AB proposal retains the old sign on every acceptance and leaves 2,591 of the previously excited starts unresolved. The one-appointment source and the repair use separate launch cohorts and the same native energy decoder. The per-owner, per-depth tests use simultaneous 99% counting allowances. [Complete repair counts](surface.json) · [Source counts](records.json) · [Checks](test_surface.py)

## Identifying the finite four-law operands

At a stage with K owners let Ω be the 2K measured pairs (address a, energy sign e), e=0 or 1. Equal address enrollment and one filter appointment give R*(a,0)=2/(3K) and R*(a,1)=1/(3K). Thus R* is normalized and faithful. The complete repaired datum is the address. On this surface each native returning repair has law R*(a,e) divided by the total reference mass 1/K of its address fibre. These are the reference, positivity, normalization and visible map of the finite common-reference repair theorem. Arbitrary initial address populations are preserved by a repair. Uniform source enrollment provides the address weights used to prepare the Gibbs family.

The independent flight and detector measurements fix E(a,e)=32+e in native energy units. Put β*=log 2 and log Z*=log(3K/2)−32 log 2. Direct substitution gives −log R*(a,e)=β*E(a,e)+log Z* on **every** measured pair. This supplies all scalar fields and equalities in the finite calibration record in native units. An empirical relation between these clock/energy units and a laboratory standard would supply additional physical information.

For any multiplier β, the canonical population Pβ is uniform over addresses at a fixed sign, with excited fraction p=1/(1+exp β). Its log density relative to R* is an affine function of the independently measured energy. Every competing distribution V on the complete address/sign space with that same mean energy consequently satisfies

\[
D(V\Vert R_*)=D(V\Vert P_\beta)+D(P_\beta\Vert R_*).
\]

This proves unique state selection at that moment and identifies the reference used by the transition on the same records. The proof allows competitors with arbitrary address/sign correlations. It applies to the discrete β values reached by consecutive half filters as well as to mathematically defined Gibbs laws.

Choose the 20-owner stage as the distinguished member. Its forty dimer readings are exactly Ω; the identification sends each reading to itself and preserves its energy. Together with the growing guarded dimer family and its coarse address maps, this gives the native-unit operands of the cofinal four-law theorem on one declared surface. Physical temperature, duration and energy in laboratory units require their own calibration. The construction's finite-time refusal and the larger all-open spectrum retain their measured meanings. [Repair and calibration fields](../../../../../Lean/Thermodynamics/FourLawAdequacySurface.lean) · [Cofinal attachment fields](../../../../../Lean/Thermodynamics/CofinalSpectralTailFamily.lean)

Replay from the RER root.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B -W error -m model.experiments.interval_echo_consensus.refinement_cooling.dimer_band.surface
.venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/refinement_cooling/dimer_band/test_surface.py
~~~
