# A material read of the sourced Maxwell action

One owner moves an excess load around a triangular face. Its closed path and a stationary excess give identical charge censuses at all five appointments. The second owner reads the signed crossings. At the second appointment after the circuit, the two histories differ in magnetic flux by −3/2 on the source face and +1/2 on each of its three neighbours. A face two dual edges away first responds at the third appointment, with relative flux +1/8. On the closed path the complete prescribed action reads −93,352 in integer units, against 432,512 for the stationary source. All 180 signed first variations vanish on each of five source histories. [Frozen protocol](PROTOCOL.md) · [Source, signed words, variations and cohort counts](records.json.gz)

PR-66 is Demonstrated for a declared finite construction at step h = 1/2. The programs implement all four terms of the registered temporal action. The earlier theorem gives its stationarity equivalences for any nonzero real step on this complex. The material realization uses addressed contacts, a supplied clock index, exact control angles and finite digital words. Full OW-29 remains Open. A freely propagating transverse electric and magnetic material field with a common limiting speed c still requires measurement. The earlier freely released optical probes have different threshold velocities. [Finite action theorem](../../../../Lean/Screen/ScaledMaxwellStability.lean) · [Free-flight comparisons](../../../candidate_laws/interval_echo/openwave/row_29_maxwell/README.md) · [Published criterion](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#L135)

## Source and receiver

The measurement uses the complete 12-port, 30-seam, 20-face chart of the certified action. Chart vertices 0, 1 and 2 map to physical ports 0, 1 and 5 of owner zero. Port 7 supplies the compensating vacancy. An addressed native π exchange moves the excess with probability one. The held open path instead begins with a π/2 contact. Its original 1,048,576 starts give 523,972 stays and 524,604 crossings. Each branch then proceeds to the third vertex. The archive contains both conditional field histories and the original branch frequencies.

The second owner exchanges each inspected source port with sensing head 14 along native seams. On the outward and return path every intermediate load is restored. Echo reads against standards 15 (load zero) and 16 (load two) return 01, 00 and 10 for source loads zero, one and two. Inverse echo restores the sharp loads. All five full port censuses enter the density, ρ = load − 1. Before-and-after observations of the particular addressed seam enter J, with orientation chosen for Δρ = −h dᵀJ. A loop gives J₀ = −2 Cᵀe₀, where e₀ selects face zero. Reversal negates J, while an idle gives zero. The independent endpoint census is identical in all three. The program audits the relevant columns of the full 99-state local law and every transposition route. [Source program](material.py)

On owner one the dimer (12,13) uses alternatives |20⟩ and |02⟩. The native seam restricts to 32.5 I − X/2. A private degree-five reference on port 13 supplies diag(0,32). Two quarter phases around a half-exchange give the Hadamard recombiner up to an overall phase. After recombination a zero standard at 17 and an echo on its seam to 13 distinguish the two occupations. All source routes avoid the dimer and this standard. Each input bit enables its own pulse. A signed product enables pulses for the active bit tuples in its factors. At digit j the combined relative phase is exp(2πi z/2ʲ⁺¹), where z is the signed polynomial in the input words. The earlier j returned bits set an inverse phase, leaving +1 or −1. A Hadamard and occupation read give the next bit. Induction proves the full integer decoder for every input whose absolute term budget is below 2³⁹. Of 29,021 signed-word reads the largest absolute budget is 28,688,016 and the largest return magnitude is 10,713,232. The serial controller stores the returned words and schedules the next contact. [Native phase and digit reader](rails.py)

## Window and stationarity

The recurrences are Aₙ₊₁ = Aₙ − h Eₙ, Bₙ₊₁ = Bₙ − h C Eₙ and Eₙ₊₁ = Eₙ + h(CᵀBₙ₊₁ − Jₙ). The first four updates use φ = 0. One more potential endpoint closes the five-interval window. Voltages, electric and magnetic fields, charges and currents are archived in units 1/L, with L = 512. Writing a star for an integer input and e* = −2ΔA* − dφ*, the action meter reads

```text
F = 4 L² S
  = Σₙ [ ||e*ₙ||² − ||CA*ₙ₊₁||² + 2〈J*ₙ,A*ₙ₊₁〉 + 2〈ρ*ₙ,φ*ₙ〉 ] .
```

The action reader starts from the potential and material source words. It constructs each incident local square, both source couplings and their signed sum through the dimer. An independent integer calculation recomputes the full window from its definition. Every interior A coordinate has trials with displacements ±1/L, as does every φ coordinate. Unchanged terms cancel exactly. For a local square w v², the difference is 4wvδ, so these symmetric readings determine the directional derivatives exactly for all inputs in the declared word domain. The companion Lean file checks the per-step normalization against the existing action, including the forward magnetic and current time placements, and the square variation identity. [Action program](action.py) · [Kernel-checked normalization](ActionBridge.lean) · [Finite variational theorem](../../../../Lean/Screen/ScaledMaxwellStability.lean)

| Source history | F | Nonzero interior A variations, out of 120 | Nonzero φ variations, out of 60 |
| --- | ---: | ---: | ---: |
| Idle | 432,512 | 0 | 0 |
| Closed loop | −93,352 | 0 | 0 |
| Reversed loop | 150,632 | 0 | 0 |
| Open, first stay | 367,456 | 0 | 0 |
| Open, first crossing | 227,472 | 0 | 0 |
| Same crossing, time-dependent endpoint-fixed gauge | 227,472 | 0 | 0 |

Gauss, continuity, Faraday, Ampère, source work and the defining E/B identities have zero integer residual on all five histories and the gauge trial. Changing the held potential at one time gives seven nonzero A and four nonzero φ variations. Removing the kinetic term gives 65 and 10. Removing the magnetic term gives 65 and zero. Reversing the seam-current coupling gives two and zero. Removing the port-charge term gives zero and ten. Reversing the current during evolution also violates Gauss by two charge units in the held crossing. Dropping magnetic evolution violates Ampère and the work balance. The compact gauge changes φ across the time window, so its charge coupling is exercised outside temporal gauge. [Acquisition and control comparisons](measure.py)

## Flux records

For the loop minus idle, ΔB₀ and ΔB₁ vanish. The first nonzero face vector is ΔB₂ = −CCᵀe₀/2. With the source now stationary, the operator gives ΔB₃ = (2I − CCᵀ/4)ΔB₂. The certificate checks both predictions independently for all twenty faces. Current sign, timing and magnetic propagation thus have distinguishable consequences even with identical port censuses.

Each pair below gives the bright X and Y counts on independent cohorts of 1,048,576 original starts per quadrature. The appended two-rail read uses returned digits from the loop and inverse pulses from idle at the same face, with π/2 phase per unit flux. The changed contact halves the private strength during the data dwell only. The positive durations, preparation and analyzer are held fixed.

| Step | Face | Distance | Native bright X, Y | Changed contact bright X, Y |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 0 | 0 | 153,150 · 153,529 | 323,256 · 1,008,539 |
| 2 | 1 | 1 | 894,646 · 894,640 | 1,008,615 · 725,078 |
| 2 | 3 | 2 | 1,048,576 · 523,615 | 1,048,576 · 524,102 |
| 3 | 3 | 2 | 1,038,363 · 627,804 | 2,561 · 473,326 |

The 720 archived settings cover the first face of each of the six dual distances, all five times, the four non-idle sources, both receiver strengths and all three field schedules. Every quadrature frequency and the half-exchange split satisfy a simultaneous 99% Hoeffding half-width 0.002449. The largest absolute deviation across the 1,440 quadrature frequencies is 0.001851 after rounding upward. The archive retains the binary64 sampling probabilities to reproduce every original seed draw across floating-point libraries. Each replay recomputes and checks the probabilities from the physical pulses within 2×10⁻¹³. The addressed order and h supply time; the experiment compares simulated observer records in that schedule. [Deposited quadrature and source data](records.json.gz)

From the repository root, with one numerical worker:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
nice -n 15 .venv/bin/python -m model.experiments.interval_echo_gauge.temporal_action.measure
nice -n 15 .venv/bin/python -m pytest -q model/experiments/interval_echo_gauge/temporal_action
(cd Lean && nice -n 15 lake env lean -j1 ../model/experiments/interval_echo_gauge/temporal_action/ActionBridge.lean)
```
