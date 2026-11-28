# What the two resolved energy lines imply for a turn comparison

IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.

The admitted band on each asymmetric three- and four-fan route of the 320-cell board contains exactly two nondegenerate energy lines for either charge sign. The finite coarse clock passes both coherently. Independent comparisons of coherent and fine-energy-dephased receiver clicks resolve beat ratios of 0.58551 and 0.36144 against the neutral-ring standard. These count-based readings establish a resolvable splitting. The spectral count and disjoint residual intervals identify two one-dimensional lines for the argument below. [Finite clock and rank](COHERENT_FILTER.md) · [Independent beat records](../fan_recoil/BEAT.md)

## The fixed-guide order implication

Test each proposed physical turn on copies of each admitted energy input. Require the same band pass and the same fine energy on every such trial. In the coherent order comparison the fine meter is absent. Each turn must act as a fixed unitary on the guide band, with external control and ancillary degrees of freedom returned to the same ready record by each individual action. The two arms use the same phase reference. The claim concerns both orders of these same two operations.

Write $e_0,e_1$ for the two calibrated energy directions and $A,B$ for the turns. The separate energy outcomes are orthogonal and each retained line has rank one. A repeat of the same energy with probability one implies $Ae_i=a_ie_i$ and $Be_i=b_ie_i$ for phases $a_i,b_i$. Each order therefore gives the same amplitude $a_ib_i e_i$ on both lines. Any accepted superposition gives the same outgoing record in both arms. Equivalently $AB(BA)^{-1}=+I$ on the band. An ideal calibrated order interferometer puts the entire admitted overlap in the positive output. The spin-½ operators for orthogonal π turns, $-i\sigma_x$ and $-i\sigma_y$, anticommute. They cannot both act on this band while satisfying the separately tested energy conditions. [Checked linear-map implication](EnergyOrder.lean)

An ancilla that stores a nontrivial operation can change the order fringe while leaving both guide energy labels fixed. For example, independent Pauli actions on a retained ancilla have negative order. The test requires a fixed coherent guide operation and matched control resets. The hypothesis concerns energy-preserving turns in this undriven guide. A turn that transfers energy to or from a control, an action that changes the environment, and exchange of two cores require separate measurement protocols.

## Full-guide pulse check

Compare the addressed controls used in the [capture-phase analyzer](README.md) on both complete admitted energy inputs. The native isolated first-ring pulse of duration $\pi$ is an exact permutation of enabled terminal head words. The two supplied outer-rim head phases are an exact diagonal sign operation. Use the same receiving fan and keep every word of every fan. For each pulse and incoming line, compare the lower line, upper line and refusal by the same band detector. The check forms the pulse overlaps with integer trial vectors and widens every amplitude by the certified eigenvector error. Charge conjugation intertwines the entire generator and both pulse actions.

| Guide and full terminal action | Refusal lower bound for the lower energy input | Refusal lower bound for the upper energy input |
| --- | ---: | ---: |
| Three fans, native first-ring pulse | 0.11607% | 0.28054% |
| Four fans, native first-ring pulse | 0.10836% | 0.13110% |
| Three fans, supplied outer-head phases | 99.94067% | 99.91444% |
| Four fans, supplied outer-head phases | 99.96162% | 99.96399% |

Both actions also have nonzero transition probability to the other retained line. The displayed refusal bounds use the exact band projector. Replacing it by the replayed finite yes-and-reset meter changes the refusal probability by at most $2\eta < 2.937\times10^{-21}$. Thus both pulses fail the energy-preservation test on the whole band. The negative order phase of the capture analyzer still demonstrates the detector's ability to distinguish phase signs.

The replay rebuilds the 320-cell fixture, verifies the inertia at both cuts and encloses the two eigenvectors before computing the outcome bounds. [Executable](energy_order.py) · [Complete receipt](energy_order.json)

Run the Python receipt from the repository root and the Lean check from the repository's Lean directory.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/fan_rotation/energy_order.py
~~~

~~~sh
lake env lean ../model/candidate_laws/interval_echo/openwave/electron_field/fan_rotation/EnergyOrder.lean
~~~
