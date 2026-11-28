# Coherent exchange of both physical face records at the 1045 junction

IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.

The two observers at 1045 and 1028 can exchange both physical face qubits and both material origins with native isolated pulses. Both receivers keep their original physical analyzer addresses and compare against independent untouched source references. A fixed **31-pulse** schedule gives exactly the marker swap on all four ordered sharp inputs. Every load remains in a different cell throughout. The equal-exposure native interferometer reads **positive** exchange phase. The two input qubits survive through the apparatus, although they are simultaneously readable on faces attached to their corresponding material origins only at the beginning and end. [Protocol](PROTOCOL.md) · [Microscopic and observer records](results.json)

## Which records cross

Keep the physical 1045 face-0 gate after clock-bypassed flight to wait 30.8125, the completed central field and the independent Bell-tagged face at 1028. Support the six written links of the original Y, all six bonds of its two physical marker triangles, the two native triangle flips and no other evolving terms. Each addressed material term has the registered rate 4 and each ring term rate $1/2$. The respective isolated exposures are $\pi/8$ and $\pi$. Other contacts and the material position costs are shuttered during every pulse. [Original junction and preparation](../README.md)

The sharp two-origin description has 240 reachable microscopic records. A simultaneous four-input search compiled the schedule with 29 material pulses and two ring pulses. The executable checks every intermediate port occupation and complete head word for both charge signs. The result returns the same two occupied addresses with the material origins exchanged. It also interchanges the two ordered face directions, with all spectator and bridge head bits restored. The position of the spectator load is always in a cell different from both options of an enabled move, so the fractional-pulse superposition also lies within the distinct-cell sector. Omitting either one of the two addressed ring pulses from this fixed schedule leaves only two of its four sharp inputs with the intended endpoint and word.

The same source can be tested with the original twelve-pulse material exchange, which keeps the field references parked. Another compiled control swaps only the face codes in 35 pulses while returning the material origins, then applies that material exchange. This second compilation gives the same terminal amplitude as the 31-pulse route and includes same-cell occupation. Thus the endpoint exchange result survives a comparison of schedules while the paths and the available intermediate records differ.

The receivers keep both local orientation analyzers. Using the inherited two-qubit arrival moments and the separately supplied second Bell pair gives these probabilities, conditioned on the original face gate.

| Joint Bell pass | Material-only arm | Both markers and material exchanged |
| --- | ---: | ---: |
| Original references | 20.89273% | 5.22318% |
| Crossed references | 5.22318% | 20.89273% |

Each of the two formerly local reference correlations is entangled with the **opposite** receiver under the full transfer. The 1045 source reference gives the inherited negative partial-transpose eigenvalue near $-4.2172368\times10^{-4}$ per original launch at the 1028 analyzer. The separately prepared 1028 reference gives $-7.9218498\times10^{-3}$ per launch at 1045. The two original pairings after marker transfer are product marginals. The receipt retains their full Pauli moments, and a local material/endpoint read includes all earlier refusals. Phase values depend on the comparison of the same physical final records and the untouched references carry no gate actuation.

## Equal-pulse native phase comparison

The 31-pulse route has five primitive terms with odd counts. For each such term, prepare a configuration in which that term acts as the identity on all four sharp inputs. The three odd material terms are inside the 1028 triangle, and their common control moves both loads out of that triangle using the tree. For each of the two ring terms, its two local spokes send whichever sharp orientation is present to an allowed face vertex. Both resulting triangles are noncyclic. Apply the disabled term once, then undo the preparation with the same native full pulses in reverse order. The concatenation of these five returns is an exact identity on all four words and has the same pulse-count parities as the crossed arm.

Balance the remaining even differences by adjacent pairs of the same full pulse. Both final phase arms use **35 material pulses and two ring pulses**, with the exact same count for each addressed generator. The total exposure is $51\pi/8$ per arm. Both schedules preserve distinct occupied cells and all port Gauss records. No inverse-time pulse or controlled negative rate is needed for the reference.

The crossed arm is $+S$, where $S$ interchanges the two face words, and the uncrossed arm is $I$. On the unresolved independent-reference ensemble, the second Bell pair has local reduced density $I/2$. Therefore the swap overlap after tracing the references is $p/2$, with $p$ the original-launch face acceptance. The phase meter gives bright fraction $3p/4$ and dark fraction $p/4$.

| Outcome per original launch | Native pulses | Supplied opposite phase |
| --- | ---: | ---: |
| Clock, ideal source or ideal write refusal | 0 | 0 |
| Face arrival refused | 98.41563005% | 98.41563005% |
| Arrived but endpoint or field refused | 0 | 0 |
| Joint bright | 1.18827746% | 0.39609249% |
| Joint dark | 0.39609249% | 1.18827746% |

The absolute numerical allowance per event is $5.5413\times10^{-10}$. Independent equal-word controls, which resolve the first arriving marker and prepare the second with the same word, give only bright for the native experiment. The bright rates per original launch are approximately 0.76709009% and 0.81727986% in the two settings; the record retains nonmatching-word refusals. Applying a minus in the crossed arm puts those same fractions in the dark output.

## The sign supported by this class of native controls

Each complete native pulse is a nonnegative permutation of the sharp head and material configurations. It follows by composition that any such fixed program taking the four sharp input words to their exact crossed counterparts has amplitude $+S$ on the complete four-dimensional code in the calibrated sharp-word basis. Intermediate routing and the number of taken steps can depend on the input word. A common scalar minus cannot come from this class of pulses alone.

A supplied scalar $-S$ comparator changes the result for **every** input, including equal words. It is distinct from selecting the antisymmetric internal vector $(|01\rangle-|10\rangle)/\sqrt2$. That vector already yields a dark fringe with $+S$ while equal words remain bright. The imposed odd-worldline comparison multiplies the whole four-column operation by minus once the material origin permutation is odd. These controls identify the actual sign relative to the native reference. They do not identify the source with an electron.

## Attachment along the route

An intermediate local single-face read requires a common charge port across all four sharp branches, two opposite cyclic words on the same face when the selected reference bit changes, and support of that face incident to the material origin it is meant to accompany. On the 31-pulse route, counting the initial and final configurations as cuts 0 and 31, the first source marker meets these conditions at cuts **0, 28, 29, 30, 31**. The second meets them at cuts **0, 1, 22, 31**. Both meet them together only at the endpoints. This is exact collision-free transfer of the complete endpoint information through an addressed apparatus. Persistent local attachment and a free exchange phase remain separate premises under the [published OpenWave criterion](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#simplest-test-per-criterion).

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/dressed_exchange/intact_path/junction/marker_transfer/measure.py
~~~
