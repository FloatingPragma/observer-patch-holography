# Terminal reflection prevents transmission of both field amplitudes

IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.

Observers at port 116 can record charge and a local field doublet at the end of either archived longer guide. With the receiving fan initially in its vacuum, every accepted receiver record has the same reflection-even field direction. This holds at every exposure for both charge signs, including when the other fans are unobserved. Complete field-code transmission fails even though every incoming superposition can yield a joint click. [Protocol](TERMINAL_CHANNELS_PROTOCOL.md) · [Exact receipt](terminal_channels.json)

## The exact reflection

Reversing the order of the five rim ports on the fan at 116 fixes the gateway and leaves all its local rim degrees equal to five. The hub stays fixed at degree six. The head counts $(3,1,2,0,2,1)$ and their charge conjugate are preserved by the same reversal. Moving the directed heads with the ports gives a permutation of the full weighted configuration graph. The executable checks this permutation on all **800 configurations** of the three-fan route and all **3,840 configurations** of the four-fan route for each sign. It commutes separately with the generator containing local moves and the term containing bridge moves. The four-fan source at port 68 has its own independent reflection.

In each of these terminal fans the vacuum is reflection-even. The decoder in the two original odd-ring axes transforms by $-X$, where $X$ exchanges the axes. Thus the reflection-even direction is their difference and the reflection-odd direction is their sum. At the receiving gateway the reflection acts only on the local field columns and leaves all other detector words intact. For an initially even receiver the odd outcome therefore has exactly zero probability in every environment word at every time. An ideal energy filter preserves this rule because its operator is a function of the same generator.

Every successful operation element to the receiver has image in the same one-dimensional even line. For an input entangled with an untouched reference, the accepted reference and receiver output is a product of a conditional reference record with that one receiver direction. The rate of success can depend on the input. The accepted observer's local field cannot carry its two complex amplitudes.

## Arrival and a sensitive control

The following orders are the first nonzero powers of $2H$ after using the original receiver and integer-scaled source and readout amplitudes. An order identifies a nonzero Taylor coefficient. It gives exact evidence of an arrival for sufficiently short nonzero time without a fitted exposure.

| Guide and environment record | Inputs with a nonzero witness | Receiver starts in vacuum | Receiver starts in odd field line |
| --- | --- | ---: | ---: |
| Three fans, other fans unobserved | Both original odd-ring inputs | 6 | 3 |
| Three fans, uniform environment vacuum required | Both original odd-ring inputs | 9 | 6 |
| Four fans, other fans unobserved | Both source reflection lines | 8 | 5 |
| Four fans, uniform environment vacuum required | Reflection-even source | 11 | 8 |

With the other fans unobserved, two independent incoming amplitude rows occur at the first nonzero order for every setting in the table. Their determinant is nonzero at leading order, so for all sufficiently short nonzero times every normalized input has a strictly positive chance of a joint click. With three fans and the uniform environment vacuum required, independent rows occur at two different Taylor orders. No fixed input is dark for every time, although the successful amplitude at a given time has at most one input direction. With four fans, a uniform vacuum requirement selects the even direction at the source as well. If the source is odd, any successful receiver leaves odd reflection in the unobserved source fan and the uniform vacuum comparison refuses it.

Replacing just the receiving fan's vacuum with its odd field preparation switches the allowed receiver output to the odd line and leaves the even click probability exactly zero. The same decoder detects that odd line at the indicated orders. This locates the selection rule in the native generator and receiving preparation. The result applies to the specified symmetric terminals on the 20-cell board. A route with degree-distinguishable terminal rims changes the symmetry hypothesis.

Run from the repository root to compare the exact receipt. [Executable](terminal_channels.py)

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/fan_rotation/terminal_channels.py
~~~
