# Can local reads erase the five-bridge exchange record?

IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.

On the concrete 320-cell cycle, no sequence of independently local readings erases the exchange tag if every individual operation preserves the incident Gauss records and each cell's material charge, leaves the other bridges frozen, and omits at least one of the five cycle bridges from its quantum support. This includes adaptive selection and discarding of local classical outcomes. A supplied coherent reversal of the **whole** field loop has a different outcome law. That comparison requires internal head reversals as well as bridge reversals. [Five-cycle source geometry](README.md#bridge-records-on-cycles) · [Executable](bridge_eraser.py) · [Complete outcome receipt](bridge_eraser.json)

## Local criterion

Let $n_c$ be the bridge bit equal to one when the head points from cycle cell $c$ toward cell $c+1$. Each internal seam has precisely one head in its own cell. With the unused external links fixed, summing the port Gauss records in cell $c$ gives a constant plus $n_{c-1}-n_c$ minus the material charge of the cell. Compare two sharp records connected by an individual measurement operation. If the operation preserves those summed Gauss records and the material charge of every cell, then

$$
\Delta n_{c-1}=\Delta n_c
$$

at every cell. All five bit changes must be equal. If the operation omits any bridge, that bridge is unchanged and all five changes vanish. This argument applies to each conditional quantum operation, including interactions with a local meter that do not change its support or conservation properties. Such operations are block diagonal in the bridge word throughout a sequence. Adaptation by earlier classical records and tracing out unread meters preserve this property.

For a balanced comparison let $\psi_0$ and $\psi_1$ carry the uncrossed and crossed bridge records. For a detector outcome with effect $E_r$, sign sensitivity requires a nonzero matrix element $\langle\psi_0,E_r\psi_1\rangle$. Its real part relative to the chosen phase reference is the bright–dark difference for that outcome. Every effect of the local class above has zero such element. Reading the bridge words directly gives either bridge record with probability $1/2$ and each phase port with probability $1/2$ conditioned on that bridge record. There are no refusals in this ideal two-record ensemble.

The local theorem preserves each cell's charge during **each** operation, not merely global total charge or final cell charges. A post-flight material route that changes a cell charge temporarily, or opening extra unfrozen bridges, changes the tested hypotheses.

## Bare bridge erasure and the full record

The five bridges on the physical cycle have ten distinct incident ports. Reversing any subset of bare bridge heads changes the incident head counts unless the subset is empty. Reversing all five preserves the **summed** count of each cell, but changes one port by $+1$ and another by $-1$ inside each fan. Consequently no complete return with unchanged charge and the same fixed internal port head-count decoder can contain both bridge words in the same Gauss sector. The internal field must also change.

As a sensitivity check, apply a hypothetical X-basis read to one bare bridge or to all five and count what happens when the original incident Gauss records are checked afterward. For $m$ addressed bridges, projecting an individual fine X outcome back into the original incident Gauss sector leaves just $2^{-m}I$. All off-diagonal terms flip at least one incident port without its internal correction. Each X-outcome string and phase port has same-sector probability $2^{-2m-1}$ and sector-refusal probability $2^{-m-1}-2^{-2m-1}$. These probabilities are independent of the interchange phase. In the complete loop control there are also internal head records that remain distinct when only bare bridges are addressed, so discarding those internal records removes the cross term even without the sector check.

| Hypothetical X read | Same incident Gauss record | Gauss record refused |
| --- | ---: | ---: |
| One bare bridge | 1/2 | 1/2 |
| All five bare bridges | 1/32 | 31/32 |

These X reads serve as outside-sector controls. Projecting back and retaining only the rare sector passes does not turn either control into an eraser.

## A Gauss-preserving comparison on the whole loop

In each actual four-face fan, join its incoming and outgoing ports along a shortest retained internal path. Combining nine internal seams with the five physical bridges makes one simple 14-edge circuit. A head word with every circuit edge directed coherently can be reversed around the whole circuit while leaving every port head count, charge, and off-loop head unchanged. The executable checks all 16,384 reversal subsets of this circuit. The empty subset and the full reversal are the only masks that preserve every port count.

On the span of the two directed circuit words, let $S$ exchange them. An externally supplied coherent detector of $(I+S)/2$ and $(I-S)/2$ is Gauss preserving on this span. It reads both bridge and internal field differences. Its action is fixed independently of the interchange sign. Outside the span one can record a separate refusal. For equal normalized crossed and uncrossed test arms carrying these two words, let $z=+1$ or $-1$ be the supplied arm phase used to calibrate the read. If the loop outcome is $\ell=\pm1$ and the phase output is $b=\pm1$, the joint probability is $(1+b\ell z)^2/8$. Both phase marginals are balanced. The product $b\ell$ identifies the sign on every test run, with no Gauss or loop refusal on the prepared pair.

This is an exact comparison between gauge-compatible sharp loop records, with charge positions and spectators held fixed. The registered native ring terms are internal to cells. The whole-loop coherent meter requires a supplied operation spanning the five bridges and their internal paths. The records do not establish that a pair of transported four-face doublets produces precisely these matching full-loop arms or that the native controls realize the meter.

~~~sh
python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/dressed_exchange/intact_path/bridge_eraser.py
~~~
