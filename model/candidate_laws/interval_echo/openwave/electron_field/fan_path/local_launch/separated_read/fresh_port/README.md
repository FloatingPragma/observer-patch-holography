# What the 1010 receiver can recover

A sharp material arrival at vertex 1010 has **one** accessible local word
when the observer reads only that vertex and its six incident head
records. The complete cell containing 1010 has **two** distinct words at
that same material address. Those words can give a nonseparable record
with the source reference on both charge signs. The difference between
these scopes is a triangle inside cell 84 that misses vertex 1010.
This test uses the four-fan written-field comparison law on an Interval
Echo board. [Boundary and timing protocol](PROTOCOL.md) · [Records](records.json)

## Local support and the two direct comparisons

The physical cell has 12 material addresses, 30 internal links, and three
boundary links. The reached comparison-law carrier uses six of those
addresses, nine internal heads, and two boundary heads. Unused links remain
frozen. The executable independently reconstructs all 3,520 global
microscopic records. For either sign, the three active head orientations
incident to 1010 are fixed whenever the load is there. So are the remaining
three incident heads, because their links are frozen. Conditioning on arrival
then yields a fixed local report tensored with the possibly changing
conditional state of the 98 reference. Every receiver channel supported
only on that strict vertex has a measure-and-prepare form on those
accepted trials. This is an exact all-wait statement on the declared
carrier; it concerns a fixed-time read that does not drive local fields
from adjoining vertices into the aperture.

The cell observer may compare all 30 internal heads and its own
material. Either convention for the two active boundary heads gives the
same accepted support, since both are fixed whenever cell 84 has the
load. At sharp 1010, the only excess-chart words are 70 and 101. They
differ by the actual local face $(1011,1012,1017)$. The vacancy uses the
complementary ordered words. The two faces in the wider 1011 analyzer
also lie wholly inside cell 84. There are 20 material-and-head words
if the cell observer leaves the within-cell address unresolved.

After the original **27.0625** 1045 herald and a native wait of **1/2**,
**0.10240687%** of original starts place the charge at 1010. The 1010
observer with access to the entire cell gives a nonseparable result on all
of that two-word bank, with unnormalized partial-transpose minimum near
$-1.694835\times10^{-5}$. Reading the actual local face and the retained
98 reference on fresh copies gives the counted witness. The earlier 1045
membership is part of this preparation. No terminal comparison relies on
1045, any route flag or any spectator vacuum. The first-gate and
outside-1010 refusals remain in separate bins. The first gate refuses
98.54360517% of original starts. Another 1.35398796% pass the first
gate but miss 1010. The maximum CHSH value over unfiltered traceless
Pauli analyzers is [1.28657573, 1.28658000] for either sign, below the
local bound 2.

At uninterrupted clock 31.0625, sharp 1010 accepts **1.68045123%** of
original starts. Its complete owned-cell sharp-address read has a
strictly positive partial transpose, with minimum near
$+2.31622\times10^{-4}$, and hence a separable two-qubit result. This is
the same uninterrupted clock at which the old 1045 face and the source
reference have a negative partial transpose. [Old-station comparison](../README.md)

## A wider local read at the same clock

At uninterrupted time 31.0625 the same signed preparation at hub **1011**
gives **3.17713822%** of original starts. There are four words at that
address, ordered (86,117,334,365) in the excess chart. Two active local
triangles at fixed charge give an actual two-bit analyzer. Face 3 toggles
the first bit and face 0 toggles the second. The two edge sets are
disjoint. Each bit has an exclusive head for the diagonal read and a
native face rotation for mixing its two orientations. A supplied local
head phase gives the remaining analyzer axis. These are terminal
read controls after the original signed-material appointment.

The measured joint comparison uses the two cell bits and the untouched
reference at 98. Its unnormalized partial transpose is negative near
$-4.264901\times10^{-4}$ on both signs. Every two-word projection in the
1011 register, including the nonadjacent pairs, has a strictly positive
partial transpose. A scan of single-face tests would miss this
comparison. Three of the four three-word subsets retain a negative
partial transpose. The unrestricted 20-word calculation in cell 84 also has
a negative partial transpose, near $-0.00601658$ on **11.4571944%** of
all starts at that clock. Its positive operational subexperiment is
the two-face reader at 1011. The receipt does not require a whole-cell
interferometer across material addresses.

## Witness and scope

### Simultaneous 99% count intervals

| Receiver and free wait | Excess | Vacancy |
| --- | ---: | ---: |
| Early gate then native 1/2, 1010 | [−0.017088, −0.015127] | [−0.017084, −0.015122] |
| Uninterrupted 31.0625, 1010 | [+0.038205, +0.038690] | [+0.038185, +0.038670] |
| Uninterrupted 31.0625, 1011 | [−0.013513, −0.013220] | [−0.013559, −0.013266] |

The fixed witnesses, chosen after the exhaustive local scan, are the
partial transposes of positive rank-one effects on the shared cell and
reference systems. Their expectations are nonnegative on every
separable record. The count banks estimate their Pauli expansions
with 2³² simulated original starts per setting and simultaneous 99%
intervals across both signs and all three displayed cohorts. Resolving
the reference in its preparation basis dephases the comparison and
gives nonnegative witnesses. Independent sparse propagation checks the
displayed appointments. Numeric enclosures include sparse Taylor
truncation and roundoff, read contraction and eigensolve checks.

The results identify local carrier degrees of freedom that can preserve
reference information. In both positive instruments the observer needs
heads outside the strict material vertex. Their measured analyzer axes
act on head orientations, so identifying them with a particle's physical
spin rotations requires a further calibration. The finite reference
records are model observations.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.fan_path.local_launch.separated_read.fresh_port.measure
~~~
