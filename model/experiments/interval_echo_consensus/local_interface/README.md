# A local operation with a material consultation record

A native echo and a complete material reading give a normalized local operation
on the two endpoint ports and their collar. A first material record also selects
which of two enrolled seams receives the second echo and reading. All 81 public
record pairs and their conditional quantum maps are retained. The executable
[construction](measure.py), [receipt](results.json) and independent
[tests](test_local_interface.py) use the registered Interval Echo board. Protected
services, their sharp material read rule and the controller are supplied;
spacelike separation and an autonomous clock require separate attachments.

## Registers and conserved load

On the 20-cell board, seam A=(0,1) contains two ternary ports and its eleven-state
collar. Its dimension is 99. Seam B=(2,3) has the same registers with disjoint
addresses. Reordering the actual port/collar register tuple gives the split;
there is no fit to a postprocessed joint distribution. Other collars incident
to these ports belong to the complement, as do the degree-reference registers.
The basic board allocation is
`3^240 * 11^630 * 33^180`; optional fuel, public readers and idle apparatus
registers are spectators. For A alone the complementary powers are
`3^238 * 11^629 * 33^180`. The receipt checks all 9,801 words of A and B,
including dirty collars, against the split and its inverse.

The conserved port load does not give a tensor product inside a fixed global
sector. With `H_A,q` denoting local load q, its dimensions are
`11,22,33,22,11` for q=0,...,4. For total load N,

`H_N = direct_sum_q (H_A,q tensor H_complement,N-q)`.

For two enrolled seams at total four this dimension is 2,299, whereas the
unrestricted product dimension is 9,801. Coherences between the q blocks can
exist with a balancing reference; discarding them is an additional restriction.
Every deposited branch preserves the local charge, and therefore restricts to
this direct sum. This states a conservation decomposition, not an assumption
that all charge-coherent preparations are forbidden. The echo can convert
(0,2) into (1,1), so a fixed material multiset is not a valid domain for this
menu, despite its invariance under free blank-memory contact evolution.

## Complete material instrument

Let U be the native 99-state echo, with the exact public formula in the source
law. Let P_r project onto material pair r in {0,1,2}² and leave all eleven
collar states intact. The one-read branch is `A_r=P_r U`. The public alphabet
has nine values. No private collar label is read, resolved, reset or discarded
by this consultation. A material reading is the supplied sharp native service
used elsewhere in the candidate, not a derived Born response or a collar reader.

The fixed-region adaptive example has

`K_r,s = P_s U F_r P_r U`,

where `F_r=U` if the first material pair has a<b, and identity otherwise. This
is a recorded measurement followed by a native unitary, with no assertion that
U reverses the measurement. Both first and second records are retained. The
exact construction checks `sum K*_r,s K_r,s=I_99` and charge preservation on
all input columns over Q(sqrt(2)). Restricting the input to blank memory gives
complete rectangular maps from nine inputs into 99 outputs. A subsequent
operation must retain the output collar: two states with identical material
input |1,1> and different collars give next-echo material laws at total
variation distance 1/2. There is no reusable material-only nine-state channel
for arbitrary dirty collars without specifying their state or a reset.

Kraus completeness gives positivity, complete positivity with arbitrary
references, and trace preservation by the usual finite matrix identities.
The existing Lean partial-trace helper consumes exactly this Gram identity.
The full reference audit uses Schmidt rank 99 across A and B, with opposite
port loads maintaining total four. Its complete remote density matrix is
preserved to floating error below 2e-18. A separate complex spectator test
retains every spectator/reference matrix entry. These inputs audit the channel
identity; the program supplies no native preparation of every entangled audit
state. Selecting only first record (0,1) has probability 1/9 and changes the
normalized remote state by trace distance 8/9. That conditional steering is not
an unconditioned remote channel change.

## A record selects the second support

The selected-support schedule first echoes and reads A. If the actual first
record a<b, the second service echoes and reads B; otherwise it echoes and
reads A. The selected seam is therefore an output of a native material
consultation under a declared finite controller rule. The two enrolled seams
and the rule are architecture inputs. With `A_r=P_r U`, its full branches are

`L_r,s = A_r tensor A_s` if a<b,

`L_r,s = (A_s A_r) tensor I_B` otherwise.

These factored maps act on all 9,801 union inputs. For each first record,
summing the second-record Gram matrices gives `A*_r A_r tensor I_B`;
summing r gives identity on the union. This exact computation avoids a dense
9,801-dimensional matrix. It also proves reference preservation for any
complementary registers and restrictions to fixed global charge sectors.
The consultation and selected service together can touch A union B, so its
no-signalling domain is the complement of that union. On a complex total-six
input with all dirty local collars and a material reference at seam (4,5), the
unconditioned complementary density error is below 2e-16. Conditioning on the
selected seam B changes that complementary state by trace distance 0.02593.

For sharp material input A=(0,1), B=(2,0), with both collars blank, the public
event first=(0,1), second=(1,1), selected=B has probability 1/2. A control that
serves A for every record gives that event probability zero. Its full public
law differs by total variation 1/2. Thus the consulted source record changes
an observable subsequent service, beyond a relabeling of static effects.

## Access and boundary controls

Protected addressed windows isolate the enrolled service and its supplied
material read. The receipt lists every crossed seam and degree reference.
A has nine seam crossings, including the external seam (1,82), and the private
degree reference at port 0. Its complement is not a freely isolated component.
Allowing the native (1,5) contact for pi/2 changes the outside load-one
probability from zero to 1/2 on input loads (1,0). That control uses the full
99-state contact exponential. Independent collars also do not make two services
with a shared port disjoint: interchanging echoes on (0,1) and (1,5) changes
the complete three-port material law by total variation 1/2. Truly disjoint
endpoint/collar services commute, including their branch maps with references.

The source layout supplies a justified finite register split, complete local
maps and record-dependent scheduler provenance for this menu. It supplies no
unique factorization, nonconstant refinement tower, free-flight isolation,
physical distance or spacelike region, readout derivation, autonomous controller
or experimentally calibrated laboratory channel. No sampled cohort is claimed:
the receipt contains exact identities and computed observer probability laws.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_consensus.local_interface.measure
OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -O -W error -m model.experiments.interval_echo_consensus.local_interface.measure
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/local_interface
```

The explicit first-generation option `--write` uses exclusive creation.
Replays compare source bindings and discrete data exactly and floating fields
within 1e-11. The JSON is simulation evidence, with no external custody claim.
