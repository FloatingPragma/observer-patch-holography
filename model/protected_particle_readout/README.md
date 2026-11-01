# An internal observer reads the protected cube response

This packet executes a complete seventeen-bit Ramsey experiment on the
original three-dimensional cube. Its unconditional pointer probabilities
are attached, with exact source and full-Fock evolution error bounds, to
the eight interacting spectral lines already derived from the same action.
This identifies **the local observable's eight-dimensional cyclic response,
plus its elastic reference component**, as something an observer in this
specified computed universe can measure. It does not identify all1820
protected states or all119 eta-three states with this one probe.

The model clock, source encoding, blank preparation, complete repair,
pointer readout and repetition/reset services are specified simulator inputs.
No laboratory clock, SI scale, Standard Model species, vacuum particle pole,
or arbitrary-law uniqueness is inferred. The reference R is excited. The
measured modes are neutral finite response modes relative to that reference.

## Actual source experiment

The field is the same open2x2x2 cube, sixteen occupation modes, c=1/8,U=1.
The inherited source first prepares

`R=(eta+)³|vac>/sqrt2016`, with `H R=3R`,

using its5379-pulse eight-even-mode signed tree and24-pulse pair copy.
The tree's exact spectator-blank embedding retains every source export.
An initialized seventeenth pointer bit is then prepared and coupled through

`controlled F`, where `F=exp(i*pi*n_origin/2)`.

Both pointer arms remain at ideal number N6. The actual source state retains
all131072 field/pointer amplitudes, including numerical residuals outside
that ideal sector. Controlled F has an exact six-term diagonal Pauli
decomposition including its scalar phase; the emitted18-pulse word and the
one-pulse pointer preparation are executed through the original M6 source.

Both arms evolve under the original H for model time1, divided into two
half-time segments with32 Strang steps per segment. Each segment is the
unchanged337920-pulse `full_fock_field.source_word`, including every routed
fermionic SWAP and every contact factor. After the first segment apply C on
all five admitted native seams; after the second apply E on all five. The
domain remains the four x-directed dimers at ports0,1,5,7, with seams
(0,1),(0,5),(0,7),(1,5),(1,7). Every native, zero and outside record is kept.
Unexpected native branches terminate and retain their output state. The
outside branch continues without normalization. Its intermediate weights
are not added repeatedly to the final probability sum.

For each final setting phi=0 or pi/2, execute controlled F inverse, the
pointer phase and its final Y rotation. The actual six-outcome source export
retains both pointer values and all four abort labels. The field output is
not uncomputed or discarded by the execution. The observer can ignore it
when calculating the pointer marginal.

## Complete block execution, not ideal gate substitution

Executing every M6 pulse directly on131072 rows is unnecessary. Each
hopping, FSWAP or contact word acts on only two local bits; the pointer word
acts on three. `execute.compile_block` composes its actual literal M6
success matrices and **every first-export failure matrix** locally.
If U_j is the actual source pulse, S_j its proper-code corner, and F_j,r
its r-th export complement, the local output maps are

`S_m ... S_1`, and `F_j,r S_(j-1) ... S_1` for every j,r.

Tensoring these maps with all untouched field/reference factors is exactly
the complete block instrument of the original word. No ideal hopping
matrix replaces a native block. No number, eta or spectral projection is
applied to the actual state. `circuits.field_binding` checks the flattened
block tape against the original complete field tape.

The immutable `Run` retains its initial full state, local native maps,
ordered blocks, pulse offsets and all first-failure weights. Its
`failure(block,event,slot)` reconstructs the actual rectangular failure
vector after the complete preceding success prefix. The output row basis
lists surviving local work axes before exterior spectators; it may differ
by a known row permutation from a flattened global evaluation. Reference
Grams and probabilities agree. Regression tests compare all failure maps
on nonadjacent/reversed targets and multiple untouched references, including
a deliberately order-one leaking M6 control.

`experiment.repair_output(law,seam_index,record)` similarly reconstructs
every native or outside repair output from its actual field/source prefix.
The binary archive retains all per-pulse four-column failure tables and
both complete final analyzer states. Pair-code rows have two real aborts;
the other two columns are explicitly identified padding, not extra outcomes.
All zero-weight outcomes remain addressable. Source failures stop: histories
with a later second failure are not spuriously added.

Local and global floating matrix checks are diagnostics with separately
reported admissions. Very small negative computed failure weights below
1e-40 are rounded to zero only in the numerical weight table; the count,
total correction and per-entry limit are reported, and the actual Kraus
maps remain unchanged. The mass/operation-scaled failure-sum admission is
also reported. None of these numerical checks supplies the source theorem.

## Exact source-to-spectrum attachment

`bounds.py` independently consumes every native routing template, including
all fixed scalar phases. The exact dyadic trigonometric and literal-port
certificates bound the represented source controls. The inherited signed
preparation tree has all255 conditional splits checked against integer
target amplitudes. Source errors telescope on the initialized preparation
and on arbitrary inputs/references for subsequent operations.

The unchanged full-Fock Strang coefficient is629/256. Consequently the two
field segments contribute exactly

`2*(629/256)*(1/2)^3/32² = 629/1048576 < 0.0006`.

Adding every literal source error gives the retained rational complete-state
bound for each setting. Exact complete repairs are contractions on the
complete output space and are the identity on the ideal protected state;
their unexpected actual output branches are included in that distance.
For normalized complete Stinespring outputs, a vector-distance bound delta
also bounds each unconditional pointer-event probability by delta.

The independent ideal comparator evolves the original physical H on its
proved invariant N6 sector; that compression is used only for comparison.
Exact residue/root intervals and dyadic sine/cosine bounds separately give

`C(t)=1/16+(15/16) sum_j w_j exp[-it(E_j-3)]`,
`p0(phi,t)=(1+Re(exp(i phi) C(t)))/2`.

Thus the two settings read Re C and minus Im C. The exact intervals enlarged
by the complete source/field error enclose the actual mathematical native
probabilities. Floating execution is checked against these intervals with
its separately stated numerical admission. Both quadratures at time1 are
separated from their zero-time values. A complete native zero-wait control
and an independent constant-energy-shift control are retained.

The full run uses681283 native pulses per setting. Its pointer-zero
probabilities are0.994249701120 and0.501979889450, compared with the independent
original-H values0.994249466284 and0.501979716912. The complete-output distance
is about2.27e-6 numerically; the exact source/field upper bound is below
0.000599861231. The numerical norm drift is about4.45e-11 and remains a
reported floating diagnostic.

There is also an exact matched contact-free control, without executing an
additional native U=0 word. On this bipartite cube let
`Gamma=(-1)^(sum_odd_sites n_x)`. Integer sparse checks verify
`Gamma H0 Gamma=-H0`, `H0 R=0`, and `Gamma R=R`, `Gamma FR=FR`.
Thus `<FR|exp(-it H0)|FR>` is real for every time, and the ideal contact-free
pi/2 pointer probability is exactly1/2. The interacting time1 imaginary
fringe remains above1/2 after both the exact source error and the declared
sampling error. This rejects that matched contact-free model; it is not a
classification of all possible alternative interactions.

Only the differences E_j-3 are measured by this fixed-N6 two-arm protocol.
H->H+cI and a common chemical potential H->H-mu N produce the same pointer
law. The absolute H energies are assignments from the supplied action,
not a measured absolute energy zero.

## Simulated counts and the separate resolution problem

The packet also gives a reproducible finite simulated-record example:
NumPy PCG64 seed20261101,2^24 independent reset runs per setting. A record
is coarsened to final pointer zero versus **all other complete outcomes**,
including every source or repair termination. This uses the unconditional
probability; there is no survival selection. The binary coarsening does not
pretend to be an archive of all individually resolved raw records.

Hoeffding and a union bound give simultaneous sampling error1/2048 with
failure probability at most `4 exp(-8) < 4(3/8)^8 < .01`. This statement is
around the numerical probabilities used to generate these simulated counts;
the separately admitted host-evaluation discrepancy is retained when
comparing to exact native probability intervals. The reported frequencies
are checked against those intervals enlarged by sampling and host errors.
These are records in the specified simulated experiment, not laboratory data.
The retained zero counts are16680454 and8423380 out of16777216 each;
the corresponding frequencies are0.994232535362 and0.502072572708.

Two settings at one time do not reconstruct eight unknown spectral lines.
`identification.py` gives a separate prospective resolution certificate:
at times0,39,...,273 the eight-frequency Vandermonde determinant exceeds1400
and its inverse infinity norm is below2744/525. Its stated noise allocation
and77,309,411,328-shot schedule would resolve the eight known-frequency
weights under the supplied long-time systematic-error and reset premises.
That long-time schedule, its required native accuracy and floating inverse
are **not executed here**. The executed result identifies the same-source
protected response through complete observer probabilities; it does not
claim blind spectral tomography or identification of its degenerate dark
directions or the other free-energy modes.

Run from the repository root:

```sh
.venv/bin/python -W error -m pytest model/protected_particle_readout -q
.venv/bin/python -W error -m model.protected_particle_readout --check model/protected_particle_readout/results.json
```
