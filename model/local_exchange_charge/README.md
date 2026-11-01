# Local preparation and transport of the exchange charge

This construction realizes the [protected ordered cover](../defect_protection/README.md)
on local port degrees of freedom, prepares an initially localized charge
perturbation, and follows its complete repair history. It supplies an optimal
local preparation bound and a sharp-charge localization obstruction. The
transported object is a two-port coherence, not an identified physical particle.

## One local law and its domain

Each of the twelve ports has states `0,A,B`. Restrict their tensor product to
exactly one A and one B: dimension 132. The isometry W maps the unordered
pair and orientation bit of the ordered cover to the positions of these two
labels. Thus the old global orientation memory has a spatial realization.
The additional coherent label at each binary port is an explicit apparatus
choice; this does not construct a new conserved observable on bare H66.

For a seam, use the existing [two-site source service](../defect_field_bridge/FULL_DOMAIN.md):
apply its coin-selected hold or label-preserving move when precisely one endpoint is occupied; hold when
neither or both are occupied. In particular, an occupied-occupied seam does
**not** swap A and B. All raw C/E outcomes measure only the native occupation
and current data. Each raw outcome is one coherent partial-isometry Kraus
map; measuring finer label/configuration data would change the law.

Dual-rail encoding `0=00,A=01,B=10` gives sixteen two-site input addresses.
The four-state extension retains `11` and its recorded collision hold. All
sixteen inputs, both coins and sixteen raw slots per coin are retained. On
the qutrit subspace the collision slot is zero. Tensor the **two-site stencil**
with spectator identities: do not replace it by the nonadjacent global CAR
hop, which has different spectator signs. This is ordinary tensor-qutrit
locality, distinct from the earlier graded CAR representation.

`core.py` checks `K_local W = W K_cover` for every native raw slot, every
seam and all 132 columns. The on-site label swap S exchanges A and B and
fixes 0. Its product Q over all ports restricts to the deck involution
`I66 tensor X`. Each complete Kraus map commutes with Q. Consequently Q and
its two sharp sectors are conserved under any finite classical adaptive
sequence of admitted repairs. The exchange loop retains its existing Z2
character; distinguishable labels and this symmetry do not derive fermionic
statistics or spatial spin.

## A localized state with a moving, conserved contrast

Let `sigma=I132/132`, let C be a set of ports, and let P_C project onto
configurations with both labels inside C. For `|epsilon|<=1`,

```
rho_epsilon = (I132 + epsilon P_C Q)/132.
```

Its eigenvalues are `(1+epsilon)/132` and `(1-epsilon)/132` once per unordered
pair in C, and `1/132` elsewhere. It is faithful on this sector for
`|epsilon|<1`, not on the full qutrit tensor space. The partial trace over C
kills every added off-diagonal term, so the **entire exterior density matrix**
initially equals the reference. The charge expectation is
`epsilon choose(|C|,2)/66`. This is a conserved excess, not a sharp Q sector.

A deterministic local preparation exists. In each both-in-C pair block,
replace the orientation by `(I+epsilon X)/2`; retain the complement and all
preparation outcomes. This is a CPTP map, not a unitary applied to the
tracial reference. The executed witness uses the two-port core `(0,1)` and
`epsilon=1/2`. Preparation environments are marginalized when defining rho;
their records are not included in the subsequent native-tape equality claim.
The initialization, disposal/readout and coherent-control resources belong
to the supplied observer service.

Write `Q_x=P_x Q` for an unordered pair. Every native branch maps Q_x to
its transported Q_y with its actual branch weight. For uniform independent
seam selection, the complete C and E channels induce the same exact kernel:
`T_yx=1/60` for each legal single-endpoint move, and
`T_xx=1-degree(x)/60`. This is the two-particle exclusion walk on the native
port graph, derived from the instrument. It is not a chosen continuum fit.
For nonnegative initial coefficients `w_x=1[x subset C]`, `w(t)=T^t w(0)`
and `sum w(t)=choose(|C|,2)` exactly.

For any region R the reduced trace distance and charge readout are

```
D_R(t) = |epsilon| sum_{x subset R} w_x(t) /132,
Delta<Q_R>(t) = epsilon sum_{x subset R} w_x(t) /66.
```

Here Q_R is the local label-swap product projected onto exactly one A and
one B in R. It saturates the trace-distance test. Every one-port marginal
contrast is zero. The part of the global contrast absent from both R and
its exterior resides in split-pair correlations. Thus continuity is on the
**pair configuration graph**, not an additive electric-current law on ports.

The same formulas hold with the complete classical native tape retained:
each unnormalized history contrast has nonnegative coefficients, and each
has zero trace. The native tape distributions of rho and sigma agree exactly.
The compared controllers use the same native-record-only policy; they do
not choose later actions using the differing preparation archives.
Tape plus a fixed exterior generally ceases to agree once both labels can
leave the initial core. Conditioning on a selected history is not used to
normalize the claimed global conservation or locality.

For C=`(0,1)`, the exact normalized contrast partitions are:

| Attempts | Both inside | Both outside | Split |
| --- | --- | --- | --- |
| 0 | 1 | 0 | 0 |
| 1 | 13/15 | 0 | 2/15 |
| 2 | 113/150 | 7/900 | 43/180 |
| 3 | 5911/9000 | 8/375 | 2897/9000 |

The exterior is unchanged through one attempt and differs at two. The
connected, holding, bistochastic pair kernel converges to `w_x=m/66`, where
`m=choose(|C|,2)`. This determines the asymptotic regional contrast exactly.
Conserved global charge therefore does not imply a permanently localized
core or a bound pair. Attempt counts are not calibrated physical lifetimes.

## Optimality and the sharp-sector obstruction

The local construction has the largest possible charge expectation under
its preparation contract. Start with sigma and an independent ancilla,
apply any deterministic CPTP map supported inside C, and require the output
to remain in the one-A/one-B sector. The state remains diagonal in the
exterior occupation/label basis. Q connects different exterior basis states
unless both labels lie inside C, so only that block contributes. Locality
fixes its probability at `p_C=choose(|C|,2)/66`, and Q has norm one there.
Hence `|<Q>|<=p_C`, attained by the displayed preparation at epsilon=`+/-1`.
Outside access, correlated initial ancillas, postselection and sector escape
are outside this contract. `optimality.py` checks the supporting matrix
elements for all 4096 port subsets.

There is also an exact obstruction for **any** sharp Q state, regardless of
its occupation marginal. Let B be the exterior and measure its projector
onto two occupied ports with the opposite exchange character. This event
has zero probability in that sharp sector, but probability
`choose(|B|,2)/132` in sigma. Thus

```
D(rho_sharp restricted to B, sigma restricted to B)
    >= choose(|B|,2)/132.
```

The uniform sharp state `(I+Q)/132` attains the bound. For a two-port core,
the minimum is `15/44`; it cannot be repaired by redistributing occupation
probabilities. This explains precisely why the localized preparation is
not a localized sharply charged particle.

## Scope, implementation and remaining joins

The source service reuses the existing compiler; it does not add a competing
repair implementation. Its local dilation is checked on the full sixteen
input columns, including arbitrary reference entanglement. Preparation,
native repair and charge readout have explicit retained records, source
complements, resource counts and numerical admission ceilings in
[`results.json`](results.json). These ceilings are execution checks, not
outward rounding certificates or hardware error rates. A short composed
history tests actual states and coherences, not only record probabilities.
The moving execution uses three ports `(0,1,5)`, prepares on `(0,1)`, applies
C or E on `(0,5)` and reads on `(1,5)`. All 64 ambient dual-rail rows remain.
Its charge expectation on the readout region is `1/264`, half the initial
`1/132`. Each repair commits five coin/record bits through native six-outcome
readouts; the 7744 noncode outcomes are retained zeros beside 32 code slots.
The moving witness requires eleven coherent logical bits including local
work and spectators, and thirteen measurement load/export pairs in addition
to its pulse services. The supplied selector and reference preparation are
explicit inputs, not a derivation of a vacuum or a cost-free reset.
The terminal detailed readout is destructive: global charge preservation is
not asserted after it. Port locality of the effective stencil does not
establish spatial locality of the supplied processor's pulse hardware.

Label dephasing at an initially occupied core port preserves native occupation records but destroys
the contrast: protection depends on label-blind service. The global tracial
sector reference has scalar modular generator and exactly zero modular-energy
excess for these states. No mass, Lorentz spin, additive electric charge,
electron/positron identification, refinement stability or nontrivial
observer clock is supplied. A physical candidate must join this same
transported excitation to a justified reference, energy/clock response
and independent particle observables; global protection alone cannot do it.

This completes a bounded local-realization and transport component of
[CP-03/CP-05](../../docs/COMPLETION_PROGRAM.md). It neither closes those
physical gates nor adds an empirical comparison. It complements the
decaying localized coherence in `native_localization` and the global
ordered-cover charge; it does not replace their distinct contracts.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.local_exchange_charge --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/local_exchange_charge
```
