# Executed law gates on one carrier

The finite development claim concerns one twelve-port carrier, retained preparation
`G`, one integer register load and zero initial tape. The observer phase curve is
a declared encoding of its current records. This audit supplies no independent
electromagnetic field sensor and awards no OpenWave row. The pointer reads current
counters; the curve calculation is a diagnostic of that encoding.

The audit implementation is [audit.py](audit.py), its independently generated
counts are [ORACLE.json](ORACLE.json), its execution receipt is
[AUDIT_RESULT.json](AUDIT_RESULT.json), and its domain is
[CONTRACT.json](CONTRACT.json). The source and the BR6/BR9 findings were read before
coding this implementation. The separate oracle finishes before any producer
function is imported. This is an independently coded audit, not a blind or
prospectively registered acceptance run.

## Exact source and receipt pins

The frozen scientific source commit is
`19288e2b53b680e244e8ec5f4c5d45e000b2ef51`. Paths in this table are relative to its
`code/` directory. The supplementary implementation is identified by its own
content hash; it is not asserted to be a member of that frozen commit.

| Object | SHA-256 |
| --- | --- |
| `openwave_row1_protocol/law_controls.py` | `27842accc12759e0a60a84a01d94fe780eea27c58c50a75ed0e075aac4ed18c0` |
| `openwave_row1_protocol/LAW_CONTROLS.json` | `756a15690f0b720d215323757fed1ddf0c9f123979f81bdf66b9526c042054ce` |
| `openwave_row1_native/manifest.json` | `d5e692b0597dbbeafabb88e6f7d417da8f7d7d1fd39a334c08fa73ffe40e8d90` |
| `openwave_row1_native/source_packet.json` | `07ecff52d11e4f7fcdf907b17488603b0d0a5ff1644b01c583a86ddbbd0ec2ea` |
| `openwave_row1_native/native_kernel.py` | `b84e5ecdfbb8727a91571ac1b3ec6614f88ebf64b72fab0e8d35c338f08b1cff` |
| `openwave_row1_native/native_charge_winding.py` | `8526fbf0cb417172d5df320a8ed75a9dea3fe1ceaaefc4591a5f478344b8dbec` |

## Law, unit and regional domain

Orient seam `e=(a,b)` with `a<b`, and give its incidence column `B_e` entries `+1`
at `a` and `-1` at `b`. Each update satisfies

\[
n'=n-jB_e,\qquad F'=F+je.
\]

Induction from `F=0` and `n=G` gives `BF=G-n` exactly for any scalar current law,
including the half-current mutation. Continuity alone therefore proves neither
integrality nor nonfragmentation.

For the native law, balanced integer division of the two endpoint loads gives
`n'_a=floor((n_a+n_b)/2)+c*((n_a+n_b) mod 2)`. If the total source load is one,
the endpoint total is zero or one. Each event preserves a single whole unit.
Consequently every regional difference `G(A)-n(A)` belongs to `{-1,0,1}`.
The unit is one register load: the same rule sends `(2,0)` to `(1,1)` for both
coins. It does not establish an indivisible physical electric unit for arbitrary
normalizations.

The oracle enumerates simple seam words through length three and selects the
lexicographically first shortest word for each ordered pair. There are 60 pairs
at distance one, 60 at distance two and 12 at distance three. Retaining every
coin word gives `120+240+96=456` histories, 888 events, 564 holds, 324 moves and
132 arrivals. Fair history weights normalize separately for each schedule.
This scope contains one selected shortest schedule per pair, rather than every
possible shortest schedule.

The regional test explicitly enumerates all 4,096 subsets for every distinct
exact prefix charge. It agrees with the producer's singleton-integrality and
extremal-subset shortcut at each of the 5,328 steps across all four arms.
Checking only the real interval `[-1,1]` would admit fractional regional charges;
the pinned implementation also requires singleton integrality.

| Arm | Histories | Exact true-tape continuity | Closed curves / correct degree | Whole unit at every step | Integral regional domain at every step |
| --- | ---: | ---: | ---: | ---: | ---: |
| Native | 456 | 456 | 456 | 456 | 456 |
| Half-current | 456 | 456 | 0 | 0 | 0 |
| Random integer | 1,368 | 1,368 | 1,368 | 504 | 1,069 |
| First moving record omitted | 456 | 456 | 456 | 456 | 456 |

The omission arm distinguishes the true cumulative tape from held observer
memory. Exactly 228 histories lose a nonzero record. Their held memory violates
`B F_held=G-n` in all 228 cases, while all 456 held tapes remain legitimate
integer curves with degree `B F_held`. Source-load comparison detects the loss;
the tape's topology has no internal detector for it. Holds write zero current.

## Actual curves, sampling and homotopy

The auditor reconstructs oriented primal triangles and dual walks separately,
then evaluates the actual producer functions against analytic curve samples and
exact current periods. For the declared family, with `|turns|<=3`, the global
phase-derivative bound per subdivision is

\[
\frac{0.62+2\pi(3+0.23)}{32}<0.654<\pi/2.
\]

The 32 subdivisions are fixed for all admitted tapes. The density is not chosen
from the observed number of turns. The frozen package's adaptive curve function
uses 16, 16, 32 and 48 subdivisions for zero through three turns; the supplementary
function uses 32 in every case. The actual fixed function rejects turns four and
32, and rejects a density of four. A smooth 32-turn loop sampled at 32 subdivisions
has a numerical readout near zero and passes a sampled-increment test alone.
The admitted turn bound is therefore a necessary part of the measurement contract.

Half-current updates preserve continuity but fail shared endpoint gluing. Their
first nonzero half-turn has endpoint defect two. The auditor tests the actual
curve outputs on the law-generated histories, including their deformed variants;
the result is not an assertion on a typed-in half-turn constant.

For an integer dual-vertex potential `h`, a change of edge lift by
`h(head)-h(tail)` admits the family

\[
H_s(e,t)=\exp\{2\pi i[s h(\mathrm{tail})+
 (F_e+s(h(\mathrm{head})-h(\mathrm{tail})))t]\}.
\]

Its endpoint phases agree at every shared vertex, and each closed contour's
degree is unchanged. Vertex phases move during this free homotopy. Holding every
vertex phase at one would forbid changing an edge's relative winding. The
executable tests a nontrivial vertex-potential example at five intermediate
stages and compares its endpoint with the actual phase readout. It does not
recompute the integral lattice certificates in the separate BR6 audit or extend
their one-carrier conclusion to glued carriers.

## Receipt checks and reproduction

The auditor verifies exact event witnesses, every arm's independently computed
counts, finite residuals, required check names, source hashes and development
scope. A receipt with empty checks cannot pass through vacuous `all(...)` logic.
Controls reject forged passing counts, omitted checks/pins/arms/events, a NaN
residual, altered density and official-status escalation. Actual constant-zero
phase and literal-true regional implementations are also challenged.

Run from the workspace root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 openwave-tuning/next/CHAIN2/law_gate_audit/audit.py
```

The program writes only in this audit directory and checks the frozen package's
file hashes before and after execution. Reserved Fable acceptance configurations
are neither read nor run. The result concerns finite development evidence;
physical identification and official validation require their own evidence.
