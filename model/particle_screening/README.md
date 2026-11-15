# Local screening and native-edge binding tests for CP-03

This packet sharpens the particle test on the existing binary total-two
carrier. It supplies exact limits for every size of exterior cap, an improved
actively correctable sharp-`q=1` code, and a finite bond-exit calculation for
the sharp-`q=2` pair. Neither construction yet meets the localized physical
particle gate: the protected qubit still occupies an exterior shell, and the
two-port pair has no invariant native-edge bond. No species, fermionic exchange phase,
rest energy or calibrated lifetime is inferred.

The [source/core sector packet](../exchange_sector/README.md) defines the
co-moving classical core `c`, its conserved charge
`q(x,c)=sum(x[p] for p in c)`, the faithful transported diagonal reference
`sigma_c`, and the complete C/E raw instruments. The fixed-core sector
weights are `p=(2/5,5/9,2/45)` for `q=0,1,2`. The reference is a declared
comparison state, not an executed exact preparation.

## Exact minimax distance on smaller caps

For a `k`-port subset `S` of the ten exterior ports, restrict both a sharp-`q`
state `rho` and `sigma_c` to the full quantum algebra on `S`. Define

\[
 d_{q,k}=\inf_{\rho\text{ supported in }q}\;
 \max_{|S|=k}\frac12\|\rho_S-(\sigma_c)_S\|_1.
\]

In a sharp sector the exterior particle number is `n=2-q`. Average the
distribution of the number in `S` over **all** `k`-subsets. For every source
state, including coherent states, it is the hypergeometric distribution

\[
 H_{n,k}(r)=\frac{\binom nr\binom{10-n}{k-r}}{\binom{10}k}.
\]

The reference's distribution is
`H_ref,k = (2/5)H_2,k + (5/9)H_1,k + (2/45)H_0,k`.
Trace distance dominates the number-measurement total variation; convexity
under the cap average gives `d_{q,k} >= TV(H_{2-q,k},H_ref,k)`. The diagonal
reference conditioned on `q` is uniform within that sector. It has this
same number distribution on *every* cap, with no off-diagonal contribution,
so it attains equality. This establishes an exact minimax over arbitrary
coherent sharp-sector states, rather than a bound on each specified cap.

| Cap size `k` | `q=0` | `q=1` | `q=2` |
| ---: | ---: | ---: | ---: |
| 1 | 29/450 | 8/225 | 61/450 |
| 2 | 26/225 | 14/225 | 59/225 |
| 3 | 23/150 | 2/25 | 19/50 |
| 4 | 8/45 | 4/45 | 22/45 |
| 5 | 17/90 | 4/45 | 53/90 |
| 6 | 1/5 | 2/15 | 17/25 |
| 7 | 7/25 | 14/75 | 343/450 |
| 8 | 28/75 | 56/225 | 188/225 |
| 9 | 12/25 | 8/25 | 9/10 |
| 10 | 3/5 | 4/9 | 43/45 |

The last row reproduces the prior full-exterior theorem. The first row shows
that sharp `q=1` can be much closer on individual exterior ports than on the
whole exterior. A state with a single good cap is not thereby screened on
the whole cap menu.

## The all-seam code has a stronger one-port obstruction

The earlier [finite all-seam flagged-code theorem](../common_source/ENCODED_MEMORY.md#optimal-native-carrier-code)
applies to any two-dimensional exact code under the complete C/E raw menu.
It forces binary support on **nonedge source pairs** and scalar compression of
every port occupation. Fix a classical core `(a,b)` and impose sharp `q=1`.
If an exterior port `j` is joined by a native seam to either `a` or `b`, at
most one of `|a,j>` and `|b,j>` survives the nonedge restriction. Thus `n_j`
has rank at most one on the allowed support. A rank-one operator cannot
compress to a positive multiple of `I_2`, so the code's occupation at `j`
must vanish. This is an actual model specialization of the finite theorem,
not a new assumption about localization.

The thirty native seams form a degree-five graph. The exact core census is:

| Core class | Number | Exterior ports nonadjacent to both core ports | Exact `q=1` qubit |
| --- | ---: | ---: | --- |
| Native edge | 30 | 4 | Yes |
| Nonedge, distance two | 30 | 2 | Yes |
| Opposite pair | 6 | 0 | No |

For core `(0,1)`, the four safe exterior ports are `2,3,4,6`. The new code is

\[
 |0_L\rangle=(|0,2\rangle+|1,3\rangle+|0,4\rangle+|1,6\rangle)/2,
 \quad
 |1_L\rangle=(|0,3\rangle+|1,2\rangle+|0,6\rangle+|1,4\rangle)/2.
\]

For each of the 60 cores with at least two safe ports, the packet constructs
the analogous crossed rectangle (two rectangles when four ports are safe).
It checks every raw C/E input projector and the injectivity of every branch,
including zero/abort slots. Equal support counts give the exact flagged
Knill–Laflamme relation `V* K_h* K_h V = p_h I_2`. Each positive branch has
two orthogonal output columns, so an exact full-carrier unitary recovery
exists for each retained flag. One explicit extension first undoes a signed
endpoint move and then applies a Householder reflection within each of the
two disjoint logical-support blocks, mapping the normalized surviving subset
to its original uniform column. This is algebraic active correction; a new
nine-qubit source encoder and source pulse count have **not** been compiled.

For any logical state, each of the four safe exterior ports has occupation
`1/4`; the other six have zero. The transported reference has occupation
`61/450` at every exterior port. Hence this code's **worst one-port trace
distance is `61/450`**, compared with `82/225` for the prior two-port code.
The result is optimal among exact all-seam sharp-`q=1` qubits on an edge core:
every such code leaves at least six exterior ports empty. With only two safe
ports the optimum is `82/225` because their occupations sum to one; with no
safe ports no qubit exists. The unrestricted sharp-sector optimum `8/225`
cannot be reached by an exact all-seam `q=1` qubit on this carrier.

The code therefore improves small-cap screening while remaining spread over
both core and exterior. It still has full-exterior distance at least `4/9`,
requires supplied active correction after each native attempt, and does not
make a passive or geometrically localized particle.

## Co-moving transport and finite first exits

The raw C/E endpoint maps preserve `q=1` exactly when the classical core is
transported with the current. For the new code on `(0,1)`, under a uniformly
selected one of 30 seams and a fair coin, the core label changes with
probability `1/15` per native attempted repair. It becomes a nonedge with
probability `1/30`. These probabilities are independent of the unknown
logical state because the flagged effects are scalar. After record-selected
source recovery to this fixed code, preserving the same sharp sector requires
an explicit reset of the seven-bit core label whenever it changed. Retaining
a status bit and old/new seven-bit labels costs 15 archive bits per attempt
under that reset format. The exact source pulse cost for this **new** code is
unknown; the older code's compiled counts cannot be assigned to it.

There is also a simple candidate for a directly localized pair: sharp `q=2`
with the two occupied ports equal to the two core ports. The pair remains
equal to the co-moving core under every raw C/E step, so its relational charge
is protected. Starting from **any native-edge pair**, however, one uniform
native attempt produces a nonedge pair with probability `1/15` under either
law. As long as the pair remains an edge, the same hazard applies at each
attempt. Its first-exit survival is `(14/15)^t` and its mean first exit is
15 native attempted repairs. This is graph adjacency, not particle decay,
and `t` is an attempt count, not a calibrated observer clock.

Finally, the existing [complete-history agreement theorem](../observer/adaptive_record_recovery.py)
puts a strictly positive lower bound on the error of **any final decoder** of
an unknown passive binary-carrier qubit by 23 uniform C attempts or 12
uniform E attempts, with all records retained and no intermediate quantum
control. This is an information bound, not a half-life. Active correction
changes that grammar and must carry its controller, record and recovery costs.

## Status and reproduction

CP-03 remains open. The new result rules out exact transported-reference
screening and an invariant edge bond for the tested sharp sectors, while giving
an optimal small-cap active-memory control. A particle claim still needs a
localized *state and observable* that remain protected under the full record
process, with an autonomous or explicitly controlled lifetime, independently
identified mass/energy, and a coherent statistics test when relevant.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/particle_screening/test_screening.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_screening.screening --check
```

The [campaign](campaign.json) pins predecessor source files and the canonical
science-result hash. The [receipt](screening_results.json) binds the prior
packet's source import closure together with this packet's code, tests and
documentation. No physical comparison is added.
