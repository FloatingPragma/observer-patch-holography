# Which action follows from native repair?

The full C or E instrument fixes its local Liouville operator. The classical
endpoint rule alone does not. This attachment derives the exact distinction
across the existing [record-law family](../record_laws/README.md), on one
carrier. It also constructs an interacting Hamiltonian **component** of a
phase-chosen repair generator. That component and the accompanying dissipation
vary while every classical coin/current record law and the compatible fixed
space remain unchanged. No federation is needed for these results.

This is additional finite evidence for the existing
[native action construction](../native_action/README.md). It does not replace
its C/E instruments, select a new physical law or derive a Lorentzian action.

## Domain and scheduler

On a binary seam use basis `0=00`, `1=10`, `2=01`, `3=11`. This four-state
sector is invariant under native repair. The collision-containing 78-state
total-two source is used separately to check actual classical records; the
Hilbert--Schmidt action theorem below is **only** for the binary sector.
The collision channel is not unital and does not inherit its positivity proof.
The three-port quantum control uses the existing full binary tensor extension,
including its three-load sector; it is not described as the total-two
preparation. Restricting the identity to invariant number sectors is valid.

Let `Phi` be a chosen unread seam channel and `Phi*` its Hilbert--Schmidt
adjoint. One discrete attempt applies `Phi`. Independently assigning a rate-one
Poisson attempt process to each seam gives the optional generator
`L=sum_e(Phi_e-I)` in a dimensionless parameter `s` measuring expected attempts
per seam. Its exponential follows directly by averaging the discrete counts.
It is not the same experiment as replacing one attempt by `exp(Phi-I)`.
No laboratory-time or observer-tick identification is assumed.

This is a concrete finite instance of the completely positive generator
framework of [Lindblad (1976)](https://doi.org/10.1007/BF01608499). No general
semigroup theorem is needed for the calculation: the Poisson sum proves
complete positivity, and Kraus completeness rewrites `Phi-I` directly in
Lindblad form. The result here is the native record-table classification,
its commutator criterion, and the attached exact source controls.

## Exact action classification

Write the record table as `(k,g,h,m,n)`, in the conventions of the existing
classification. The coefficient `n` involves collisions and drops out of the
binary restriction. Its channel has the following orthogonal eigenspaces:

| Operator | Channel eigenvalue |
| --- | --- |
| `E00`, `E33`, `E11+E22` | `1` |
| `E11-E22` | `0` |
| `E03`, `E30` | `k`, `conj(k)` |
| `E01`, `E02`, and their adjoints | `g`, `g`, `conj(g)`, `conj(g)` |
| `E31`, `E32`, and their adjoints | `h`, `h`, `conj(h)`, `conj(h)` |
| `E12`, `E21` | `m`, `m` |

Here `Eab=|a><b|`. The table follows by inserting the allowed same-current
transitions, independently of any continuum ansatz. Its positivity constraints
imply `|k|<=1`, `|g|,|h|<=1/sqrt(2)` and `|m|<=1/2`, with `m` real.
Consequently `I-Phi` is self-adjoint and positive semidefinite **exactly when**
`k,g,h` are real. Without that additional condition it is not a positive
Hermitian action, although its symmetric part always is.

For every operator `X`, the latter has the exact quadratic form

```text
D = I-(Phi+Phi*)/2,
<X,D X> = (1/2)|X11-X22|²
          + sum_(a != b) (1-Re(lambda_ab)) |Xab|² >= 0.
```

Thus the native classical averaging fixes one dissipation term. The other
coefficients depend on record coherences. A common compatible space does not
fix them. The action is a form on operators, not an energy on state vectors.

## When is the skew part a Hamiltonian?

Put `J=(Phi-Phi*)/2`. It kills every diagonal matrix unit. If
`J(X)=-i[K,X]` for a Hermitian `K`, then `K` must commute with every diagonal
projector and hence be diagonal. The entries of `J` therefore give the
necessary and sufficient condition

```text
Im(k)+Im(h)-Im(g) = 0.
```

When it holds one may choose
`K=diag(0, Im(g), Im(g), Im(k))`. Inserting it in all sixteen matrix units
proves sufficiency. Necessity also follows from the Leibniz rule on
`E03 E31=E01`.

For a counterexample, E followed by phase `-i` on the full-seam output has
the admissible table `(i,0,0,0,0)`. Its skew part obeys

```text
J(E03 E31)-J(E03)E31-E03 J(E31) = -i E01 != 0.
```

So even the skew part of a valid repair generator need not be a Hamiltonian
commutator. This does not contradict a Lindblad representation, whose jump
terms can themselves have a Hilbert--Schmidt skew part. Taking the imaginary
part of a transfer is not a general method for deriving physical energy.

## A concrete interacting component, with an unselected coefficient

Keep C's exact Kraus labels `(coin,current)` and multiply each waiting
single-load amplitude by `exp(i theta)`. Every actual classical selective
map remains identical, including its probability, outgoing configuration and
record label. The instrument remains complete. Its table is

```text
k=1,   g=h=exp(-i theta)/2,   m=0,   n=1/2.
```

This construction therefore preserves every finite adaptive classical history
whose preparations and controls preserve the diagonal algebra and use those
records. This assertion uses the selective maps, not only the unread channel
conditions; it does not include coherent interventions that mix populations
and coherences. The receipt checks all 156 source state/coin branches for each
of four exact phases. All sixteen local quantum matrix units are compared
with an independent construction using the native endpoint function.

All these laws have the **same** compatible operator space: arbitrary
operators on the even span `{00,11}`, plus a scalar on the odd span `{10,01}`.
It is specified independently of their updates. The eigenvalue table proves
that their fixed spaces equal it, including dimension five. Their symmetric
actions instead have waiting-coherence coefficient `1-cos(theta)/2`.

Since the commutator criterion holds, the optional Poisson generator splits as

```text
L_theta(X) = -D_theta(X) - i[K_theta,X],
K_theta = -(sin(theta)/2) sum_edges P_odd,e,
P_odd,ab = N_a+N_b-2 N_a N_b.
```

Here `D_theta=sum_e[I-(Phi_theta,e+Phi_theta,e*)/2]`. Both the Hermitian
operator `K_theta` and the dissipative superoperator are derived once this
phase-dependent instrument and scheduler are supplied. The coefficient of
`N_a N_b` is `sin(theta)`, so this is an actual two-port interaction whenever
it is nonzero. The receipt checks the full equality on all 64 matrix units
of the two overlapping native seams `(0,1)` and `(0,5)`, using the same carrier.
On any finite binary graph it follows by adding these local identities.
This is the specified Hilbert--Schmidt symmetric/skew decomposition. When the
commutator criterion holds, `K_theta` is fixed up to a scalar within this
split. Other Lindblad jump representations have their usual decomposition
freedom; this calculation does not select a unique physical energy operator.

For this family, the symmetric channel `(Phi+Phi*)/2` is again unital and
completely positive, so `-D_theta` is a valid dissipative generator. Positivity
also shows that any stationary operator of the summed generator lies in every
local compatible space: the real part of its inner product with `L_theta`
is minus a sum of nonnegative forms. Conversely every such compatible
operator is fixed. Thus adding seams does not remove the phase freedom by
demanding the same global stationary space.

| Waiting phase | Symmetric waiting-coherence cost | Coefficient of `P_odd` | Quarter-phase interference probability |
| --- | ---: | ---: | ---: |
| `0` (C) | `1/2` | `0` | `3/8` |
| `pi/2` | `1` | `-1/2` | `1/8` |
| `pi` | `3/2` | `0` | `3/8` |
| `3pi/2` | `1` | `1/2` | `5/8` |

The interference column uses the existing coherent preparation/detection
protocol at phase `pi/2`; its probability is `3/8+Im(g)/2`. Such a coherent
experiment distinguishes laws that classical records cannot. Its preparation
and readout services are declared apparatus premises.

## What this resolves and what is open

There is a derivation of a local dissipative action from the **full**
C/E instruments. Their source-reduced process cannot be autonomous unitary
evolution: every admissible table sends the pure `|10><10|` to
`(E11+E22)/2`, reducing purity from one to one half. Keeping a dilation
environment can represent this process unitarily, as the existing history
action does; that representation adds its own implementation contract.

The result above goes further by deriving a Hamiltonian component for a
specific phase choice, while showing exactly why repair agreement and
classical records do not select it. It supplies no spatial dimension,
Lorentz symmetry, relativistic field action, physical energy units or lab
seconds. An independently specified coherent observation can constrain
`theta`; a physical realization must additionally identify the accessible
fields, clock and refinement regime with controlled errors. Enlarging the
observer into a federation supplies none of these missing identifications
merely by adding ports. Federations are useful for locality and scaling tests, after the local law and observation contract are explicit.

```sh
python -m model.repair_action_audit --check model/repair_action_audit/results.json
python -m pytest -q -W error model/repair_action_audit
```

The receipt uses exact arithmetic and binds its source bytes. `--write` refuses
to overwrite a retained result. No laboratory measurements enter this packet.
