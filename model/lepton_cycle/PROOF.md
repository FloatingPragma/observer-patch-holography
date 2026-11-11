# Complete native lepton cycles and their finite accuracy

## Source and preparation

The source consists of two independent cells, each with a blank double load
at port 2 and the sixty ordered addresses of the public seam program.
`instrument.source_data` calls the actual `Walk` completion through
`lepton_register_round.rounds`. The native event gives fifty blank, five
positive-current and five negative-current addresses per cell. The forward
and inverse native ticks preserve the full source record. The analyzed
pullback returns to the blank double at the original port; the address
register carries its nonreturn information.

A fresh pair has address state `S=u tensor u`, where every component of
`u` is `1/sqrt(60)`. An explicit real Householder unitary maps the first
address basis vector to `u`: for normalized `v=(e_0-u)/||e_0-u||`,
`H=I-2vv*`. This supplies an executable preparation map on the address
space. Fresh blank double loads and the Householder service are declared
inputs. Neither their autonomous production nor their elementary native
pulse compilation is asserted.

The two return tests are the orthogonal projectors
`P_x tensor P_y`, `(I-P_x) tensor P_y`,
`P_x tensor (I-P_y)`, `(I-P_x) tensor (I-P_y)`, with `P=|u><u|`.
Every result has its own two-bit record, including zero-weight outcomes.
For a source amplitude matrix `V_r` with one column for each register mode,
the outcome map multiplies register entry `rho_kl` by
`G_r(k,l)=sum_a V_r(a,k) conjugate(V_r(a,l))`.
These are Kraus maps with diagonal Kraus operators indexed by the source
coordinate `a`. Their Gram form proves complete positivity, including
arbitrary untouched quantum references.

After the tests, SWAP the old pair with a fresh prepared pair. The old pair
is exported together with its record and every correlation with the
register or reference. The active pair is exactly `S` for the next round.
This is a unitary exchange with supplied fresh storage, not an erasure of
an unknown state. `reset_branch` represents that factorization without
allocating the square of the 3,600-dimensional source space. The complete
Kraus sum is trace preserving, so the untouched reference marginal is
unchanged. These are the concrete input objects required by the finite
Kraus/no-signalling framework in
[`SlotLocalitySurface.lean`](../../Lean/EventAlgebra/SlotLocalitySurface.lean).
No additional Lean proof or physical locality attachment is claimed.

## Complete oriented-current round

Let `J=diag(-1,0,1)`, `U_h=exp(-ihJ)` and `p=2/9`. The native forward,
exclusive-read kick and inverse tick give

```
S tensor psi + X tensor (U_h-I) psi,
<S,X> = ||X||² = p.
```

The source's four projected squared norms of `X` are exactly
`4/81`, `11/324`, `125/1296`, `55/1296`. Thus the complete register
instrument has Kraus operators

```
A_yy = (1-p)I+p U_h,
A_xn = sqrt(11/324)(U_h-I),
A_yn = sqrt(125/1296)(U_h-I),
A_nn = sqrt(55/1296)(U_h-I).
```

The first operator includes interference with the unmarked standing source.
The other three weights sum to `p(1-p)=14/81`. Their effects sum to identity;
discarding the outcome and exported source gives exactly

```
Phi_h(rho) = (1-p)rho + p U_h rho U_h*.
```

The reset is essential to this composition. Native source execution without
reset changes the exclusive count in the four return sectors to
`2/9,7/9,7/9,2/9`. Counting only successful returns does not execute the
complete repeated channel.

## Complete mass-read round

The same carrier has a three-character register. Its analyzer is the equal
superposition of the three normalized, mutually orthogonal source flags
(blank X/written Y, positive X/blank Y, negative X/blank Y).
Writing `q_k=theta+2 pi k/3`, native compression gives

```
lambda_k = sqrt(5/108) + 2 sqrt(5/216) cos(q_k),
e_k = lambda_k²,       E=diag(e_k),       tr E=5/18.
```

On character `k`, let `U_k` be the actual native tick with its character
phase, `phi` the analyzer, and `X_k=U_k*|phi><phi|U_k S`.
Then `<S,X_k>=e_k`. The analyzer gives equal weight to each character phase,
so the root-of-unity identity `1+w+w²=0` proves

```
<X_l,X_k> = delta_lk e_k.
```

The code checks this identity on the actual 60-by-60 source arrays. A kick
`z=exp(-ih)` on the analyzer, followed by the inverse tick, gives the column
`S+(z-1)X_k`. All four return projectors are executed separately. The
returned-source branch has operator `I+(z-1)E`. The sum of the three
nonreturn Gram matrices is

```
B = diag(e)-e e^T.
```

For any vector `v`, weighted Cauchy-Schwarz gives
`|sum_k e_k v_k|² <= (sum e_k) sum e_k |v_k|²`; therefore `B` is positive.
At `theta=2/9` every `e_k>0` and `sum e_k<1`, so `B` has rank three.
In particular, replacing all failure environments by one scalar Kraus
operator would change the off-diagonal records.

The exact complete mass channel preserves diagonal entries and has

```
Psi_h(rho)_kl = [1+(z-1)e_k+(conjugate(z)-1)e_l] rho_kl,  k != l.
```

Its four outcome maps are the actual source Gram maps. This supplies a
complete finite instrument for the declared squared-eigenvalue rates.
The native-pair register's read effect has trace `2/9`, with the same shape;
the double-load rates therefore differ by the common factor `5/4`.
Their ratios coincide. These are distinct declared preparations; equality of
their rate ratios does not identify their complete source laws. No physical
pole-mass or GeV identification follows
from that common shape or from this clocked apparatus.

## Repeated cycles, references and controlled limits

Because every source pair is reset by a fresh-state SWAP, the same complete
channel acts at each cycle. Its `N`-cycle entrywise multiplier is the `N`th
power of its one-cycle multiplier. Explicit finite histories retain all
`4^N` words; summing their unnormalized outcome maps reproduces that power.
For `h=t/N`, the weak channel is an average of unitaries with total angle
`h B_N`, where `B_N` is binomial with mean `Np` and variance `Np(1-p)`.
Its limit is conjugation by `exp(-itpJ)`.

For a Hermitian generator of spectral diameter `D`, the commutator
superoperator has diamond norm at most `D`. Taylor's integral remainder
for a unitary channel is at most `D² s²/2` after its linear term. Centering
the binomial angle cancels that linear term. With `D(J)=2`,

```
(1/2)||Phi_(t/N)^N - Ad(exp(-itpJ))||_diamond
    <= p(1-p)t²/N = 14 t²/(81N).
```

For the mass source, the full pulled-back analyzer is a projector `A` on
source and register, with diameter at most one. Its compression on `S` is
`E`, with diameter at most `tr E=5/18`. The first-order terms of the reduced
source channel and `Ad(exp(-ihE))` agree. Applying the two Taylor remainder
bounds, then the finite-channel telescoping inequality, gives

```
(1/2)||Psi_(t/N)^N - Ad(exp(-itE))||_diamond
    <= [1+(5/18)²]t²/(4N) = 349 t²/(1296N).
```

These inequalities cover every register state and untouched reference.
The exact identities and numerical native attachment are distinct from the
analytic norm argument; floating-point sample agreement is not used to
prove the bounds. At `t=1`, errors below `10^-6` require, conservatively,
172,840 ideal weak cycles or 269,291 ideal mass cycles. If each implemented
cycle has complete-channel error `epsilon` and fresh-pair trace error `eta`,
the error grows by at most `N(epsilon+eta)`. The preparation bound includes
independence from the continuing register and reference: a source-marginal
error alone does not bound correlated preparation. The receipt gives a finite design
splitting a target equally between truncation and supplied service accuracy.
Hardware accuracy is not inferred from these allowances.

Every diagonal register density is stationary for both finite channels;
the uniform density is a faithful example. For `h=1/N` and the declared
nonzero rates, only diagonal entries are fixed. There is no unique reference
selection. The register/last-record marginal can be stationary under that
diagonal input, while the exported archive continues to grow. Neither this
family nor its weak limit establishes a physical Gibbs state or vacuum.
At finite `N` the coherences contract; the exact channel is not unitary.
Even its phase per current difference is not common to gaps one and two.

## Identifiability and residuals

Write the common register shape as `f_k=1+c cos(theta+2 pi k/3)`, with
`c=2b/a`. Ratios of the squared eigenvalues are independent of the common
normalization. On the positive-eigenvalue chamber,

```
Koide = 1/3+c²/6.
```

At `theta=2/9`, `c=sqrt(2)`, the Jacobian of
`(log(mu/e),log(tau/e))` with respect to `(theta,c)` is approximately

```
[[56.1438827453, 32.6113343121],
 [51.2267935781, 34.4543952961]],
```

with determinant `263.829438446`. The two shape parameters are locally
identifiable from the two ratios; the common scale is not. This rank
calculation does not license choosing those parameters using target masses
and then calling the same comparison a prediction.

The exclusive phase is `p_F+p_H-2p_FH`. Fixed marginals `1/12,1/6` allow
`p_FH` anywhere in `[0,1/12]`, giving phases `[1/12,1/4]`. Independence of
the fresh cells selects `p_FH=1/72` and phase `2/9`. Thus the reset and
preparation map supply a premise that marginal counts alone do not supply.

The comparison-only rounded reference `206.768283` gives a muon/electron
residual of `+9.83213044 ppm`. The linear phase displacement diagnostic is
`-1.75122945e-7`; it is not applied. Finite-round diagnostic phases and ratio
biases are reported for fixed round counts. Selecting a round count to
cancel the measured residual would be calibration. These calculations
supply neither an uncertainty for the physical theory nor a mass-pole
interpretation of the register rates.
