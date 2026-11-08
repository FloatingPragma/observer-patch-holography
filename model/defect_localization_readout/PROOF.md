# Local creation, a conserved symmetry sector and a measurable lifetime

Use the two-flavor centered source and its protected single-occupancy
carrier from `protected_spin_limit`. The source Hamiltonian is
`K=T+UD-(U/2)N`; in the signed hidden-spin frame its hopping is scalar
in the two internal states. The complete C/E extension is exactly
`defect_field_bridge/full_domain.py`, including collision holds, both
coins and every raw record. These are the specified field-extension
interfaces, not an interpretation of the bare one-bit carrier.

## 1. A symmetry preserved in every recorded repair branch

Let `R=Gamma(I_sites tensor X)` be the second-quantized common internal
Pauli X. On an onsite Fock basis `(empty,up,down,double)` it fixes the empty
state, exchanges up and down, and multiplies double by minus one. Consequently
`R*=R`, `R^2=I`. The minus sign is the fermion determinant, not an optional
phase convention on double occupation.

Both scalar hopping and `UD-(U/2)N` commute with R. The local staggered
selection term proportional to `S_x^1` also commutes with R. Each full
repair Kraus operator obeys

```
[M_(coin,raw),R]=0.                                        (1)
```

A hold is a function of the two endpoint occupation counts. A recorded
move on its one-particle/empty domain is the sum of identical-flavor CAR
hops; simultaneous internal X at the endpoints preserves that sum. The
record depends only on the counts and current. This proves (1) on all
Fock sectors, including separated endpoints and arbitrary spectator
references. The executable check verifies every full-domain column on
two, three and four sites, every seam, both laws and all 32 coin/raw slots.
It also contracts R with the actual full-Fock centered dimer Hamiltonian.

On the protected carrier P there is one fermion at each site and

```
M_(coin,raw) P = delta_(raw,0) P/sqrt(2).                   (2)
```

The code verifies all 12,288 protected cube columns. Thus any finite
complete repair history, even with seams and laws adapted to earlier
records, leaves its spin state unchanged and supplies independent fair
coins. No branch is discarded. On the full source (1) preserves the
eigenvalue of R even when Hamiltonian evolution leaves P. Ideal repairs
alone therefore protect the encoded spin information indefinitely.

This symmetry is not occupation excess `Q=N-|V|`: both created and
reference spin states have Q=0. It is not a winding number. Coherent
local spin operations are admitted controls that can change R, so there
is no invariance under the entire control algebra. An infinite product
of X is not assumed to exist in the infinite-volume representation.

An extension to a charge-fluctuating local input must retain failure. The
spin-odd branch `Z_s=n_up-n_down` and `F=I-Z_s^2` obey
`Z_s^*Z_s+F^*F=I`, `R Z_s=-Z_s R`, and `RF=FR`. Its failure is precisely
the empty/double subspace; it cannot be discarded or replaced by a
deterministic parity flip. For the inherited low-band isometry W with
`||W-P||<=2tau/U`, its failure probability is at most `4(tau/U)^2` on
every normalized encoded input. This follows from `FP=0`, without a
thermal or ground-state assumption. The deterministic experiment below
uses P. Conditioning the general instrument on success can change
exterior correlations, so (3) is not asserted for that conditioned input.

## 2. Local creation on the actual prepared cube reference

Take the open eight-site cube, J=1 and B=1/4, with the actual selection
Hamiltonian

```
H_B=sum_<xy>(S_x.S_y-1/4)-(1/4)sum_x epsilon_x S_x^1.
```

`native_ordered_ground` proves its unique ground g, complete finite gap
and finite native preparation. In its Marshall coordinates H_B commutes
with the bit-complement permutation. Positivity of the unique Perron
vector therefore gives Rg=g. The Marshall transformation has four odd
sites, so it commutes with this global parity; the assertion is unchanged
in the encoded physical spin coordinates.

Apply the one-site logical unitary `U=Z_0`. It preserves P and satisfies
`R U P=-U R P`. The state `U g` is therefore an exactly odd response,
orthogonal to the even reference, on this same Hilbert space. For every
exterior observable A supported away from site 0,

```
<Ug,A Ug>=<g,A g>.                                         (3)
```

The even onsite CAR extension of U commutes with the exterior algebra;
equivalently its phases cancel in the partial trace. This is exact initial
localization relative to g, including all exterior correlations. It does
not assert that the excited state factorizes from the exterior.

The reference X pointer at site 0 has effect `(I+X_0)/2`. Since
`U X_0 U=-X_0`, its probability contrast is `<g,X_0 g>`. Exact rational
evaluation of the retained ground proposal, with the independently
verified vector error, gives

```
pointer probability contrast      [0.670950141860,0.670950141862],
<Ug,H_B Ug>-<g,H_B g>              [1.825861195059,1.825861195062],
<Ug,H_0 Ug>-<g,H_0 g>              [1.658123659594,1.658123659597].
```

Here `H_0=H_B+(1/4)sum epsilon_x S_x^1`, evaluated on the same input g.
The finite-field gap is greater than 1.092284677676. Since U g is
orthogonal to g, all its H_B spectral weight lies above that gap.
The H_0 comparison does not call g an H_0 ground or stationary state.
Neither selection energy nor energy contrast is a relativistic rest mass.

These intervals use `||Z H_B Z-H_B||<=3J+B` and
`||Z H_0 Z-H_0||<=3J`. For normalized vectors separated by epsilon, the
expectation error of A is at most `2||A||epsilon`. All 256 proposed
eigenvectors enter the parent ground certificate; its eigensolver is not
an assumed spectral oracle.

## 3. Actual complete local creation and readout

The native encoding has `Zbar=Z_up`. The existing phase-exact two-pulse
word `I(-pi/2),Z(pi/2)` on the up register implements this logical Z.
Both ideal operations preserve single occupancy. The finite program
inserts these pulses immediately before the prepared ground's complete
25-pulse X-pointer word. The reference program has 1,090 pulses and the
created program 1,092. Both execute on the full 17-bit, 131,072-row work
register through the actual six-level source.

The source coefficient bounds and exact angle bounds give complete
isometry errors below `1/10^8` for each program. Every first-failure label,
both binary pointer values and reserved output slots are retained. Failure
does not lead to normalization of the surviving histogram. The predicted
postmeasurement target exchanges the two pointer branches and applies
Z to the field, so the finite check also tests back-action.

The full preparations and local detectors are executed. The complete
recorded lifetime estimate below is a separate analytic construction;
its ten-repair tape is not claimed to have run. An approximate state or
compiled intervention contributes its full half-trace error to that bound.
Coherent memory, load/export and refinable source control are declared
interfaces. Bare repair records are identical on g and Ug by (2) and
cannot replace the spin pointer. Repairs inserted between constituent
unprotected control pulses are outside the stated schedule.

## 4. Quantified signal and spatial spreading on the same clock

For the Heisenberg bond `h_xy=J(S_x.S_y-1/4)`, direct Pauli algebra gives
`||[h_xy,X_x]||=J`. Thus for a vertex of degree d,
`||[H_0,X_x]||<=dJ`. Let both g and Ug evolve with the same H_0; no
stationarity assumption is made. Each pointer probability changes at
rate at most dJ/2, so their signed probability contrast obeys

```
d_X(t)>=d_X(0)-dJ|t|.                                     (4)
```

Arbitrary complete repairs inserted in this ideal spin history leave
(4) unchanged. To compare with the actual centered CAR dynamics and
its complete recorded repairs, use the all-input/reference estimate
proved in `protected_spin_limit`, for m free intervals and total duration T:

```
D=8m tau/U+6T tau^3/U^2,       tau=c|E|, U>=4tau,
||rho_source,cq-rho_spin,cq||_1<=min(2,D).                  (5)
```

It applies separately to g and Ug, on the same initial encoding and
clock. A binary event error is at most half the unhalved trace error.
If each complete experimental protocol has an additional half-trace
error at most eta, subtraction of the two probabilities yields

```
d_X,source(T)>=d_X(0)-dJT-D-2eta.                          (6)
```

No sum over accepted branches replaces the full record in (5). On the
cube trajectory `c=12^4`, `U=12^8`, J=1, take any at most ten repairs,
at most eleven free intervals, `T<=1/16`, and `eta<=1/10^8`. Since
`tau/U=12^-3` and `tau^3/U^2=1/12`, (6), using the conservative
`d_X(0)>67/100`, gives

```
d_X,source(T) > 540437473/1350000000 > 2/5.                (7)
```

This is a finite measurable signal lifetime on the specified source
clock, for every allowed seam and adapted C/E choice. It is not a
particle lifetime or a conservation law for fixed spatial support.

For an exterior contraction A at graph distance r from the creation
site, the same source's spin commutator-chain estimate gives

```
|<Ug,alpha_t(A)Ug>-<g,alpha_t(A)g>|
 <= sum_(n>=r)(24J|t|)^n/n!
 <= 2^-r (exp(48J|t|)-1).                                 (8)
```

This follows by commuting the initial local unitary through the evolved
observable, using the fixed-spacing spin bound in `protected_spin_limit`.
For an exterior effect the complete CAR comparison adds at most D for
the two protocols, and their implementation errors add at most 2eta.
For a general norm-one observable the corresponding trace errors are
unhalved. Equation (8) describes controlled spreading, rather than
permanent localization.

## 5. What this identifies

The construction provides one operationally created local excitation,
an exact finite symmetry label preserved by all repair outcomes and
the chosen Hamiltonian, positive finite response energy and an explicit
same-clock detector lifetime. It consumes the selected finite state and
declared coherent preparation and detection services. Source creation of
a topological defect, an infinite-volume charged sector, a relativistic
dispersion component and a measured particle assignment do not follow.
Its physical particle mass table is empty.
