# Finite response and complete observations

The results concern declared native instruments and explicitly supplied finite
field attachments. A density transfer eigenvalue, an amplitude frequency and a
physical rest mass are distinct objects. The executable certificate uses exact
rational or algebraic identities for its mathematical checks and labels its
floating controls separately.

## Comparing capacities without comparing every history

Let `A_y` and `B_y` be nonnegative matrices for every retained record `y` of two
finite classical instruments. Their sums are column stochastic. Let `C` be a
nonnegative column-stochastic map from the larger state space to the smaller
one. Suppose, separately for every outcome,

```
C A_y = B_y C,                 C p = q.
```

For a word `w=(y1,...,yn)`, induction gives

```
C A_yn ... A_y1 p = B_yn ... B_y1 q.
```

Multiplication by the all-ones row proves equality of the full joint record
probabilities at every finite length. Any identical estimator applied to these
records, including a specified spectrum estimator, therefore has the same
distribution. Deterministic systems are included as zero-one transition
matrices. Determinism alone supplies neither the intertwining nor the compatible
initialization or record map.

The condition can be checked on a finite list of matrices instead of an
exponential list of histories. The code retains the outcome alphabet and
checks normalization before admitting it. The same induction applies to
quantum instruments if `C` is a trace-preserving quantum channel and the
outcome maps intertwine as superoperators. The rational executable certificate
implements the classical case; it does not assert a quantum realization for
the particle source.

For an approximate comparison put

```
epsilon0 = ||C p-q||_1,
delta = max_j sum_(y,i) |(C A_y-B_y C)_ij|.
```

The combined record-and-state maps are contractions in the `l1` norm. At each
step their comparison defect is at most `delta` on a normalized nonnegative
state. A telescoping sum, retaining every record, bounds the final joint
state-record distance by `epsilon0+n delta`. Marginalization contracts that
distance, so the record total-variation distance is at most
`min(1,(epsilon0+n delta)/2)`. The exact case has zero error for every `n`.
This is a sufficient comparison criterion, not a claim that a capacity map
for the physical source has been selected.

Equality of unread channels is insufficient. The one-state instruments
`A_0=1,A_1=0` and `B_0=0,B_1=1` have the same unread channel and deterministic
evolution. Their one-step records are disjoint. The certificate also checks a
two-state alternating record process against a four-state refinement whose
hidden bit flips, including all 127 words of length at most six.

## Native neutral memory in four dimensions

`model/common_source/native_family_record_memory.py` enumerates all 128 native
record labels. Their C and E probabilities are `1/12` and `1/6`, with joint
probability `1/72`. The source projectors at memory addresses 2 and 3 have
rank one and overlap `Tr(P2 P3)=1/5`. Orthonormal coordinates on their bright
plane give

```
P = diag(1,0),              Q = [[1,2],[2,4]]/5.
```

Both annihilate the orthogonal family line. Under the declared neutral current
gate, the E probability is zero and the C probability is
`p=sin(theta sqrt(9/8))^2/12`. The recurrent transition is
`[[1-p,1],[p,0]]`, with stationary column `(1,p)/(1+p)`. Its average bright
projector `(P+pQ)/(1+p)` has trace one and determinant
`4p/(5(1+p)^2)`. Thus its complete static spectrum follows from a quadratic;
a numerical three-family search is unnecessary.

For the explicitly unshifted amplitude extension put `z=exp(-ih)`,
`U=I+(z-1)P`, `V=I+(z-1)Q`. Its recurrent bright transfer is

```
F = [[(1-p)U, U], [pV, 0]].
```

Direct symbolic determinant evaluation gives

```
f(lambda) = lambda^4 -(1-p)(1+z)lambda^3
 + ((1-p)^2 z-p(2z+(z-1)^2/5))lambda^2
 + p(1-p)z(1+z)lambda + p^2 z^2.
```

At `h=0` this is `(lambda-1)^2(lambda+p)^2`. The retained memory roots are
part of the response. On the dark line the recurrent roots are `1,-p`.
The other ten memory inputs of the native controller are retained in its full
instrument: every step staying outside addresses 2 and 3 requires a C event,
so their transient amplitude block has induced block-sum norm at most `p`.
The recurrent calculation therefore conceals no other roots of modulus above
`1/2` in the domain below.

For `0<=p<=1/12`, `0<h<=1/2`, set `x0=h/2` and
`lambda=exp(-ix0)x`. The quartic reduces using `y=x-p/x` to

```
y=(1-p)cos(x0) +/- i sqrt(C0)sin(x0),
C0=(1-p)^2+4p/5.
```

The slow quadratic root is `xs=(y+sqrt(y^2+4p))/2`, with square root of positive
real part; the fast root is `-p/xs`. For the positive imaginary `y`, write
`xs=r exp(i alpha)`. Then

```
(r-p/r)cos(alpha)=(1-p)cos(x0),
(r+p/r)sin(alpha)=sqrt(C0)sin(x0).
```

The bounds `r >= (1-p)cos(x0) >= 341/384 > 1/2` and
`p/r <= 32/341 < 1/2` separate the two branches. Write
`v=sin(x0)^2`, `u=sin(alpha)^2`. Elimination gives

```
4p u^2-B0 u+C0 v=0,       B0=(1+p)^2+(4p/5)v.
```

On `0<=u<=1` its derivative is at most `-73/144`; at `u=0` it is nonnegative,
and at `u=v` it is `(16p/5)v(v-1)<=0`. Hence `0<=u<=v`,
`0<=alpha<=x0`. Also `r<=1`: if `r>1`, then `r-p/r>1-p`, forcing
`cos(alpha)<cos(x0)` and contradicting `alpha<=x0`. For `p>0,h>0` these
inequalities are strict. The equal-modulus slow roots have positive phase
coefficients `omega=1/2-alpha/h` and `1-omega`.

Since `v<=1/16`, `C0>=653/720`, `B0<=53/45` and
`B0-4pu >= (1-p)^2>0`,

```
u >= (653/848)v,
alpha/x0 >= (95/96)sqrt(653/848),
omega <= a = (1-(95/96)sqrt(653/848))/2
            = 0.06580877602554946...
```

The sine estimate uses `sin(x0)/x0 >= 1-x0^2/6 >= 95/96`. Under the added
normal-order attachment `(m1,m2,m3)=M(0,omega,1-omega)`, `M>0`,

```
(m2^2-m1^2)/(m3^2-m2^2) <= a^2/(1-2a)
                         = 0.0049871977631631886...
```

For the inverted assignment `m3=0,m1=M omega,m2=M(1-omega)`, the adjacent
ratio `(m2^2-m1^2)/abs(m3^2-m2^2)` is at least
`1-[a/(1-a)]^2`. No arbitrary-pulse or shifted-reference bound is asserted.
All current angles enter only through `p<=1/12`; the feedback-pulse bound is
an independent hypothesis. This is a conditional diagnostic exclusion, not a
physical neutrino mass theorem.

## State, detector and field reference

Let `J=[[0,1],[-1,0]]` on the bright plane and zero on the dark line. For either
projector, `Pj J+J Pj=J`. With X-eigenvectors `|+>,|->`, define
`A=|+><-| tensor J`, `O=A+A*`, and bright/dark projectors `PB,PD` on the
six-dimensional response. Then `O^2=PB` and the three effects
`(PB+O)/2,(PB-O)/2,PD` are a complete projective readout with values `+1,-1,0`.
Each feedback pulse conjugates `A` to `exp(-ih)A`. Sequential measurement,
including the first measurement's backaction and all nine joint outcomes,
gives respectively

```
rho=PB/4:       C(n)=cos(nh),
rho=PD/2:       C(n)=0,
rho=I6/6:      C(n)=(2/3)cos(nh).
```

All three states are stationary with the stated stationary memory. Thus a
transfer pole can have nonzero or zero recorded residue depending on the
admitted preparation. Stationarity alone does not select that preparation.

Replacing every pulse `Uj` by `exp(-i delta h)Uj` changes no original density
history, including an arbitrary spectator. With an additional reference
coordinate and block unitary `diag(1,exp(-i delta h)Uj)`, however,
reference-to-response coherences acquire that phase. Absolute amplitude
frequencies are therefore not selected by the original density instrument.
The extension supplies no Lorentz-scalar mass or physical vacuum.

For `G(z)=a(I-zF)^(-1)b`, a pole is visible only if its residue survives `a,b`.
An invertible coordinate change of the full response cannot move its poles:
`Gnew=L G R` has the same pole order and leading-residue rank. Analytic factors
invertible at a pole have the same local property. A change of detector or
reference is a different physical attachment.

## Conserved passive zero and scalar projection

The native passive Schur calculation supplies a family frame `B`,
`B^T B=G I`, and harmonic extension `U` satisfying `B^T U=G I`. Therefore

```
Z=U^T U/G = I+(U-B)^T(U-B)/G >= I.
```

The native reference stationarity identity is
`sigma^T(I-T)sigma=(157/90000) Wcross`, where `Wcross` sums nonnegative
core-to-outer seam weights. Stationarity excludes all eight such seams.
Every remaining recorded branch conserves core charge. The centered vector
`c=(5,5,-1,...,-1)` has nonzero family coordinate
`x=B^T c/G=(0,3sqrt(5)/5,0)`. Each admitted seam annihilates `c`.
Harmonic elimination therefore gives `S0 x=0`; the maximal 22-seam mixture
attains rank two. Positive normalization is an invertible congruence and
preserves rank. It cannot generate the third charged family coefficient.
This argument covers the nonnegative stationary-compatible seam cone, not
active repair or a replacement reference.

With the declared spatial symbol `alpha_i=Z2 tensor sigma_i tensor I`, a
spin-invariant Hermitian term has auxiliary factors `I,X,Y,Z`. Exactly its
`X,Y` parts anticommute with all spatial `alpha_i`; its scalar projection is
`(B-alpha_1 B alpha_1)/2`. A family potential `I4 tensor A` has zero scalar
projection. For `H_h(k)=alpha.k+X2 tensor I2 tensor M+h I4 tensor A`, direct
expansion and spin traces give

```
Tr[(H_h(k)^2-|k|^2 I)^2]-Tr[H_h(0)^4]
    =16 h^2 |k|^2 Tr(A^2).
```

A nonzero coherent family potential is therefore not a collection of constant
relativistic masses for this truncated Hamiltonian. Higher-order kinetic or
scalar terms in a different complete effective action are outside this test.

## Finite translations and quantum links

For finite matrices, cyclic trace gives

```
i Tr(H^2)=Tr(K[P,H])-Tr(([K,P]-iH)H).
```

If `H` is Hermitian and `[P,H]=0`, `[K,P]=iH`, then `Tr(H^2)=0`, hence `H=0`.
No Hermiticity premise on `K,P` is needed. The Hilbert-Schmidt inequality also
gives `||H||_F^2 <= ||K||_F ||[P,H]||_F+||[K,P]-iH||_F ||H||_F`.
This constrains an exact global finite representation. It does not constrain
a controlled low-energy sector or an infinite limit without their specified
approximation maps and errors.

An operator-valued one-particle projector on matter times a quantum link need
not define integer fermion occupation. For
`P=|Phi><Phi|`, `|Phi>=(|0,0>+|1,1>)/sqrt(2)`, its entries are
`Pij=|i><j|/2`. The naive lift `N=sum c_i* Pij c_j` on the complete
four-state matter Fock space times two-state link has eigenvalues
`0 (multiplicity 5), 1 (multiplicity 1), 1/2 (multiplicity 2)`.
The certificate computes this exactly using the full CAR matrices.
Gauge covariance and one-particle spectral isolation do not by themselves
supply a canonical fermion factor with quantum-link coefficients.

## Connections to the finite theorem stack

The [native passive package](../common_source/full_repair_response/README.md)
supplies the harmonic Gram identity and seam-cone conservation certificate.
Reusing those objects replaces charge-by-charge kinetic searches. The neutral
bright/dark decomposition replaces large eigenvalue searches by a
four-dimensional calculation while keeping transient-memory bounds.

`Lean/SourceLaw/RepairAdmissibility.lean` proves its time-coupling locality
selection under a Euclidean history-register premise, including admissibility
at every sufficiently small time step. It does not attach the neutral family
metric to a Lorentzian kinetic action. `Lean/QFT/SourceHistoryGNSDynamics.lean`
constructs a selected finite history generator and its represented dynamics;
its stated scope supplies no physical time-slice theorem. The comparison
criterion above identifies the further joint-record maps needed to transfer
such a construction between capacities. No additional Lean formalization of
this package is asserted.

A scalar fixed algebra of a particular autonomous channel excludes nontrivial
exactly invariant observables there. It does not exclude transient responses,
metastable modes, a changed seam set, record-assisted recovery, or another
occupation domain. In particular the full-occupation native model has the
empty/full coherence block; a scalar result in the fixed-two-occupation
sector does not cover that block. These scope distinctions prevent using a
finite fixed-algebra calculation as a universal particle exclusion.
