# Finite regional entropy of the literal prepared source

## Source, algebra and retained outcomes

Use the open four-site Weyl/contact plaquette, with eight signed
canonical anticommutation modes, speed `c=1/8`, contact coefficient `U=1`
and inverse temperature `beta=1/32`. The full field has dimension 256;
the region is the complete matrix algebra on modes `(0,1,2,3)`, of
dimension 16. No occupation sector is removed.

Let `Sigma` be the literal mathematical source preparation's success
field marginal, after summing its 256 classical column records. Write
`sigma=Tr_out(Sigma)` using the signed Fock factorization. These blocks
are subnormalized: `Tr Sigma=Tr sigma=p<=1`. The source retains every
classical record and abort outcome; summing records to define a marginal
does not assert that the recorded source is pure or that a conditional
column is thermal.

The [literal control certificate](../certified_reference/README.md)
proves a complete preparation error `e<17/10^13`. The independent
[polynomial design theorem](../plaquette_reference/preparation.py)
supplies a rational density error `d` and spectral floors `m0` regionally
and `m0/16` globally. Put

```
d_exec=(2+e)e,
m=m0-d_exec,       m_full=m0/16-d_exec,
E=R8+(d+d_exec)/m,
```

where `R8` is the rational eighth-order regional logarithm remainder.
The exact checker recomputes and compares every value with the literal
source certificate, checks that certificate's source hashes, and proves

```
sigma >= m I,      Sigma >= m_full I,
m>1/20,           16 m_full>1/20,          E<2e-9.       (1)
```

For the source's regional polynomial `P8` and its Gibbs normalization
scalar `a`, the inherited operator theorem is

```
||-log sigma-P8-a I|| <= E.                             (2)
```

The scalar is kept in the absolute comparison. Its contribution vanishes
only for trace-zero variations. The simulator law, conditional scheduler,
coherent load/export and record memory are the specified source services.
Their physical calibration is not a consequence of (1).

## Every small regional variation has a positive full-field marginal lift

Let `X=X*`, `Tr X=0` and `r=||X||<=1/40`, where the unlabelled norm is
the operator norm. For every `|t|<=1`,

```
sigma+tX >= (m-r)I > I/40.                             (3)
```

Thus this is an actual two-sided positive family at fixed trace `p`.
Let `V` be the unitary signed Fock factorization into the region and its
sixteen-dimensional complement. Define

```
Sigma_t=Sigma+t V* (X tensor I16/16) V.
```

Its signed regional trace is `sigma+tX`, its trace is `p`, and

```
Sigma_t >= (m_full-r/16)I > I/640.                      (4)
```

The stronger global floor in (1) is checked explicitly; positivity of a
regional marginal alone would not imply (4). The aggregate success/abort
state `Sigma_t direct_sum A_abort` is normalized when the original
aggregate output was normalized and keeps its entire abort block fixed.
No entropy of a success-conditioned state is used.

Equations (3)--(4) concern field marginals and the aggregate abort flag.
They do not specify a positive lift preserving every original
record-conditioned field block, or execute a native intervention. The
original column states can have kernels even though their sum is faithful.
An arbitrary perturbation of that sum cannot be assigned to each column
without an additional construction. All statements below are regional
entropy statements, not joint field-and-record entropy statements.

## Entropy identity and noncommuting quadratic control

For a positive finite matrix `A`, define `S(A)=-Tr(A log A)`, including
subnormalized matrices. Let `K=-log sigma` and

```
D(sigma+X || sigma)=Tr((sigma+X)(log(sigma+X)-log sigma)).
```

Equal traces give the exact identity

```
S(sigma+X)-S(sigma)=Tr(X K)-D(sigma+X || sigma).         (5)
```

For `A>0` the resolvent representation of the logarithm gives

```
Dlog_A[X]=integral_0^infinity (A+sI)^(-1) X (A+sI)^(-1) ds.
```

This follows by differentiating
`log A=integral_0^infinity[(1+s)^(-1)I-(A+sI)^(-1)]ds`.
A positive lower spectral bound makes the derivative integrals locally
uniformly convergent; no entrywise logarithm or commutation is assumed.
The entropy Hessian is the negative of the positive quadratic form

```
g_A(X,X)=Tr(X Dlog_A[X]).
```

In an eigenbasis of `A`, its coefficient on `|X_ij|^2` is
`integral_0^infinity[(lambda_i+s)(lambda_j+s)]^(-1)ds`.
If `l I<=A<=M I`, each coefficient is between `1/M` and `1/l`.
Every matrix `sigma+tX`, `0<=t<=1`, has trace `p<=1`, so its largest
eigenvalue is at most one and its least eigenvalue is at least `m-r`.
Twice integrating the entropy derivative proves

```
D(sigma+X || sigma)=integral_0^1 (1-t) g_(sigma+tX)(X,X) dt,
(1/2)||X||_HS^2 <= D <= ||X||_HS^2/[2(m-r)] <=20||X||_HS^2. (6)
```

In particular, positivity of relative entropy follows here from the
integral formula, including at fixed subnormalized trace. It is strictly
positive for every nonzero admitted `X`. The same fixed abort complement
contributes zero to entropy differences; its entropy is not silently
discarded from a varying joint state.

Combining (2), trace duality and (5)--(6) gives the finite polynomial
first-law estimate

```
|S(sigma+X)-S(sigma)-Tr(X P8)|
 <= E ||X||_1 + ||X||_HS^2/[2(m-r)].                    (7)
```

Thus the polynomial's linear error and the actual nonlinear entropy
remainder are distinct and uniform over the complete variation ball.
This is more than a derivative evaluated at one state or one direction.

## Cubic control and a symmetric finite estimator

For `f(t)=S(sigma+t Delta)`, differentiation of the same resolvents gives

```
f'''(t)=2 integral_0^infinity
 Tr(Delta R Delta R Delta R) ds,
R=(sigma+t Delta+sI)^(-1).
```

The integrand is real since it is the trace of
`(R^(1/2) Delta R^(1/2))^3`. Schatten norm inequalities give
`|Tr(Delta R Delta R Delta R)| <= ||Delta|| ||Delta||_HS^2 ||R||^3`.
Consequently on `|t|<=h`, with `l=m-h||Delta||>0`,

```
|f'''(t)| <= ||Delta|| ||Delta||_HS^2/l^2.              (8)
```

Taylor's integral remainder yields both

```
|D(sigma+h Delta || sigma)-(h^2/2)g_sigma(Delta,Delta)|
 <= h^3 ||Delta|| ||Delta||_HS^2/(6l^2),

|[S(sigma+h Delta)-S(sigma-h Delta)]/(2h)-Tr(Delta P8)|
 <= E||Delta||_1+h^2||Delta|| ||Delta||_HS^2/(6l^2),     (9)
```

for `h>0`; at zero the estimator means the derivative. The cancellation
in the symmetric estimator is controlled by (8), not inferred from an
unquantified first law. The calculation does not bound sampling noise or
the cost of measuring the two entropies.

## All 255 directions on the same state

For `0<=i<15` take `D_i=E_ii-E_15,15`. For every `0<=i<j<16`, take
`R_ij=E_ij+E_ji` and `J_ij=i(E_ij-E_ji)`. These 15+120+120 Hermitian
trace-zero matrices each have operator norm one, trace norm two, and
Hilbert--Schmidt norm squared two. Their diagonal Gram block is
`I15+11^T`, of determinant 16. Their other Gram block is `2 I240` and
the cross blocks vanish. Therefore the full Gram determinant is
`16*2^240>0`: they span every Hermitian trace-zero regional variation,
including number-changing and noncommuting ones.

The full matrix algebra is a mathematical domain. Fermionic parity
restrictions are explicit: the parity-preserving directions consist of
all 15 diagonal directions and both off-diagonal directions for the 56
pairs of equal occupation parity, giving dimension 127. They span the
trace-zero Hermitian part of `M8 direct_sum M8`. If number preservation is
also imposed, the 27 pairs with equal particle count give dimension 69,
the trace-zero part of `M1 direct_sum M4 direct_sum M6 direct_sum M4
direct_sum M1`. Every bound holds on these subspaces without change. The
other 128 directions are not asserted to be physical even interventions.
In the signed factorization, even `X` gives an even added full-field
operator in (4). An `X` commuting with regional number gives an added
operator commuting with full number. These statements do not impose
parity on the literal reference or fix the changed mean particle number.

For real coordinates `theta` with `||theta||_1<=1/40`, the triangle
inequality places `X=sum theta_j Delta_j` inside (3). This defines one
common finite-dimensional positive variation family. There is no claim
that all 255 coordinates are independent spacetime stress components.
For any individual basis direction, (9) gives the simple bound

```
|symmetric entropy estimator-Tr(Delta_j P8)|
 < 4e-9+(1600/3)h^2,       0<h<=1/40.                 (10)
```

The receipt also reports the sharper exact rational bounds from (1),
rather than only the simplified constants in (6) and (10).

## Gravity scope

The same prepared interacting finite reference has a faithful regional
algebra, a nondegenerate entropy Hessian, a full two-sided variation
space, a positive full-field marginal lift, and controlled finite
first-law remainders. The modular parameter is dimensionless and belongs
to that reference. No spatial or small-ball refinement is constructed.
The region is a source site region, not an identified geometric diamond.
This proves no equality with a conserved geometric stress charge, no
edge-area normalization, no positive Newton constant and no joint
fixed-proper-volume generalized-entropy stationarity. Those identifications
cannot be inferred from the positivity of (6).
