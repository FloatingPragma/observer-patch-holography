# Independent algebra review

These are conditional circuit identities, not new physical observations or
an implementation certificate. They assume the declared coherent native
dilation/inverse and memory gates. The actual finite record, precision bound,
calibration procedure, complete costs and common-source execution remain
separate obligations.

Let Π project onto the known-success workspace and let U be a Hermitian
involution with success block P=ΠUΠ. The current-derived reversal is a
permutation/involution on the admitted one-unit position/coin domain.
Conjugating it by any declared choice unitary preserves those properties.
Thus the finite weighted source U_P is an appropriate U, including its
inactive choice labels.

Write A=(I+U_P)/2 and B=(I−U_P)/2. The complete contrast circuit is

```
X H diag(I,U_P) H = [[B,A],[A,B]].
```

It is Hermitian and squares to the identity because A²+B²=I and AB=0.
Its success code block is (I−P)/2. The construction does not discard the
other ancilla or choice outputs.

For S_z=I+(z−1)Π, |z|=1, the full amplified circuit U S_z U S_z U has
successful block

```
F_z(P) = (2z−1)P + (z−1)²P³.
```

The preliminary π circuit therefore has block 4P³−3P and is itself a
Hermitian involution on the complete dilation space. Denote its block by
T, or T=P when no preliminary boost is used. Applying the complete forward
instrument, selecting zero, then the complete adjoint instrument and
selecting zero gives F_z(T)†F_z(T). For α=2(1−Re z), direct multiplication
using z conjugate(z)=1 gives

```
K(T) = I − (I−T²)(I−αT²)².
```

This identity holds for the actual finite rational phase; it does not require
exact phase matching. Exact matching α=1/s² makes K(s)=1. Nearest rational
tangent phase generally does not, so exact-dark assertions based on K(s)=1
must remain ideal statements. K is positive because it is F†F, and is a
contraction because it is a selected Kraus product of complete instruments.
Any final negative closed-path interference comes from the ordered product
of these source-dependent maps, not a final supplied minus sign.

Completeness must be applied to the two sequential instruments, not just K.
If {F_j} is the forward instrument and {G_j} the adjoint instrument, with
G_0=F_0†, the full first-failure instrument is

```
{F_j : j≠0}, {G_j F_0 : j≠0}, {F_0†F_0}.
```

Its effects sum to I. Repeating filters or composing arms preserves this
only when all first-failure branches, shared workspace flags and conditional
probabilities are retained. Normalizing after a selected branch does not
replace its success cost.

For the five-of-eight selector H S_0(z) H S_good(z) H, its first-column
amplitudes are

```
active   (5z²+6z−3)/(8√8), on each of five labels;
inactive (5z²−2z+5)/(8√8), on each of three labels.
```

Their exact total weights are a=5|5z²+6z−3|²/512 and
b=3|5z²−2z+5|²/512, with a+b=1. The source reversal is identity on the
inactive labels, so the actual mean block is aP+bI. The selector's phase
and its declared memory gates therefore contribute to the calibration input;
the native current columns alone do not determine a. At ideal Re z=1/5,
b=0; at a finite nearby phase, b remains in the complete instrument.

For the always-swap star mean, its restriction S to the whole optical plane
is traceless with S²=uI. The corresponding finite filter is the exact positive
scalar

```
k = 1 − (1−u)(1−αu)²,
```

replacing u by u(4u−3)² when the preliminary boost is used. If the moving and
stationary arm products are aI and bI with a,b>0, the common subsequent
re-encoding map acts identically on the same actual input density. Hence

```
p± = (a±b)² / [2(a²+b²)].
```

The minus port is exactly dark only if a=b. At finite frame-dependent
calibration phases, retain its actual small nonzero value; both laws must
use the same parity-selected port. This scalar argument works on any actual
optical mixed state and its reference correlations, without replacing it by
an ideal eigenray.
