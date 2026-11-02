# Joint native observations of finite spin energy and spatial symmetry

The selected spin ground has a complete joint energy and spatial-symmetry
instrument. It supplies an operational spectral readout of the actual state
prepared in `native_ordered_ground`, including every clock bin and native
failure outcome. The open eight-site cube admits **combined graph-reflection
characters**, not a small-momentum dispersion curve. Ordinary combined
translation momenta occur in the separate theorem for even periodic cubes.
Neither construction assumes a sharp mass shell.

## The exact symmetry and the prepared reference

Let

```
H_B = J sum_<xy> (S_x.S_y-1/4) - B sum_x eps_x S_x^1,
J>0, B>0.
```

On the open two-cube, R_j flips coordinate j of every vertex. It exchanges
the two sublattices. Let Z_all be the product of all physical Pauli Zs and
put g_j=R_j Z_all. Each g_j is a unitary involution, the g_j commute, and
each commutes with H_B: the reflection and Z_all each reverse the selected
field, while preserving exchange. In the Marshall frame U=product_(eps=-1) Z,
the exact identity is U g_j U=R_j. The source's real matrix in this frame
has negative hopping and transverse-field entries. The previously proved
Perron argument gives a unique positive ground. Permuting its coordinates
therefore fixes it, and consequently g_j Omega=Omega in the physical frame.
These graph actions permute the protected spin tensor factors; they do not
replace the binary cubic CAR rotation representation in FA95.

For dimension d, G=(Z/2Z)^d has projectors

```
Q_k = |G|^-1 sum_t (-1)^(k.t) g_t.
```

`algebra.py` verifies the full configuration permutations, multiplication,
Marshall signs and Hamiltonian commutation exactly; no rounded eigenvector
chooses the symmetry of the ground.

On an even periodic cubic torus of side L, replace R_j by the actual
one-site translation T_j. Then g_j=T_j Z_all commute, g_j^L=I, and the same
Marshall calculation gives U g_j U=T_j. Their characters lie on the usual
finite momentum grid 2*pi*n/L. This is a separate periodic graph and a
separate finite Hamiltonian. The open two-cube data below are not silently
assigned that periodic Hamiltonian or interpreted as p tending to zero.

## A complete positive energy–symmetry instrument

Start with a unit-norm local Pauli probe A applied to the ground. In the
retained example A is encoded Z at vertex zero. This preparation is unitary:
its norm is one, without a small-probability spin-filter postselection. Set
K=H_B-E_0>=0. For M clock values and time step tau>0 define

```
B_j(K) = M^-1 sum_(r=0)^(M-1) exp(2*pi*i*j*r/M) exp(-i*r*tau*K),
F_j(K) = B_j(K)^* B_j(K).
```

The actual joint effects on the probed state are Q_k F_j(K). They are
positive, commute, and sum over every k,j to I. To implement them, prepare
a uniform group register, perform its controlled g_t, and Fourier transform
that register. Prepare a uniform time register, perform controlled
exp(-i*r*tau*K), and apply the positive-sign Fourier transform. The joint
post-measurement vector is Q_k B_j(K) A Omega. Summing squared norms proves
normalization and the complete recorded instrument, not merely a formal
collection of expectation values. The two clocks commute in the ideal
instrument. A retained reference or subsequent operation does not alter
this dilation proof.

All energies are positive source excitations of the selected H_B reference.
The clock is periodic in tau*E with period 2*pi. If the circular phase
distance from a bin center is at least delta in (0,pi], its probability
kernel is at most 1/[M^2 sin(delta/2)^2]. A union of b such bins has leakage
at most b times this bound. This is a resolution/leakage bound, not an
assertion that a finite bin is an eigenvalue or a rest mass. The complete
FA98 spectrum certifies that tau=1/2 puts every cube excitation within one
phase period; all four bins are retained even when their centers exceed
the spectrum. An arbitrary larger tau must retain its possible aliases.

## Certified cube probabilities without a degenerate-eigenvector oracle

The unchanged FA98 integer witness contains a full matrix V and diagonal D.
Its original exact certificate gives

```
||V^*V-I||<=d<1,   ||Q-V||<=d,   ||H_B Q-Q D||<=r,
```

where Q is V's polar unitary; r is the certificate's complete polar radius.
Let e0 approximate E0 with error e and let the normalized first column
approximate the ground with vector error epsilon. Duhamel gives
`||exp(-itH_B)-Q exp(-itD)Q^*||<=|t|r`. The Fourier coefficients of F_j
are `(M-|n|)/M^2`, with total absolute sum one and
`sum |n| (M-|n|)/M^2=(M^2-1)/(3M)`. Thus

```
||F_j(H_B-E0)-V F_j(D-e0)V^*||
 <= tau*(M^2-1)/(3M)*(r+e) + d*(2+d).
```

The normalized ground error adds at most 2*epsilon to any joint probability.
For each k, the computed expression is the positive sum

```
sum_i F_j(D_i-e0) |v_i^* Q_k A v_0/||v_0||||^2.
```

Every overlap squared is evaluated as an exact rational from integer
vectors and permutations. Cosines use rational pi enclosures, exact range
reduction, a bounded dyadic rounding of the reduced argument, and a
degree-32 Taylor polynomial with its degree-34 remainder. The argument
rounding error is added by the global unit Lipschitz bound for cosine.
The displayed operator bound supplies intervals for the true probability.
No eigenspace separation, choice of basis in a degenerate level, or
eigensolver output is accepted as an independent proof. This calculates
the complete 8-by-4 table from the original source and original witness.

## Native compilation and an executed dimer control

The spin encoding, twelve six-level port generators, complement exports,
and source-read precision algorithm are those of `native_spin_observer`.
For any Pauli P and a separate clock qubit c,

```
controlled exp(i theta P)
 = exp(i theta P/2) exp(-i theta Z_c P/2).
```

This includes scalar P=I terms; dropping them changes the relative clock
phase. A logical site SWAP is exactly exp(-i*pi*h_xy/J), including the
exchange scalar. Hence graph permutations and their controlled versions
are products of native controlled exchange words. Controlled Z_all is a
controlled Pauli, and the finite clock Fourier phases have explicit Pauli
expansions. The selected field uses controlled encoded-X rotations.

For a general finite graph, let C be half the sum of norms of commutators
between the chosen exchange and field terms. First-order splitting of
all controlled binary clock intervals with maximum step Delta has
operator error at most C*Delta*(M-1)*tau. One explicit cubic bound is
`C<=A J^2+E J B`, where A counts overlapping unordered edge pairs and E
counts edges. All controls are blocks of the same unitary, so this bound
does not acquire a clock-dimension factor. With N emitted pulses, dyadic
port precision p adds at most 29*N*2^-p. Approximate ground preparation,
the certified origin error `(M-1)*tau*e`, Fourier-angle errors and these
simulation errors add in the complete dilation norm. Twice their sum
bounds the unhalved trace error for all recorded outcomes and references.
This is a finite constructive circuit, with the stated precision interface;
a fixed binary64 source is not extrapolated to zero error.

`experiment.py` executes the complete two-site J=1, B=1/4 instance, with
one combined-reflection bit and two energy-clock bits. Its ground is
proportional to `(r,1,-1,-r)`, r=sqrt(5)-2, in the physical spin frame.
The local encoded-Z probe has exact excitation energy
`Delta=1/2+sqrt(5)/4`. Its reflection probabilities are
`(1-2/sqrt(5))/2` and `(1+2/sqrt(5))/2`. Both reflection sectors have the
same energy here; this is a useful control against inventing a dispersion
from different character labels. Binary clock intervals one and two use
32 and 64 splitting steps. The uniform arbitrary-input splitting bound is
3/128. For this particular prepared state the splitting error is exactly
zero: its transverse probe lies in the span of `(1,0,0,1)` and `(0,1,1,0)`,
and both the exchange and selected-field terms annihilate this span. The
combined reflection preserves it. Every controlled split factor therefore
has its exact source action on every ideal prefix. The initial ground
error and all literal native errors are still added; no exactness is
claimed for an arbitrary input outside this common kernel.

The execution evaluates the real six-level matrices on every one of the
128 work-register rows. It retains every first-failure probability and
the complete successful output vector; no branch is success-renormalized.
The cube table is certified source spectral computation, not an assertion
that its much larger native program was executed. Complete C/E repairs
may be inserted between completed spin gate words: their ideal action on
single occupancy is identity with fair recorded coins. Their actual
off-code action and their own full pulse/error costs must be included;
they are not inserted during unprotected partial pulse sequences.

This joins a known finite reference to an actual finite spectral instrument.
It does not supply the ground-specific cancellation needed for a sharp
scaling shell, an effective box-size rate for the infinite ordered state,
Lorentz symmetry, a physical mass calibration, or Einstein dynamics.

## Finite transverse packets with every rejected outcome retained

The same instrument applies to a prescribed finite packet
`A_f=sum_x f_x S_x^a`. Set `lambda=sum_x |f_x|/2>0`. Prepare a site register
with amplitudes `sqrt(|f_x|/(2lambda))`, apply the controlled unitary
`exp(i arg f_x) sigma_x^a`, and undo that preparation. The blank-register
Kraus operator is exactly `A_f/lambda`; all other register outcomes complete
the instrument and are recorded. The coefficients are a specified finite
list of computable complex numbers with certified approximants, explicitly
given nonzero support, and a certified positive lower bound on lambda.
The argument of each supported coefficient is consequently computable;
zero coefficients require no phase operation. Uncomputable coefficients
are outside this algorithmic statement. Native tree preparation and controlled
Pauli words implement the construction at arbitrary declared dyadic source
precision. If the coefficient-state preparation has isometry error
epsilon_p and the controlled phases have maximum error epsilon_phi,
preparation followed by its inverse and SELECT adds at most
`2epsilon_p+epsilon_phi` in complete-dilation norm, before the separately
budgeted native pulse errors. These errors enter delta below.

After the joint clocks, the unconditional probability of packet success
and outcome `(k,j)` is
`||Q_k B_j(K) A_f Omega||^2/lambda^2`. Summing successful clock outcomes
gives `s=||A_f Omega||^2/lambda^2`. No normalization by s is hidden in the
raw readout. If a complete native approximation has unhalved trace error
at most delta, every recorded event has probability error at most delta/2.
If an independent argument supplies s>=s0>0 and delta<s0, the approximate
success probability is positive. The conditional successful state has
unhalved trace error at most `2delta/s0`: restrict to the success
block and account separately for its changed normalization. A state with
zero or uncertified success weight receives no conditional guarantee.
This exposes the small-weight cost relevant to Goldstone packet work,
while the complete joint probabilities require no such lower bound.
