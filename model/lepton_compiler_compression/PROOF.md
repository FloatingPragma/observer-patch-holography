# Exact active-mode and archive reduction

## Fixed source and its native span

Use exactly the finite neutral compiler, phase, physical CAR algebra,
kinetic law and coarse vacuum/nonvacuum records of
[lepton_mass_compiler](../lepton_mass_compiler/PROOF.md). The source has
90 physical CAR modes per site. Its address analyzer acts on 3,600
addresses and three characters. The fixed neutral scalar background is
one; changing it or compiling all Higgs orientations at once is outside
this reduction.

In character k let p be the blank address, v_k=U_k†p, and
lambda_k=<p,v_k>. The literal analyzer gives lambda_k>0 and
e_k=lambda_k² strictly between zero and one. Define

```
a_k=(v_k−lambda_k p)/sqrt(1−e_k),    B_k=(p,a_k).
```

Then B_k†B_k=I, <p,a_k>=0 and the **full**, uncompressed address projector
has the exact factorization

```
W_k=v_k v_k†=B_k w_k B_k†,
w_k=[[e_k, sqrt(e_k(1−e_k))],
     [sqrt(e_k(1−e_k)),1−e_k]].
```

The executable obtains v_k by applying the actual analyzer adjoint to
the blank columns. It checks the factorization on all 10,800 address
columns and directly through the original U†QU action. The exact
rank-one and orthogonality identities establish the theorem; the
binary64 execution checks its actual native attachment. The imported
exact rational intervals certify the strict positive chamber.

## Complete CAR and kinetic intertwiner

Each character has two spin-oriented charged pairing blocks. Each block
contains one physical left and right mode and needs one ancillary left
and right mode. Thus six blocks need twelve ancillary modes. All ninety
physical modes remain, including the 78 physical modes absent from this
fixed neutral pairing. The total is 102 modes per site, compared with
324,000 in the full address realization.

Use B_k on the left and conjugate(B_k) on the conjugate right field. The
exterior lift Gamma(T) of this one-particle isometry maps the complete
compressed Fock space into the original Fock space. T fixes every
physical mode. The pair kernel transforms as T Delta T^T. Consequently
the oriented original block is B_k w_k B_k†, with the reverse transpose
required by the native CAR pairing. Using B_k on both legs fails when
the ancillary frame is complex; the executed control detects this error.

All original address directions perpendicular to the selected spans have
zero pair coupling. The same Weyl hopping acts identically on each
address copy. The selected spans are independent of position and spin,
so hopping preserves them, including its spin mixing. The orthogonal
spectator address modes begin in vacuum and remain in vacuum under this
number-preserving kinetic law. This is an invariant subspace statement
for the complete Hamiltonian, not a vacuum matrix-element approximation:

```
H_original Gamma(T)=Gamma(T) H_compressed.
```

It implies the corresponding exact propagator identity for every time,
all physical occupations, every finite spatial volume and arbitrary
untouched quantum references. Neither a low-particle sector nor a
one-particle replacement of the physical source is used.

The number twelve is minimal within the stated category of
number-preserving CAR mode embeddings that fix all physical modes and
reproduce the first-order coupling. The physical-to-ancillary block has
twelve independent nonzero singular directions: two charged fields,
two spins and three characters, with e_k(1−e_k)>0. A smaller ancillary
one-particle space cannot have that rank. This does not assert a
universal minimal dimension among arbitrary encoded Stinespring spaces.

## Original records and quantum archives

Gamma(T) factors into the identity on physical modes and a
number-preserving ancillary isometry, with the usual fermionic ordering
signs. It maps ancillary vacuum to ancillary vacuum and the orthogonal
sector into the original nonvacuum sector. Therefore the production
vacuum/nonvacuum projectors intertwine. Export the compressed ancillary
state after each round and retain its flag. The tensor product of the
ancillary isometries reconstructs the full original quantum archive,
including every fine-address amplitude and all within-flag coherences.
Unused directions are fixed vacuum under the inverse coordinate map.
They are not traced out as unknown environmental states.

The independent control uses a six-mode raw Fock space with a complex
ancillary rotation and its four-mode active realization. It checks all
sixteen embedded occupation columns, the Hamiltonians, their actual
exponentials, every one-round Kraus operator and all three-round
history operators. The raw archive has 4,096 fine coordinates and the
compressed archive 64. The eight committed classical histories agree
under the archive isometry. A maximally entangled physical/reference
preparation checks the complete law; equality of every input operator
establishes the arbitrary-reference statement. The fine coordinates are
archive coordinates, not additional production measurements. Their
nonzero off-diagonal density entries are retained.

## Exact leakage norm

Let P embed the physical system with fresh ancillary vacuum and let Q=I−P.
For one character and spin-oriented pair, put r²=e(1−e). The sixteen-state
CAR algebra gives

```
P V² P−(P V P)²
  =(1−e²)I−e(1−e)(n_left+n_right).
```

Its four eigenvalues are 1−e², 1−e, 1−e and (1−e)². In particular, the
fully occupied physical pair still admits ancillary pair creation. The
largest eigenvalue is 1−e², attained on the physical pair vacuum.

Different pair blocks and sites emit into mutually orthogonal nonvacuum
ancillary sectors. Their leakage forms therefore add, rather than
requiring the square of a sum of full interaction norms. For V sites,

```
B²:=||Q H P||²=||Q V P||²=V(6−2 sum_k e_k²)<6V.
```

The equality is attained when all six physical pair blocks per site are
empty; arbitrary states, including a filled negative band and entangled
references, are bounded by the same value. The kinetic Hamiltonian
preserves P, so it contributes no off-diagonal block.

## Complete-output error at finite round count

Write A=PHP, D=QHQ, F(t)=P exp(−itH)P and G(t)=Q exp(−itH)P.
The coupled block equations and unitarity give

```
G(t)=−i integral_0^t exp(−i(t−s)D) QHP F(s) ds,
||G(t)||≤B|t|,
||F(t)−exp(−itA)||≤B²t²/2.
```

This finite-dimensional Duhamel estimate is independent of ||A|| and
||D||. In N fresh rounds of duration h=t/N, the all-vacuum branch is
F(h)^N. Contraction telescoping bounds its operator distance from
exp(−itA) by delta=B²t²/(2N). For every physical/reference density matrix,
the corresponding unnormalized all-vacuum density operators differ in
trace norm by at most 2delta. The total probability of all other
classical histories is at most 2delta, either by trace loss or by the
one-round bound ||G(h)||². These histories are positive and occupy
orthogonal committed record blocks. Keeping their complete quantum
archives therefore yields

```
||actual complete output−ideal output with vacuum archive||_1
    ≤min(2,4delta)=min(2,2B²t²/N)≤min(2,12Vt²/N).
```

The last strict comparison is interpreted before clipping; both clipped
bounds can equal two. The nonvacuum probability is at most
min(1,B²t²/N). A per-round complete half-diamond implementation error
eta adds 2N eta to the unhalved output bound. This theorem uses the
committed classical flags. A wholly coherent unmeasured flag dilation
does not inherit the same trace-norm rate.

The existing source mass, four-outcome Ramsey observable and ratio
denominator guard are unchanged. Replacing only the compiler resource
bound permits fewer rounds for the same probability-law tolerance.
The receipt computes conservative integer round counts for tolerance
10^−12 and keeps every exported ancillary mode and flag in the resource
account. These large campaigns are not executed. Finite-shot noise,
hardware error, empirical-scale uncertainty and physical-model error
remain separate; no observed particle identification is promoted.
