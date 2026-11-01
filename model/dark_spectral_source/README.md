# Native preparation and Ramsey readout of the dark spectral sector

The previous [repaired-spectrum result](../repair_spectral_sectors/README.md)
identified an exact 15-multiplicity, spin-one, eta-two dark sector on the
unchanged six-site Weyl/contact chain. Its states preserve nonconstant
energy coherences under every complete existing block repair. The missing
source preparation and detector now have concrete native words.

This package prepares one sparse dark state with **292 native pulses**;
its spectral measure has strictly positive weight on all **13 distinct
dark energies**. A controlled preparation and complete two-phase Ramsey
readout use one additional pointer bit. Every pulse acts through the
existing twelve-port, six-level source and its complete proper-code
export. No target statevector or desired unitary is a processor primitive.

This advances CP-03/#1139 and CP-07/#1143 on a fixed finite source. It
does not establish physical charge coherence, particle identity, a
localized three-dimensional defect or a continuum mass. The source's
coherent qubit encoding, initialized registers, load/export, memory,
readout and phase controls remain explicit specified services. The
vacuum/N=4 Ramsey superposition respects fermion parity but uses a supplied
number-phase reference; an electromagnetic charge reference or a physical
superselection rule is not derived or circumvented.

## Sparse deterministic preparation

Use the Hubbard onsite frame `d_x=(i sigma_x)^x c_x`, x=0,...,5. In this
paragraph sigma_x means the usual Pauli X matrix. Define

```
eta_+ = sum_x (-1)^x d_(x up)^* d_(x down)^*,
eta_z = (N-6)/2,
|D01> = eta_+ |0 up,1 up>/2.
```

The normalized state has only four occupation amplitudes. In little-endian
up/down CAR order they are `(53,+1/2)`, `(197,-1/2)`, `(773,+1/2)`, and
`(3077,-1/2)`. The original-frame target is `frame^dagger |D01>`.
The unchanged exact dark-column operator supplies an independent equality.

The circuit starts with up bits at sites 0,1,2. Three balanced two-bit
rotations split the site-2 bit to sites 2,3,4,5. Four exact copy gates put
the same occupancy in the corresponding down bits; Z phases at sites 3
and 5 provide the alternating eta signs. The original inverse onsite
frame is compiled last. There is no charge measurement, postselection,
generic 4096-amplitude tree, supplied spectral eigenvector or uncomputed
ancilla. SymPy checks the two-bit rotation and copy identities exactly;
a separate sparse symbolic calculation verifies all four signed terms.

Every operation is a Pauli exponential with a rational multiple of pi.
Its scalar phase is included before adding the Ramsey control. In
particular, `exp(i*pi*X/2)=iX` must not be silently substituted for X.
For any Pauli P the controlled word uses the exact commuting identity

```
C exp(i theta P) = exp(i theta P/2) exp(-i theta P Z_pointer/2).
```

This acts as identity for pointer zero and as the full original
exponential for pointer one, including scalar phases. The complete
controlled preparation, including the initial balanced pointer, has 807
native pulses. Its initial input is the 13-bit blank. Two independent
pointer input columns also check the coherent map and reference Gram.

## A complete detector with a positive spectral measure

Write U for the preparation circuit. The Ramsey sequence prepares

```
(|0> |vacuum> + |1> U|vacuum>)/sqrt(2),
```

evolves the field with its **original H**, inserts complete E/C/E repairs,
applies controlled U*, a pointer phase phi, and the inverse balanced
pointer rotation. The actual incoherent six-level export gives two
pointer outcomes and four abort outcomes. The uncomputed field is retained
in each pointer branch. Neither vacuum return nor pointer success is
postselected.

Both vacuum and the entire dark sector are annihilated by all native
repair Kraus operators; Q=I-P is identity on their direct sum. Consequently
their cross terms also survive the ideal complete repair, and the ideal
pointer probabilities are

```
P0(phi,t) = [1+Re(exp(i phi) C_H(t))]/2,  P1=1-P0,
C_H(t) = <D01|exp(-itH)|D01>.
```

The reference is explicit. Previous dark energies used K=H-N/2. On N=4,
H=K+2, while the vacuum H energy is zero. Hence

```
E^H_pq = 1+[cos(p*pi/7)+cos(q*pi/7)]/8,
w_pq = (16/49) sin(p*pi/7)^2 sin(q*pi/7)^2
                 [cos(q*pi/7)-cos(p*pi/7)]^2,
C_H(t) = sum_(1<=p<q<=6) w_pq exp(-it E^H_pq).
```

Every weight is strictly positive: sine is positive and cosine strictly
decreasing on (0,pi). Exterior-power Parseval gives sum w=1. Exact integer
Hamiltonian moments give mean H energy 1 and variance **1/256**, so the
response contains more than a common number-sector phase. The three
pairs p+q=7 share energy 1; the previous exact characteristic polynomial
proves that the remaining twelve are distinct. The receipt keeps all 15
mode-pair contributions and all 13 exact groups, without thresholding.
An independent two-by-two determinant checks the characteristic function:

```
C_H(t) = exp(-it) det([exp(-it A_chain/16)]_{0,1;0,1}).
```

This gives a source-accessible spectral measure. Two phases at one time
do not empirically resolve thirteen frequencies or fit a particle mass.
The exact spectral formula is derived from the same Hamiltonian; the
native experiment tests its complete detector probabilities. Simple pair
annihilation from the earlier eta-singlet ground still cannot reach this
eta-two sector. The new source preparation is a different explicit
operation, not a reversal of that selection rule.

## Complete source and field accounting

The executed history uses t=1/2 and two original-H source segments of
duration 1/4, each with two Strang steps and 3360 native pulses. Between
them it attempts E on blocks (0,1), then C on (1,2) and E on (0,2).
Each Ramsey analyzer has 809 pulses, so one setting has **8336 native
pulses**, 8336 coherent loads and 8336 coherent exports, plus its final
incoherent pointer export.

The declared complete protocol stops on the first source abort or first
native repair record. The continuing Q branch is never normalized. The
receipt retains all 40 possible native repair records, every first source
export failure, and all six final pointer records. Source abort vectors
remain reconstructible from their immutable input, tape and full six-level
matrices. Four-column failure tables include explicitly identified padding
for pair-code events; padding is not counted as an additional outcome.
The ideal protocol takes only Q and the two proper pointer outputs.

The earlier exact-source checker independently verifies the native gate
routing and all fixed pi angles. Literal binary64 port coefficients are
interpreted as exact rationals in the retained Q(sqrt(5)) source response;
outward dyadic pi bounds and anti-Hermitian Duhamel estimates include
proper-code leakage. Complete export isometries do not enlarge these
operator-error bounds. This supplies a mathematical source-specification
bound, not a claim about the accuracy of experimental hardware or the
host matrix-exponential implementation.

The two Strang errors use the unchanged complete-Fock support-count bound.
Source and splitting errors compose through the complete repair/output
isometries, for arbitrary untouched references. The total proved bound
is below **0.007215** per setting; its source-only contribution is below
`1.048e-12`. The executable floating diagnostic is about **2.753e-7** and
is kept separate from that proof. Numerical pointer-zero probabilities
are about `0.938577077511` and `0.739595749433`; the analytic values are
`0.938577079690` and `0.739595750623`. All probability is retained to the
stated numerical admission, without normalization.

Countercontrols remove a required scalar phase, dephase the pointer,
inject an actual source-code leakage rotation, alter a native target or
angle, and omit field evolution. They distinguish a complete coherent
source/readout experiment from a target statevector or successful branch.
The inherited active90/dark15 classification, repair energy disturbance,
3D alignment criterion and transverse-dimer obstruction remain unchanged.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.dark_spectral_source --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/dark_spectral_source
```
