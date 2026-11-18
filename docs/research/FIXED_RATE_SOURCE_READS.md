# Fixed-rate noise and protected source reads

The [fermionic source construction](FERMIONIC_SOURCE_CLOCKS.md) supplies an
all-sector clock and actual read instruments, but its sufficient unprotected
noise budget decreases as the lattice is refined. Here that restriction is
removed for the named operational outputs: **a fixed finite positive rate of
the same independent physical dephasing admits vanishing clock/read error**.
The price is error-correcting redundancy and a paid increase in local control
speed. Whole-device microscopic fidelity is a different requirement and can
still tend to zero. Both statements are proved below.

This is a construction inside the same specified proper-code source and its
cofinal CCG capability. It does not select that source, a noise rate, masses,
or a laboratory energy standard from A1--A3. Each noisy finite member
approximates the exact CCG reference; it does not satisfy exact read agreement
despite nonzero errors. Classical records remain reliable in the parent's
declared quantum-noise model.

## 1. Standard results used, and what is new here

We use the constructive concatenated-code accuracy theorem of
[Aliferis, Gottesman and Preskill (AGP), sections 2--7](https://arxiv.org/abs/quant-ph/0504218):
fixed finite universal gadgets can suppress local stochastic faults
quadratically at each encoding level. Their extended-rectangle argument
handles overlap between consecutive recovery circuits. Their section 7.3
supplies verified non-Clifford state preparation and feed-forward; a
transversal Clifford circuit alone would not suffice. We use finite bounded
ancilla attempts, rather than their optimized numerical threshold estimate.
Finite universal synthesis is the constructive
[Solovay--Kitaev algorithm](https://arxiv.org/abs/quant-ph/0505030).
These are imported mathematical results, not new OPH premises or discoveries.

The OPH work is the translation from continuous native interface noise to
those hypotheses, a checked bounded recovery circuit, the causal source
schedule, the fixed-rate refinement and the distinction between raw and
decoded accounting. We do **not** transplant AGP's numerical threshold into
our differently scheduled hardware. Constants belonging to a fixed universal
gadget library remain symbolic. The finite executions below establish their
specified circuit identities, not a numerical threshold for a full machine.

## 2. Noise during a pulse is a local circuit fault

Let an active interface of dimension d have generator

```
L_t(rho) = -i[H_t,rho] + lambda (Delta(rho)-rho),             (1)
```

where Delta is its fixed computational-basis dephasing channel. A Poisson
expansion realizes (1) by ideal driven evolution interrupted by independent
Delta jumps at rate lambda. For an interval of length h, the no-jump branch
has weight exp(-lambda h) and is precisely its ideal driven channel. The
sum of all other branches is a CP trace-preserving channel after dividing
by its **known total weight**, so

```
T_h = exp(-lambda h) U_h + (1-exp(-lambda h)) N_h.            (2)
```

This is a channel decomposition, not selection of successful laboratory
records. The experiment retains the full channel T_h. The calculation does
not commute dephasing past a noncommuting drive. The evidence checks the full
six-dimensional driven channel and the CP/TP remainder using an independent
matrix exponential, including a noncommuting pulse.

The source's qubit and pair codes are coordinate subspaces. H, T and CZ
preserve their respective proper codes throughout their pulses; Delta does
too. H-CZ-H gives CNOT. Resets, reads, load/store and pack/separate use the
existing complete instruments. Thus these noisy gates remain channels on
one or two qubits, with no unhandled leakage. The complementary outcomes
continue to exist on general physical inputs. All extra code qubits have
their own counted buffers/processors; parallel operations never share one.

Partition the complete scheduled hardware history into elementary gate,
preparation, measurement, flight and idle locations. A native decomposition
of a one- or two-qubit circuit location has bounded size. Group its active
interfaces and time intervals; a Poisson hit is a fault on at most those
qubits. Idle intervals on every live qubit are included. Disjoint groups
have independent jump indicators. Padding makes the exposure bounds fixed
even when classical outcomes select a different ancilla. Consequently, for
any specified set S of locations,

```
Pr(all locations in S are faulty) <= p^|S|,
p <= C_native lambda a/c,                                  (3)
```

when a/c bounds a coarse flight and a local control slot. C_native is a
finite constant for the fixed native gate interface, independent of q.
Conditional on a fault the operation can be arbitrary; it need not be a
Pauli channel. Quantum interfaces at rest, active processors, visitors and
helpers all enter the grouping. An instantaneous remote gate is not used.

## 3. A complete, bounded recovery circuit

Use the [[7,1,3]] CSS code with binary check columns 1,...,7. Its three
row masks are 85,102,120. The six stabilizer generators are X(h_j), Z(h_j).
The logical basis is the uniform row-span state and its all-bit complement.
There are 64 orthogonal syndrome subspaces of dimension two. For syndrome
s=(s_X,s_Z), choose C_s with one Z at column s_X and one X at column s_Z,
omitting a factor when its three-bit column is zero. If J is the code
isometry, the complete decoder has Kraus maps

```
K_s = J* C_s*,           sum_s K_s* K_s = I_128.             (4)
```

The code projection is implicit in J*: K_s vanishes on all the other
syndrome spaces. Equation (4) includes all 64 outcomes and is defined on
every input. For every single-position Pauli E, each K_s E J is a scalar
multiple of I_2. This proves correction of coherent linear combinations
and spectator entanglement, not just correction of a sampled error label.
The independent checker constructs stabilizer projectors and verifies these
full operator identities against the emitted sparse Kraus maps.

The physical recovery does **not** assume that this ideal decoder is a free
operation. For each weight-four stabilizer prepare two separate four-qubit
cats. Verify each cat's three adjacent ZZ checks with a reset verifier.
Retain both candidates' check outcomes; use the first accepted candidate.
Couple its four qubits transversally to the four data positions and read
the cat in X. Both candidates are prepared even when the first succeeds.
If neither succeeds, retain a local failure flag and use the second cat as
the fixed fallback. There is no unbounded retry or missing output channel.

Read all six stabilizers in four rounds. Select the first pair of consecutive
equal six-bit syndromes; if no pair exists, flag local failure and use the
zero correction. Later rounds and their records remain in the schedule even
if agreement occurred earlier. Finally
apply the selected Pauli correction with padded slots for all seven X and
all seven Z positions. Every unused qubit has an explicit idle location in
each slot. There are 16 qubits, 1,262 primitive slots, 18,402 idle locations
and **19,664 fault locations** in this particular abstract Clifford schedule.
The native source decomposition charges its additional pulses and events.

The ideal reference instrument is checked independently of Pauli frames.
Contract the actual cat preparation, verification, controlled-data gates
and X reads on the full four-qubit active data space. For stabilizer S,
first-cat read string u and unused-cat environment index v, its Kraus maps
are `(I + (-1)^popcount(u) S)/8` for v=0 or 15, and zero otherwise.
All 256 maps sum to a trace-preserving instrument; grouping by read parity
gives the two projectors `(I +/- S)/2`. Each of the six rejected ideal
verification maps is zero on every data input. This is an operator check,
not a sampled logical state or conditioning on a nonzero success probability.
The other three data qubits are arbitrary spectators. Omitting cat
preparation, a data coupling or a read, or using the wrong coupling basis,
fails this check even without relying on an error-propagation census.

A single preparation fault producing a harmful multi-qubit cat X error is
detected before the cat touches data. A global cat X is harmless modulo the
measured stabilizer; a single remaining cat X affects one data position.
A verification fault can either reject the candidate or affect at most one
relevant position. If the first candidate is rejected by the only fault,
the second is clean. Phase faults can change a syndrome result; repeated
rounds prevent an erroneous recovery from becoming a logical error. A data
fault can cause one transition between true syndromes and at most one mixed
round. Four rounds suffice for consecutive agreement, or an earlier agreed
syndrome leaves just that one later data error.

The exact fault census covers all **65,328 nonidentity one-location Pauli
faults**, including all fifteen possibilities on each two-qubit location.
Every result differs from the intended output by at most one Pauli error
modulo a stabilizer; no such case aborts. Forward symbolic propagation and
independent backward propagation give identical responses. A separate
adaptive executor checks actual selection of the second cat. Removing cat
verification produces 192 failing cases; one syndrome round produces 1,565.
Thus ideal recovery has not been used to bless a non-fault-tolerant circuit.

The syndrome-selection rule is invariant under adding the same initial
syndrome to all four rounds. Its final syndrome defect therefore does not
depend on that initial syndrome: if the injected propagated Pauli has
syndrome e and the selected fault-induced read is s, an initial syndrome t
leaves output syndrome `t XOR e XOR (t XOR s) = e XOR s`.
The clean-input census puts this in a one-error syndrome subspace. An
arbitrary state is a superposition over all 64 such initial subspaces;
complete stabilizer measurement and linearity give the same valid-output
property. Together with (4), this establishes the
usual recovery properties on arbitrary inputs, on one-error inputs without
new faults, and on clean inputs with one fault. Pauli expansion extends the
argument to arbitrary one-location fault operations. Transversal H, physical
S-dagger (logical S in this basis), and CNOT preserve one-error propagation
between blocks; full encoded-basis identities are checked as well. Universal
protection also uses the verified non-Clifford gadgets cited in section 1; the Clifford
census is explicitly not an execution of those entire gadgets.

For those fixed universal gadgets, cap each verification supply at two
fresh attempts and retain exhaustion as a **local** failure. Continue with
a fixed, possibly erroneous local output; a failed supply never acts outside
its gadget's wires. Ideal verification accepts;
exhaustion requires at least two faulty attempts. Legitimate ideal positive
and negative measurement outcomes receive their prescribed feed-forward;
they are not rejected as preparation failures. For the repeated nondestructive
logical measurements, use four padded rounds with the prescribed recovery
between rounds and the first adjacent agreement. With at most one fault,
there is at most one disturbed round and one change of the underlying
eigenvalue, so four rounds suffice. No agreement is a local failure with a
fixed fallback; it is never an unbounded retry.
In the non-Clifford preparation, the measured observable is the Hermitian
`Q=T X T*`, with eigenvalues +/-1. With our convention `S=diag(1,i)`, this
is `exp(-i pi/4) S X`, not the unphased product S X. Its controlled phase is
part of the finite native circuit. A negative eigenvalue is corrected by Z;
both teleportation read outcomes are retained with the required S correction.
The nominal two-outcome injection channel and these phase conventions are
checked directly. This does not replace the cited noisy encoded-gadget proof.
This retains one-fault correctness with a larger finite gadget constant.
Lower-level failure flags are diagnostic records, not commands to abort the
whole experiment: an upper code can correct a faulty lower gadget. Only a
failure declared by the outermost protected gadget becomes an application
abort. Aborting on every lower-level flag would invalidate the claimed
power of p. Bad extended rectangles include their own exhaustion and
undetected logical faults, and higher levels treat them as local faults.
The established concatenation argument then supplies a finite A such that

```
p_k <= (A p)^(2^k)/A,       A p < 1.                        (5)
```

A depends on the complete fixed universal library, including idle slots
and bounded ancilla supplies. **It is not inferred from the recovery count
19,664 alone.** Our proof needs its finiteness, not an invented numerical
threshold. Overlapping recovery rectangles are handled by the cited
extended-rectangle theorem, not by pretending they fail independently.

## 4. Protection without a diverging propagation slowdown

Fix physical cube length L=qa, masses, read kernels and a native observation
horizon. The parent's fully padded preparation, walk, read and report
program has O(q^4) primitive quantum locations, including waits, and O(q^3)
live qubits. Approximate each requested one-qubit unitary to operator-norm
error at most q^-16 using the finite native universal set. CNOT and the
standard preparations/reads need no synthesis approximation. The sum of
the gate errors is O(q^-12), also for instruments with feedback and arbitrary
spectators, by channel contractivity and telescoping.

Solovay--Kitaev gives O(log^4(2q)) gates per requested rotation. Allowing
classical identifiers, padding and local bookkeeping, a safe envelope for
the protected program's level-zero location count is

```
M(q) <= C_M q^4 log^6(2q).                                  (6)
```

Choose **four concatenation levels**. The code blocks and the complete
level-four universal gadgets have constant finite size, independent of q.
Software ancillas are prepared and verified with the same noisy gates;
they are not granted error-free. Known source preparation starts with noisy
standard resets and is included in the protected program.

There is a scheduling detail that matters. Serially transporting every
encoded constituent over a macroscopic link would slow propagation as
redundancy increases. Instead move the whole finite block in parallel, with
one counted source and destination interface per constituent. Each moves
once, resets its sender, overwrites its receiver and records its result.
Recovery is done locally at the endpoints. A block flight is a transversal
identity gadget. The recursively nested leading recoveries can be executed
at departure and trailing recoveries at arrival; there are no extra trips
over the physical link. Every constituent still travels at speed at most c.

Within each owner allocate a finite cluster of processors and buffers for
its protected gadgets. Let B_q be an upper bound on the actual native and
classical elementary operations of its longest local compiled control slot.
The fixed four-level library and (6) give B_q=O(log^6(2q)). Determine B_q from
the finite compiled instruction table, before quantum preparation. Choose

```
eta_q = a/(1024 c B_q),
Omega_q = pi/eta_q,
cluster radius r_q <= c eta_q/128.                          (7)
```

All these quantities are finite and positive at each q. Gates within a
cluster are serialized if they share a processor. Their total native pulse,
event, control-message and internal-flight duration is bounded by a/c;
the factor 1024 leaves room for the bounded native decomposition and both
endpoint recoveries. An external block flight has the old link length plus
at most 2r_q, and retains its finite propagation time. Pad each complete
source slot to the same conservative constant multiple of its original
duration. Idle blocks run their local maintenance in parallel. This produces
a fixed q-independent wall-time factor K_protected, finite quantum exposure
for each q at fixed L and horizon, and O(q^3) quantum interfaces with a larger
constant. A conservative event/identifier ledger is O(q^4 log^7(2q)).

This uses the parent's cofinal availability of small positive service times
and finite per-owner inventories. It **does** require faster local controls
than the particular choice Omega=pi c/a of the parent: a polylogarithmic
factor is paid. Neither unbounded fixed-device strength nor a decreasing
physical noise rate is assumed. The fixed classical program is distributed
in finite charged prehistory, as in the parent; dynamic syndrome processing
and all quantum waits are charged during the experiment.
Identifiers and the fixed control table are preloaded in that prehistory.
Local recovery decisions have constant-size inputs; retained records are
written to separate counted central slots in parallel. No q-dependent global
query is hidden in a local gadget or in the constant A.

## 5. Fixed-rate operational refinement theorem

For every fixed finite lambda>0 and c,L>0, (3) has p=O(q^-1). The library
constant A is fixed, so A p<1 holds on a nonempty cofinal tail. Equations
(5)--(6), with k=4, give

```
Pr(a bad protected rectangle anywhere)
    <= M(q) p_4 = O(q^-12 log^6(2q)).                       (8)
```

Include the synthesis error. For each specified source program, the joint
state of its decoded logical output and its **named application records**
differs from the ideal program by trace distance

```
epsilon_q = O(q^-12 log^6(2q)).                             (9)
```

The same conclusion is a channel bound, including arbitrary spectators, for
inputs already in the protected code. The actual clock experiment, including
its gauge vacuum, packet preparation, flight history and detector, starts
from protected standard preparations, so it needs no ideal unknown-state
injection. Every abort is a distinguished output and contributes to (8).
No division by an acceptance probability occurs in (8) or (9). The abort
flag here is the outermost application abort just defined; lower-level
failure flags remain in the auxiliary record and may be corrected above.

The parent's resolved number instruments, complete weighted detector and
non-Gaussian finite interaction are source programs and therefore inherit
(9). In particular each named event probability changes by at most
epsilon_q. Differences of two event probabilities change by at most
2 epsilon_q. The ideal positive clock swing consequently survives at fixed
physical noise on a sufficiently fine member of the family. Macrostep time
is a common constant rescaling; the same clock/velocity dilation ratio
survives. This does not assign a laboratory value to that rescaling.

For the parent's finite positive accounting operator E_q, including its
invalid-fermion-code outcome, 0<=E_q<=E_max(q) I with E_max(q)=O(q^4).
Let D_q be the complete decoding channel, extended to send a flagged abort
to an orthogonal failure symbol with accounting value E_max(q). Define the
fixed physical accounting effect of the protected representation by

```
E_q^protected = D_q* (E_q direct-sum E_max(q)).               (10)
```

It is positive, bounded by E_max(q), and defined on every noisy input.
To extend it to the full physical interface space, let P_carrier be the
tensor product of native qubit-carrier projectors. Apply the complete
syndrome decoder to `P_carrier rho P_carrier`, and send the complementary
weight `Tr((I-P_carrier) rho)` to the failure symbol. These two terms sum
to unit trace on every physical input, including carrier leakage.
P_carrier is not the much smaller error-correcting code projector: all
correctable syndrome spaces inside the carrier remain decoded normally.
This supplies an accounting effect on general physical inputs even though
the specified dephasing and native pulses preserve the carrier.
Unitality of the adjoint and (9) give

```
|actual accounting - ideal accounting|
    <= E_max(q) epsilon_q = O(q^-8 log^6(2q)) -> 0.           (11)
```

This is the same logical accounting observable in an error-correcting
representation, with every syndrome included. It is not a physical
Hamiltonian and does not identify drive work with Floquet accounting.
Three levels would give only an O(log^6 q) bound after this full spectral
range factor. That fails to establish accounting convergence by this
estimate; it does not prove that three levels physically fail or that four
are optimal.

## 6. Why microscopic fidelity is the wrong exit condition here

The old encoded-vacuum obstruction remains valid for the unprotected
encoding. There is also a direct obstruction for protected hardware. In
each lowest-level seven-qubit block choose one physical Z whose syndrome
column is nonzero. Conditional on all other independent phase errors,
there is at most one value of this chosen error bit that makes that block's
syndrome zero. If its error probability is p_Z<=1/2, the probability that
all N such blocks are syndrome-free is at most (1-p_Z)^N. Hence fidelity
with any state in the exact protected code is no larger than that quantity.
For one a/c wait at fixed lambda, p_Z=(1-exp(-lambda a/c))/2=Theta(q^-1),
while N=Theta(q^3): the raw fidelity tends to zero.

This does not contradict (9). A correctable syndrome changes the physical
state while leaving its decoded logical state unchanged. The complete
seven-qubit calculation gives, for independent phase flips of probability p,

```
raw entangled fidelity = (1-p)^7 + 7 p^4(1-p)^3,
decoded phase error = 21 p^2(1-p)^5 + 7 p^3(1-p)^4
                    + 28 p^4(1-p)^3 + 7 p^6(1-p) + p^7.     (12)
```

All 128 error patterns and all decoder branches are included. The memory
polynomial illustrates the distinction; it is not substituted for the
noisy-gadget proof of (8). The evidence also reconstructs D*E and tests its
expectation on faulty encoded inputs that a raw code-violation penalty would
assign the maximal value.

Auxiliary syndrome and verification records remain in the source history.
They generally distinguish noisy from noiseless hardware. The semantic
comparison in (9) retains the requested source-read/clock records and abort flag,
while ignoring auxiliary error diagnoses; it is not equality of every
microscopic transcript. Ignoring a field of an unconditional distribution
is different from conditioning on its successful value.

Finally, an arbitrary **unencoded** unknown input cannot be made arbitrarily
accurate after an initial noisy encoding step: a logical phase error before
redundancy is added is indistinguishable from a different input. The evidence
contrasts phase noise before encoding with a correctable physical phase
error after encoding. The theorem covers protected inputs and the complete
known-preparation experiments, rather than smuggling a perfect raw encoder
into the source.
More generally, phase noise of probability p<=1/2 reduces the trace distance
between the raw states |+> and |-> from one to 1-2p. Any subsequent encoding
and recovery are contractive. The triangle inequality therefore forces a
worst-case output error of at least p for this pair, whatever encoder is used.

## 7. What this closes

Decreasing physical dephasing is not necessary for the **constructed
operational refinement** of the named fermionic source programs. Preparation,
recovery, non-Clifford computation, geometry, causal transport and physical
time all have an implementation in the same cofinal source capability.
Corrected outputs and their positive accounting converge at fixed rate;
raw fidelity, raw syndrome penalties and all-record indistinguishability do
not need to converge and in general do not.
The new syndrome and control events also change raw event counts. No equality
of microscopic event rates or entropy coefficients with the unprotected
representation follows from the named-output bound.

The code in [m1_fixed_noise](../../code/m1_fixed_noise/README.md) provides the
complete decoder, continuous M6 checks, bounded recovery, independent fault
census, adaptive controls and hostile tests. Its receipt is compact and binds
the proof, parent evidence, verifier and scoped claims. There is no asserted
numeric full-machine threshold or executed astronomical fault-tolerant
clock. Those would be different deliverables from this all-size existence
and operational-stability theorem.
