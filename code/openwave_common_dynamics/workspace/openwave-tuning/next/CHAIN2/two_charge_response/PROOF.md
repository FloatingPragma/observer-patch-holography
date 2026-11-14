# Two charges on the shared native field

The construction is the fixed full-Fock program
\[
 V_m=U_E(\alpha)U_B(\beta)\Gamma(S)_F
       \Gamma(R\exp(-ihmR)),\qquad h=1/8.
\]
There are 2,048 matter modes, 1,152 physical matter modes, 288 shared
four-bit electric registers, and three source work bits. The physical
graph has 96 ports. The field modulus is 16. The member uses
\(\alpha=.17\), \(\beta_0=.41\), and zero coupling on the other 207
faces. The literal binary64 coefficients have the ancestry specified in
the frozen native field and mass programs. The mass values 0 and 1 are
architecture inputs. No observed mass or fine-structure constant is input.

`INPUT_SNAPSHOT.json` binds the exact 299-file input snapshot. Its many-body
manifest digest is
`230a908ce84bab53de63dc6db6b99adc57b0bafa336bb30c611af1f6af7b27bf`.
The field manifest digest is
`aa4a18d89c1a8893ee86eb5e9151bc75e8a0b3a527640cee66aa25382ca55e76`;
the particle manifest digest is
`678983fa40a2a97b8eb399d6b5ad4ab6ace77a5f441baf448cb0c00bc8a6f9a9`.
The separately hashed import annex supplies static Python dependencies
absent from that manifest closure. It changes none of the 299 input files.

## Exterior square and the complete source

Let \(u\) be uniform on the twelve physical event modes at one port,
with zero entries on its four padding modes. Then
\(R=2|u\rangle\langle u|-I_{16}\). On every number sector,
\[
 \Gamma(R)=(-1)^N(I-2n_u),\qquad
 \Gamma(e^{-i\mu R})=e^{i\mu N}
       [I+(e^{-2i\mu}-1)n_u],\quad\mu=hm.
\]
Choose a one-particle basis with first vector \(u\). In that basis
\(R\) has one eigenvalue +1 and fifteen eigenvalues -1. A fermionic
occupation vector therefore has phase \((-1)^{N-n_u}\), proving the
first identity. The sum of its occupied one-particle eigenvalues is
\(2n_u-N\), proving the second. Both identities include vacuum,
padding, spectators, and double occupation in distinct modes.

An antisymmetric coefficient matrix \(A\) transforms as \(KAK^T\),
where \(K=R\exp(-i\mu R)\). For two different port blocks the stored
matrix is the ordered coefficient block; for coincident ports its strict
upper triangle stores the two distinct occupied modes. The code expands
the latter into \(A-A^T\), applies \(K\), and stores the strict upper
triangle. It never deletes same-port, distinct-mode outcomes.

For \(N=2\), the reflection identity is also
\(A\mapsto A-2(PA+AP)\), where \(P=|u\rangle\langle u|\).
The missing quadratic term is zero because \(PAP=0\) for antisymmetric
\(A\). All 120 local two-mode basis vectors check this identity.
`native_contract.py` compares the literal emitted reflection and mass
words with their frozen generators and executes their actual native
M6-derived pulses on a reference-column array containing vacuum, a
padded one-particle state, a complex superposition of all 120 two-mode
pairs, and the fully occupied sixteen-mode state. This replay retains
the global phases.

Transport applies the frozen 2,048-mode permutation to both occupied
modes. Sorting their images contributes the canonical anticommutation
sign. The shared current is the sum of their two oriented edge currents,
including cancellation when opposite occupants exchange an edge. It
updates one shared field. `mode_field.mode_transport` supplies the
independent routing comparison, including padded modes. The transport
preserves \(n+B F\) modulo 16. The local coin preserves each port load;
the face vector satisfies \(BA_0=0\). Consequently the whole update
preserves the declared Gauss sector.

The electric phase is \(\exp[-i\alpha\sum_e\bar F_e^2/2]\), with
\(\bar F_e\in[-8,7]\). The plaquette step acts on all sixteen members
of each complete cyclic face orbit. Every orbit coefficient is retained.
No classical pair potential, field assigned independently to each charge,
or amplitude cutoff appears in this execution.

## Preparation and an identifiable probe

The first orbital is uniform in the twelve event modes at port 0. The
second orbital is supported at another port \(y\), so the two are
orthogonal. The initial field is \(F=P_0+P_y\pmod {16}\), with all
cycle additions zero, and the background is fixed at \(G=2\delta_0\).
This is a definite finite field-string preparation. It is distinct from
a uniform affine-field ensemble. Each nonzero electric basis bit has an
explicit native flip in the preparation program.

For the identifiable-probe experiment the second orbital is
\(d_y=(|y,0\rangle-|y,2\rangle)/\sqrt2\). Both event slots are held
slots and its local sum is zero. Thus the source obeys, on the entire
Fock-field Hilbert space and with arbitrary spectator references,
\[
 V_m a_{d_y}^{\dagger}
 =[-e^{ihm}]a_{d_y}^{\dagger}V_m,
 \qquad [V_m,n_{d_y}]=0.
\]
The coin gives the displayed phase. The transport and its field-current
controls commute with this held-mode creation operator. Field-only
electric and plaquette gates commute with it as well. This proves the
identity without assuming any field state or field-particle product.

In the preparation with exactly one occupied \(d_y\), the other charge
lies in its orthogonal complement. The probe's measured port law is the
total expected port occupation minus \(\delta_y\). This definition
uses an invariant occupied orbital; it does not attach persistent labels
to identical fermions. `reduced_dark` executes the resulting exact
invariant representation while retaining the shared field. Direct
generic exterior-square runs at two steps check its isometry and the
invariant occupation for both mass settings and both positions.

The identity proves compact charge-occupation stability for arbitrary
step counts. It does not make an arbitrary field state a Floquet
eigenvector. An initial dark-particle product evolves into the same
matter orbital times its evolving field factor; a full product
eigenvector also requires an eigenvector of the field Floquet operator.
Relative to the vacuum, the retained dark creation phase is
\(-e^{ihm}=e^{i(\pi+hm)}\). With convention \(V=e^{-ihH}\), its
quasienergy difference is \(-(\pi+hm)/h\) modulo \(2\pi/h\).
The phase branch at \(\pi\) is not identified with a low-frequency
physical rest energy. A physical clock, unit conversion, species
dictionary, spin assignment, and empirical mass anchor are separate
obligations.

## The field-string counterexample and two mobile charges

The same commutation identity proves a useful limitation. If the shared
initial field is held fixed, adding or removing the occupied dark
orbital leaves the probe evolution identical, apart from the retained
global creation phase. Moving \(y\) while holding \(G\) fixed changes
the prepared field from \(P_0+P_y\) to \(P_0+P_{y'}\). The nonzero
probe response in this comparison is therefore a susceptibility to
the charged field preparation. It is insufficient evidence for a
dynamically sourced two-mobile-charge force.

The electric-off control has no such position response because its
field translation commutes with the remaining source. The fully
coupling-off control has the same property. The separate-field control
executes two 288-edge registers and records their tensor-product law;
it deliberately changes the source. Coupling interventions at fixed
\(y\) use the identical complete initial state and Gauss background.

The two-mobile experiment prepares the first uniform orbital and
\(R|y,1\rangle\) for the second. Its first coin makes the latter's
first hop deterministic. Subsequent steps use the unrestricted common
source on both charges. This finite preparation reduces support growth
without cutting any outcome. Both coupling settings use the same
initial field, occupations, reference, and source work. Their complete
unordered-pair, electric-field, and current laws diagnose a source
coupling response. Overlapping identical particles have no separately
measured persistent first/second labels in that experiment.

## Records and errors

The full native detector reads all 2,048 matter occupation bits and all
1,152 electric bits. Pair positions and events are deterministic coarse
records. The last step's oriented current is recovered by applying the
inverse transport permutation to the terminal occupied modes. Field
gates commute with this current record. Hence the terminal occupation
and field law also determines the complete joint last-current law,
without an intermediate measurement or discarded trajectory.

The seven-bit one-particle address-parity detector is inapplicable to
two charges. The detector used here retains zero, multiple, padded,
timeout, source/export, and pointer-complement outcomes. The full-Fock
native source operator bounds apply to this number sector, and the
native preparation and full readout costs are added explicitly.
Reset failures use worst-conditional bounds and a union bound, with
no product assumption for correlated initial reset states. The bound
concerns the terminal observable marginal and the original quantum
reference. Fine reset tapes are retained by the instrument; no
tape-reference independence is asserted.

The sparse decimal law uses binary64 arithmetic. Norm, Gauss,
continuity, source-pulse, and invariant-sector residuals are reported.
The native operator certificate compares the exact ideal circuit with
the complete native instrument; it is not an interval certificate for
every printed decimal. No success probability is divided out.

The source carrier contains the full 2,048-mode Fock space and the
full \(16^{288}\) field space. The executed states have explicitly
declared invariant supports. The three-step and block limits are
resource refusals, never a prescription to prune or renormalize.
At modulus 16 the reported electric records are modular source records.
No one-particle integer-field cutoff bound is transferred to two
charges. There is no Coulomb \(1/r\) result, fine-structure prediction,
Maxwell refinement, physical mass, observed species, or stable-matter
claim.
