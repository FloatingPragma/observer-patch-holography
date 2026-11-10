# Actual native charges on a transported flux composite

The mobile link object and the native occupation register carry different
group actions. A protected face position is not an internal color or weak
component. The construction below embeds the fifteen-state representation
into actual native occupations and supplies a finite charge-preserving
transport for a specified bound composite. It uses the same link domain,
native update and current matrices as the owning source packages.

## 1. Charges of the native-empty flux object

The [native occupation law](../species/PROOF.md) has global internal
commutant `U(3) x U(2)`. For a five-by-five Hermitian current T its native
operator is

\[
 Q(T)=\sum_{p,r,a,j,k}c^\dagger_{p,j,r,a}T_{jk}c_{p,k,r,a}.
\]

It acts on native wave occupations. It annihilates the native vacuum,
independently of the literal link word or routing memory. Hence the
native-empty flux route is a gauge singlet under these actual currents.
Its five face labels do not carry their defining five-dimensional action.

The retained fifteen-state representation R has hypercharges
`-2/3,1/6,1,-1/2,1/3`, with multiplicities `3,6,1,2,3`. Its hypercharge
matrix Y_R is invertible. A current-preserving linear map B from any
native-empty flux code into R would obey `Y_R B=B*0=0`; therefore B=0.
This excludes that particular charged identification. It does not exclude
gauge-singlet fields, a charged native dressing, or a changed current
dictionary with its own proof.

There is also a direct obstruction to assigning an added color/weak action
to the five face positions in the native-coupled route. Restrict to one
native wave in `(port5,internal3,repair0,address0)` and use the exact encoded
generator `K=I-Re(U_F)`. For every input face, execute U_F and its inverse
on both Gauss partners. Ninety resulting native basis outputs give fifteen
distinct five-by-five face coefficient matrices B_n. A pure face current
T must commute with every B_n. The executable solves these exact equations:

\[
 \operatorname{rank}\bigl(T\mapsto([B_n,T])_n\bigr)=24.
\]

There are twenty-five matrix entries, so only scalar pure face currents
commute. Scalars commute with the entire K, making the conclusion exact.
As simpler controls, all ten ways to mark three of the five faces with
charge two and the other two with charge minus three have a nonzero
commutator on the same source input. These are actual source operator
tests. They do not exclude a current involving both native and face data.

## 2. An exact native occupation module

For each internal component j, define the normalized native orbital

\[
 d_j^\dagger=\frac{c^\dagger_{2,j,0,0}-c^\dagger_{2,j,0,1}}{\sqrt2}.
\]

Addresses zero and one are the two orientations of seam `(0,1)`. Their
amplitudes sum to zero, so the address reflection gives minus the orbital.
The seam reflector acts only on ports zero and one, leaving port two
untouched. Port two is absent from every controlled interior-transport
endpoint. Consequently, for all twelve allowed link words g,

\[
 W_g d_j^\dagger|0\rangle=-d_j^\dagger|0\rangle.
\]

The full occupation update is the exterior lift of W_g followed by the
actual endpoint-number controls. On products of these orbitals the controls
are identity. Thus U_F acts as `(-1)^n` on degree n. In particular the
nonempty even subspace

\[
 R_{\rm dark}=\Lambda^2\operatorname{span}(d_j)
             \oplus\Lambda^4\operatorname{span}(d_j)
\]

has dimension fifteen and satisfies `U_F=I`, hence `K=I-Re(U_F)=0`.
The full even module also contains its degree-zero vacuum singlet and has
dimension sixteen. Selecting the fifteen nonempty components is a stated
preparation/domain choice; the idle law does not exclude the singlet or
identify it with a physical neutrino.

The certificate executes all sixty one-body source columns and all 180
degree-two/four columns. For the latter it applies the actual one-body
columns to the exterior polynomial and then executes the occupation
controls. This equals the defined full Fock update without expanding and
discarding cancelling bright basis terms. No projected-matrix agreement
stands in for a full operator identity.

The dark orbitals obey the actual primitive current rule

\[
 [Q(T),d_j^\dagger]=\sum_kT_{kj}d_k^\dagger.
\]

Exterior differentiation therefore gives the owned fifteen-state matrices
`rho_R(T)`. On a positive Gauss link orbit let B map each selected subset
to its literal dark-occupation wedge. Direct CAR evaluation verifies

\[
 Q(T_A)B=B\rho_R(T_A)
\]

for all twelve currents and all fifteen columns. Native occupations are
at port two, so the interior Gauss action has no doublet sign there; both
link partners remain in every encoded state. This is an actual native
current attachment, with no empirical masses or trial action coefficients.

The representation weights and anomaly cancellations are those checked in
[the spin/current attachment](../physical_attachment/PROOF.md) and the
[finite gauge attachment](../../standard_model_source_map/GAUGE_ATTACHMENT.md).
Their exact table is one all-left-handed Standard Model family at the
representation level. This statement specifies internal weights, not the
spacetime chirality of a propagating field. The twelve tested currents span
`su(3) + su(2) + u(1)_Y`, the trace-free subgroup of the native `u(3)+u(2)`.
The additional central degree current is also an actual conserved native
operator; fifteen further column checks give charges two on the ten
exterior-two states and four on the five exterior-four states. If promoted
to a gauged all-left Weyl symmetry, its linear and cubic anomaly sums would
be forty and four hundred. Thus the anomaly-free Standard Model subgroup
must not be confused with gauging every native central symmetry. This
representation-level anomaly test does not invalidate the finite global
source current or make its degree charge a physical fermion number.

## 3. A complete bound-composite transport

Use an odd graded flux object and the actual even native cloud. Their
product has odd total grading. At a fixed face/routing component, one
source slot has the code

\[
 \mathcal B=\mathbb C|\mathrm{flat},\mathrm{native\ vacuum}\rangle
             \oplus R_{\rm bound},\qquad\dim\mathcal B=16.
\]

A bound basis vector consists of the literal paired-flux Gauss orbit,
occupied routing/flux flag and one dark native wedge. The same link word
is retained; no duplicate link is introduced by the native cloud. Native
Q acts as `0 direct_sum rho_R(T)` on this code. The odd flux grading is a
declared architecture, not a conclusion from an even native wedge.

For two independent source slots retain all 256 states. Define

\[
 S|a,b\rangle=(-1)^{p(a)p(b)}|b,a\rangle,
 \qquad p(0)=0,\quad p(a)=1\ (a\ne0).
\]

S is a Hermitian involution. It transports the native cloud together with
its flux object. It commutes with every total native current `Q_L+Q_R`;
the executable checks all 256 columns for all twelve currents. Therefore
`E=(1-i)I/2+(1+i)S/2` is a complete unitary with the same conservation law.
An initial left bound cloud gives position probabilities one half left
and one half right. All 256 basis detector outcomes are retained. Two
occupied slots have exchange sign minus one; for two identical internal
labels E has phase minus i. Their even native clouds add no exchange sign.

The source lift through the code isometry J is
`J E J*+(I-J J*)`. It is unitary on the entire physical support and
commutes with the admitted interior Gauss constraints and total native
currents. The range of J is invariant under those self-adjoint currents,
so the orthogonal complement is also invariant. The identity extension
does not omit a source sector. The exact charge continuity operator is
`I_T=Q_L-E*Q_L E`; the right charge changes by the opposite amount.

This is a specified bound transport. The native idle update does not
prepare the bound code, bind the constituents energetically, or move the
native cloud between spatial cells. A positive penalty
`Delta(I-J J*)` can be declared to penalize its complement, but Delta and
that penalty are additional dynamics. K vanishes on the bound dark module,
so the native weak coupling supplies no hidden binding energy there.

The two-slot construction has an occupation constraint. Its local matrix
`|a><0|` does not satisfy the full fifteen-mode CAR on the sixteen-state
slot, as an explicit anticommutator control verifies. Odd exchange grading
is not used to claim an unconstrained composite Fock algebra.

## 4. A finite whole-object CAR extension

There is a legitimate separate extension that retains actual native
currents. Let K denote a finite one-object space containing position,
routing spin, flux face and the native R module, with isometry B into a
complete source object. Put D equal to the dimension of K and provide D
independent whole-object source slots with vacuum labels. For `n<=D`, set

\[
 E_n(v_1\wedge\cdots\wedge v_n)=\frac1{\sqrt{n!}}
 \sum_{\pi\in S_n}\operatorname{sgn}(\pi)
 Bv_{\pi(1)}\otimes\cdots\otimes Bv_{\pi(n)}
 \otimes|0\rangle^{\otimes(D-n)}.
\]

On orthonormal exterior basis vectors this is an isometry; different n
have orthogonal occupied-slot supports. Summing actual native currents
over the physical object slots commutes with antisymmetrization. Hence

\[
 Q_{\rm total}(T)E_n=E_n\,d\Gamma(\rho_K(T)).
\]

The direct sum E supplies the full finite exterior Fock carrier. Operators
`E c_i E*` obey CAR with identity `EE*` on its code, and complete unitaries
extend as identity on its orthogonal complement. No scalar product or
native current is changed. The executable checks the two-object internal
instance: dimension 105 maps into two fifteen-state source clouds of
tensor dimension 225, with all 1,260 current-column identities verified.
The general formula gives every finite occupation sector.

This encoding uses whole objects, rather than storing one fixed native
cloud vector in each nominal internal CAR slot. That latter proposal
fails: a native current changes the cloud inside its slot, whereas the
target CAR current moves occupation between internal slots. They have
different actions unless another controller and its intertwiner are supplied.
The whole-object encoding is a declared finite architecture. Geometric
locality of arbitrary gates lifted through this particle-slot code does
not follow from the isometry. It supplies no interacting field limit.
Every native cloud is fixed by U_F, so the encoded cloud products also
have U_F identity and native K zero.

## 5. Chirality and mass freedom

Tensoring the four-component Dirac routing space with this actual R module
gives both spacetime chiralities in the same gauge representation. The
beta mass coin mixes them while conserving the currents. This is a
vectorlike Dirac spectrum. An all-left-Weyl convention for the internal
weight table does not remove the opposite spacetime chirality.

The exact commutant of the twelve R currents has five scalar irreducible
blocks, of dimensions `3,2,6,3,1` when ordered by hypercharge. The code
checks each block's commutant rank. Therefore a Dirac mass matrix can be

\[
 M=\sum_{r=1}^5\mu_r P_r
\]

with five independent real coefficients, without violating any of these
current equations. Their values are not fixed by the representation.
Using the unit face-count coin instead gives rest-frequency magnitude two
for every corresponding internal component, and ratios one. Changing its
declared count-to-coin coefficient from one to two gives magnitude four
with exactly the same topology and current action.

For an exclusively left-Weyl R field, all hypercharge pair sums are nonzero,
so a constant symmetric gauge-invariant Weyl mass matrix is zero. The
owned Higgs tensors supply the allowed up, down and charged-lepton channels
with independent coefficients. Their source-field realization requires the
scalar state, a local chiral kinetic/gauge law and the same current action.
In particular K=0 on this native cloud cannot produce those interactions.
A different declared coupling or a controlled departure from this invariant
module is needed. Chirality, mass terms and graded fields cannot be exchanged
between the different constructions by relabeling their weights.

If a source rest frequency omega is associated with a physical duration
tau seconds per source time unit, its mass energy in GeV is
`hbar_(GeV seconds) * omega / tau`. This is an allowed unit attachment.
The conversion constant alone does not select tau or a dimensionless mass
ratio. A specified reference clock can supply that attachment; matching an
experimental mass used to select it is a calibration.

## Finite theorem interfaces

| Existing result | Reused conclusion | Source-specific boundary |
|---|---|---|
| `Screen/ExteriorSelection` and `WeylYukawaConventions` | Exact exterior weights and charge-conjugation conventions | The dark code supplies their actual native current module; spatial chirality is separate |
| `Screen/FermionSectorAssembly` | Charge-compatible algebraic pairings on its declared bundle | Its coefficients and derivative-free functionals do not select this source's kinetic law or masses |
| `Screen/ChargedCurrentDictionary` | Typed current conventions | The 180 current-column checks provide the literal native intertwiner for the selected cloud |
| `standard_model_source_map/FERMION_ATTACHMENT.md` | Full exterior CAR and transferred one-body currents | Whole-object antisymmetrization supplies an admitted carrier; a hard-core slot alone does not |
| `standard_model_source_map/GAUGE_ATTACHMENT.md` | Exact finite quantum-link Ward identities | Attaching those dynamical links changes the source; their nilpotent link operators are not identity backgrounds for a free Dirac theorem |
| `Geometry/MassShellKinematics` | A supplied positive mass and timelike frame satisfy its shell | The theorem takes the mass as input; it does not choose any mu_r or the physical clock scale |
| `QFT/InheritanceMatrix` | Explicit domain limits of its finite-net inheritance claims | It supplies no automatic spin-statistics, vacuum or scattering transfer to a different source |

The finite result is a genuine charge-carrying native source module and a
complete transport for a stated composite. A Standard Model identification
requires a common chiral interacting field carrier, its gauge/scalar state,
physical mass observables and measured predictions. The current embedding
alone selects neither those dynamics nor any measured particle name.
