# Quantum Higgs screening and fermion equations of motion

## Source and domain

The carrier is the full native thirty-mode fermion Fock space at each cell,
the five-state scalar memory from `chiral_higgs/model.py`, and the quantum
links of the declared source Hamiltonian. Every occupation from zero through
thirty is retained. The native internal order, all twelve currents and the
three symmetric Yukawa tensors are consumed from
`physical_attachment/chiral.py`. No Higgs expectation value, measured mass,
or empirical mass ratio enters these identities.

Write `c_{alpha,i}` for the annihilator with mode `15 alpha+i`, and let
`epsilon=[[0,1],[-1,0]]`. The scalar operators are

\[
 \phi_k=|0\rangle\langle1+k|+|3+k\rangle\langle0|,
 \qquad k=0,1.
\]

For channel `a`, set `Phi^a_k=phi_k` for up and
`Phi^a_k=phi_k^dagger` for down and charged leptons. The exact local operator is

\[
 H_Y=\sum_{a,k,i,j} y_a C^{a,k}_{ij}
       c_{0,i}c_{1,j}\Phi^a_k+\mathrm{h.c.}
\]

The receipt checks the equality of this polynomial to every word of the
actual periodic-family Yukawa compiler. Coefficients `y_u,y_d,y_e` are
independent real architecture parameters in the executable checks. The same
creation equation holds with complex coefficients in the pair term; its
adjoint contains their conjugates.

## Screened channels with inherited signs

Define

\[
 J^a_{i,\beta}=\sum_{k,j}C^{a,k}_{ij}\Phi^a_k c_{\beta,j},
 \qquad O^a_{i,\alpha}=-\sum_\beta\epsilon_{\alpha\beta}J^a_{i,\beta}.
\]

The table suppresses the common spin component. The indices refer to the
actual native exterior basis, rather than an independently relabeled color
basis.

| Channel | Singlet index | Exact `J` |
| --- | ---: | --- |
| up | 0 | `phi0 c8 - phi1 c7` |
| up | 1 | `-phi0 c6 + phi1 c5` |
| up | 4 | `phi0 c3 - phi1 c2` |
| down | 12 | `phi0† c7 + phi1† c8` |
| down | 13 | `-phi0† c5 - phi1† c6` |
| down | 14 | `phi0† c2 + phi1† c3` |
| charged lepton | 9 | `phi0† c10 + phi1† c11` |

For all twelve native currents `G_A`, the complete operator identity is

\[
 [G_A,O^a_{i,\alpha}]=\sum_j(T_A)_{ji}O^a_{j,\alpha}.
\]

Thus each composite has the current representation of the corresponding
conjugate-singlet creator. All three weak-current commutators vanish for
these rows. The hypercharges of the up, down and charged-lepton composites
are respectively `-2/3`, `1/3` and `1`. Up and down retain their color
representation; the charged-lepton operator retains hypercharge. An
apparatus or link dressing must compensate these residual local charges
before it defines a gauge-invariant preparation.

The neutral contraction

\[
 J_{\nu,\alpha}=\phi_0 c_{\alpha,11}-\phi_1c_{\alpha,10}
\]

commutes with all twelve gauge currents. The native fifteen-dimensional
family contains no neutral conjugate-singlet creator and no associated
Yukawa coefficient. This is an exact current statement. It establishes
neither a massless neutrino excitation nor a propagating neutrino pole.

For `S_a=dGamma(sigma_a/2)`, the executable verifies

\[
 [S_a,O_{i,\alpha}]=\sum_\beta(\sigma_a/2)_{\beta\alpha}O_{i,\beta}.
\]

The uncontracted neutral `J` has the annihilator spin transformation. The
antisymmetric tensor obeys `B^T epsilon B=epsilon` for `B in SL(2,C)`.
These are local tensor identities. A physical little-group representation,
boost symmetry of the cutoff Hamiltonian, and a relativistic dispersion law
require separate dynamical results.

## Full-domain equations of motion

The canonical anticommutators give

\[
 [c_a c_b,c_p^\dagger]=\delta_{bp}c_a-\delta_{ap}c_b,
 \qquad
 [H_Y,c_{\alpha,i}^\dagger]=\sum_a y_a O^a_{i,\alpha}.
\]

`model.py` compares both sides in the independent basis of normally ordered
fermion monomials times scalar matrix units. Consequently the thirty checked
equalities hold on all `5*2^30` local states, including fully occupied
inputs. They are transferred to the complete graded tensor carrier at every
cell. The twelve Ward identities `[G_A,H_Y]=0` are verified in the same
basis. The implementation also compares 495 independently evaluated columns
with the inherited full native compiler, including all 320 scalar/lepton
occupation columns and occupied quark controls.

For every path contribution `c_p^dagger A_pq c_q`, where `A_pq` is the actual
quantum-link transporter and spin coefficient,

\[
 [c_p^\dagger A_{pq}c_q,c_r^\dagger]
       =\delta_{qr}c_p^\dagger A_{pq}.
\]

Scalar and link penalties commute with the elementary fermion creator.
Thus for the full three-dimensional periodic family,

\[
 [H,c_{x,\alpha,i}^\dagger]
    =[H_{\rm hopping},c_{x,\alpha,i}^\dagger]
      +\sum_a y_a O^a_{x,i,\alpha}.
\]

The Hamiltonian equation is `dO/dt=i[H,O]` in the declared clock with
`hbar=1`. The hopping term retains every oriented path, link matrix unit,
spin coefficient and antiperiodic boundary sign of the source. The identity
does not set a link to its expectation value.

The composite's own equation follows from the exact Leibniz rule:

\[
 [H,\Phi c]=[H,\Phi]c+\Phi[H,c],\qquad
 [H,\Phi]=\omega[N_H,\Phi]+[H_Y,\Phi].
\]

The second scalar commutator contains the actual finite-memory matrices.
For all seven charged/color rows, it produces nonzero cubic fermion terms;
the receipt gives explicit witnesses. The neutral composite also has a
nonzero equation of motion with cubic terms. A closed free Dirac equation
and a massless neutral pole do not follow by discarding these terms.

The composites are fermion odd but fail the canonical anticommutation
relations. For example, the diagonal values of `{J_e,J_e^dagger}` on empty
fermions and scalar labels `0,1,2,3,4` are `2,1,1,0,0`; on filled fermions
they are `2,0,0,1,1`. The receipt retains analogous witnesses for up, down
and neutral channels. No uniform rescaling turns these operators into
canonical fermion creators on the complete carrier.

## Local parity and separated correlations

The exact dressed local parity from
[`field_limit`](../field_limit/PROOF.md) is

\[
 C_x=(-1)^{N_{\mathrm{matter},x}}
       \prod_{e\ \mathrm{incident\ to}\ x}(-1)^{N_{\mathrm{link},e}},
 \qquad [H,C_x]=0.
\]

Every composite above contains one local fermion and changes no link. Thus
`C_x J_y C_x=(-1)^(x=y) J_y`, including the neutral gauge singlet. For a
reference with `[rho,C_x]=0`, conjugating the correlation by `C_x` gives

\[
 \operatorname{Tr}(\rho J_x(t)J_y^\dagger)
 =-\operatorname{Tr}(\rho J_x(t)J_y^\dagger)=0,\qquad x\ne y.
\]

The hypotheses hold for a Gibbs reference and a unique ground projector.
The result is a selection rule for these local operators on this source.
It does not exclude transport correlations with retained link dressings.
A transporting operator must account for the endpoint star characters and
actual gauge-current contractions of its link string; a physical dressed
neutrino propagation theorem is outside this calculation. The same-cell
spectral overlap below obeys this selection rule.

## Same-H periodic spectral bridge

The spectral calculation uses the separate two-cell periodic lepton control
with two six-state oriented links, two five-state scalars, all twelve active
CAR modes, `y_e=omega=kappa=1` and `y_u=y_d=0`. It consumes the complete
neutral Gauss ground of `periodic_spectrum`, with `L=2`, `2Sz=0` and
translation involution `+1`. This is a different source from the full
three-dimensional family. The operator theorem transfers; its computed
ground and spectrum do not.

At cell one and spin zero let `T=c_ec^dagger` on the add branch and `T=-c_ec`
on the remove branch, matching the inherited reservoir convention. Let
`J=[H_Y,T]` and `K=[H_hopping,T]`. The code constructs these two maps
independently from scalar/CAR operations and quantum-link hopping columns,
then projects each into the exact charged Gauss basis. With the source and
target coordinate Hamiltonians,

\[
 H_t T-T H_0=K+J.
\]

This exact matrix equality covers the entire 251-dimensional source block
and each 330-dimensional charged target block. No approximate eigenvector
enters its verification. For exact eigenvectors it gives

\[
 (E_n-E_0)\langle n|T|0\rangle
 =\langle n|K|0\rangle+\langle n|J|0\rangle.
\]

Every target eigenstate is retained through the complete isolated spectral
projectors and their ranks. Residue sums and cross residues within each
projector are invariant under changes of basis inside a degenerate eigenspace;
individual degenerate eigenvectors are not physical alternatives here. At the first transfer,
approximately `0.97009612837` in source energy units, the elementary residue
is approximately `0.44499394005`, the Higgs-composite residue is
`0.42308412196`, and the hopping residue is `0.00001101405` in both branches.
The residues are squared operator matrix elements. Their interference obeys
the amplitude equation; the individual squared values are not additive.
They are not probabilities for an implemented composite preparation.

The complete residue sums are approximately `0.5`, `0.89476914311` and
`0.02888543770`. The complete numerical amplitude identity has error below
`5e-15`. These decimal values are binary64 diagnostics. The first composite
residue also has an outward rational enclosure in the receipt, using the
inherited complete energy enclosures, rational approximate eigenvectors and
projector error bounds. For both branches the interval is
`[56785389/134217728,28392695/67108864]`, whose lower endpoint is positive. For the error
estimate, the two exact norm-one scalar/CAR summands give `||J||<=2`;
`delta_ground` bounds pure-state trace distance, and
`8 delta_ground+4 delta_projector` conservatively bounds the expectation error.

This proves that the elementary and quantum-Higgs composite operators couple
to the same isolated finite-source charged energy line. A physical rest
mass, a continuum pole and a conversion to GeV are outside this result.

## Relation to gauge-invariant Higgs constructions

Gauge-invariant fermion composites and their relation to elementary spectral
channels have an established field-theory context. Afferrante, Maas,
Sondenheimer and Törek test the Fröhlich-Morchio-Strocchi mechanism in a
vectorial lattice model in [*Testing the mechanism of lepton compositness*,
SciPost Physics 10, 062 (2021)](https://arxiv.org/abs/2011.02301). Their
spectral conclusions concern that model. Here the exact finite-operator
commutator and same-H residue identity provide the attachment; no
weak-field expansion or imported pole-mass equality is used.

Only one native generation occurs in this carrier. Three observed
generations, quark confinement, massive neutrino mixing and empirical mass
identities require further source and observation structure. None is
inferred from the current representation or the count of composites.

## Reproduction and source custody

The receipt binds every package Python/prose file, the transitive local
Python import closure, native representation inputs, the consumed periodic
spectrum receipt and its pinned sources, and the local-parity proof.
`--write` regenerates the receipt; `--check` validates its pins and performs
an independent read-only numerical and exact-algebra replay. The consumed
periodic receipt's own input hashes must match before either operation.

Certified interval endpoints, operator identities and custody hashes compare
exactly during replay. Binary64 cluster diagnostics use absolute tolerance
`1e-10`; a mutation test rejects changes to certified bounds and hashes.
The first certified positive residue is independent of that diagnostic
tolerance.
