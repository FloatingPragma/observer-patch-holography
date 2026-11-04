# Conditional Standard Model representation and mass action

This theorem composes the finite gauge and matter classification, the
same-current exterior/Higgs construction, and a declared canonically
normalized field action. It determines the charged field content, invariant
mass channels and quadratic dispersion in terms of supplied couplings. It
does not require those couplings to be selected before making these
conditional conclusions. The final section states the additional maps
needed to transfer them to propagating defects of a microscopic source.

## Premises and conclusion

The finite premises are:

1. A `SMStructurePremiseData` instance: the barycentric face weighting,
   compact-classification dimension inputs, and a nonempty chiral
   anomaly-free mask in the committed ten-component exterior grammar.
   To use the executed matter construction, supply the existing
   `native_sm_current` law, its declared spin lift, and the nonvacuum even
   exterior sector. Its exact current intertwiner identifies the abstract
   charge table with this concrete representation. The conjugate mask is
   a second classified choice, not an additional generation.
2. A `SMThreeFamilyPremiseData` instance and its family-factor attachment:
   one complete faithful adjacency band, dimension between three and five,
   and the stated strict seam-cost selection rule. The band acts on the
   family factor and the gauge currents act as the identity there.
3. One complex scalar doublet with the charges of the colorless
   `q=6Y=3` row. Its canonical kinetic term and positive quartic potential
   are supplied. The three charged Yukawa coefficients are arbitrary
   complex matrices `Yu,Yd,Ye in Mat3(C)`. No value, texture or hierarchy
   is assumed for them. The charged quadratic fermion action uses the
   complete same-current Higgs invariant lines, physical CAR statistics,
   the declared Weyl principal symbol, and a canonical positive kinetic
   metric. The flat-space action is a local gauge- and Lorentz-invariant
   polynomial in these fields and their derivatives, with operator dimension
   at most four, modulo integration by parts. No additional light fields
   are included in this minimal target.

Under these premises, the selected finite representation has the five
one-generation multiplets below and three copies of each. Its four
charge-weighted anomaly forms vanish. The three Higgs invariant lines
give charged mass matrices `Ma=w Ya`, where `w` is the neutral Higgs
norm. Consequently the charged action masses are `w sigma_f(Ya)` and
their quadratic frequencies obey `E^2=v_signal^2 |p|^2+m_a,f^2`.
The weak current in the quark mass basis has mixing matrix
`V=UuL^dagger UdL`. The one-doublet canonical bosonic action has one
unbroken electromagnetic direction, eight unbroken color directions,
three broken electroweak directions, and one physical radial scalar
direction. Its tree mass coefficients are derived below.

This is a conditional representation and field-action theorem. A
continuum interacting QFT, a physical vacuum, and a defect-to-field map
are stronger premises than the finite ones just listed. Source selection
of the listed inputs is not part of this implication.

## Gauge type, charges and family count

`OPH.SMStructureComposition.smStructureComposition_receipt` combines the
finite dimension classification with the measure-dependent coarse-family
selection and the anomaly-free mask. Its Lie-type reading consumes the
compact classification input; the arithmetic conclusion is center
dimension one and simple-ideal dimensions three and eight. Its selected
coarse label does not identify a point in a continuous bracket family.
The executable native current supplies a concrete bracket under its
declared complex-versus-real control architecture. It is lifted to the
five-dimensional representation before taking exterior powers.

In `native_matter_attachment.representation`, all 66 pairs of the twelve
source currents are checked under the spin and exterior lifts. Restriction
to `Lambda^2(C5) direct_sum Lambda^4(C5)` gives dimension 15. Joint
primitive-charge and color/weak Casimir eigenrecords yield:

| All-left multiplet | Color representation | Weak representation | `q=6Y` | Complex components |
|---|---|---|---:|---:|
| `Q` | `3` | `2` | `1` | 6 |
| `u^c` | `bar3` | `1` | `-4` | 3 |
| `d^c` | `bar3` | `1` | `2` | 3 |
| `L` | `1` | `2` | `-3` | 2 |
| `e^c` | `1` | `1` | `6` | 1 |

These names abbreviate the derived finite quantum numbers. Reading them
as laboratory quarks and leptons requires the physical current and Spin
dictionary. The table is computed from the same currents, rather than
being supplied to the source record decoder as species labels.

`OPH.ExteriorSelection.selection_unique` checks all 1024 subsets of the
declared ten-component grammar and leaves exactly the two conjugate
parity masks. `masked_generation_count`, `masked_anomaly_forms` and
`masked_z6_kernel` in `SMStructurePremiseData` bind the count, anomalies
and central kernel to the same selected mask. The pure color cubic
anomaly also cancels as `2-1-1=0`; the weak global parity count is four
doublets per generation, as in `witten_automatic`.

The faithful action has the effective quotient `(SU3 x SU2 x U1_q)/Z6`.
The same local weights allow a physical global group with any subgroup
`Gamma<=Z6` in the denominator: the trivial, order-two, order-three and
order-six choices. Selection requires a same-source line/loop attachment. The
auxiliary exterior parity alone is not physical fermion parity.

The family-band conclusion is not inserted as the number three:
`SMThreeFamilyPremiseData.band_identified` derives the rank-three band
from the stated complete-band, faithfulness, window and strict-cost
premises. `three_family_receipt` then gives 45 internal components and
the three-copy anomaly and central-kernel identities on the same mask.
The executable action requires an isometry from that band to its `C3`
family slot; choosing an orthonormal family basis does not choose a
Yukawa matrix. Without the band premise this paragraph is a theorem
for a supplied family multiplicity, not a selection of three.

## Invariant action and same-current Higgs lines

Give H and a derivative dimension one, each Weyl field dimension 3/2,
and curvature dimension two. The local invariant action grammar then has
the following exhaustive reduction. Gauge invariance organizes derivative
terms into covariant derivatives and curvatures. At dimension at most four,
the gauge terms are invariant quadratic curvature forms and their duals.
Distinct simple ideals have no mixed invariant bilinear form; there is one
abelian factor. Positive quadratic forms give the three canonical gauge
kinetic terms with free couplings after normalization. The dual forms carry
free theta coefficients and can depend on the global/topological domain.

The scalar invariants are polynomials in `rho=H^dagger H`, giving
`Lambda0+mu2 rho+lambda rho^2` and the scalar kinetic term. The would-be
antisymmetric contraction of two identical commuting Higgs fields vanishes.
The fermion derivative terms are positive family kinetic matrices within
each irreducible multiplet; whitening them fixes the canonical convention.
There is no gauge-invariant fermion bilinear without H in the selected
chiral representation. Four-fermion terms have dimension six, and a
two-fermion/two-Higgs term has dimension five. A formal two-fermion/curvature
dipole has dimension five; gauge-invariant dipoles in this chiral content
need an additional Higgs field and have dimension six. The remaining
two-fermion/one-Higgs terms are exactly the invariant lines computed below.
Consequently this operator grammar yields the canonical gauge, Higgs and
Yukawa action up to its free couplings, theta terms and a constant. The
broken-potential choice `mu2=-lambda vSM^2<0` is a supplied phase premise;
local covariance does not select its sign or prepare its quantum vacuum.
Higher operators require relaxing the dimension bound, rather than being
excluded by the finite charge table.

For each possible symmetric matter-pair charge origin with total
absolute charge three, `native_matter_attachment.invariant_higgs_lines`
solves the full linear invariance equations

```
Ga^T Ck + Ck Ga + sum_j (ha)jk Cj = 0
```

for every one of the twelve matter currents `Ga` and the corresponding
Higgs current `ha`. Each of the three allowed origins has a
one-dimensional solution space. They are `(-4,1)`, `(1,2)` and
`(-3,6)`, giving the up, down and charged-lepton channels. The first
uses H; the last two use its conjugate in the all-left convention.
This is an exhaustion result within the stated matter-pair/Higgs grammar.

`action_intertwiner` constructs the 15-row permutation between the
source-derived tensors and the existing field action, retaining each
tensor's exact scale and sign. `kinetic_attachment` checks both Higgs
components with 27 independent complex symbols, including every
off-diagonal family entry. Its 90-row map is
`I_spin tensor(P15 tensor I_family)`. It is unitary, commutes with all
three Weyl kinetic matrices, and preserves the CAR pairing. A frame
change therefore cannot silently change the kinetic normalization or
the mass singular values.

Because the gauge action on the family factor is the identity, each
one-dimensional invariant line tensors with `Mat3(C)`. Thus the full
charged Yukawa freedom is exactly three complex 3-by-3 matrices before
family basis changes. In the physical four-component convention these
are `Yu,Yd,Ye`. In the all-left source convention the directed block
is `conjugate(Ya)`, and the reversed block is its transpose, not its
adjoint. `OPH.WeylYukawaConventions` checks the relative conjugation
charges and the weak epsilon contraction. Its mixed-convention controls
are charged, so replacing `u_R` by `u^c` in only half of the notation
does not preserve the theorem.

Choose the neutral scalar `H=(0,w)`, `w>0`. For any channel, a
singular-value decomposition `Ya=UaL diag(ya,f) UaR^dagger` gives
`ma,f=w ya,f`, with `ya,f>=0`. Quark masses repeat in the three color
components. Left-handed quark basis changes transform the same charged
weak current to `UuL^dagger UdL`; its unitarity follows from those two
unitary basis changes. No simultaneous diagonalization of `Yu` and
`Yd` is assumed. Degenerate or zero singular values are retained and
do not define unique family eigenvectors.

The symmetric 45-row all-left mass matrix K has a Dirac pair of
singular values for each charged family/color entry and three neutral
zero directions in the minimal representation. The physical pairing
coefficient is `epsilon_spin tensor conjugate(K)`, an antisymmetric
90-row matrix. For `h(p)=-v_signal sigma.p tensor I45`, the complete
Nambu symbol is

```
B(p) = [[h(p), Delta], [Delta^dagger, -h(-p)^T]].
```

The identity `sigma_j epsilon=-epsilon sigma_j^T` cancels its squared
off-diagonal blocks. The remaining blocks have eigenvalues
`v_signal^2 |p|^2+sigma(K)^2`. This proves the quadratic relativistic
dispersion without fitting a propagation speed separately for each
species. The 90 physical annihilator components and 180 Nambu rows are
not 90 or 180 species; Nambu doubling also includes redundant conjugate
coordinates. Minimal neutral content consists of three left Weyl
fields with zero mass in this charged Yukawa grammar.

An optional symmetric dimension-five coefficient `C5/Lambda` gives
`Mnu=w^2 C5/Lambda=vSM^2 C5/(2 Lambda)` under the explicitly stated
operator convention `C5(L.H)(L.H)/(2 Lambda)+h.c.`. Takagi singular
values then give Majorana action masses. This is a separate supplied
effective operator, not a consequence of the charged three-line scan
and not an assertion of right-handed neutrinos in the 15-component table.

## Canonical bosonic normalization

Use canonical scalar kinetic energy `|D_mu H|^2`, canonical real-vector
kinetic terms, and

```
V(H)=lambda (H^dagger H-vSM^2/2)^2,
lambda>0, vSM>0, w=vSM/sqrt2.
```

`OPH.ElectroweakBreakingComposition.minimizer_iff` proves that the
minimizers are exactly the shell `|H|^2=w^2`. Its
`stabilizer_exactly_one_line` and `goldstone_exactly_three_parameters`
give one unbroken electroweak direction and three gauge-orbit tangent
directions. They are algebraic statements about the declared potential;
the selection of this phase of an interacting quantum theory is a
separate vacuum premise.

The Lean `gaugeMassForm` is the coefficient `|D H|^2`, not a real-vector
mass squared by itself. Its exact expression is

```
w^2/4 * [gW^2(W1^2+W2^2)+(gY B-gW W3)^2].
```

Comparing it with the canonical `1/2 A^T M^2 A` multiplies the form
by two. Substituting `w=vSM/sqrt2` therefore gives

```
mW^2=gW^2 vSM^2/4,
Mneutral^2=(vSM^2/4)*[[gW^2,-gW gY],[-gW gY,gY^2]],
mZ^2=(gW^2+gY^2)vSM^2/4,  mPhoton^2=0.
```

The neutral eigenvectors are proportional to `(gW,-gY)` and `(gY,gW)`.
Color is unbroken because H is a color singlet. With canonical real
radial coordinate `H=(0,(vSM+h)/sqrt2)`, the scalar Hessian gives
`mH^2=2 lambda vSM^2`. The three zero scalar Hessian directions are
gauge-orbit directions, not three additional physical massless scalars.
The generator-dependent gauge couplings and Yukawas remain supplied.

This resolves the normalization difference between the finite Lean
shell radius, the executable neutral Higgs norm and the conventional
`vSM`. If noncanonical positive kinetic matrices are used, they must
first be whitened; mass eigenvalues are then generalized kinetic/mass
eigenvalues, not singular values of an unnormalized coefficient array.
`calibrated_matter` is an explicit empirical instance of the same tensor
construction, through `spectrum.calibrated_inputs()` and `ActionInputs`.
Its conventions are `gY=sqrt(3/5) g1_GUT` and
`lambda_SM=lambda_paper/4`. Its running Yukawas, Higgs norm and mixing
parameters are calibrations. Its tree-comparison potential differs from
the loop-matched potential by a stated comparison prescription. Those
running/tree coefficients are not exact unstable or confined pole masses.

## Local covariance and the stronger continuum premises

`OPH.LocalGaugeJetAction.covariant_first_jet` retains the derivative of
the local gauge parameter and proves its cancellation against the
inhomogeneous connection term. `curvature_second_jet` proves curvature
covariance, and `quadratic_trace_ward` proves the cyclic-trace quadratic
identity. Together with the current-intertwined Higgs tensors they
verify the algebraic gauge-covariance ingredients of the canonical
local action. A source-selected connection or a controlled chiral gauge
regulator is not obtained by these ring identities.

For an interacting continuum Standard Model interpretation, additionally
supply a local regulator/refinement with these same currents and physical
Spin/CAR fields; a positive unitary continuum theory with the required
Ward identities and renormalization; and an appropriate vacuum/phase on
its observable algebra. These are independently testable analytic
premises, rather than a conclusion inferred from the finite matrices.
They must retain all gauge and matter degrees of freedom and specify
running-to-pole matching where masses are compared with observations.
Gauge anomaly arithmetic is necessary but not a construction of that
continuum theory. Gravity and the Einstein equation are separate from
this conditional matter-action result.

An optional gravity attachment can consume the matched regional hypotheses
of [matched_einstein_limit](../matched_einstein_limit/PROOF.md) on this same
total matter action, stress, reference, metric and clock. Its three regional
errors must be `o(ell^4)` and its area coefficient converge to a positive
spacetime constant. That theorem then gives the Einstein equation. A curved
action also admits `xi R H^dagger H` and curvature terms; their treatment
must be specified in its stress and cannot be inferred from this flat-space
operator classification.

## Transport to a microscopic defect sector

The massive-family/defect-family association is an additional conditional
premise. For each massive target family, supply a protected source class,
its allowed finite-energy deformations and an embedding of the relevant
excitation sector. Massless fields are excluded from that association.
Protection of a total topological class does not imply that every
excitation in it is stable: collision, fusion and decay maps must preserve
the declared total class and reproduce the target interactions.

The minimal Higgs shell itself is `S3`: it is the unit sphere in `C2`.
Its `pi1` and `pi2` vanish, so this shell supplies no topologically protected
string or point-monopole class. Its `pi3=Z` labels normalized textures;
that homotopy statement alone supplies no finite-energy stable soliton.
The defect-family premise therefore requires additional source topology
and dynamics. It is not a consequence of the rank-three family band or
a rule equating three generations with three windings.

There is a further constraint if the protecting label is additive and
massless neutrinos carry zero. With all nine quark charged-current entries
nonzero, conservation at `ui-dj+Wminus=0` and `Wminus-ef=0` forces
`up=D-E`, `down=D`, `lepton=Wminus=E` for every generation. Neutral
conjugate-pair vertices give `Z=Higgs=0`. The integer relation group is
therefore `Z^2`; a distinct additive winding for each generation is
incompatible with these interaction premises. If the three-generation
electroweak `(udd nu)^3` process is also supplied, its relation is
`9D-3E=0`. The universal group becomes `Z direct_sum Z3`; integer-valued
charges satisfy `E=3D`, the electric-charge direction, while torsion labels
must be retained. `defect_constraints.py` verifies these integer relations
and their Smith normal forms. This is a constraint on additive labels,
not a classification of nonabelian fusion sectors or a source instanton
construction.

To identify one of these fields with a source defect, provide a carrier
map J from the target field sector into a source sector which preserves
the inner product/CAR structure, intertwines the derived currents and
the Higgs tensors, and binds the target clock and reference to the source
clock and reference. It must also intertwine spatial propagation and
the action; a matching degeneracy or topological index is insufficient.
Topological class spaces or their sector projectors require their own map;
a defect winding label does not equal an electric or color charge by
definition. The catalogue's origin-pinned radial zero doublet and charged
poles have no such map to the 45-component target in this theorem.

There is a useful quantitative consequence when this premise is provided.
For finite self-adjoint generators Hs,Ht, an isometry J and
`||Hs J-J Ht||<=epsilon`, every normalized target eigenvector of energy E
maps to a normalized source vector with residual at most epsilon. Hence
`dist(E,spectrum(Hs))<=epsilon`. More precisely, its source spectral measure
outside `[E-eta,E+eta]`, for `eta>0`, has weight at most `epsilon^2/eta^2`: integrate
`(lambda-E)^2` against that measure and use its lower bound `eta^2` outside
the interval. `matching.check_intertwiner` checks the entire rectangular
residual, including leakage out of the image, and uses the certified bound
`epsilon^2=Tr[(HsJ-JHt)^dagger(HsJ-JHt)]`. It also checks supplied conserved
charge maps `QsJ=JQt`, not just the compression `J^dagger HsJ=Ht`.
For `Im z!=0`, the resolvent identity gives

```
||(z-Hs)^-1 J-J(z-Ht)^-1|| <= epsilon/|Im z|^2,
||J^dagger(z-Hs)^-1 J-(z-Ht)^-1|| <= epsilon/|Im z|^2.
```

Indeed the first difference equals
`(z-Hs)^-1 (Hs J-J Ht) (z-Ht)^-1`; each resolvent has norm at most
`1/|Im z|`. The spectral-distance conclusion follows by the spectral
theorem applied to the residual. When epsilon is zero, the image of J is
reducing for the self-adjoint Hs, so `1_B(Hs)J=J1_B(Ht)` for every Borel
set B. This transfers spectral measures and their probabilities exactly.
For a k-dimensional target eigenspace and `eta>epsilon`, the bound
`||(1-P_interval)J psi||<=epsilon/eta<1` holds for every normalized vector
in that eigenspace. Thus `P_interval J` is injective there and the source
interval has dimension at least k. The bound does not determine exact level multiplicities or
exclude additional source states. Isolated-pole residues additionally
require spectral separation and a state/observable intertwiner; the
resolvent estimate alone does not assert stable particles on the real
axis. Infinite-volume unbounded generators require their corresponding
domain/resolvent convergence statements. A many-particle field attachment
must also transport the full Fock representation, collision and
pair-changing terms where present; a two-quasiparticle calculation does
not supply those maps.

Thus the finite stack supplies concrete charge, Higgs and normalized
mass-action targets for a defect realization, with an explicit error
criterion for transporting spectra. The missing defect and continuum
maps are named premises. The values of the source coefficients can be
selected or calibrated separately without changing this conditional chain.
