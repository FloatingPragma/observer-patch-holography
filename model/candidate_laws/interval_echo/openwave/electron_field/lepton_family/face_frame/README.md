# Three charged lines seen from independently calibrated observer frames

Two observers can agree that the same preparation contains three distinct charged mass lines even when their local coordinates differ. The following field construction combines the exact rank-three port band, the two charged exterior components and a specified mass term. Its three positive rest energies have ratios **1 : 206.77031597 : 3477.47283710**. Local calibration predicts all 24 held electroweak detector intervals. On registered Interval Echo, three complex direction memories on one oriented face give the two winding senses different measured return probabilities. Their free spectra on the cell cannot give the lepton hierarchy. [Field and native comparison](measure.py) · [Joint charge and spin calculation](field.py) · [Two-observer readout](observers.py)

## Which particles does the gauge theorem describe?

Write a component as (colour, weak isospin representation, hypercharge), with electric charge Q = T₃ + Y. One orientation of the exterior construction gives five chiral multiplets. Dimensions count internal states before spacetime spin is attached.

| Multiplet | Colour | Weak | Y | Components in one generation |
| --- | --- | --- | ---: | ---: |
| Left quark doublet | 3 | 2 | 1/6 | 6 |
| Conjugate up quark | conjugate 3 | 1 | −2/3 | 3 |
| Conjugate down quark | conjugate 3 | 1 | 1/3 | 3 |
| Left lepton doublet | 1 | 2 | −1/2 | 2 |
| Conjugate charged lepton | 1 | 1 | 1 | 1 |

On this table the common central subgroup Z₆ acts trivially, giving the faithful group (SU(3) × SU(2) × U(1))/Z₆. In the separate defect-flux dictionary, the charged and neutral members of the left lepton doublet share exterior flux class 3 modulo 6, and the conjugate singlet has class 0. Every generation repeats these classes. A flux read alone therefore leaves the generation undetermined. [Defect-flux and observer dictionary](../../../../../../../code/defect_particle_classes/PROOF.md)

The charged left component is occupation mask 23, its neutral partner is 15, and its conjugate singlet is 24. The weak raising operator transfers 23 to 15 and annihilates 24. All eight colour generators annihilate the three states. The electric charges are −1, 0 and +1 respectively. The right component of the negatively charged Dirac field uses the conjugate of the singlet. This dictionary describes the fermion representations of the stated exterior and gauge construction. A physical protected defect must give the corresponding currents, including the null controls, on that same excitation. [Exterior theorem](../../../../../../../Lean/Screen/ExteriorSelection.lean) · [Charged-current theorem](../../../../../../../Lean/Screen/ChargedCurrentDictionary.lean) · [Multiplicity and weak counting](../../../../../../../code/a5_closure/chiral_multiplicity/README.md)

## A persistent face and the three mass flags

Take the physical twelve-port projector P with adjacency eigenvalue +√5. The implementation verifies its integer tables against the Lean file and binds them to the cell's actual port ordering. For each outward face f, its proper 120° turn R defines three rank-one projectors on the band, Pⱼ(f) = P Σₙ ω⁻ʲⁿRⁿ/3 with n = 0, 1, 2. The mass operator for that face is M(f) = Σⱼ mⱼPⱼ(f), where mⱼ = [1+√2 cos(2/9+2πj/3)]². The phase convention and use of this operator as a mass are assumptions of the field construction. The gauge currents act on the exterior factor. Every proper cell rotation carries the entire projector and mass construction to the rotated face, as checked for all 60 rotations on all 20 faces. [Exact band theorem](../../../../../../../Lean/Screen/A5FamilyBand.lean) · [Koide theorem](../../../../../../../Lean/ObserverPatchHolography/KoideCirculant.lean) · [Spectral and field records](results.json)

An observer need not know which face was prepared. If the same face persists through the history, each of the three energy flags repeats with certainty. Redrawing the face independently between two reads makes all three later flags equiprobable. Averaging the generators first gives a single operator 2P, whereas averaging the records of the twenty persistent frames preserves three lines. This difference is experimentally testable. With the Dirac kinetic term, common units and spin supplied, each line gives two positive-energy spin states of each charge. Squaring the finite one-particle Hamiltonian gives Eⱼ(p)² = |p|² + mⱼ². Joint reads check spin, charge and colour at rest and for three momentum directions. This proves existence of three free charged species **in the declared field**. Their laboratory names and weak decays require evidence from the substrate.

The clock compares an evolved and a stored path of the same charge. Its two output quadratures measure the sum of three phases with positive weights. The determinant of the three-by-three matrix of successive readings, with entries sᵢ₊ⱼ, is abc [(x−y)(x−z)(y−z)]² for phase factors x,y,z and weights a,b,c. It is zero for every stationary two-line clock. The five readings can be taken at a common integer stride. Two independent port preparations yield determinant lower bounds greater than 0.78 at stride 15 with simultaneous 99% coverage over both sources and all sixteen tested strides. The bound encloses errors by multilinearity and the product of the three row norms. [Clock proof](ClockLines.lean)

A six-reading fit poorly resolves the light line. A fresh equally spaced series fits 256 complex readings and holds 32 in reserve. Its frequencies use no mass-dependent settings. Three frequencies predict the held complex amplitudes at root-mean-square errors below 0.00083 for each source. The fitted two-tone comparisons miss by more than 0.57. Rank three is sufficient in this model, and the measured determinant excludes every two-line stationary clock under the stated instrument. The frequencies carry the field's chosen unit. The exact geometric ratios come from the operator. The finite counts have their own accuracy. [Clock counts and held errors](results.json)

## Can two face choices account for quark mixing?

Suppose the left-handed up- and down-quark mass eigenstates each use the three characters of one outward face. Take nondegenerate masses and a charged weak current that acts as the identity on the family factor. Then the CKM matrix is the overlap of the two face bases. Its CP-odd Jarlskog invariant vanishes for any two faces and any choice of the six masses.

In the real rank-three band a face turn fixes its unit axis $n$. Choose a real unit vector $y$ perpendicular to both face axes and put $x_i=y\times n_i$ for each face. Its three normalized characters can be written as $n_i$ and $(x_i\pm i y)/\sqrt{2}$. All overlaps between these two bases are real because $y$ is perpendicular to both $n_i$ and $x_i$. Rephasing either basis leaves each quartet $V_{ik}V_{jl}\overline{V_{il}V_{jk}}$ unchanged. Parallel axes admit the same construction.

For signed axis cosine $c$, the entry magnitudes are $|c|$, $\sqrt{(1-c^2)/2}$ and $(1\pm c)/2$. The icosahedron's face normals give $|c|\in\{1,1/3,\sqrt{5}/3\}$. Among all eight distinct entry magnitudes, the closest to the kaon leptonic determination $|V_{us}|=0.2250\pm0.0004$ is $(3-\sqrt5)/6=0.127322$. Thus even arbitrary permutations of the three species cannot produce that entry. The executable checks all 400 ordered face pairs against the exact menu and all 3,600 quartets. This exclusion applies to these two single-face mass bases and the stated current. [Mixing calculation](mixing.py) · [PDG 2024, section 12.2.2](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-ckm-matrix.pdf#page=3)

## What do Alice and Bob actually compare?

Alice sends a preparation and reports her source setting. Bob reads local charge, a weak-transition flag and an energy channel at his device. Their numerical components can differ. With local changes of basis Gₐ and Gᵦ, the link transporting a state becomes Gᵦ U Gₐ† and Bob's detector effect becomes Gᵦ D Gᵦ†. The scalar event probability is then identical. A family coordinate change acts on a different factor from the gauge transformation. Relabeling a family basis preserves the energy projectors and their probabilities when they are transformed together.

The program lifts two noncommuting weak axes, two colour frames and two hypercharge phases to the full exterior algebra. For a charged left state, reading the old numerical mask at Bob's site gives probability 0.33836. Calibrating the detector to Bob's local Higgs and charge axes gives probability one. The colour check includes a transforming triplet. Six known doublet preparations, read on three axes, calibrate Alice-to-Bob weak transport without using any charged-lepton mass. These counts predict **24/24** held intervals for four new inputs, three new analyzers and both conjugate sectors at simultaneous 99% coverage. Treating their coordinates as identical fits **0/24**. In a separate joint experiment each observer resolves the same energy flag and weak click or refusal for two charged rest preparations at both signs. Route interference compares a supplied electromagnetic loop against an equal-charge reference and reads opposite phases on opposite charges. The change of frame itself contributes zero loop flux. [All calibration and joint counts](observers.json)

Bob also keeps a material record of each individual outcome. His field instrument draws one of six outcomes, combining three energy labels with the presence or absence of a weak response. A supplied encoder writes that classical result once into balanced ternary loads. A second device reads his loads and sends a copy to Alice along five paths of 22 or 23 native bonds. Alice decodes all 4,194,304 reports across four preparations. Blocking the paths leaves Bob’s original intact and gives Alice a refusal on every start. Reversing the face convention without translating the energy labels produces 385,880 to 387,093 disagreements per setting. The experiment tests communication of an orthogonal detector record between owners. The source of its six possible outcomes is the declared field instrument. [Shared event experiment](shared_report.py) · [Joint records and blocked controls](shared_report.json)

## Testing the proposed frame in the native substrate

On one prepared cell the external seams are shuttered and private degree references have zero contrast. One vacancy or its independently stock-balanced conjugate can move on all thirty internal edges. Twenty-seven are blank. On the edges of one oriented face, three even, unit-contrast direction memories form a material record. The source is an ideal projection onto the port band and one of its three face characters with a specified head state. Preparation, a local phase reference for each complex head, and joint return projections are supplied apparatus. The release uses the projected 99-state Interval Echo columns with no mass operator or family-dependent rate. All 96 positions and head words enter the calculation. The two conjugate stock preparations give the same generator in their matched direction conventions.

With three real ready directions, the two opposite cyclic inputs have identical spectral measures. With three identically phased direction bits (|0⟩+i|1⟩)/√2, their spectral weights have total variation distance **0.73415**. After one native time unit, the two inputs return to the original port band and material head state with probabilities **0.70073 and 0.54008**. Conjugating the head preparation interchanges the results. This is a physical intervention. Independently rephasing the two observers' port coordinates preserves all six calibrated probabilities, including refusals. A static change of coordinates supplies no gauge force.

All three initial family vectors couple out of their prepared subspace under free evolution. The exact integer row bounds put every matched-host formation energy between **8 and 45** at the registered contrast. The largest allowed ratio is 5.625. Halving the contrast on the same alphabet and geometry changes the bounds to [4,25]. This test distinguishes a measurable family-sensitive pattern from three freely propagating lepton rest energies. [Native chart and energy records](results.json) · [Native phase-frame comparison](observers.json) · [Related material clock](../../../../../../experiments/interval_echo_lepton_register/phase_feedback/README.md)

From the repository root run the following commands with one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.face_frame
python -m "$p.mixing"
python -m "$p.measure" --write
python -m "$p.observers" --write
python -m "$p.shared_report" --write
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/electron_field/lepton_family/face_frame/ClockLines.lean)
```
