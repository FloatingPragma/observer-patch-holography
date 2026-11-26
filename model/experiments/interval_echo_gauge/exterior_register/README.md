# A five-mode exterior algebra read at two stations

Two observers on the twenty-cell Interval Echo federation can implement the declared five-mode matter grammar on their material registers. Both observers resolve every one of 32 occupations. Each identifies the same ten nontrivial components and reports exactly two admissible subsets among all 1,024 choices. Coherent insertions on the source and a phase read at the receiver also test the relative sign of exterior multiplication. All 80 paired trials give opposite deterministic readings for signed and ordinary insertion. [Fixed protocol](PROTOCOL.md) · [Complete receiver words, refusals and ballots](results.json)

## The ten rails

Cells 13 and 0 each have five disjoint internal dimers. One vacancy per source dimer encodes a bit. A vacancy on its first rail means absent, and a vacancy on the marked rail means occupied. Each marked rail has a physical degree-five reference. Five remote doubles compensate the five vacancies. The sender can prepare every sharp code word from the same stock using addressed full exchanges. The receiver starts with ten ordinary filled rails.

Each observer compares each rail successively with load-zero and load-two standards through an echo, an agreement read and the inverse echo. The pair of clicks is 01 for a vacancy, 00 for a single and 10 for a double. Dimer words 0100 and 0001 encode absence and presence. Other four-bit words refuse. An observer can use the same comparisons as a nondestructive projection onto either value of a specified dimer. The source binding checks the relevant columns of the complete 99-state seam operators, including both orders of the endpoint and reference. [Material gates, routes and records](material.py) · [Sharp echo and route certificate](../../../candidate_laws/interval_echo/openwave/electron_field/lepton_family/face_frame/shared_report.py)

Ten endpoint transpositions carry all ten rails. There are 438 scheduled contacts, including 140 across cell boundaries. At each contact the native blank-memory propagator at duration π exchanges endpoint loads with unit amplitude. Thus the verified permutation of the 240 port labels also transports all coherent superpositions and entanglement with an untouched reference. An actual read at the sender is included on the 32 sharp calibration flights. The phase-sensitive flights use fresh coherent inputs. With external exchange suppressed, the identical sequence of appointments acts within cells. All 32 blocked runs give twenty zero clicks at the receiver, an invalid register. The receiver uses the same standards in both arms.

## Exterior insertion and block rotations

Number the five modes from 0 to 4. The addition instrument at j first reads that j is absent. A full dimer exchange flips j, and a reference dwell π/8 gives a relative minus sign on the occupied branch of each preceding dimer. The removed common phase is fixed by the empty-mode response. Acting on a subset S, the success operator gives zero if j belongs to S. Otherwise it gives

$$
a_j^\dagger |S\rangle
= (-1)^{|\{i\in S:i<j\}|}\,|S\cup\{j\}\rangle.
$$

Removal uses the inverse operations and the presence read. These formulas give the canonical anticommutation relations on all 32 dimensions. Two successful additions in opposite orders differ by one prefix sign. At a single index, adding after removing projects onto presence and removing after adding projects onto absence. The four combinations at distinct indices cancel pairwise in their respective anticommutators. All matrix columns are checked against the native dimer rotations and zero branches. [Operator and interference checks](coherent.py)

For each i < j prepare i in the plus state and j absent. Repeat for each of the eight configurations of the three spectators. After insertion, the receiver undoes the plus preparation on i and reads every rail. In all 80 cases it reads i present with exterior insertion and i absent when the prefix pulses are omitted. Every other receiver bit agrees across arms. In 80 separate attempts to add at an occupied j the source records the presence word 0001 and refuses. The receipt lists the unit-probability raw outcomes of these exact finite instruments. All other terminal words have probability zero.

Take V to be the one-occupation subspace and identify each ordered basis wedge with its corresponding register word. The declared colour and weak subspaces are the first three and last two modes. Their block actions are executable. For neighbouring modes, native endpoint transpositions bring the two marked rails to an actual adjacent pair of degree-five ports and restore the transit registers. A refocused contact on this pair gives controlled-Z up to one scalar. A dimer basis change makes controlled-X in either direction. [Refocusing identity and positive-time pulse words](../../interval_echo_consensus/signed_effects/RECORD_TRANSFER.md)

With mode a first and b second, controlled-X from a to b maps the single-occupation states to the b-occupied subspace. Rotate a conditionally there, then undo the first gate. The conditional rotation uses two opposite half rotations of a separated by controlled-X from b to a and a closing controlled-X. For b vacant the rotations cancel. For b occupied, conjugation by X negates the rotation angle, so the rotations add. An isolated native dimer gives rotations about X with a continuous dwell, and private quarter phases convert them to real plane rotations. The resulting two-mode unitary fixes the vacuum and top wedge and acts by the defining real rotation on the two single occupations. The proof applies to every real angle. All four columns are compared at π/5 and −π/7, with physical connector paths for each of the three adjacent within-block pairs.

An addressed phase acts on an occupied mode, hence on precisely those wedges containing that mode. Adjacent real rotations and such phases generate each complex unitary block by elimination of successive matrix entries. Restricting the determinants to one gives the special-unitary actions. A one-parameter diagonal word with relative phases −2θ on the colour modes and 3θ on the weak modes acts with integer weight q = −2c + 3w for c colour and w weak occupations. These are declared charge units and mode assignments. The pulse algebra realizes their exterior representation.

## Independent ballots on all component subsets

At each station the 32 measured words group by (c,w). A group supports the component Λᶜℂ³ ⊗ Λʷℂ². The all-empty and all-full words are the two neutral trivial lines and are retained in the physical census. Removing these lines leaves ten components. Each observer computes the dimension and charge of each component from its own local records. Complementation maps (c,w) to (3−c,2−w).

For colour and weak dimensions d and e, a row contributes (deq, eq if d = 3 else 0, dq if e = 2 else 0, deq³). These integer weights use the normalization of the finite selection theorem. A subset passes if it is nonempty, contains no conjugate pair and has zero total in each of the four columns. Each observer applies the same declared decision rule. The complete ballots coincide. Each nontrivial component has a distinct integer charge, which fixes its index in the theorem independently of the anomaly test. Station indices 0 through 9 map to theorem indices 0, 1, 2, 3, 4, 7, 6, 5, 9, 8. This permutation preserves both degree-parity masks. In the theorem's indexing, masks 227 and 796 contain respectively the odd and even nontrivial degrees. Each contains fifteen occupation words. Complementation interchanges the two sectors. Toggling any one selected or unselected component of either surviving mask fails the rule in all twenty challenges. [Census and integer scan](selection.py) · [Independent component and ballot comparison](test_exterior_register.py) · [Finite selection theorem](../../../../Lean/Screen/ExteriorSelection.lean)

This realizes the declared finite grammar of PR-59 on a programmed apparatus of the main candidate. The two block labels, primitive charges, projection convention and subset filter are inputs. Claims about dynamically selected matter require a physical Spin and chirality attachment, protected conserved defects, and a common gauge field and current. The measured sign concerns an encoded exterior insertion into specified dimers. Physical exchange of freely moving defects requires a separate experiment.

From the repository root, with one numerical thread:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
nice -n 15 .venv/bin/python -m model.experiments.interval_echo_gauge.exterior_register.measure
nice -n 15 .venv/bin/python -m pytest -q model/experiments/interval_echo_gauge/exterior_register
(cd Lean && nice -n 15 lake env lean -j1 Screen/ExteriorSelection.lean)
```
