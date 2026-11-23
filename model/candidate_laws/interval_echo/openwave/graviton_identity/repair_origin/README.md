# The original graviton claim and the finite repair experiment

An earlier OPH paper explicitly claimed the inevitability of a graviton. Its argument passed from geometric modular flow and fixed-cap generalized-entropy stationarity to a semiclassical Einstein equation. It then asserted that the quantum of the resulting dynamical metric must exist. That final step identified neither an operator and positive spectral weight in the finite repair law nor a preparation and receiver for the proposed state. The paper made the analogous step for the photon after selecting the electromagnetic gauge factor. [Original “Inevitability of photon and graviton” section](https://github.com/FloatingPragma/observer-patch-holography/blob/8051af7bf7b270b48d53bb9837e34a6fff8eb8a6/paper/tex_fragments/PAPER.tex#L4639-L4681)

The accompanying mass-bundle producer wrote a literal zero in GeV for both photon and graviton and labelled it a structural zero. That number therefore came from the branch interpretation. An observer never estimated it from a repair correlation function. [Historical producer and value provenance](https://github.com/FloatingPragma/observer-patch-holography/blob/8051af7bf7b270b48d53bb9837e34a6fff8eb8a6/code/particles/scripts/build_exact_nonhadron_mass_bundle.py#L330-L365)

The July 2026 revision changed that exact passage. The linearized pure Einstein–Hilbert action on Minkowski space gives two transverse-traceless null modes. The revision separately requires a physical metric-field Hilbert space, positive massless one-particle spectral weight and the asymptotic-state hypotheses. The archived string-presentation theorem includes a positive-norm spin-two state and a positive pole among its explicit inputs. These statements locate the additional construction required to identify a graviton in material records. [Revised theorem](https://github.com/FloatingPragma/observer-patch-holography/blob/2c52cb56ad7ef5682f59721afc37ea99869c4e2b/paper/tex_fragments/PAPER.tex#L4774-L4812) · [Historical carrier acceptance criteria](https://github.com/FloatingPragma/observer-patch-holography/blob/2c52cb56ad7ef5682f59721afc37ea99869c4e2b/code/particles/CARRIER_MODE_ACCEPTANCE.md) · [Conditional quantum completion](../../../../../../paper/archive/observer_patch_holography_as_string_vacuum_selector.tex)

## A spin-two orbital representation in the finite material

The registered one-load Hamiltonian has an automatic five-dimensional eigenspace. For a port with unit direction $n=(x,y,z)$ take the five amplitudes proportional to $x^2-y^2$, $2xy$, $2xz$, $2yz$ and $3z^2-1$. They are the components of a symmetric traceless spatial tensor. In the native port labelling we checked all sixty proper rotations. On the identity, the fifteen half-turns, twenty third-turns and twenty-four fifth-turns their characters are respectively 5, 1, −1 and 0. The squared character average is one, so this is an irreducible finite representation. Only one scalar material amplitude per port enters the generator. [Protocol](DESIGN.md) · [Native Hamiltonian](../../../candidate.py) · [Calculation](measure.py)

On a single intact cell the level has energy 51. The separate uniform one-load preparation has energy 48, both in native units. To read the difference, a native pulse sequence prepares a superposition of the uniform reference and each quadrupole with four relative phases. The observer measures all twelve terminal ports at seventeen exposures. Differences of counts for phases separated by π yield the cosine and sine of the evolving relative phase. The finite counts estimate the five angular frequencies as **2.99963, 3.00074, 3.00255, 3.00040 and 3.00007**. The exact common difference is **3**. Disabling hopping leaves all five fitted frequencies within 0.0029 of zero. A finite-volume orbital transition at this energy measures a different quantity from the rest mass of an unconfined particle. [Complete preparations, twelve-outcome counts and controls](results.json)

## What crosses the cell boundary?

Each source cell on these boards meets three intercell seams. A quadrupole supported entirely on that cell remains an exact energy-51 eigenstate whenever its amplitude vanishes at the three gateway ports. For a unit-norm quadrupole amplitude $v$, extending by zero into the other cells gives
$\|(H-51)v\|^2=\tfrac12\sum_{p\;\mathrm{gateway}} |v_p|^2$.
Each external seam contributes a diagonal increment $+1/2$ at the occupied end and a hopping term $-1/2$ at the other. Antipodal ports have identical quadrupole rows. For distinct axes the row Gram matrix has diagonal $5/12$ and off-diagonal $-1/12$. Each set of at most three distinct rows is independent. Thus the dimension of the trapped quadrupole space is exactly five minus the number of distinct gateway axes. [Full-board check on every cell](measure.py)

| Native board | Cells with two trapped modes | Cells with three trapped modes | Adjacent-cell arrivals from a bright quadrupole at exposure 8 | Arrivals from a trapped quadrupole |
| --- | ---: | ---: | ---: | ---: |
| 20 cells | 20 | 0 | 10,133 / 16,384 | 0 / 16,384 |
| 80 cells | 56 | 24 | 10,670 / 16,384 | 0 / 16,384 |
| 320 cells | 216 | 104 | 10,730 / 16,384 | 0 / 16,384 |

The material detector sits in the cell directly across the first bridge of cell zero. Failed arrivals count as quiet. The bright preparation is a quadrupole with amplitude on this bridge. The other preparation lies in the exactly trapped subspace, so its remote arrival probability is zero at every exposure. The experiment establishes a genuine finite spin-two transformation law, an observed orbital frequency and a distinction between stationary and escaping material preparations. The moving pointer reads the material load. Identifying a gravitational mediator requires independent tidal response to a varying stress, the two helicities of a massless propagating field and a positive quantum pole on that same dynamics. [Tidal and material detector controls](../README.md) · [Complete intercell histories](results.json)

From the repository root, use a fresh output path and compare the complete receipt.

```sh
.venv/bin/python -B -W error model/candidate_laws/interval_echo/openwave/graviton_identity/repair_origin/measure.py --output /tmp/oph-quadrupole.json
cmp /tmp/oph-quadrupole.json model/candidate_laws/interval_echo/openwave/graviton_identity/repair_origin/results.json
```

## What survives coarse-graining?

Distribute a single load uniformly over every descendant of one base face and open all contacts. At a stop, each cell can undo a prepared local quadrupole and read its sharp return. Five such ensembles span the local orbital quadrupole. Summing their return frequencies over cells retains the quiet outcomes. The upper bound below applies at every stop. The peak readings use six free-flight appointments. [Instrument and all count records](coarse_response.json) · [Acquisition and replay](coarse_response.py)

| Cells | Descendants in initial region | Native upper bound | Largest native sampled reading | Largest sampled reading with slower internal exchange |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 1 | 4.167% | 1.100% | 13.068% |
| 80 | 4 | 2.083% | 0.664% | 7.281% |
| 320 | 16 | 1.042% | 0.343% | 3.905% |

Let $Q$ project onto the local five-dimensional quadrupole in every cell. The free one-load generator above its common offset is the sum of the internal and crossing half-Laplacians. The internal kinetic operator has eigenvalue three on $Q$ and is positive on its orthogonal complement. The crossing term is positive. Thus $H-48I\geq 3Q$, and the orbital return probability is at most one third of the conserved excess energy. A normalized collective excitation in a spectral interval $[0,\epsilon]$ above the uniform reference has orbital weight at most $\epsilon/3$. By taking the adjoint of the product of projections, the same bound controls the low-energy spectral weight of any normalized preparation in the orbital quadrupole space. The statement uses native exchange units and these normalized one-load operators.

The uniform face covers $4^l$ cells and meets $3\,2^l$ crossing seams. Only these contribute to the initial excess energy, measured with an antisymmetric two-port analyzer. Their total click probability is $2^{-l}/8$. Reducing the internal kinetic coefficient to one sixteenth keeps the same graph and still transmits the load. Its sampled orbital weight exceeds the native bound on each board. This constrains a proposed soft orbital graviton readout. A rescaled tidal operator, a many-load collective state or a distinct metric field would require its own spectral measurement.

Replay from the repository root with `python -m model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response`.
