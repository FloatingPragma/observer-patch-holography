# The lowest two native energy lines carry a wait-dependent receiver record

On the same registered 20-cell source, a single admitted energy gives a stationary and separable receiver record. Admitting the lowest two energies under one binary flag retains a local reference correlation that changes with free flight, for both signs of material charge. The receiver measures the charge at port 1 and its written coin, with the source seam and all other records discarded. A joint reference and receiver witness changes sign after a wait of 32. [Design](DESIGN.md) · [Replay](measure.py) · [Receipt](results.json) · [Parent preparation and coin analyzers](../README.md)

## Same flight and read for every window

Use the 630-contact registered generator at \(\kappa=8\) with the original degree references. The written absolute-one collars are at 74–82 and 0–1. The source Bell qubit enters the first coin at port 82 after the two original collar commissioning reads. Both collars are prepared in their even reversal sectors. The binary energy filter acts on the same physical generator for all four reversal blocks. No reversal label or individual energy is detected as part of the successful flight. The charge and coin detector is the fixed port-1 joint read on the receiver’s marked seam.

The cutoff spans the bottom of the spectrum through the selected gap. The predeclared rank targets are 1, 2, 4, 8, 16 and 32. Verified inertia at both endpoints of each gap counts every state across the four reversal sectors, so the chosen cutoffs retain whole eigenspaces. The receipt checks all waits \(0,1,2,4,8,16,32\) and a calibrated near-half-beat wait. The higher-window probabilities and witness eigenvalues in the receipt are numerical comparisons. The rank-one exclusion and the rank-two witness have mathematical certificates.

| Full rank under the cutoff | 1 | 2 | 4 | 8 | 16 | 32 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Vacancy rank in the populated sector | 1 | 2 | 4 | 4 | 4 | 8 |
| Excess rank in the populated sector | 1 | 2 | 3 | 4 | 4 | 8 |

For vacancies the first two energies are approximately 36.9823076436 and 37.0025463767. For excesses they are 53.0157702720 and 53.4362513509. These are native common-zero frequencies relative to the blank unit-filled board. The interval above each pair has certified whole-space rank two. The vacant odd sectors above these cutoffs contribute exactly zero events on this prepared bank.

## Witness on the same denominator

Let \(F\) be the binary success operator of the ideal filter, \(S\) the two-column embeddings of the commissioned coin, and \(C\) the four-row sharp charge and receiver-coin aperture. At wait \(t\) the capture matrix is \(C e^{-itH}FS\). The sum over the sender coin is taken at the density matrix level, with an untouched two-state source reference. The unnormalized reference and receiver density uses \(1/8\) in each outer product, incorporating the \(1/4\) commissioning success and the Bell normalization. Per sign, commissioning refusal, filter refusal, admitted but missing port, and admitted port sum to one. The energy acceptance is independent of the subsequent wait.

| Rank-two record at wait 32, per original start | Vacancy | Excess |
| --- | ---: | ---: |
| Commissioning refusal | 0.75 | 0.75 |
| Coarse energy refusal | 0.2460798905194 | 0.1903389664709 |
| In the band, port miss | 0.0038375419383 | 0.0317248433931 |
| Charge at port 1 with receiver coin read | 0.0000825675423 | 0.0279361901360 |
| Same witness at zero wait, enclosure | \([8.8272\times10^{-13},1.35554\times10^{-12}]\) | \([2.08376\times10^{-10},2.20113\times10^{-10}]\) |
| Witness at wait 32, enclosure | \([-1.68664\times10^{-12},-1.10670\times10^{-12}]\) | \([-1.55879\times10^{-10},-1.41831\times10^{-10}]\) |

The witness in each sign bank is a fixed normalized integer vector for the partially transposed joint density. Its expectation is a linear combination of reference and receiver Pauli outcomes, with failures assigned zero. The negative upper endpoint establishes entanglement on the full outcome branch. The positive lower endpoint at zero wait for that same witness establishes a changed physical record under free exposure. The receipt certifies a wider separation at the near-half-beat waits \(19869/128\) for the vacancy and \(7651/1024\) for the excess. Their relative spectral phases are within \(1/1024\) radians of \(\pi\).

The zero-wait coarse projection has its own negative witness, with upper endpoints below \(-1.1602\times10^{-12}\) and \(-1.4298\times10^{-10}\). The energy filter is spatially nonlocal and can furnish an entangled zero-wait record before material has had a free-flight interval. The result is the measured change between two uses of the same filter, source and detector. In a fine-energy comparator, project onto each of the two simple lines and forget the line label. Each term factors the material from the reference. Their sum is separable, and its phases cancel for every elapsed wait. The one-line window is a single instance of the same exact obstruction. This makes rank two the minimum in the tested nested sequence.

## Finite bounds and instrument

Each populated eigenvector is rounded to a dyadic frame. Integer multiplication by \(2H\) encloses its spectral residual and its Gram error. The certified empty interval separates the selected spectrum from its complement. The spectral theorem bounds the complementary frame component by residual divided by separation. Orthonormalizing that frame gives the projector error used in the calculation. Integrating the residual along exact free evolution bounds the phase-evolved frame. Rational phase Taylor enclosures and Gaussian-integer witness expectations give the stated intervals after restriction to the receiver and exterior trace. The fixed Pauli bank would need experimental uncertainty below the witness scale to reproduce this calculation with sampled outcomes.

A finite coherent filter can realize the binary read within the stated margins. Use controlled evolution under the same native generator at elementary exposure \(1/64\), 22 phase bits per register, and 51 coherent repetitions. Compute the majority admission bit, uncompute the fine work, and retain success only when the qubits return to their initial values. The accepted amplitude for an eigenvalue is the majority admission probability. Inertia and a Gershgorin enclosure verify a phase margin and exclude aliasing. A phase-estimation tail bound and binomial estimate give success-operator errors below \(3.43\times10^{-25}\) and \(8.86\times10^{-24}\). These also bound a possible shift of each original-start witness, and all signs survive. The controlled clock costs a total forward and inverse exposure of \(213909453/32\) native units in this circuit. The controlled generator, coherent auxiliary clock, postselection and protected analyzers are supplied. The circuit is bounded analytically. The reported records are computed from the ideal projector with explicit error allowances.

The freely changing record is a coin encoded on two fixed written collars. The calibrated read supplies no pair of perpendicular physical rotations of that captured material, and the two energies are separately resolved by the external calibration. No spin degeneracy, intrinsic spin-½, exchange sign, or identification with the physical electron follows from this comparison.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.theorem_interfaces.material_spin_band.low_window.measure
~~~
