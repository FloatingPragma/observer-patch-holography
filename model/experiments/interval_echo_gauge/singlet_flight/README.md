# Colour during orbital transport

On the registered one-load current instrument the first three-dimensional orbital band carries the fundamental colour action. A terminal energy measurement and an independently compiled colour-centre interferometer both resolve this band in an escaping packet. In the leading arrival from a state initially supported on one complete cell, exactly one quarter of the probability is in this band. A receiver can still herald the complementary colour outcome and recover orbital quantum information. The single-load colour action is trivial on the complement. Its use for a vacancy or an excess on a uniform material background is a specified identification of the orbital meter. [Protocol](PROTOCOL.md) · [All synthetic counts](results.json)

The simulated band instrument applies the specified isolated-cell projector P as a terminal orbital measurement. The one-load energy on a standalone degree-five cell realizes this projector. On a filled federation, cutting bridge hopping while keeping the degree-reference potential leaves a nonuniform energy on the three gateways. An energy-only implementation there would also have to guard that offset and isolate the interior Laplacian. The signed-sector runs assume the calibrated P analyser. The dark-port probability alone cannot distinguish a triplet from its conjugate representation. The closed-bridge check concerns probability transfer alone; it does not establish this analyser. The frozen protocol and raw counts identify the observable originally computed.

## A local test that holds for every launch

Write A for the cell adjacency and J for its antipode permutation. The spectral projector onto the first triplet of the isolated cell is P = (I − J)(I + A/√5)/4. Each of the twelve sharp ports has squared projection 1/4. The first nonzero power of the native Hamiltonian to put amplitude in a neighbouring cell acts on a vector still confined to the source cell. Each neighbour has exactly one bridge to the source. The arriving vector at this order is proportional to that bridge's sharp port. Consequently every localized launch that escapes has the same leading colour fraction. At finite exposure distinct paths interfere. On 320 cells at exposure one, the projected-gateway launch gives triplet fractions among arrivals of 0.169, 0.216 and 0.287 for an isolated load, vacancy and excess respectively. The archive includes the uniform-cell preparations and the exchange-only comparison. [Calculation](measure.py)

## Constraint on all-times pointwise singlets

If P annihilates the amplitude in every cell at every time, it also annihilates the time derivative. The interior Laplacian commutes with P. For each of a cell's three gateway ports p, let q be its opposite bridge endpoint. The remaining tangent equation sums the three columns P|p⟩ with coefficients αψ(p) − ψ(q). Here α = 1 in the isolated-load sector and α = 1 − 4κs for a unit-filled background with signed defect s. The column Gram matrix has diagonal 1/4, antipodal entries −1/4, and all other off-diagonal entries ±1/(4√5). Without an antipodal pair it is strictly diagonally dominant. Thus each coefficient vanishes. With an antipodal pair the two coefficients agree, and the third vanishes. A zero fixes a real proportionality between endpoints, forcing the imaginary hopping current to vanish.

| Cells | Bridges | Bridges forced by one endpoint | Largest connected set left for signed flow | Diameter in unexcluded bridges |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 30 | 30 | 1 | 0 |
| 80 | 120 | 108 | 2 | 1 |
| 320 | 480 | 432 | 4 | 3 |

For an isolated load, the coefficient at the far end is the negative of the near-end coefficient. Propagating forced zeros through the antipodal equalities exhausts every bridge on all three boards. Every pointwise-singlet trajectory in that sector has constant cell populations. For the vacancy and excess the table gives upper bounds on the possible connected flow. The signed-sector calculation does not establish flow on the unexcluded bridges. The argument pertains to a pointwise local singlet in these sectors. A collective excitation carrying a different representation would require its own test. [Finite-board incidence enumeration](measure.py) · [Gateway algebra](../../../../Lean/SourceLaw/ColourSingletGateway.lean)

## Singlet readings at both endpoints

A second experiment transports two coherent orbital modes of a single signed defect. At both endpoints the local apparatus selects the complement of P. Six orbital preparations and three receiver settings reconstruct the channel and hence the output it would produce with a retained maximally entangled reference. All twelve sign, distance and board combinations resolve a negative eigenvalue of that reference and receiver matrix after partial transposition. On the 320-cell board, two cell steps away, 0.576% of original vacancy preparations and 0.655% of excess preparations pass both filters and reach the two-mode decoder. The measured magnitudes of the negative eigenvalue are 0.0001192 and 0.0002619. The simultaneous resolution threshold is 0.00000692. These are unconditional channel weights, including failed preparations. Even a coherent singlet reading at launch and receipt does not establish a singlet throughout the flight. The same port modes give an ordinary orbital phase under a full spatial turn. [Joint protocol](ORIENTATION_PROTOCOL.md) · [Independent count panels](orientation.json) · [Rotation and exchange](../../interval_echo_angular/charged_orientation/README.md)

Run the two modules from the repository root with an unused destination for each receipt.

```sh
python -m model.experiments.interval_echo_gauge.singlet_flight.measure --output /tmp/colour-flight.json
python -m model.experiments.interval_echo_gauge.singlet_flight.orientation --output /tmp/colour-orientation.json
```
