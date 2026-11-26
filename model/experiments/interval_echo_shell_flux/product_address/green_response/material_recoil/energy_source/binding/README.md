# Does binding change the apparent pull at fixed material content?

Two observers at the receiving face measure less inward displacement when the same two loads are prepared next to each other. Their independently read binding energies are 17, 32.5 and 65 native energy units for ordinary–ordinary, ordinary–double and double–double loads. With the ordinary–double pair, for example, the two probes measure reductions of [0.003769, 0.004674] and [0.003672, 0.004578] receiver links per preparation. The source clock predicts [0.004230, 0.004234]. These are finite observer records from a switched material source, an energy-sensitive release apparatus and three stocked guides. [Frozen protocol](PROTOCOL.md) · [All terminal counts](records.json) · [Simultaneous intervals in integer millionths](results.json)

## Two independent measurements of binding

Fix a resident of load $p$ and enroll a second load $q$ on two exchange-coupled terminals, one adjacent to the resident and one separated from it. These words specify contacts used by the apparatus. The source has two loads in both preparations. Every occupied port has the same degree. Keep the native diagonals on all incident seams and every onsite completion. Switch on the mobile exchange only. Write $e_n,e_f$ for the two sharp configuration energies above a blank apparatus. The two-dimensional Hamiltonian has these diagonals and off-diagonal entries $-1/2$. It acts during an unobserved dwell of 0 or 0.13 exchange-clock units.

Moving the second load next to the resident replaces its two vacuum contact terms by the shared one. The diagonal on a contact is $d(p,q)=\kappa(p-q)^2+\tfrac12{\bf1}_{p\ne q}$. At the active stiffness $\kappa=8$ the binding difference is

$$B=e_f-e_n=d(p,0)+d(0,q)-d(p,q)=16pq+1-\tfrac12{\bf1}_{p\ne q}.$$

The spectrometer mixes the two sharp positions with an isolated half exchange, parks hopping while leaving the full diagonal active for $1/32$ clock unit, and undoes the mixer. The return probability is $[1+\cos(B/32)]/2$. All three gaps lie on its first branch. The independent clock bounds are [16.992996, 17.005283], [32.497001, 32.504326] and [64.996073, 65.003022]. This differential reading is invariant under a common additive shift of both energies. The energy difference compares two preparations in the same particle-number sector. It does not measure a bound-state eigenvalue or inertial mass.

A separate observer measures two incompatible source observables on different copies. One reads the sharp occupations and computes the native sharp diagonal minus half a unit. The other transfers the two mobile positions by isolated swaps to the calibrated antisymmetric exchange reader and parks the resident clear of the phase operation. Swaps preserve relative phase at the active stiffness. The board is connected, so a sequence of adjacent transpositions realizes the required permutation. Every spectator term remote from the mobile trajectory contributes a common phase. In the near–far basis the odd effect is $P_-=(I-X)/2$. Thus

$$\operatorname{diag}(e_n-1/2,e_f-1/2)+P_-=H.$$

Of 1,024 equally likely tickets two use the odd reader. On the other 1,022 tickets the occupied branch's score $2(e_i-1/2)$ decides which tickets click. Every score is a nonnegative integer within the available range. The exact mark rate is $\langle H\rangle/512$ for **every** two-position density matrix and at every free-flight time, by energy conservation. In particular the even and odd preparations have identical sharp probabilities at launch, and the full meter resolves precisely one unit of exchange energy between them. Silencing the odd tickets removes that contrast. The mark probability uses the declared blank-energy convention. The interaction contrast measured by the Ramsey clock is independent of that convention. [Native reader and replay](measure.py) · [Algebra in Lean](../EnergySource.lean) · [Parent odd instrument](../instrument.py)

## What a distant probe can predict

Each click starts one independent guide release. Its first outward face has mean signed displacement $-1/15$ link per pulse for both materials at the registered matched shutter, as measured in the independent receiver experiment. On all starts the near-minus-far displacement is therefore $B/7680$. Its sign is positive because both experiences are inward and the paired source produces the smaller pull. The three parameter-free differences in these apparatus units are 0.002214, 0.004232 and 0.008464. All six held receiver intervals contain their respective predictions. The same pulse normalization multiplies the large-guide-distance inverse-square response of the parent construction. This experiment measures one finite face. [Receiver calibration, counts and guide-chart domain](../../common_receiver/README.md)

There are $2^{31}$ independent synthetic preparations per source setting. All 24 local banks contain both clicked and silent outcomes, three Ramsey banks count both return alternatives, and twelve receiver banks keep both arms and all material outcomes, including null responses. The two receiver species can share a guide history. The audit uses a union bound and requires no independence between their conclusions. Half the familywise 1% budget covers the three source clocks and 48 full or parity-silenced mark means. The other half covers the twelve conditional receiver means using the parent bound for independent guide walks. Bounds include dyadic event-probability rounding and apply to the specified commands. No orbit, lens or astronomical coupling enters. The prepared two-position source, the conversion from a mark to a guided pulse and the receiver length unit are laboratory choices. Inferring the Newton coupling of an uninstrumented freely moving body requires an independent inertial mass, proper distance and the source's physical emission law.

From the repository root run the archived audit with the project Python environment.

```sh
p=model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.energy_source.binding
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m "$p.measure"
python -m pytest -q --pyargs "$p"
```
