# Can the arriving charge be converted into a separately detectable neutral partner?

**A vacancy crossing a cell boundary can exchange one unit with a prepared material donor and emerge as a neutral vacancy–excess pair.** The receiver acts on the arriving one-defect sector. It resolves the reverse transition from an independently prepared neutral pair and leaves a positive defect arriving at the same address untouched. The driven apparatus preserves an arbitrary coherent vacancy state on an entire twelve-port cell. A second bench performs the conversion under a constant Hamiltonian after one encoding pulse. These results realize a regional charge ladder. A physical electroweak assignment requires an independently attached current and chirality on the same propagating mass branches. [Protocol](PROTOCOL.md) · [Driven records](results.json) · [Constant-contact records](static.json)

## What the owners record

Use the twenty-cell federation, with one load per flight port. The source at 82 releases a vacancy or excess. The entry owner reads signed load at port 1 after a quarter exchange unit. All incoming flight contacts outside the parked workstation are open. The original private degree references contribute to the flight Hamiltonian. Addressed exchanges along the route 1–8–3 relays the captured load to the workstation. For a vacancy the complete ideal capture probability is 0.00240853. For an excess it is 0.00244222. Misses are included in every end-to-end rate.

Write the work word as (incoming occupation at 3, target at 2, donor at 10). The measured regional charge on the first two ports is their load sum minus two. The negative word is (0,1,2), its independently prepared neutral partner is (0,2,1), and the positive control is (2,1,2). The latter uses the same physical arrival port. Two supplied compensating vacancies balance the double reference and donor during preparation. The signed echo against a zero reference gives a sharp load read on all three possible occupations. The donor is also read, so an undisturbed pair of units cannot be called an outgoing neutral particle.

| Input and read | Native result at the full driven setting |
| --- | --- |
| Free vacancy arrival, then charged-to-neutral conversion | 10,088 neutral pairs from 4,194,304 original launches |
| Independently prepared neutral pair, then inverse conversion | 4,194,304 charged outputs from 4,194,304 preparations |
| Free excess arrival at the same entry | All 10,287 captured positives retain their port and donor word |
| Filled sham | All 4,194,304 inputs retain the filled target and its donor |
| Negative vacancy at the other matter port | All 10,064 captured negatives retain their port and donor word |
| Depleted donor | All 10,038 captured negatives retain their port and donor word |

On a separate bank of $2^{34}$ simulated starts, open a five-port native path after conversion from port 2 through 6, 3 and 48 to 49. The third owner reads all nine ordered words on 48 and 49. Both neutral words, (2,0) and (0,2), count. At time 16 this owner registers **89 neutral arrivals** from the charged source, zero for zero exposure of the converter, zero from the positive source and zero from the depleted donor. The independently prepared neutral state gives 35,314 neutral reads across the same number of starts without an incoming flight requirement. The calculated probability given conversion is $2.08077\times10^{-6}$, or approximately $5.01\times10^{-9}$ per original charged launch. The receiver observes redistribution of material supplied at the source and donor. The vacancy retains its location through the conversion service. [Complete banks](results.json) · [Independent 243-state receiver calculation](test_measure.py)

## Why the charged and neutral outputs are connected

The stored symmetric extreme-contrast collar on a seam with endpoint loads $u,v$ has the energy $32(1-\mathbf1_{u=v=1})$. Its preparation from ordinary material and its action on all nine load words follow from the full 99-state echo. Prepare these lines on two adjacent seams. The incoming occupation $x$ is on one end, the target $a$ is their common end, and the guard has occupation two. Apply the native echo between $x$ and the guard. If $x=0$ their output is (1,1). If $x=1$ it is a superposition of the two orderings of (1,2). If $x=2$ it is (2,2). Both coherent orders have exactly $2-x$ unit loads. Expose both stored contacts for half the phase in their energy units, undo the echo, and apply the opposite half exposure to the input–target contact. The direct correction subtracts $\mathbf1_{x=1}$. Since $2-x-\mathbf1_{x=1}=2\mathbf1_{x=0}$, the restored actuator acts, up to a common phase, as

$$
C(\phi)\,|x,a\rangle=\exp\big(i\phi\,\mathbf1_{x=0}\mathbf1_{a=1}\big)|x,a\rangle.
$$

The start and end contacts form a balanced interferometer on the two donor orientations. Its transition probability for a vacancy is $\sin^2(\phi/2)$, and its effect on the positive control is the identity up to a common phase. Two full operations return the original negative word with probability one. The oriented donor-resolved branch raises regional charge by one. These facts recover the finite raising, lowering and singlet-null rules of the specified code. On halving the native contrast at fixed exposures the full negative conversion probability becomes 1/2 and the filled sham transfers the donor. No physical weak coupling, helicity selection, or particle rest energy enters the contact schedule. [Physical pulse calculation](measure.py) · [Count and charge theorem](../../../../Lean/SourceLaw/ArrivalChargeLadder.lean) · [Physical storage of a collar bit](../../interval_echo_consensus/signed_effects/COLLAR_STORAGE.md) · [Exterior current dictionary](../../../../Lean/Screen/ChargedCurrentDictionary.lean)

## Preserving the internal state across conversion

Take all twelve rails 72–83 in owner 6. For each rail, use blank native swaps to bring its load to working port 3, apply the conditional phase and return every material port to its starting location. All paths avoid the guard, target, donor and compensators. Apply the donor beam splitter once on either side of the twelve visits. Exactly one address is vacant in every basis state of the one-vacancy code. Every such state gives the same phase to the same donor, with all guard records restored. Thus the complete conversion amplitude on the twelve rails is $iI_{12}$ at the full setting. On any embedded two-dimensional input $V$, the success effect for an aperture $A$ is $V^\dagger P_A V$; for a fully covered isometry it is the identity. This preserves a reference entangled with the vacancy's position. An interferometer on neighboring ports 72 and 73 recovers all 4,194,304 inputs, compared with 2,096,374 when the address is read in between.

The neutral output occupies the twelve rails **and the common target port in the adjacent owner**. Its double is supplied by the donor. Moving the marked vacancy outside the aperture removes its conversion. Any two vacancy spin or orbital states wholly inside the aperture receive the same conversion strength. This gives a specific control for interpreting an apparatus-selected doublet as a chiral weak response. The function <code>aperture_map</code> returns the complex forward, reverse, survivor and singlet blocks in the caller's port order, together with the actual restored relay paths. [Circuit and transfer maps](measure.py) · [Two-arm control](results.json)

## Conversion with fixed contacts

The static bench uses target 2 and donor 11, both of degree five, and a filled comparator at 10. Write symmetric extreme collars on the incoming–target and comparator–donor seams. After encoding the arriving occupation against the double guard, keep just those two contacts and the ordinary target–donor contact open. Keep the physical private references throughout. The captured negative gives encoded incoming load one, so the two stored energies add to a constant on the two target–donor orientations (1,2) and (2,1). Their private terms also add to a constant. The off-diagonal native exchange is $-1/2$. A constant free exposure $t$ converts the negative and neutral in either direction with probability $\sin^2(t/2)$.

For an incoming positive, the encoded incoming load is two. The energy difference of the two donor orientations is 32. Its donor exchange probability is therefore

$$
P_{+}(t)=\frac{1}{1025}\sin^2\left(\frac{\sqrt{1025}\,t}{2}\right).
$$

This off-resonant transfer has positive regional charge at both ends. It is a sensitive check of the complete singlet record. At $t=\pi$ the simulation gives full negative conversion and full inverse conversion in $2^{30}$ sharp preparations each. The positive input gives 592 donor exchanges in $2^{30}$, against probability $5.87296\times10^{-7}$. Flipping only the comparator memory parity shifts its energy to 33 and reduces the negative conversion probability to 0.316564. A filled input at the same exposure has 49.96% unreturned coding memory or guard, with the rejection bin retained. An independent incoming-charge herald therefore matters. This constant-contact result uses the encoded source, parked contacts, a populated charged reservoir and a fixed-duration read. For an arriving vacancy this experiment measures resonant material conversion at a fixed appointment. An autonomous decay width and a source-produced electroweak field require different records. [Native Hamiltonian](static.py) · [Complete calibrated inputs](static.json)

From the repository root run

```sh
OPENBLAS_NUM_THREADS=1 python -m model.experiments.interval_echo_gauge.arrival_ladder.measure
OPENBLAS_NUM_THREADS=1 python -m model.experiments.interval_echo_gauge.arrival_ladder.static
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/experiments/interval_echo_gauge/arrival_ladder/test_measure.py
```
