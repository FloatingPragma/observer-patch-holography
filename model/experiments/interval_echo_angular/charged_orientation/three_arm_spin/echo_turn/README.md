# Turning the received three-arm charge

An observer can coherently erase the three direction heads of an arriving charge, turn its material record on the receiving icosahedron, and write the record back. Three successive 120° turns give the ordinary positive return on the measured qubit. The same interferometer resolves a deliberately applied minus phase. The result holds both for free arrival at time 32 and for arrival through the energy filter centered on 16.43120649. [Source, energy selection and detector](../README.md) · [Pulse program](pulse.py) · [Simulated records](results.json)

## What actually turns

The charge arrives at physical port 160. The receiver checks the vacancy there, leaves its three direction heads unread and stops the flight for addressed operations. The emitting owner's three heads are discarded in the final read. The lab supplies three vacancy–excess probe pairs, six reversible staging slots, one ordinary private reference at port 137 and access through the equal-load shutter on 152–145. The original remote compensator stays at 36. The three probe pairs belong to the instrument. Running the experiment requires their sharp stocks to be held until the appointment. [Federation geometry and stock routing](material.py)

On each marked edge call the preparation address $a$ and the opposite end $g$. A full native magnitude-one pulse acts on a probe vacancy at $a$. For incoming direction zero the vacancy moves to $g$ with amplitude $+1$. For incoming direction one it stays at $a$ with amplitude $-1$. In both cases the head ends pointing toward $g$. The unknown bit now resides in the vacancy position. This is a reversible transfer of an arbitrary superposition. All three head bits, including their entanglement with an untouched reference, are transferred this way. No bit value is read. [Registered 99-state contact columns](pulse.py)

While the three probes are parked, the original vacancy visits the known written heads. A partial native contact and a routed private phase align the two written terms with the range of the echo write. The actual addressed echo erases the head. A reversed blank splitter returns the original vacancy to its entry port. Each step has an inverse attainable with positive dwell time because the isolated native contact spectrum is integral in these units. Both physical endpoint orderings give unit reset probability. Substituting the identity for the echo leaves zero amplitude in the blank sector. The 18 equal-load internal shutters also open at their stock loads. At the three gateway-preserving orientations the receiving owner can clear and restore its seven private shutters. [Direct matrix-exponential checks](test_turn.py)

The probes return to six real receiver ports, one vacancy on each marked edge's pair of endpoints. A proper icosahedral port permutation now physically transports the original vacancy and these six load coordinates. Each elementary step uses a full native exchange on an adjacent blank bond. After the turn, the receiver runs the inverse preparation on the rotated edges. The three probe pairs, all intermediate ports and the phase reference return to their supplied states. For the 120° step the spatial permutation alone takes 16 adjacent exchanges. Rotating the same records without the echo would leave the written magnitudes on their old edges. [Flight-only conservation result](../ROTATION_DOMAIN.md)

## Which phase is measured

Write $E_D$ for transfer to material vacancy positions in orientation $D$, and $P_R$ for the compiled permutation through the spatial rotation $R$. The operation between two mounted head registers is $E_{RD}^{\dagger} P_R E_D$. Each marked bit contributes its measured pulse phase. The eight resulting ternary load words are distinct. The physical permutation takes each to the corresponding load word of the next mount. Thus, on the occupied subspace, successive intermediate encodings cancel. A closed full-turn sequence gives $E_D^{\dagger} P_1 E_D=I$ including its phase relative to the idle arm. This identity extends to any reference entangled with the eight-word register. It also explains why changing an encoding convention cannot generate a spinor sign. Numerical replay of all 60 proper rotations and their powers has maximum operator error $3.2\times10^{-13}$. [Turn composition and readout](measure.py)

Two nonidentity 120° rotations preserve the three physical gateway addresses as a set. The measured cycle uses one of these three times. In the 15 unordered pairs of perpendicular proper half-turns, the compared orders also give $+I$, to $3.9\times10^{-13}$ in operator norm. Each of these half-turns changes which internal vertices meet the external gateways. Those measurements restore the internal shutter pattern in each new frame, hold the external connections isolated and restore the private pattern on returning to the original frame. The distant arms remain in place throughout. These operations use addressed site exchanges representing the discrete rotations of the receiving cell. They give no continuous, energy-preserving rotation of the freely propagating pole.

A phase-coherent comparison of two pulse histories is supplied. Its bright and dark amplitudes are half the sum and half the difference of the histories. In the independent phase control the three $g$ rails visit the ordinary private reference, each for $\pi/8$, while the unknown word is in material vacancy positions. On the accepted one-hot subspace exactly one of those rails has unit load. The registered private cost $8n^2$ supplies one minus phase. An ideal spin-half comparator also gives the dark result for the two half-turn orders through $(-i\sigma_x)(-i\sigma_y)=-(-i\sigma_y)(-i\sigma_x)$. The rail phase is a compiled native intervention. The spin-half comparator is a mathematical control.

| Receiver input | Full turn, bright / dark | Full turn and applied rail phase, bright / dark | Accepted probability per original start | Conditional minimum eigenvalue after partial transpose |
| --- | ---: | ---: | ---: | ---: |
| Free release at time 32 | 42,514 / 0 | 0 / 42,462 | $1.55221\times10^{-7}$ | −0.464109 |
| Finite energy filter | 1,372 / 0 | 0 / 1,345 | $4.86937\times10^{-9}$ | −0.495552 |

Each simulated bank has $2^{38}$ original starts, including both local one-third apertures, missing charge and any clock refusal. On the ordinary turn the simultaneous 99% upper bound on the accepted dark fraction is 0.000141 for free arrival and 0.00436 for the filtered arrival, paired within each preparation with its material phase control. The second source uses the inherited coherent global delay service at 16.43120649. Reconstructing the complete remote amplitudes before the direction read gives the same transferred reference state after the full turn. This demonstrates coherent handling and ordinary orbital return for this apparatus. An intrinsic spin assignment would require a physical rotation acting on the intact state that determines the reference phase independently of this unwrapping and remounting circuit. [Full receipt](results.json) · [Independent original forward and reverse measurements](../RECIPROCAL.md)

From the repository root run the following commands.

```sh
p=model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.echo_turn
d=model/experiments/interval_echo_angular/charged_orientation/three_arm_spin/echo_turn
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m "$p.measure" --output "$d/results.json"
python -m pytest -q "$d"
```
