# Can an observer read an arriving quantum's energy without starting a current?

Yes, with a coherent two-path energy comparison. For each of 8,249,494 one-quantum arrival records, a supplied path switch let the wing excitation spend the waiting interval either in the native receiving wing or in idle storage. The switch was undone before the final reads. Its pointer resolved **31.80150, 32.51524 and 33.21471** radians per unit model exchange time against an independently prepared ordered ground. All 8,249,494 terminal one-quantum reads retained the initial flag. The phase belongs to the difference between the two histories. No local current read was made before waiting. [Frozen protocol](ROUTING_PROTOCOL.md) · [All pointer, flag and energy counts](routing_counts.npz) · [Count-only results](routing_results.json)

## How the energy origin is defined

Use the original four-dimer ring collision and the same broad left-wing flag as in the [free-end experiment](README.md). In each ordered ring phase, prepare a separate matching three-dimer memory in a specified energy eigenstate of its matching pin. A balanced control initially in its positive X state chooses between leaving the registers in place and swapping their entire load configurations. During the dwell, apply native wing evolution U to the active register and identity V to the storage register. Undo the controlled swap; read X or Y on independent copies of the control and record the terminal number flag. Both ring phases, both initial and final numbers, both signs and unordered misses are in the archive. The procedure needs neither a phase in the arriving state nor a sinusoid programmed at its excitation energy.

For the two switch values the post-routing unitaries on target and memory are A = U⊗V and B = V⊗U. With identity storage, a reference eigenstate of energy Eᵣ and target energy weights pⱼ conditional on the herald, the measured complex mean is

$$
\langle X\rangle+i\langle Y\rangle
=\sum_j p_j\exp[i(E_j-E_r)t].
$$

This follows by tracing the two path amplitudes against each other on the same registers; it holds for an initially mixed target and also after tracing out its unread partner. Both branches commute with the target number flag after unrouting. Postselection on a pointer sign can change *which* mode is likely, but cannot move an event out of its initial occupation sector. Adding any common scalar to the active Hamiltonian gives the same phase in both histories and cancels exactly. The implementation also checks a shift by +7 on each dense native wing. With identical evolution U in both registers, A=B and the pointer stays in its initial positive X state. [Implementation](routing.py) · [Direct real-space controlled-swap and terminal-number checks](test_routing.py)

| Source and reference | Measured switch record |
| --- | --- |
| Wall source, ground memory | +31.80150, +32.51524, +33.21471; predicted held complex RMS error 0.0156 |
| Same wall, middle-mode memory | −0.71494, −0.00019, +0.69891; held error 0.0147 |
| Ordered sham, ground memory | Mean X = 0.999997; largest distance of X+iY from 1 was 0.00644 |
| Same wall, identical evolution in both registers | X = 1 on every accepted X read |
| Same wall, dephased switch | Mean X = 0.000670 and mean Y = 0.000255 |

Each comparison used 268,435,456 starts over the two pointer settings. The sham table includes every ordered event. Other pointer rows condition on the actual single-quantum herald (about 3.07%). An independent energy measurement on 262,144 launches classified 7,941 single-quantum arrivals as 2,051, 4,010 and 1,880 events in ascending energy. These proportions, 0.2583, 0.5050 and 0.2367, can be compared with the ground-memory pointer weights 0.2708, 0.4983 and 0.2309 (individual deviations below the simultaneous 99% Hoeffding radius 0.0201 of the direct energy counts). The pointer used the first 256 appointments for blind signed-frequency discovery and reserved the next 256 for testing. The old absolute wing energy calibration was consulted only after the fit.

## What the apparatus supplies

In the installed collision geometry the ring is clamped and the two wing generators then act on separate wings. Each such product evolution, with operations local to those separated registers, lacks a switch-controlled coherent transfer between them. The experiment explicitly supplies a phase-coherent whole-wing router, an independently ground- or mode-prepared memory with idle storage, and a quadrature pointer. It does not compile a native three-dimer controlled swap or prove that another board with additional couplers could never implement it. The simulation audits the complete supplied switch by evaluating its real-space control, target and memory amplitudes directly. Statistical sampling uses the equivalent reduced instrument.

There is a finite obstruction on this particular clamped fixture. On both recorded ring phases the full native generator restricted to that phase is exactly the sum of the two independent wing generators and a scalar (a direct equality on all 64 rows and columns). Within each side allow arbitrary signed-echo readouts, transports using its open contacts, and shared classical outcomes that choose later operations. With the independent memory separated from the target, each complete record still has product Kraus operators across that cut; averaging or forgetting records gives a mixture. Feed the intended switch a balanced control, an excited target mode and a ground memory. Its two output branches place orthogonal states in the memory and have orthogonal states on the other side. Their Schmidt weights are equal. Every separable output of the installed operations has fidelity at most one half with this intended pre-dwell switch output, also after a successful classical herald. This bound follows by Cauchy–Schwarz on the two branch overlaps, then convexity. Reopening a coherent connecting route or supplying pre-shared quantum correlations would change the apparatus. [Full-generator equality test](test_routing.py) · [Checked fidelity inequality](FreeEnd.lean)

The result establishes relative energy spectroscopy on the collision daughter. It resolves formation costs of the finite material's three standing modes on the same herald. It gives neither an independent inertial mass nor a free particle's autonomous local rest clock; OW-03 stays partial. Passing the stored switch through the same flight removes every phase, demonstrating exactly what physical difference the pointer measures.

From the repository root, run the following in the scientific Python environment. Add `--reacquire` to repeat the collision and sampling using the installed numerical stack.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_03_clock.free_end.routing
```
