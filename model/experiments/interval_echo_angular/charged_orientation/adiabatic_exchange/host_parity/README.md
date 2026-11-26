# Two mobile cores acquire a material-dependent exchange sign

Two identical unit loads cross a Y junction by different routes and reach the same two exit ports with exchanged assignments. On the parity-odd material, one exchange gives a relative phase of −1 and two exchanges give +1. On the parity-even material one exchange gives +1. A native contact prepares either material from the same sharp loads. The first simulated bank of 16,777,216 launches gave 8,388,532 odd-source dark exits and 8,388,684 even-source bright exits for one exchange. Every launch entered the count. [Frozen protocol](PROTOCOL.md) · [Records](results.json)

The held charged-core transfer uses two **actual vacancies** on a passive four-port ring. At native time 64, odd-source launches give zero hits on the opposite vacancy pattern and the even-source hits exceed 7.30 million per bank of 16,777,216 unfiltered starts. A single vacancy winding through any three-host input has dark probability at most 3/4. A host-only permutation accounts for the full two-vacancy phase. [Charged-vacancy transfer](VACANCIES.md)

| Observer comparison | Even material, dark probability | Odd material, dark probability |
| --- | ---: | ---: |
| Two three-hop routes to the same two destinations | 0 | 1 |
| Closed exchange versus retracing the route | 0 | 1 |
| Two closed exchanges versus two retracings | 0 | 0 |
| Swap the two host loads while both cores stay still | 0 | 1 |
| Exchange versus the matched host swap | 0 | 0 |
| Change only the host preparation, repeat the two routes | 1 | 0 |

A balanced timing arm probes one blank contact at its swap time. The two energy lines of an initial sharp empty–double pair occur at 32 and 33 native energy units. Both outcomes have probability one half. The path comparison begins after this source read. Fresh copies read the locations of the two unit loads at each appointment. Both physical paths end with units on the centre and the spare leaf. An independent coherence cohort bypasses the intermediate reads. Without the energy read the two route outputs split equally. The ancillary timing arm, comparison arm, coherent shutter control during the addressed exchanges and terminal analyzer are laboratory services. [Source, routes and complete position laws](experiment.py)

A second detector opens all three contacts together and reads the two unit positions after a wait. The cores then propagate without direction commands. At the native contrast and time 128 the even source reached the opposite two-port configuration in **9,355 and 9,410 of 16,777,216** trials. The odd source reached it **zero times** in both banks. Reflection forbids that event at every wait. At half contrast and time 64 the even counts rise to 63,728 and 63,726. The deposited counts are independent pseudorandom samples of the finite model's probabilities. The exact free-flight selection rule holds on passively prepared 20-, 80- and 320-owner federations with every seam and private reference open. Thirteen zero-contrast contacts isolate each four-port junction. The external boards have more owners but the enrolled material has the same four ports. [Held setup](PASSIVE_PROTOCOL.md) · [Stock preparation and release](passive.py)

The exchange sign belongs to the four-port material sector. An exchange of its two physical hosts gives the same phase while the cores remain parked. A phase applied to one host before launch reverses the result with identical core trajectories. Enabling another spare leaf allows two routes for exactly the same two-core exchange. With one extra sharp empty host, the two resulting dark probabilities are **1 and 1/4**. Their mean is at most **3/4 for any joint three-host preparation**. Thus this experiment establishes a native, locally reproducible negative exchange phase in a specified medium. Intrinsic spin-½ statistics of a portable electron core and a physical 2π rotation remain open. [Proof and bound](PROOF.md)

From the RER root, set the four math libraries to a single thread and run the following commands.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.experiments.interval_echo_angular.charged_orientation.adiabatic_exchange.host_parity
python3 -B -m "$p.measure"
python3 -B -m pytest -q model/experiments/interval_echo_angular/charged_orientation/adiabatic_exchange/host_parity
```

The free records order the six unit-location pairs cyclically, beginning at the source. The joint source and interferometer records order their four exits as even-bright, even-dark, odd-bright and odd-dark. The held records are repeat calculations on physical mounts. The two shot banks are independent simulations of the registered detector law.
