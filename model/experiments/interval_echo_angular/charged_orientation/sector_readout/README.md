# When do two charge readings distinguish different particles?

The two-loop material can carry a quantum signal to the far owner in 32 of 64 preparations at the tested time. Resolving one positive formation energy first prevents an unknown spin qubit from reaching that owner in **every** one of the 64. The search includes all six written links in both contrast parities, a blank bridge, both receiver winding records and the entire receiving owner. A known spin-half control also loses its qubit under energy resolution when an external field splits its spin levels. The same control exhibits the spinor rotation sign. Losing coherence under an energy measurement therefore diagnoses the tested encoding and background. Establishing particle spin requires a separate physical rotation measurement. [Protocol](PROTOCOL.md) · [Exact certificate](certificate.py) · [All outcomes](results.json)

## What was searched?

A source owner and a receiver owner each contain a charged three-port loop. The two gateways share a blank bridge. The matched host and all six signed material addresses come from the [native apparatus](../two_winding_loops/README.md). On a unit-contrast transport seam the odd-parity formation term is exactly 9/8 times the even-parity term. The code checks this relation against the original local matrices and examines all 64 assignments on the 16 reachable charge and head records for each. The parity settings are preparation choices. Electron, muon and tau masses play no role in selecting the space or the instruments.

At the fixed release time 16, the calculated two-state channel at the single remote gateway has a negative partially transposed Choi eigenvalue below minus $10^{-10}$ in 32 settings. One of these channels has already been confirmed by finite-count local tomography and an independent rational witness. Across all 64 spectra there are **464 simple and 24 double positive energy lines**, counting a line again on each material in which it occurs. Twelve preparations have positive doublets. These counts describe modes of the prepared material.

For the sharp-energy test the sending material charge may occupy any of the three source ports, with *any* coherent input on all eight compatible head and position records, including the two advance settings of the distant winding memory. The observer at the receiving end collects all three receiving ports and their local heads. Both possible head histories left at the source are traced. The detector also retains the failed arrival. The following channel identity applies even to an arbitrary new choice of source encoder and local qubit decoder in this domain.

## Exact obstruction on a sharp line

Let $H_8$ be the integer matrix $8H$, with square-free minimal annihilating polynomial $m$ obtained from the characteristic polynomial. For an eigenvalue $x$ of $H_8$, the polynomial divided difference
$$
A(x)=\left.\frac{m(y)-m(x)}{y-x}\right|_{y=H_8}=m'(x)P_x
$$
supplies the full spectral projector $P_x$ up to a nonzero scalar. On a simple line the projection prepares a fixed state. For a double line, restrict the columns of $A(x)$ to the eight input records and its rows to the four local receiving records at each source-history value $b$. These two real matrices $K_b$ give the accepted map up to a common scalar.

For every pair of output rows $r,s$ and input columns $i,j$, the program checks that the polynomial
$$
\sum_b [K_b(r,i)K_b(s,j)-K_b(r,j)K_b(s,i)]
$$
vanishes modulo the irreducible factor containing $x$. The 2,016 integer polynomial identities cover every positive repeated root. They say that the Choi matrix of each accepted sharp-line map equals its partial transpose. In particular it is positive under partial transpose. Local source encodings and receiver decodings preserve this property, including upon heralded success. Any decoded two-qubit Choi state is then separable, so this mapping cannot transfer an unknown spin qubit between owners. The result leaves broadband coherent transport, different head magnitudes and additional quantum references shared with the receiver outside the tested source open. [Integer matrix and polynomial checks](certificate.py)

## What can the charge meters themselves tell apart?

A stable label unaffected by free evolution and every charge-position measurement must have a projector commuting with $H$ and every measured position projector. Solving the linear equations exactly gives only multiples of the identity on each of the 64 reachable components. Using only the coarser two-owner meters gives a four-dimensional commutant for the all-even bench and a two-dimensional commutant for the bench with one odd source edge. The finer observer can distinguish and couple orbital records that the coarser observer groups together. This criterion tests exact conservation under these specified scalar instruments. A physical species that interacts or mixes under an added probe needs its own dynamical description.

On the one-odd-seam source a single material position measurement with its outcome forgotten transfers population between positive lines. For initially equiprobable positive energy eigenstates, every partition of the eight lines into three nonempty groups has at least 0.0029716 unconditional probability of moving from one positive group to a different one. The 966 possible partitions are enumerated. A further 0.1736714 goes to states whose formation energy is negative relative to the same filled host. These are spectral transition calculations; they do not assume that an observer has built an ideal global energy meter.

A separate direct experiment uses the independently prepared bare and once-wound vacancy words from the same gateway. After 16 native units, record nothing, the owner only, or the exact charge address. After four more free units count the final address. For the bare input the receiving gateway reports **199,764**, **279,682** and **250,653**, respectively, per 1,048,576 original starts. For the wound input it reports **338,608**, **377,223** and **290,476**. All intermediate and terminal outcomes are retained. The first detector is a physical intervention with different back action at its two resolutions. Both detectors record the same conserved vacancy across the bridge. [Free evolution and counts](measure.py)

## A control with known spin and flavor

In the comparison Hamiltonian, a particle travels on a three-site orbital path with a known spin-half factor. Spin-independent charge projectors yield commutant dimension four at zero applied field and two with a spin splitting. Forgetting a perfect energy measurement gives the same far-site acceptance $2/9$ in both cases. The accepted Choi partial transpose has minimum $-1/9$ without the splitting and zero with it. Independently calibrated perpendicular spin half-turns give the negative order sign in both cases. Three distinct imposed flavor copies give commutant dimension twelve and three explicit invariant flavor projectors. The controls illustrate how orbital lines, a split spin and distinct imposed flavors appear in the same analysis.

From the repository root run with one numerical thread.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.experiments.interval_echo_angular.charged_orientation.sector_readout.measure
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/sector_readout
```
