# Charge and colour of a bound triangular face

Three ordinary Interval Echo loads are detected together on a triangular face with 98.11% probability after one free exchange time on a cell. An occupied tip gives one positive material unit against an empty reference. A collective colour-centre turn of three calibrated one-load receivers gives a dark interference port on 42.26% and 41.98% of all launches from two different faces. A supplied colour singlet stays bright. The distinction concerns this reversible extraction and the registered one-load colour action. [Design](DESIGN.md) · [All detector probabilities](results.json)

## What an observer can read

Start from three sharp occupied ports on a face of an otherwise empty cell. All thirty contacts act during flight. On fresh preparations three seam flags identify a face. A signed echo then compares an occupied port with a reference of load zero, one or two. For each occupied tip it reads +1, 0 or −1 respectively. A separate contact interferometer also distinguishes occupied and empty ports. Against an empty probe, the imaginary-quadrature bright probabilities are 0.129971 for an occupied port and 0.5 for an empty port. At the first launch face one initially occupied port has unconditional occupation 0.991464 after flight. Each of the three detected loads has the same calibrated response. The experiment measures signed material contrast. Its correspondence to electromagnetic charge remains open.

For the collective test, empty neighbours receive the home ports by reversible swaps. Each address has a fixed destination. A count at each receiver checks whether it holds exactly one load while preserving coherence between addresses in that cell. On that outcome one arm applies the same guarded colour word in each receiver. Both arms undo the routing and interfere. The table counts routing and source failures in the original denominator. The face and colour experiments use fresh copies.

| κ = 8 preparation | Face probability given source | One load in each receiver per launch | Centre-dark per launch | Centre-dark given three receivers |
| --- | ---: | ---: | ---: | ---: |
| Sharp face (0,1,5) | 100% | 100% | 42.19% | 42.19% |
| Same face after one free time | 98.11% | 98.52% | 42.26% | 42.89% |
| Ground filter from that face | 98.59% | 2.45% | 0.84% | 34.22% |
| Sharp face (0,1,7) | 100% | 100% | 42.19% | 42.19% |
| Second face after one free time | 98.11% | 98.48% | 41.98% | 42.63% |
| Ground filter from the second face | 98.59% | 1.97% | 0.68% | 34.61% |

The ground filter accepts 4.93% at either launch face. With exchange alone the ground face weight is 9.09% and the face weight after sharp release is 20.03% at time one. The sham and the turn followed by its inverse give dark probability below the receipt's precision. Three identical triplet orbitals also stay almost bright under the centre but go dark with probability above 0.99999 under a half-angle colour turn. The antisymmetric three-triplet control stays bright under both. This second control matters because the centre acts trivially on any tensor of three triplets. [Registered one-load instrument](../../../../../experiments/interval_echo_gauge/singlet_flight/README.md)

## Can another routing reveal a singlet?

There are 5,775 ways to distribute twelve port addresses among three unlabelled receivers of four addresses each. Between two and fourteen of the twenty triangular faces put one load in each receiver. On **every** partition and **every coherent combination of its accepted face addresses**, the ideal collective centre has dark probability at least **63/256**, or 24.61%. The calibrated finite pulse changes a three-receiver probability by at most 0.000118 using its measured operator error. This bound applies to the face-supported accepted input. The table above gives the actual full-state results, including configurations off the faces. [Enumeration and bounds](routings.json) · [Follow-up design](ROUTINGS.md)

Let $P$ be the one-load triplet projector. The ideal word multiplies a triplet by a cube root of unity and acts as the identity on its complement. The collective dark probability is $3/4$ on configurations with one or two triplets. It vanishes with zero or three. The corresponding effect is $3(\sum_i P_i-\sum_{i<j}P_iP_j)/4$. On sharp accepted face addresses its diagonal is $27/64$. For addresses differing in one place, the off-diagonal is $3P_{ab}/8$. For two places it is $-3P_{ab}P_{cd}/4$. All other off-diagonals vanish. A distinct-port entry of $P$ has magnitude $1/4$ for antipodes and $1/(4\sqrt{5})<1/8$ otherwise. Subtracting a rational upper bound on the sum of absolute off-diagonal entries in each row from its diagonal gives the stated lower bound. Independent proper rotations of receiver frames preserve the centre projector.

The inverse experiment supplies a three-triplet antisymmetric singlet to the readers and asks where its loads go when the swaps are undone. In the common port frame, each accepted face has squared overlap $\det(P_{abc})/6=(1+1/\sqrt{5})/960$, where $P_{abc}$ is the submatrix on its three ports. Even the best routing returns the three loads to a home-cell face with only **2.11% probability**. With arbitrary relative triplet frames or different sharp destination addresses, each projected sharp port has norm $1/2$. The squared determinant divided by six is then at most $1/384$. The fourteen-face maximum gives **7/192**, or 3.65%. These probabilities include material left in other cells in the denominator. A constituent spread coherently over several ports has a different overlap with the triplet subspace.

## Does the full ground state contain the three-triplet singlet?

With all three readers aligned to the common port frame, the singlet made from three colour triplets has **exactly zero overlap with the one-cell three-load ground state** under every balanced routing. This includes the amplitudes on all 220 configurations. The sharp releases break the ground state's symmetry. At time one either launch has up to 0.199% three-triplet singlet weight per preparation, depending on the routing. Keeping the ground state and rotating just the first reader gives a maximum of 0.0852% at κ = 8 for the tested rotation. The zero therefore depends on the common colour frame. [Full-state test](GROUND.md) · [All pair-fixing witnesses and control maxima](ground.json)

The hopping graph is connected and every hopping amplitude is negative, so the ground state is unique and can be chosen positive. Every incidence symmetry preserves the Hamiltonian and fixes this state. Write $v_a$ for port $a$ projected into the registered triplet and $\psi_{\{a,b,c\}}$ for the ground amplitude on three distinct sites. Define the alternating amplitude $T_{abc}=\psi_{\{a,b,c\}}\det(v_a,v_b,v_c)/\sqrt{6}$ and set it to zero on repeated sites. For each of the 66 pairs $a,b$, an incidence reflection fixes both ports and reverses triplet orientation. The enumeration checks that its square is the identity and its exact trace on the triplet is one, so its determinant there is $-1$. Summing over the third site pairs opposite contributions and gives $\sum_c T_{abc}=0$.

Now assign all sites to receiver sets $A,B,C$. Sum this identity over $a$ in $A$ and $b$ in $B$. Terms with $c$ in $A$ cancel on interchanging $a,c$, and terms with $c$ in $B$ cancel on interchanging $b,c$. The remainder is the singlet amplitude with one load in each receiver. The argument concerns the antisymmetric three-triplet singlet. Components in which all three loads are individually invariant can contribute to the total invariant weight recorded by the apparatus.

## Reproduce

Run from the RER root. The tests compare the binding and flight with an independent three-load generator and diagonalize the compressed centre effect for all 5,775 partitions.

~~~sh
OPENBLAS_NUM_THREADS=1 python -B -m pytest -q -W error model/candidate_laws/interval_echo/openwave/quark_winding/bound_face
~~~
