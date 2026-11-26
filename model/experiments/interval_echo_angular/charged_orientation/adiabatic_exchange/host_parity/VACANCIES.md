# Two charged vacancies exchange on a passive square

Two empty ports in a ternary material each carry a local vacancy of signed load −1 relative to a unit port. Native swaps transport the vacancies from opposite corners of a four-port ring onto the other two opposite corners. Clockwise and counterclockwise routes exchange the destination assignments. Their relative phase is −1 if the two occupied ports were prepared on their odd native contact line, and +1 on the even line. A closed exchange squared is bright on either line. [Frozen transfer](VACANCY_PROTOCOL.md) · [Complete results](vacancies.json)

| Test | Even host, dark probability | Odd host, dark probability |
| --- | ---: | ---: |
| Two vacancies on opposite arcs | 0 | 1 |
| Two repetitions of the closed exchange | 0 | 0 |
| Adjacent hosts swapped with both vacancies stationary | 0 | 1 |
| Compare vacancy exchange with the matched host permutation | 0 | 0 |
| One vacancy winds through the same parity pair and one extra unit | 1/4 | 3/4 |
| Three windings of that single vacancy | 0 | 0 |

The source is a sharp 0102 on the cyclic ports. A native swap brings the two hosts together. A balanced timing measurement of their contact distinguishes native energies 8 and 9 at the registered contrast. Undoing the routing pulse places the two hosts back on opposite corners. The two source readings occur with equal probability. In flight quiet-written bonds around the four sites isolate the ring. Each unwanted chord or outward bond is written on equal sharp loads before release. All external loads remain at one. All four ring sites carry the same native private degree term. The mounting procedure found ports 0,5,4,11 in the deposited 20-, 80- and 320-owner federations. Eleven physical collars keep each mount closed under the unrestricted generator. The coherent arm and addressed scheduling of the two exchange routes are separate from the unrestricted flight.

The native contrast on the square satisfies $\sum_{ab}(n_a-n_b)^2=2\sum_a n_a^2-2\sum_{ab}n_an_b$. Only two sites are occupied. Their product term is either zero or the same product $1\cdot2$ when adjacent. Thus the generator preserves the parity under interchanging the two host levels despite their unequal individual contrasts against a vacancy. A spatial reflection in the axis through the initial hosts fixes the source even when its **host** parity is odd. When the vacancies arrive at the complementary diagonal this reflection exchanges the two host slots. Odd material forbids that destination at every free-flight time. The even material has nonzero transfer at every registered wait. In each simulated bank every sharp launch counts, including both source outputs.

| Native free hold | Even-source opposite-diagonal hits, bank 1 | Bank 2 | Odd-source hits, both banks |
| ---: | ---: | ---: | ---: |
| 1 | 7,184 | 7,246 | 0 |
| 8 | 478,545 | 478,109 | 0 |
| 32 | 5,685,113 | 5,681,974 | 0 |
| 64 | 7,301,325 | 7,302,344 | 0 |
| 128 | 3,731,868 | 3,729,270 | 0 |

Each column uses 16,777,216 total starts per bank and wait. At time 64 the model conditional destination probability is 0.87034024 on even material. The source-conditioned distinction survives the separate half-contrast run. Probing all terminal spoke differences on three consecutive blank links and using the known load sum identifies the two vacancy positions. The deposited free counts use the order even-target, even-other, odd-target, odd-other. [Native implementation](vacancies.py)

At the chosen timing a single vacancy winding once moves the three other hosts by a three-cycle $T$. Its eigenphases are $0,2\pi/3,-2\pi/3$ since $T^3=I$. Consequently the dark read $(1-\mathrm{Re}\langle T\rangle)/2$ is at most $3/4$ on *any* three-host input. A sharp ordinary host gives zero. Filling the second vacancy with one extra unit and repeating the original host-pair preparation gives $1/4$ and $3/4$ in the two branches. Three windings give zero in both. A bare one-body π winding cannot explain the negative two-vacancy phase on this apparatus. A *two-host permutation* can explain it exactly. This identifies a relative charge-vacancy exchange phase in a specified passive material and leaves portable electron statistics open. [Held winding comparison](ONE_VACANCY.md)

Reproduce with the one-thread environment in the [junction replay](README.md) and run `python3 -B -m model.experiments.interval_echo_angular.charged_orientation.adiabatic_exchange.host_parity.vacancies`.
