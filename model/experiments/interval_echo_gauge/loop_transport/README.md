# A vacancy distinguishes two closed paths

Ordinary blank-memory Interval Echo transports a material mark noncommutatively. Take a unit-filled cell with one vacancy and one excess. Two triangular vacancy paths with a common root and spoke permute the excess in different orders. On the 20-owner board the two finite-exposure sequences produced different terminal excess counts, with every failed vacancy checkpoint retained. A separate disjoint-cell experiment on the same sector ground resolved no neighboring response. [Path protocol](DESIGN.md) · [Original records](results.json) · [Placement protocol](PLACEMENT.md) · [Held receiver records](placement.json)

## What travels around each loop?

Write $|x,y\rangle$ for the excess at $x$ and vacancy at $y$. The exact generator block from vacancy $y$ to a neighbor $z$ is $-P_{yz}/2$. The map $P_{yz}$ sends an excess at $z$ to $y$ and fixes the other excess positions. This follows from the native seam exchange. All other contacts and the diagonal contrast cost have zero matrix element between these two vacancy fibers. Thus the first derivative of one observed crossing is $iP_{yz}/2$ in the free exposure parameter. The code checks the complete off-diagonal block against the native two-defect Hamiltonian.

Take neighboring faces $(r,a,b)$ and $(r,b,c)$. Returning the vacancy to $r$ along the first face exchanges marks $a,b$. The second face exchanges $b,c$. Their group commutator is a three-cycle. On the two-dimensional plane of excess amplitudes with sum zero, its normalized trace is $-1/2$. On the constant vector it is one. With only a vacancy and no excess the transport fiber is scalar and the commutator is one. The two triangle operations are reflections in the marked plane and generate the finite permutation group $S_3$. There are no independent link variables in this preparation. The same infinitesimal permutations occur with the contrast cost set to zero.

At finite exposure the observers monitor the vacancy location after each of six ordinary holds and read the final excess. A projection records one specified vacancy checkpoint. Other contacts act during every hold. The resulting success operators include interfering paths and loss at the checkpoints. The exact group commutator above refers to the leading transport blocks. Finite-wait counts are a separate order test. On the 240-port board the initial excess was at the shared spoke (port 5).

| Order at exposure 1/2 per checkpoint | Success probability | Excess reads at ports 1, 5, 11 and elsewhere |
|---|---:|---|
| First triangle then second | $3.0112\times10^{-10}$ | 1, 19, 1,379, 27 |
| Second triangle then first | $3.1665\times10^{-9}$ | 32, 3,560, 10,070, 303 |

Each row uses a multinomial draw for $2^{42}$ simulated starts from the full history probabilities. The records contain all six first-failure counts. On the first route at exposure 1/4 there were no successful 240-port trials. Isolated-cell records with and without the contrast cost use the same program. Interchanging the path order also changes the finite-time dynamics at zero contrast cost.

## Vacuum and disjoint receiver

The triangular exchange has a negative-eigenvalue read $A_{r,ab}$. With the vacancy at $r$, it projects the excess onto $(|a\rangle-|b\rangle)/\sqrt2$. At any other vacancy position it gives zero. The unperturbed state for this measurement is the numerical lowest state of the complete ordered excess-vacancy sector, of dimensions 132 and 57,360 on the cell and first federation. The 99-state native seam generator preserves its blank-memory subspace during these free holds. The sparse-sector ground uses the independently moving charges only as a solver preconditioner. The full interacting Hamiltonian supplies the residual, the Perron energy interval, the flights and the correlation. The connected Euclidean readout is calculated from the sector spectral measure. The real-time detector uses fresh simulated binary clicks.

Rotation symmetry makes the unique isolated-cell ground amplitudes equal on $a,b$ at a fixed adjacent root. Its negative-loop variance vanishes exactly. On the federation the first enrolled triangle gave a numerical variance of $5.98\times10^{-13}$. A phase flip on this home projector followed by ordinary flight produced zero clicks at the disjoint receiver in $2^{24}$ copies per setting, including the sham. The follow-up enumerated 1,200 rooted triangles and compared 13,500 ordered eligible pairs across a gateway. It selected the pair with the largest smaller ground variance, $6.12\times10^{-11}$. Across the six held times there was one pulse-arm click and no sham clicks, with $2^{30}$ trials per arm per time. No comparison crossed the prespecified familywise threshold. At exposure 1 the selected home correlation was $1.98\times10^{-14}$, down from $6.45\times10^{-11}$ at zero. This finite correlation does not identify a mass.

For any ordered-sector state, $0\le A_{r,ab}\le P_r$, where $P_r$ is the vacancy-at-$r$ projector. Each binary variance is at most $\langle P_r\rangle$ and $\sum_rP_r=I$. For one choice of triangle per root the mean original-denominator variance is at most $1/N$ on an $N$-port screen. In a sector vacuum its connected Euclidean correlation at any nonnegative wait is at most that variance. This rules out a translation-homogeneous nonzero local variance floor for this specific bounded click in the single-vacancy refinement, regardless of clock scaling. More generally, the sum of the root vacancy projectors is the vacancy number operator. For any binary read supported on a root vacancy the mean variance is at most the expected vacancy density. Any subextensive number of vacancies still loses this signal at typical roots. Localizing at a special root, increasing defect density, changing the observable normalization or finding an independent field register changes the hypothesis. [Finite variance proof](../../../../Lean/ObserverPatchHolography/LoopVacuumBounds.lean)

From the RER root, use fresh output paths.

```sh
.venv/bin/python -B -m model.experiments.interval_echo_gauge.loop_transport.measure --output /tmp/vacancy-paths.json
.venv/bin/python -B -m model.experiments.interval_echo_gauge.loop_transport.placement --output /tmp/vacancy-receivers.json
```
