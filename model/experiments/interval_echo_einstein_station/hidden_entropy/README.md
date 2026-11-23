# A stationary repair state with a measurable entropy slope

One load is privately placed in any descendant of the first base face. Each available cell holds two orthogonal orbital modes which vanish on all three gateways and both have energy 51. Vary the probability $p$ of preparing the first mode. Cell receivers undo four locally compiled waveforms and read their sharp returns. Pooling those records measures the occupation imbalance and both possible coherences. The full regional entropy is $\log M + h(p)$ nats on $M$ equiprobable source cells, with $h$ the binary entropy. The decoded two-state entropy is $h(p)$ and has the same first variation. A meter that discards the internal bit loses this variation. [Protocol](DESIGN.md) · [Source and instrument](measure.py) · [Complete simulated records](results.json)

| Board | Descendant cells | Observed slope about $p=1/2$ | Observed slope about $p=3/4$ | External clicks for either dark launch | External clicks for the bright launch |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 1 | +0.00038 | −1.10780 | 0 | 117,977 |
| 80 | 4 | +0.00189 | −1.10225 | 0 | 83,236 |
| 320 | 16 | −0.00059 | −1.11284 | 0 | 61,744 |

Each column uses fresh ensembles of 262,144 copies at the final appointment. Slopes are symmetric finite differences with half-step 1/16, in nats per unit mixing probability. The infinitesimal slopes are zero and $-\log 3$. Increasing the internal kinetic weight at gateway endpoints gives 127–591 exterior clicks per dark setting on the same boards.

Each dark vector is an eigenvector of the isolated cell operator. Its zero at every gateway makes the exterior half-Laplacian vanish on it, so it is an eigenvector of the full free generator with the same energy. Every mixture above is therefore stationary throughout free flight. Across the parent cut the exterior and gateway ports stay in the same vacuum with blank seams. In any adaptive sequence of passive measurements that leaves this vacuum and the open contacts intact and injects no load, each outcome acts as a scalar on the vacuum factor. Induction on the records makes their entire probability law independent of $p$. Consequently an area or volume functional of those records has zero derivative along this family. At $p=3/4$ the generalized entropy slope for every such area and every finite area coefficient equals $-\log 3$.

At $p=1/2$ the same family has a stationary entropy. Free native evolution preserves every value of $p$, so it cannot drive these preparations specifically to the equal mixture. This result concerns the selected one-load algebra and passive exterior readers. A physical vacuum selector, an active probe coupled to interior matter, and areas or volumes supported on additional records would change the experiment.

Replay with `python -m model.experiments.interval_echo_einstein_station.hidden_entropy.measure`.
