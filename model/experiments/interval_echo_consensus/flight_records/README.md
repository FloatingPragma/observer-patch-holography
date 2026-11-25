# What can an observer retain when flight and a public repair exchange order?

At fixed total load two, a complete priority sweep can leave a detectable record in either schedule. Prepare one doubled port with all collars blank and the sweep writes exactly once. Prepare two singly occupied ports and the same sweep leaves all collars blank. An observer who collects the local occupied-seam records can answer “Was any mark written?” in both schedules. On the one-double sector alone, a stronger restriction holds. Every fixed terminal measurement that gives the same probabilities on every doubled input when flight and the public sweep exchange order has a scalar input effect on each tested screen. This includes measurements of coherent sign within a winning seam. [Protocol](PROTOCOL.md) · [Observer counts and board certificates](results.json)

## Exact restriction on the doubled sector

Let \(F(x)\) be the highest-priority seam incident to port \(x\). A complete uninterrupted sweep \(S\) publishes the winning address. On the blank doubled input space, the address fibers \(F^{-1}(e)\) consist of one or two endpoints of \(e\). Within a two-port fiber the native write takes the two orientations to orthogonal written sign states. A terminal analyzer can mix these two states. The publication of the winner separates distinct fibers. Consequently \(S^*(O)\) belongs to the algebra \(\mathcal A\) of matrices block diagonal in \(F\), for every terminal effect \(O\). [Native sweep columns](../FLIGHT_CONFLUENCE.md) · [Strict sign and blank-binary column checks](../strict_pair.py)

Write \(U_t=\exp(-itH)\) for flight on blank doubles, \(H=L/2+24\kappa I\), with \(L\) the port graph Laplacian including the degree on its diagonal. The scalar includes the vacancy-reference contacts. Let \(T_t\) denote flight after the sweep. For one fixed terminal effect \(O\), equality of the two orders for every input requires

\[
 U_t^* E U_t=S^*(T_t^*(O)),\qquad E=S^*(O).
\]

The right hand side belongs to \(\mathcal A\). If this equality holds for all sufficiently small positive exposures, taking derivatives at zero gives \(E,[H,E],[H,[H,E]]\in\mathcal A\). These are necessary tests on a common effect. An analyzer chosen separately at each exposure, agreement at just one dwell, and a coherent service that never publishes the winner have different hypotheses.

To evaluate the first derivative, consider a two-port block \(\{x,y\}\). The graph certificate gives a singleton \(z\) adjacent to \(x\) and outside the neighborhood of \(y\). For \(E\in\mathcal A\), the off-block entry \([H,E]_{zy}\) is \(H_{zx}E_{xy}\). Its vanishing forces \(E_{xy}=0\). The transposed entry eliminates \(E_{yx}\). Thus \(E\) is diagonal. On each edge crossing distinct winners, \([H,E]_{uv}=H_{uv}(E_{vv}-E_{uu})\) forces equal diagonal values. Those edges connect precisely the ports of each individual cell. Conversely any cellwise constant diagonal passes this first test. Apply the same argument to \([H,E]\) using the second derivative. This commutator must be diagonal and already has zero diagonal, so it vanishes. Equality of the entries of \(E\) propagates along every edge of the connected board. Scalar effects pass at every order. With an untouched reference the argument applies coefficient by coefficient, giving the identity on doubled positions tensored with a reference effect.

| Cells | Singleton winners | Two-port winners | Dimensions after zero, one, two derivative tests |
| ---: | ---: | ---: | --- |
| 20 | 164 | 38 | 316, 20, 1 |
| 80 | 720 | 120 | 1,200, 80, 1 |
| 320 | 2,880 | 480 | 4,800, 320, 1 |

These are dimensions of the complex linear spaces of matrices satisfying the listed necessary conditions. The Hermitian parts have the same real dimensions. The incidence calculation checks the singleton witnesses, the first-derivative components and connectivity on each board. A second implementation constructs both commutators from matrix units and performs elimination modulo 101 on the 20-cell graph. Its ranks bound the rational ranks from below. The explicit cellwise and scalar solutions attain the complementary dimensions. [Incidence check](check.py) · [Independent commutator rank check](rank_check.py)

## A phase readout on one receiving seam

On the 20-cell screen prepare equal-amplitude superpositions of doubled ports 82 and 0 with relative sign \(+\) or \(-\). The preparation uses two native exchanges and a vacancy-reference phase. Port 82 shares its winner with port 1. Port 0 has a different winner. After one quarter exposure and the public sweep, the receiver undoes the write on seam (1,82) when that seam wins, mixes its two endpoint origins by a half exchange, and reads which origin it detects. In 262,144 simulated starts per phase it registered port 1 on **92,442** positive-phase starts and **61,106** negative-phase starts. The other two outcomes are the second origin and all other winners, including every unselected seam in the denominator. With zero exchange during flight the respective port-1 counts were 65,391 and 65,450, with expected probability \(1/4\) for both. If the public sweep goes first, its two occupied input components have distinct winners and the complete output density operators for the two phases coincide. Every subsequent observer reading then has equal probabilities. [All count banks](results.json)

## The two-load normal form

The two sectors have the same total load. Blank free exchange permutes occupation numbers and contrast is diagonal. Thus any flight before the sweep preserves the distinction between a doubled port and two singly occupied ports. A strict write maps the double to two unit loads with one occupied collar. All later free flights and strict attempts preserve that occupied address. The native strict instrument is quiet on every blank binary word. A sharp reading of whether any collar was written therefore distinguishes the two sectors even if a free interval is split around the complete sweep. For superpositions the two terminal outcomes have the Born weights of the corresponding input sectors, by linearity of their orthogonal outcome effects.

The resulting coarse repair relation has three classes: blank binary, blank double, and marked binary. Its only strict step maps blank double to marked binary. Assign rank one to blank double and zero to both normal classes. On a blank input the protected boundary bit is the indicator that some port holds two units, a property fixed by the source occupation multiset. On a marked output it is the indicator of any occupied collar. These two definitions use distinct instruments. Input preparations or sharp occupancy readings on calibration copies identify the former. The observer reads the latter from the published marks. The source bit never enters the sweep. It identifies the two normal classes, and every strict step decreases the rank. The Lean file proves completeness of the normal forms, the local diamond, descent, weak normalization and equality of normal endpoints from equal boundary bits for this three-state relation. The graph and native seam arguments above supply its physical interpretation. The registered hardware, initial blanks, fixed total load, published priority and complete sweep bound the result. A sweep that can be interrupted indefinitely, full-address confluence and arbitrary preparations require separate hypotheses. [Coarse relation and Lean proof](NormalForm.lean)

From the RER root the Python commands check the sealed count banks and graph certificates. In the RER Lean workspace, build the imported endpoint theorem before type-checking the local proof.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -W error -m model.experiments.interval_echo_consensus.flight_records.check
OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -W error -m model.experiments.interval_echo_consensus.flight_records.rank_check
cd Lean
lake build ObservableNormalForms.ObserverConfluence
lake env lean ../model/experiments/interval_echo_consensus/flight_records/NormalForm.lean
```
