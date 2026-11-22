# Private record walk

Two unit loads and their comparison signs move on the [fixed OPH carrier](../../executable_model/README.md). A material move exchanges a load with a vacancy across a seam. An enabled sign move reverses one private bit. Each distinct move is one undirected edge of the configuration graph. The free generator is half its combinatorial Laplacian. The clock unit is the splitting of one isolated edge. There is one time scale and no input for a pair energy. [Microscopic neighbors](../../experiments/interval_echo_consensus/sign_graph.py) · [Candidate interface](candidate.py)

At a contact, compare the transported signed words on as many different incident bystander seams as one instrument can address everywhere. One bit belongs to each question at each contact. Unary comparison bits bring every port up to the largest single-load count. These choices specify the aperture, address order and one-load reference. All switches are open during free release. Only preparation, survey and terminal readout use addressed interventions. The terminal patch reads one of its thirty internal seams or quiet. In the sign experiment it also reads the comparison bits at the launch contact on the same copy.

## Predicting pair dynamics with a single load

On the registered 20-cell board, the probability of recovering both loads on an internal seam of cell 13 exceeds **94.099% throughout exposures 0 to 37 for every state of the eight launch signs**. The other 6,472 signs are prepared odd. Setting *every* sign even gives 20.291% at exposure 37. Thus the interacting background matters, while the launch contact needs no parity selection. The all-odd preparation exceeds 97.803% throughout the window. Starting the eight signs in sharp zero exceeds 97.120%.

| Exposure | Home pair, sharp launch signs | At least one changed sign | Home pair and at least one changed sign |
|---:|---:|---:|---:|
| 1 | 0.985553 | 0.078542 | 0.071671 |
| 7 | 0.983486 | 0.267414 | 0.262394 |
| 37 | 0.978566 | 0.868164 | 0.850549 |

Hold a sharp single load at port $a$, initialize the sign pointers to zero and apply a comparison pulse of duration $\pi$ in exchange units. An enabled bit evolves by $(I-X)/2$ and flips with certainty. An idle bit keeps its value. The native signed write supplies the enabling predicate. Swapping unequal endpoint loads changes each selected bystander's signed word. Equal endpoint loads leave it alone. Interrogating 00 and 11 tests this last statement separately. A rule enabling a flip at either occupied endpoint gives identical single-load records and a different 11 record. These are computed finite instrument laws. [Survey and local checks](probe.py)

Let $F_{aj}$ be the bit read at address $j$ for one load at $a$, and $S_a$ the set of flipped bits. On the distinct binary pair $(a,b)$, the native comparison enables exactly $S_a\mathbin{\triangle} S_b$. In an odd sign state each enabled move costs one. Thus the two-load cost equals $|S_a|+|S_b|-2(FF^T)_{ab}$. The two single-load terms are constant on the prepared carrier. The remaining term, including its sign, is predicted on every pair by the overlaps of the single-load records. On the three federations the overlap is eight on a seam and zero on every nonseam. This gives the complete prepared material generator $K-16P$, with $K$ material exchange and $P$ adjacency. The implementation constructs the Gram matrix before reading any pair outcome. A changed law using the union of the two sets gives $K-8P$. Preparing every sign even gives exchange alone. Their terminal laws use the same detector.

The ordinary preparation puts every sign in its odd state. Both outcomes of a local parity measurement can be used to prepare this state by feedback. Its supplied apparatus and finite log are described in the [preparation](../../experiments/interval_echo_consensus/README.md#what-fixes-the-binary-interaction). In the second preparation the launch contact's eight signs start in their sharp zero states, with every other sign odd. Every one of their 256 parity components participates in flight. There is no selection on a preparation result or a final word.

## Joint readout of the material and its signs

For a particular parity pattern with $m$ odd launch signs, subtract the **same** two-load rest reference as in the all-odd case. The material generator is $H_m=K-16P-(8-m)D$, where $D$ selects configurations with exactly one occupied launch endpoint. Hence nine flights give all 256 amplitudes. For a specific sharp output word with $r$ ones, its unnormalized material amplitude is

$$
 2^{-8}\sum_{m=0}^8 \left[\sum_j (-1)^j {r\choose j}{8-r\choose m-j}\right]e^{-itH_m}\psi_0 .
$$

Square after summation. Multiplication by ${8\choose r}$ gives the recorded weight bin. Individual words of equal weight are equiprobable. Quiet sums the probabilities of all configurations outside the thirty home seams, including two loads on separated ports in the home cell. For a reader ignoring the signs, any initial sign density matrix gives a convex mixture of the nine material probabilities, with weights given by its parity populations. This includes initially correlated launch signs with an untouched reference and a fixed sharp material pair. The window census retains every sector.

At launch both endpoints are occupied, so a comparison sign cannot flip until a load moves. Nine material exchanges can split the two occupied endpoints. For each split and each launch sign there is one path of length two to a one-flip record. Its amplitude starts at $-t^2/8$. The probability of any changed word starts at $9\cdot8t^4/64$. The either-occupied rule enables all eight flips immediately and gives $8t^2/4$ to leading order. An inert sign register attached to a material Hamiltonian gives zero. The sign readout distinguishes these mechanisms without interpreting an eigenvalue as a mass.

The Taylor replay uses the same energy origin for every sector. Independent removal of sector phases changes the measured word probabilities. The normalized Walsh transform is unitary, so its full-state error is bounded by the largest sector state error. Its nine-term dot products contribute at most $9\gamma_{20}$ in norm, with $\gamma_n=nu/(1-nu)$ and binary64 unit roundoff $u$. Summing positive probabilities on at most 28,680 material configurations per word weight and then nine weights costs at most $\gamma_{28690}$ times the total norm squared. The replay allows an additional $2\times10^{-10}$ for these operations. The replay takes the home probability and its first derivative at both ends of each interval of length $h=1/16$. For $A=H-\langle H\rangle$, exact integer powers of $2A$ on the launch bound $M=\sum_{k=0}^4 {4\choose k}\|A^k\psi_0\|\|A^{4-k}\psi_0\|$. Differentiating the home probability four times bounds $|p^{(4)}|$ by $M$. The cubic Hermite interpolation error is at most $Mh^4/384$. Subdividing each cubic three times and taking its smallest Bernstein coefficient bounds its minimum from below. The endpoint derivative error is bounded by $2\|H\|$ times the probability error allowance because the derivative measures $i[H,P_{\rm home}]$. This also bounds mixtures by linearity. Exposure labels are apparatus settings. A uniform window bound applies at every reading in that interval. [Propagation and complete outcomes](release.py) · [Numerical record](results.json)

## Provenance and replay

The independent checker reads the hash-verified carrier and builds its exchange graph from allowed single hops. It evolves directly in the sharp-word basis grouped by Hamming weight. Between normalized weights $r$ and $r+1$ the off-diagonal is $-\sqrt{(r+1)(8-r)}D/2$. This gives 258,120 coupled coordinates. It compares every one of the 279 joint bins at each appointment. On a four-port cycle the tests also enumerate every microscopic sharp configuration with all eight private bits explicit. [Independent calculation](check_words.py) · [Literal graph test](test_walk.py)

This is the [private-record graph](https://github.com/FloatingPragma/observer-patch-holography/commit/f41a6669e) on the binary material sector. The design choices include private addresses, maximal portable aperture, one common edge rate and preparation of the spectator signs. The historical Interval Echo design tested strengths 0, 0.5, 1, 2, 4, 8, 16, 64 and 128 with persistence in view. Its results at 8 and 128 and the odd-bank transfer were known when the sharp-sign protocol was specified. The returned material energy is derived for this graph. The historical pin supplies its coefficient. OpenWave scores that pin and has awarded partial credit on row 4. [Frozen inputs and gates](protocol.json) · [Design sweep](https://github.com/FloatingPragma/observer-patch-holography/blob/972d81b9b6068938e52457b59497766499222b5d/model/polynomial_sv/SPEC.md) · [OpenWave decision](https://github.com/openwave-labs/openwave/pull/745#pullrequestreview-6005365373)

```sh
python -B -O -W error -m model run --law private_record_walk --check
python -B -m pytest -q -W error model/candidate_laws/private_record_walk
python -B -O -W error model/candidate_laws/private_record_walk/check_words.py
```
