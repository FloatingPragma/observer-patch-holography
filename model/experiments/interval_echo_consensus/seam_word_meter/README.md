# When a seam averages what both owners read

An owner first resolves the load at one end of an addressed seam. It restores the port and opens the seam for half an exchange. For every input on a protected blank collar, both endpoint observers then find the arithmetic mean of their two incoming **three-bin distributions**, after pooling the intermediate outcomes. A pair of coherent controls starts at the same mean load of ½ on each side. Omitting the first read sends the first-port mean to 0 for one phase and 1 for the opposite phase. Performing it returns ½ in both cases. These are probabilities calculated from Interval Echo. [Frozen protocol](PROTOCOL.md) · [All pointer and terminal probabilities](records.json.gz)

## The two readers

The working cell has twelve ternary ports. Its 30 internal seams supply the addressed repair alphabet. Stations on owners 1 and 3 have separate sharp zero and double standards. A sequence of native π exchanges transposes a queried port into the station, restores every transit port, then reverses the transposition after sensing. Each reference query uses the strict echo, a balanced-pair pointer, and the inverse echo. The resulting two flags are 01 for zero, 00 for one, and 10 for two. The return 11 has zero effect. The code checks the entire restored column, including phase and all memory entries, for each reference and load. The complete 99-state generator preserves the blank subspace throughout the working and routing pulses. A new offer can therefore revisit the same seam. Directly applying the strict echo twice to one collar on a sharp 02 input instead gives first-port loads 0, 1, 0.

Write S for the endpoint swap and Pᵣ for the projector on load r at the first port. In the blank block the contrast term contributes unit phase at π/2. The working pulse is U = (1−i)I/2 + (1+i)S/2. For a sharp bin effect D at either endpoint, expansion gives a diagonal part (D+SDS)/2 and a term proportional to DS−SD. The latter couples only an ordered pair ab with ba. If a differs from b, the first-port read separates the two columns. If a equals b, their commutator vanishes. Hence

$$
\sum_r P_r U^\dagger D U P_r = (D+SDS)/2 .
$$

The identity also holds on tensoring the effects with any reference observable. The operator on a spectator port stays fixed. Giving the three bins arbitrary real scores now gives endpoint agreement, conservation of their sum, and the scalar pair-averaging map. Reading the opposite end first gives the same six endpoint effects. The order of the two pointer flags is retained in the raw record. [Native apparatus](apparatus.py) · [Exact complex-rational matrix check](test_word_meter.py) · [Formal score and word identities](RecordMeans.lean)

A fixed finite calibration makes the load mean a bounded real-valued workspace. Every finite vector of real coordinates fits some common affine chart: choose lower and upper reference values in advance of its run. A data port paired with an auxiliary port in the 20/02 sector interpolates continuously between the two through a native partial exchange. A subsequent preparation check routes twelve independent stock ports to a common writing seam. Two common charts, [−3,5] and [−5,7], prepare and recover the same twelve mixed-sign coordinates. Their sharp work registers have 1,024 and 4,096 allowed joint words, respectively. The same offset and gain apply to every port throughout a scheduled program. [Preparation and held chart test](chart.py) Affine calibration commutes with every product of pair means. An individual sharp pointer returns a ternary bin. Conditioning the workspace on a particular intermediate flag can therefore change its value. Programs that feed such conditional values or nonlinear decisions into later offers require additional contracts.

## Arbitrary finite calendars

A separate dimer on owner 2 starts with one unit at its first port. For a requested crossing probability p, a contact lasting 2 arcsin(√p) gives that probability. A ternary read distinguishes a crossing from a stay. A full return exchange after a crossing restores the dimer on either branch. The independent program keeps the selected symbol in the ordinary observer record. The archive includes endpoint controls p = 0 and 1 and an irrational setting √2/3.

For an arbitrary distribution on length-k words over the 30 seams, add the weights of leaves sharing each prefix. At any positive-mass prefix, split the possible next symbols into two groups and use the dimer with the ratio of the second group's mass to the parent mass. Repeat until one symbol remains. A zero-weight group is never emitted, and an admitted symbol always makes an offer, including on equal inputs. Multiplying the conditional probabilities along a path returns its specified leaf weight. The empty word has unit mass. The construction allows correlations among positions and uses an externally supplied finite record/controller and calibrated contact angles. The iterative controller tracks exact rational or explicitly specified algebraic weights. The formal split applies to ideal real weights and contact angles. The numerical dimer diagnostic reports unresolved probability splits explicitly when they fall below its precision. A 1,401-symbol control retains a positive branch of weight 10⁻⁴⁰⁰ in the exact law.

## All measured words

Both orders of every pair of the 30 seams, including all 30 repetitions, were replayed from each of the twelve one-load source ports. The 10,800 calculated terminal distributions match two pair means exactly in quarter-launch units. All intermediate read combinations are archived. A load-two result has zero probability throughout this one-unit bank. The universal local operator identity also covers the other ternary sectors and arbitrary input coherences.

The 30 distinct maps recovered from the repeated-offer terminal tables have uniform mean I−L/60, with L the native cell graph Laplacian. Its slowest nonconstant level is 0.95393446629 (multiplicity three). The other nonconstant levels are 0.9 (five) and 0.87939886704 (three). The centered response has the rank-three slow space used in the finite Gram construction.

In the held four-offer mixture the first two touching seams are 01 and 05, and the disjoint seam is 23. For source port 0 the complete receiving probabilities at ports 0, 1 and 5 are 0.30625, 0.35416667 and 0.33958333. The other nine ports are dark. Changing the pulse to π/3 raises the first reading to 0.430078125. Keeping π/2 while omitting the input read lowers it to 0.08125. Randomizing each position independently with its original marginal changes the total variation of the terminal record by 1159/27000. A uniform-amplitude quarter-phase input is uniform after the nominal calendar. Omitting the pre-read produces nonuniform outputs even there. [Program and complete schedule weights](measure.py)

The scalar observer process supplies PR-146 and PR-147. For PR-180, a second apparatus lets two stations read the exact midpoint of sharp material words and use it to choose their next seam. An independently prepared pooled apparatus makes sign decisions from fresh cohorts with simultaneous error coverage and can abstain near a threshold. [Sharp readout and branch theorem](finite_program/README.md) · [Pooled comparison](finite_program/pooled_feedback/README.md)

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -W error -m model.experiments.interval_echo_consensus.seam_word_meter.measure
python -m pytest -q model/experiments/interval_echo_consensus/seam_word_meter
(cd Lean && lake env lean -j1 ../model/experiments/interval_echo_consensus/seam_word_meter/RecordMeans.lean)
```
