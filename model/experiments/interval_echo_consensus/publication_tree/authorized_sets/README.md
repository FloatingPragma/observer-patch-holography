# Which readers can distinguish the 24 oriented records?

Twelve separate stations use the collective reader family `universal-twelve-port-publicness/v1` fixed by the finite capacity packet. Each station can read one bit at its own address. On 20-, 80- and 320-owner Interval Echo federations, the collective recovered every one of the 24 oriented records after each of the 40 reported continuations. Receiving eleven named readings preserves all 24 under a separate diagnostic decoder. A code that must work for *every* choice of ten stations can retain only two records. All settings, reader omissions and controls are fixed in the [protocol](PROTOCOL.md). The [simulation archive](records.json) gives exact sharp detector outcomes for 2,880 combinations, without statistical sampling.

## From the port packet to local readings

The [fixed packet](../../../../../code/capacity_readback/runtime/source_derived_public_checkpoint_packet.json) names twelve ports and two orientations. Its nonempty publicness family consists of one group containing all twelve. In the material implementation, each port names a distinct enrolled owner on the board. The enrollment order specifies the coordinate dictionary. It makes no identification of physical distance between two named stations with an edge length of the canonical icosahedron.

Each owner reads one native equality echo against a parked vacant rail. One high digit identifies the write port; one low digit identifies the check port. These are supplied preparations. For a port permutation the material controller exchanges both rails of each named transposition through actual adjacent contacts. A forward path followed by its return, omitting the final repeated step, transposes its endpoints and restores every intermediate address. A reversal exchanges the local rails at every station. Each physical pair pulse and sharp pointer read inherits the [99-state check](../../../../candidate_laws/interval_echo/clock.py). The apparatus uses no simulator history or time derivative to identify a label. [Native preparation and readers](apparatus.py)

The decoder reads the exceptional bit and its address. Its 960 answers per board agree with the independently deposited joint-kernel rows. A separate internal check verifies that every continuation returns the other rails to their prepared values. The deposited 1,600 composition entries agree with compositions of the measured tables on each board. Induction on the length of a known setting sequence then transports the records through every finite continuation word. With inter-owner exchanges parked during the held compound continuation, every one of the 24 full detector profiles fails its prediction. This control retains the raw local pointer readings if its routes disturb a reference. [Comparator, support calculation and acquisition](measure.py) · [Lean relations for record separation and composition](Access.lean)

## Why losing two readers costs more than one

Two distinct inputs of the same orientation differ at precisely two addresses. Opposite orientations differ at twelve addresses for the same named port and at ten for different ports. Thus every eleven-reader subset separates all 24 inputs. For every pair of same-orientation inputs there is a ten-reader group seeing the same ten values on both, namely a group omitting both exceptional addresses. Each orientation class forms a clique in the compound confusability graph, so a code with two orientation classes has at most two mutually distinguishable words. The two opposite orientations at one chosen port complement each other at every address. They attain the upper bound down to a single reader.

The certificate evaluates all 4,095 nonempty named subsets at identity, then computes the largest possible common-read support for every pair of inputs over every known continuation. At a *particular* ten-reader subset, 22 different read profiles occur at identity. For the promise to recover a common code through *every* ten-reader subset, the maximum code size is two. The group fixed in the publicness rule always requests all twelve named positions; a missing position leads to a primary refusal even when the separate eleven-reader decoder could infer it. [Exact clique covers and decoder witnesses](records.json)

## What one station can learn by communicating

On the 20-owner board, one addressed station's own bit is insufficient to identify all 24 inputs. A second apparatus sends one tagged bit in each of twelve rounds over local register handoffs. Its terminal reader retains the vector of twelve local echo reads across the rounds. All 288 terminal profiles (24 inputs, twelve receivers) give the individually named atoms of the fixed packet. The held compound continuation is a bijection, so these 24 inputs also exhaust the relay's codeword alphabet. Omitting just the first tree handoff of round one makes seven enrolled readers receive the distinct missing tag. Across the 24 inputs there are 168 refusals and 120 complete reads. The relay gives a physical operation that distinguishes a local bit from a collectively available word.

The setup realizes the specified publicness policy and finite readout on addressed sharp registers. It supplies the 24-input encoding, twelve selected addresses, requested continuations and routing schedule. The capacity calculation is for that prepared alphabet. A value for the public capacity of a cosmological horizon requires an independent physical source, accessible continuations and a geometric identification. [Fixed-cutoff theorem and horizon requirements](../../../../../paper/tex_fragments/OBSERVERS_N_CLOSURE_PRO5.tex)

From the repository root:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.experiments.interval_echo_consensus.publication_tree.authorized_sets
.venv/bin/python -B -m "$pkg.measure" --write
.venv/bin/python -B -m pytest -q --pyargs "$pkg.test_checkpoint"
(cd Lean && lake env lean -j1 ../model/experiments/interval_echo_consensus/publication_tree/authorized_sets/Access.lean)
~~~
