# Reading a prepared-lock decision at two observers

Four Interval Echo owners exchange signed votes for a binary value. Two additional observers can decide independently from the certificate they receive. At quorum three, a leader that offers both values to different signers obtains a decision for zero at both readers. At quorum two, one equivocating signer and two separate honest signers produce valid conflicting certificates. The two readers then decide zero and one. With one silent signer, quorum three succeeds in the next view and quorum four stalls before a proposal. [Frozen protocol](PROTOCOL.md) · [Signed messages and complete pointer outcomes](records.json.gz)

| Separate enrollment | Quorum | Signed packet transfers | Decisions of readers 4 and 5 |
| --- | ---: | ---: | --- |
| Omitted orphan lock, then a decision and view change | 3 | 63 | 1, 1 |
| Competing older and newer locks in one justification | 3 | 58 | 1, 1 |
| Equivocating leader | 3 | 29 | 0, 0 |
| Silent leader, then a live view | 3 | 29 | 1, 1 |
| Silent leader at the higher threshold | 4 | 6 | neither |
| Equivocating leader at the lower threshold | 2 | 30 | 0, 1 |

## Which observations authorize a signature?

Owner 0 leads the first view. Leadership then rotates through the four validators. A view-change report carries the sender's highest assembled lock, if any. The leader checks the signatures, the view, the number of distinct enrolled signers and the certificates they report. It selects the value at the highest reported view. It may choose either bit when a full quorum reports no lock. Each recipient repeats the checks on its **received** proposal. An honest participant remembers the first value it supported in its current view across prepare, acceptance and commit. Entering the next view closes its earlier signing slot.

A raw prepared certificate contains a quorum of distinct prepare votes. A provisional acceptance signs the digest of that particular prepared certificate. A lock certificate contains both the prepared certificate and the acceptance quorum. Before an honest signer issues a commit, it verifies the complete lock certificate, stores it in its own classical record and writes its value to a persistent local dimer. A read with the native pointer checks the write. The balanced dimer selects the corresponding certificate from the owner's record during the next view change. A decision certificate carries that lock together with a quorum of commits signed over its digest. Both decision readers authenticate and decode the entire received packet. In the orphan trial, reader 4 forwards its received certificate to reader 5. All four durable dimer values also receive an independent routed read by owner 5. [Observer guards and validators](certificates.py) · [Acquisition and controls](measure.py)

In the orphan trial, three provisional acceptances leave every durable dimer untouched. Only owner 1 initially receives the assembled lock and commits to zero. The next justification comes from owners 0, 2 and 3. Owner 1 signs a new prepare for one and retains its old durable lock through acceptance. Receipt of the assembled new lock changes the material selector before its commit. All delayed lower-view prepare, acceptance and commit requests are refused. In another trial, the reports carry a view-zero lock on zero, a view-one lock on one, and one empty report. A fresh candidate of zero gives a proposal for one. A leader-signed request for the older value is refused. The registered checks also retain rejected duplicate and short quorums, a raw prepare presented as a lock, a modified signature, a different valid assembled parent for the same value, and a changed decision value.

## Native channel and authentication boundary

Each byte is serialized into eight balanced, sharp material pairs on fresh appointments. The sender holds (0,2) or (2,0). Addressed nearest-neighbour π exchanges transport the two rails to the destination. The receiver queries each rail against load two using the strict echo and inverse. Its two pointer returns are 10 and 01 respectively. The path is undone and every transit register is restored. The signed digest, signature and complete nested evidence travel through the same channel as the value. With the inter-owner transit exchanges parked on the held route, both inputs give 01. The resulting constant byte stream fails to parse and the reader refuses.

The simulated channel uses the exact reduction of every sharp bit slot to its two checked source columns. The receiver uses its measured two-rail code and canonical packet parser. Across six enrollments the receipt represents 215 packet transfers and 1,137,624 bit slots. It stores each unique packet, the two full raw response pairs on each route, the wire order, local signing refusals, dimer reads and final decisions. These deterministic numbers count represented logical slots. The underlying 99-state operators and the route permutations are evaluated once per calibrated setting. In the persistent selector the two material rails hold the value; the enrolled local classical tape retains the complete signed certificates, including view and issuer. Ed25519 keys are generated independently for every enrollment. The deposit publishes only the 24 public keys and 130 signed-message entries. Replay verifies those original signatures and neither regenerates nor publishes their private keys. Signer authority and private-key isolation are supplied assumptions. [Channel implementation](transport.py) · [Independent raw receiver and signature checks](test_quorum.py)

## Theorem and claim scope

The [finite-set proof](Quorum.lean) gives a nonfaulty witness in the intersection of any two size-at-least-q quorums whenever 2q ≥ n+f+1. It also proves the complement of the faulty set has at least q members when q+f ≤ n. For four validators at threshold three all 64 ordered quorum-pair and single-fault cases have an honest witness. At threshold two, 48 of 144 cases have none. The cross-view induction and progress hypotheses belong to [Appendix B](../../../../paper/appendix_B_bft_qecc_extensions.tex); the [earlier finite combinatorial checks](../../../../code/consensus/verify_issue_517_proof_obligations.py) cover additional sizes and orphan-lock schedules. The current addition implements the precise certificate and durable-lock decisions using the material source and its observer records.

PR-122 is demonstrated for this declared protocol. The silent-signer experiment supplies an addressed schedule and a timeout intervention. The arbitrary-view liveness theorem still assumes partial synchrony and the required quorum connectivity. For example, the registered route from owner 0 to 3 traverses owner 4, outside the set {0,1,3}. This experiment supplies neither a physical delay bound nor connectivity within every quorum. PR-120 and PR-144 also ask for semantic read-set closure or a single atomic label per conflicting repair component. Certificate agreement has no transition that performs that aggregate repair.

From the repository root, replay the public deposit with one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
nice -n 15 .venv/bin/python -B -O -W error -m model.experiments.interval_echo_consensus.prepared_quorum.measure
nice -n 15 .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/prepared_quorum
(cd Lean && nice -n 15 lake env lean -j1 ../model/experiments/interval_echo_consensus/prepared_quorum/Quorum.lean)
```
