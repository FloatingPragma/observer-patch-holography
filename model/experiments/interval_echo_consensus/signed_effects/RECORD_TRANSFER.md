# A measurement record stored in material dimers

The two populated source rails can transfer their address to a material dimer through an explicit native contact/reference word. Two prepared dimers retain the first and second measurement records coherently until one final occupation read. Their branch maps reproduce the two-read instrument, including reference-entangled inputs. No intermediate source projection or abstract controlled copy is used. The final sharp occupation instrument is supplied. [Complete records and words](record_transfer.json) · [Native execution](record_transfer.py) · [Independent tests](test_record_transfer.py)

## Domain and physical resources

The source occupies port 0 or 1. Pointer dimers occupy park/cross ports (6,10) and (4,5). Each holds one mark of the same sign as the source: either three vacancies or three excesses in otherwise unit-filled material. The input pointers are at ports 6 and 4. The identical pulse settings act in both signed preparations. This does not construct a fixed ancilla that works on an unknown source sign. Matching the pointer sign is a preparation resource.

The five contacts used are (0,1), (6,10), (0,10), (4,5) and (0,5), all native seams of the twelve-port cell. Every other contact is parked. The degree-reference terms and protected durations are the specified control menu. All twenty placements of three marks on these six ports are simulated, including states outside the one-mark-per-dimer code. The protocol replaces the pointer on the two populated source rails. It does not replace an exhaustive source pointer on twelve arbitrary positions.

The existing echo-read-undo instruments and reversible reader words implement declared sharp predicates and copy apparatus. The row-30 construction also uses controlled phase on its prepared connectors. Here the binary blank-memory contacts and port references compile the material copy itself. No supplied reader XOR or arbitrary unitary executes that step.

## Refocused native interaction

On the binary endpoint alphabet for either matched sign, a blank seam has

\[
h=(I-\mathrm{SWAP})/2+8(n_a-n_b)^2.
\]

A degree-reference pulse of duration \(\pi/8\) at \(a\) gives \(Z_a=(-1)^{n_a}\). Its relative phase on a marked endpoint is minus one for both signs. On the four endpoint words, the two toggled generators commute and

\[
h+Z_a h Z_a=17\,\mathrm{diag}(0,1,1,0).
\]

Apply contact \(\pi/34\), reference at \(a\) for \(\pi/8\), contact \(\pi/34\), and the same reference again. Exchange cancels exactly, and the unequal endpoint words receive phase \(-i\). Reference pulses \(\pi/16\) on \(a\) and \(b\) each give relative marked phase \(i\). Their combination is controlled-Z on marked endpoints, up to one common minus sign. The native 99-state blank-memory columns and all twenty physical placements are checked before using this reduction. This gate uses the declared strength eight; it derives no coupling or physical time scale.

On a pointer dimer, apply reference \(3\pi/16\) on its cross rail, contact \(\pi/2\), then the same reference. Its encoded action is Hadamard up to one scalar. Conjugating the controlled-Z by this word and its positive-wait native inverse gives minus controlled-NOT on every source/pointer code word. The full operation is unitary on all twenty states. In particular, a sharp source port is unchanged while the pointer moves to its cross rail exactly when the source is at port 0. The source phase is retained when the pointer is kept coherent.

This copies an orthogonal address record. It entangles a superposed source with its pointer, without copying its arbitrary quantum state.

## Two stored records and final read

For each pointer, rotate the source by the specified contact \(\pi/2\), reference \(\pi/20\), contact \(3\pi/2\) word. Execute the material copy, then undo the source rotation. After the second copy, read only the occupations of the four pointer ports. The two minus signs from the material copy cancel. With \(E_{s,y}=U_s^*P_yU_s\), each decoded branch on the source is exactly

\[
K_{s,y_1,y_2}=E_{s,y_2}E_{s,y_1}.
\]

All repeatable branches with unequal labels vanish. In the control arm one native source swap is inserted between the acquisitions. The branch is \(E_{s,y_2} X E_{s,y_1}\). The first effect is unchanged, while the probability of unequal labels is \((3-\sqrt5)/8\) for vacancy and \((3+\sqrt5)/8\) for excess, for every admitted normalized source input. The final source state is retained in both cases.

Each pointer's occupation word is 0 (empty), 1 (park), 2 (cross) or 3 (both). The final read retains all sixteen paired words, including malformed pointer configurations. Four source preparations, two signs and two arms give sixteen synthetic banks of 8,192 starts, totaling 131,072. Their complete final branch maps include every complement before any decoding. No middle-of-run state is replaced or sampled. The simulator uses matrix propagation and rank-one tree sampling, with no claim of literal native-event call counts.

## Coherence, controls and error

Independent exact algebra checks refocusing and its shared marked phases. Every encoded source/pointer word is checked in the full physical sector. A source entangled with an untouched two-level reference gives the same complete labeled branch states as the direct projective instrument; the unread reference marginal is unchanged.

Keeping the first pointer unread and reversing its actual native word restores a complex source superposition. A complete pointer read with its outcome ignored, followed by the same inverse, instead leaves the source dephased. Its return probability to the original equal superposition is one half. This distinguishes a coherent stored record from an unread classical record.

Omitting the cross contacts makes the two sharp source inputs give identical pointer laws, whereas their native record laws have total variation one. Omitting the parity refocusing pulses gives resolved population outside the one-mark-per-dimer code. These controls test record transfer and nondemolition on the same source and readout.

For preparation half-trace error \(\epsilon_p\) and absolute duration error at most \(\delta\) on each native wait, the complete two-record classical/quantum output has half-trace error at most \(\min(1,\epsilon_p+793\delta)\). Each copy uses four contact pulses of norm at most nine and eight reference pulses of norm at most 32. Each source rotation/undo adds four contacts and two references. Two acquisitions plus the optional source swap give total generator norm budget \(2(4\cdot9+8\cdot32+4\cdot9+2\cdot32)+9=793\), across at most 37 waits. Telescoping and contractivity require no division by a herald probability. This is a conditional error contract, without physical pulse certification.

The final four-port occupation read and its projective state-update law, matched-sign pointer initialization, protected isolation, physical fidelity and external custody are supplied or separate obligations. The result supplies a finite native record-transfer implementation under those hypotheses. It does not retire PR-03, PR-64 or PR-65.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_consensus.signed_effects.record_transfer
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/signed_effects/test_record_transfer.py
```

Replay preserves the receipt. Explicit `--write` creates a missing receipt exclusively. Source hashes bind producer-side software execution and supply no external authentication.
