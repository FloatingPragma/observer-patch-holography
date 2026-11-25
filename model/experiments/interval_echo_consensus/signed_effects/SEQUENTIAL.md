# A native measurement with its state update and second reading

One common contact/reference word measures a different coherent basis for a vacancy and an excess. An exhaustive sharp address read followed by the inverse word leaves the measured eigenstate on the same material rails. Repeating the instrument returns the same label. Appending a native swap preserves all first-read effects and changes both the residual state and the second-read statistics. This realizes a finite instrument on the specified protected menu, including its state update. Quantum sharp-pointer semantics and pulse controls are supplied. [Replay and complete banks](sequential.json) · [Native implementation](sequential.py) · [Independent map tests](test_sequential.py)

## The same physical word for both signs

The input is one vacancy or excess on ports 0 and 1 of a unit-filled twelve-port cell with blank memories. Every other contact is parked. The native 99-state operator restricts to

\[
H_{01}=\begin{pmatrix}17/2&-1/2\\-1/2&17/2\end{pmatrix},
\qquad D_s=8I+q_s|0\rangle\langle0|,
\qquad(q_-,q_+)=(-8,24).
\]

The first operator vanishes on the other ten positions. The second has spectator eigenvalue eight. The chronological word is a contact for \(\pi/2\), the reference at port 0 for \(\pi/20\), and a contact for \(3\pi/2\). Write its unitary as \(U_s\). Both generators have integer spectra. The chronological positive-wait inverse uses contact \(\pi/2\), reference \(39\pi/20\), and contact \(3\pi/2\). All settings are identical for both signs.

The observer reads the complete sharp address \(y\in\{0,\ldots,11\}\) between these words. Its residual material state is propagated through the inverse, without replacing it by a prepared source. With \(P_y=|y\rangle\langle y|\), the actual branch operator is

\[
K_{s,y}=U_s^*P_yU_s=E_{s,y},\qquad
\Phi_{s,y}(\rho)=E_{s,y}\rho E_{s,y}.
\]

The twelve orthogonal effects sum to identity. The maps are completely positive, their sum is trace preserving, and \(K_{s,z}K_{s,y}=\delta_{zy}K_{s,y}\). These identities hold on the full twelve-position source and after adjoining an arbitrary untouched reference. The read is a sharp position measurement. No coarse projector preserving an unread within-block phase is assumed.

Ignoring a common phase, on the two active rails

\[
U_s=\begin{pmatrix}\cos(\phi_s/2)&\sin(\phi_s/2)\\
-\sin(\phi_s/2)&\cos(\phi_s/2)\end{pmatrix},\qquad
(\phi_-,\phi_+)=(-2\pi/5,6\pi/5).
\]

Both measured bases have nonzero off-diagonal effects. This is an explicit finite pulse construction; the Lie-controllability existence theorem is unnecessary for its implementation.

## Equal effects with a different state update

In the control arm append a contact for \(\pi\) after each inverse word. This is exactly the swap \(X\) on the active rails, with identity on every other position. Its branch is \(L_{s,y}=XK_{s,y}\), so \(L_{s,y}^*L_{s,y}=E_{s,y}\) on every input, including reference-entangled inputs. Its residual state is \(X E_{s,y}\rho E_{s,y} X\). The following measurement uses that state.

For every normalized input supported on the two active rails, the probability of different first and second labels is

| Instrument | Vacancy | Excess |
| --- | ---: | ---: |
| Read and undo | 0 | 0 |
| Read, undo and native swap | \((3-\sqrt5)/8=0.0954915\ldots\) | \((3+\sqrt5)/8=0.6545085\ldots\) |

The control changes the operation within the same protected menu. First-read statistics cannot distinguish these instruments. The complete ordered label pairs can. A claim that effects alone determine the poststate is falsified by this pair of implemented instruments. A nonzero repeat-label mismatch would reject the ideal read-and-undo claim; the control provides a representable mismatch on the same reader.

## Records, error and scope

Four input settings are port 0, port 1, their equal real superposition, and their equal superposition with relative phase \(i\). Each sign and arm has 16,384 synthetic starts per input, totaling 262,144. Every bank retains all 144 ordered address bins, including the zero-support complement. The receipt also retains the subnormalized first-outcome matrices on rails 0 and 1; all other matrix entries and the other ten first-outcome states are zero for these inputs. Those matrices are audit data. The observer record is the ordered address pair. Sampling uses the analytic rank-one binary tree and is checked against propagation of the full instruments; it represents no count of literal native-event executions.

Independent exact two-level algebra checks the effects and both disagreement probabilities. Matrix-unit checks compare every first-read effect and every repeat branch on all twelve positions. An untouched two-level reference checks the unread reference marginal and distinguishes the two labeled output channels. Changed integer bins or source hashes fail replay. Floating map serialization may vary by at most \(10^{-11}\); this is a replay tolerance, with no claim of physical measurement precision.

If preparation has half trace error at most \(\epsilon_p\) and every native pulse duration has absolute error at most \(\delta\), the complete two-use classical/quantum output has half trace error at most

\[
\min(1,\epsilon_p+218\delta).
\]

One read-and-undo has four contact pulses of operator norm at most nine and two reference pulses of norm at most 32. The optional swap adds one contact. Telescoping the unitary errors, contractivity of the complete sharp read and two instrument uses give \(2(4\cdot9+2\cdot32+9)\delta=218\delta\). No outcome probability divides this bound. A physical pulse-error or trusted-pointer certificate is not supplied.

The bank is a software construction and replay, without a prospective data claim or external authentication. It supplies a complete finite state-update and sequential readout on this declared menu. The complete event grammar, physical preparations and pointer fidelity, source-selected controls and producer/run custody are separate PR-03, PR-64 and PR-65 obligations. No physical particle identity or OpenWave award follows.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_consensus.signed_effects.sequential
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/signed_effects
```

The replay never writes evidence. Explicit `--write` creates a missing receipt exclusively and refuses to replace an existing one.
