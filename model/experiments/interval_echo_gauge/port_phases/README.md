# Twelve physical port phases and the five-level gauge code

Two stations on the registered 20-cell Interval Echo federation can write and read twelve independently chosen continuous phases. A second experiment connects those controls to the five-level code of the gauge paper. On an initially unit-filled cell, a single excess realizes the specified flux-two triplet and a vacancy realizes the flux-seven doublet under a programmed apparatus. One common material detector measures membership of the two bands without resolving an admitted mode. The same two ordered kicks give different terminal records when the order is reversed. These are synthetic records of the candidate's declared quantum instrument. [Gauge specialist, “The global form from port phases”](../../../../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf)

## How the two stations measure a port phase

The native degree-five reference gives a phase proportional to \(n^2\) on ternary load \(n\). An addressed blank contact with a unit-load guard has generator \((I-S)/2+8(n-m)^2\). In each unequal block accessible from \((n,1)\), the data parity \(e^{-i\pi n^2}\) reverses the hop. Two equal dwells separated and followed by that parity give \(e^{-it(8+1/2)(n-1)^2}\) and return the guard to one. The equal word is fixed throughout. The identity \(n^2-(n-1)^2=2n-1\) thus gives \(e^{-i\theta n}\) up to a common scalar, for any angle. Periodicity supplies positive waits for inverse phases. Native duration-π swaps route each port to this service and restore every transit port. Applying the construction to all twelve ports gives the tensor product of twelve independent phases, including for a data state entangled with other rails. [Native reductions and physical routes](apparatus.py)

An emitter splits unequal ternary loads on two rails. It sends one rail to the cell and one along a reference. A second owner receives the two rails and measures the interference in two quadratures. Its signed head reads each terminal load against a transported unit witness. Both roles are exchanged on fresh launches. The receiver uses measured 01 rates to predict 12 and 02 rates. A π/2 change on each individual port tests addressing on all other ports. Each bank retains all 27 load/guard words. The separate free-flight preparation records the vacancy on the original edge or outside it. [Frozen settings](PROTOCOL.md) · [Complete banks](records.json) · [Counts-only analysis](measure.py) · [Numerical receipt](summary.json)

| Readout | Maximum absolute frequency discrepancy |
| --- | ---: |
| Native writer, held 12/02 predictions | 0.00510 |
| Direct square writer, same 01 calibration | 0.99265 |
| Exchange-only law, same waits | 0.45465 |
| Native writer, all single-port changes and untouched ports | 0.00512 |

The 728 banks contain 95,420,416 original starts. Every native guard returned to one. The simultaneous 99% radius per raw bin is 0.00761, propagating to 0.02936 for the held predictions. The exchange-only law produced nonunit guards in up to 10.09% of starts. In free flight, an active π kick changed the vacancy return probability by 0.05127 at time 1/16. With exchange parked the measured change was −0.00074. Global phase conventions and physical port writes can therefore be distinguished experimentally.

## The charge in the actual echo records

A written seam head records a coin \(c\) and an oriented old contrast \(d\). Its signed delivered load is
\[
J(c,d)=\big[d+(1-2c)(d\bmod 2)\big]/2.
\]
On every write branch it equals the load lost by the first endpoint and gained by the second; set \(J=0\) on a blank. The native echo therefore intertwines \(Q_a=n_a+J\) and \(Q_b=n_b-J\) across blank and written heads. Exact integer entry checks cover its complete 99-dimensional local space. Summing oriented head contributions at each vertex gives commuting integer-valued \(Q_p\) preserved by any sequence of echoes. A diagonal detector of the written records commutes with them. An arbitrary spectator reference carries its own identity action. [Lean proof of the head identity](../../../../Lean/SourceLaw/IntervalEchoHeadCharge.lean) · [Complete echo certificate](certificate.py)

The memory flip reverses the contrast, and \(J(c,d)-J(c,-d)=d\). Under a local phase convention, its transition consequently acquires \(e^{\pm i d(\theta_a-\theta_b)}\). With sources, transfers and effects consistently rephased, the intermediate factors cancel on every fixed outcome history. The same history has exactly the same Born weight including the final effect. Holding the free Hamiltonian fixed gives a different test: the nonzero blank transition \(10\leftrightarrow01\) forces equal active phase angles on the two ends of each edge. A connected board permits only their common angle for this charge-linear family. The measured free-flight difference exhibits this distinction directly.

## How the gauge paper's code is implemented

1. **Declare the magnetic bands and gains.** Use the oriented uniform fluxes \(k=2,7\) from the gauge paper, with per-face character \(e^{2\pi i k/20}\). Their respective ground bands have ranks 3 and 2. At the ports the primitive load charges, relative to the background, are \(+1,-1\). The programmable observable \(f(n)=-n^2/2+7n/2\) instead has gains \(f(2)-f(1)=2\) and \(f(0)-f(1)=-3\). The two physical reference terms can set any pair of signed relative phases. Neither flux number nor assigned gain is learned by the compiler.
2. **Synthesize the magnetic transfers.** Conjugate each isolated native contact by the pair of phases specified on *that edge*. The sum of the resulting actual local generators equals both target Hamiltonians on the two signed sectors. A symmetric unit-time product at 8, 16, 32 subdivisions has operator errors 1.0283, 0.3368, 0.0896 or smaller in both sectors. Finite-dimensional product convergence gives arbitrarily small error under finer protected schedules. A single relabelling of the twelve vertex phases instead telescopes to unit holonomy on every triangle; the unmodified zero-flux leading band has rank one.
3. **Build the material band detector.** A common frame maps the admitted modes into three possible pointer addresses. One additional initially unit rail accommodates the different ranks. Two contact dwells with an intervening sign-selective parity cancel exchange for one sign and rotate the other. Native relative phases and these rotations implement QR elimination on the two sectors within the same pulse sequence. A neutral pointer toggles at addresses 1, 3 and 5. Its single binary read, followed by the inverse frame, preserves superpositions within the admitted set. The complete enrolled logical frame contains 1,138 operations. The code expands it to 11,392 timed native waits and checks all 240 possible position rows on both signs. The material flag costs another 105 native waits. Both physical branch operators on all input columns agree with \(P_k\) and \(I-P_k\) to operator error below \(5\times10^{-13}\). All phase scalars removed in the reduction are common to both signs. The declared controls include the sharp pointer read, protected timing, prepared unit guard and neutral pointer, and preparation and reset between reads. [Frozen code test](CODE_PROTOCOL.md) · [Deposited frame](frame.json) · [Native pulse compiler, expansion and branch checks](code_frame.py)
4. **Measure the compressed kicks.** A raw port pulse can leave the band. For \(h=(P_k)_{pp}=r/12\), the normalized bright vector \(P_k|p\rangle/\sqrt h\) loses \(4h(1-h)\sin^2(q\theta/2)\) probability. A sign kick on a block loses 3/4 or 5/9 on its bright state. With repeated physical admission the surviving branch on the code converges to \(e^{-i\theta qP_k |p\rangle\langle p|P_k}\). This follows by the first-order expansion of one finite-dimensional step, including an \(O(1/n^2)\) remainder. For \(\theta=\pi/3\), 64 steps leave worst branch-operator error 0.01066 and bright-state survival at least 0.9788. The twelve Hermitian compressed matrices have total trace zero and their numerical Gram matrix has smallest eigenvalue 0.10. This reproduces the dimension of the paper's twelve generators.

## What independent receivers can distinguish

Prepare one sharp signed load at port 0, 3, 5 or 8. Use the same material filter for the initial admission and after each of 16 pieces on each of two ports, with both orders and two held angles. A full native address-pointer write and a separate signed head at the spare rail measure the final sharp address. Each bank retains all 33 possible first rejection stages, twelve final addresses, the spare address and the malformed pointer bin. Every reported probability uses the original enrollment. The receiver and the pointer use the declared Interval Echo outcome law. [Prospective readout protocol](ON_CODE_READOUT.md) · [All 32 banks](code_records.json) · [Acquisition and alternative predictions](code_readout.py)

| Prediction | Maximum absolute raw frequency error |
| --- | ---: |
| Published code with configured gains | 0.00212 |
| Primitive unit load gains | 0.01667 |
| Exchanged signed gains | 0.05090 |
| Zero gain | 0.02834 |
| Reversed port order | 0.01535 |
| Fine address measurement at the filter | 0.03556 |

The simultaneous 99% radius over 1,504 frequencies is 0.00491 for 8,388,608 fresh starts. The gain and code were fixed before this readout; a held source and detector can discriminate their physical action. To identify these gauge currents in freely propagating Interval Echo material, observers would still have to identify the relevant fluxes, the charge responses and a common preserved code through independent source and transport records. This controlled experiment realizes the gauge construction on a declared apparatus domain.

From the RER root, use one numerical thread to reproduce all results:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python -m model.experiments.interval_echo_gauge.port_phases.measure
.venv/bin/python -m model.experiments.interval_echo_gauge.port_phases.code_readout
.venv/bin/python -m model.experiments.interval_echo_gauge.port_phases.certificate
.venv/bin/python -m pytest -q model/experiments/interval_echo_gauge/port_phases
(cd Lean && lake env lean -j1 SourceLaw/IntervalEchoHeadCharge.lean)
~~~
