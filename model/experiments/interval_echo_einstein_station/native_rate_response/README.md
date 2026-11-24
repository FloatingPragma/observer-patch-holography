# Can a seam memory slow a material clock?

An isolated native seam transfers a double load between two ports. A memory can reduce the transfer probability. A second preparation distinguishes an intermittently open contact from a clock running at a single reduced rate. At one half of the blank-memory weight, the two terminal probabilities are 0.25 and 0.75. Slowing the complete blank seam to half speed gives 0.146447 and 0.853553 at the same appointment.

## Two preparations and one detector

Write \(L=|2,0\rangle\), \(R=|0,2\rangle\), and let \(p\) be the initial probability of a blank memory. Launch either \(L\) or \((L+iR)/\sqrt2\). The latter is the output of a blank reference seam acting on \(L\) for \(\pi/2\) in exchange units, up to a global phase. Either load preparation can be transferred to the test ports by complete blank-seam swaps. At the appointment, read both port loads and count the outcome \(R\). Each launch uses a fresh copy of the same memory preparation. Exposures are specified by the apparatus.

On the isolated contact the native integer matrices give \(\langle R|h|L\rangle=-P_{\rm blank}/2\), the two diagonal memory blocks agree, and there is zero leakage to other load words. The blank projector commutes with both blocks. Hence the respective yes probabilities are

\[
 P_L(t)=p\sin^2(t/2),\qquad P_Q(t)=\tfrac12(1+p\sin t).
\]

Tracing out memory gives these expressions for every memory density matrix, including a memory entangled with an unread spectator. The changed law scales both terms of the blank contact by \(r\) on the same alphabet and route. Its probabilities are \(\sin^2(rt/2)\) and \((1+\sin rt)/2\). The phase preparation measures \(2P'_Q(0)=p\), requiring \(r=p\). Transfer from the sharp preparation measures \(2P''_L(0)=p\), requiring \(r^2=p\). Both conditions hold precisely at zero and full rate. The discrepancy \(p(1-p)\) is first order in the departure from a completely blank or completely written seam. The full 99-state propagation checks five blank weights, three appointments and a coherent memory with all ten written amplitudes present. The reference construction gives the [integer block identities and written-sector flight law](../native_tidal_response/ROUTING.md).

**The derivative test also applies on a simple native graph.** Prepare either endpoint word as specified above with any common environment density matrix \(\rho\), including arbitrary other loads and correlations among all memories. Every site term and every other edge has zero matrix element from \(L\) to \(R\), since each acts on at most one of the two marked ports. Thus \(H_{RL}=-P_{\rm blank}/2\) tensored with the identity on the remaining environment. The sharp preparation gives \(P''_L(0)=2\operatorname{Tr}(H_{RL}\rho H_{LR})=p/2\). For the phase preparation the diagonal block has real expectation and drops out of the first derivative, giving \(P'_Q(0)=p/2\). Tests apply both overlapping native Hamiltonians on a three-port chain with entangled memories and an occupied spectator. At finite times the two isolated-seam formulas need the other contacts closed. The derivative conditions concern sharp microscopic load probes at release. They place no condition of this form on long-wavelength modes in a collective state.

## A probe that reads correlation energy

For comparison, append a two-level probe prepared along \(+x\) and read along \(+y\). In common exchange units its Hamiltonian together with the seam is

\[
 H=I\otimes h_0+\tfrac{\Omega}{2}\sigma_x\otimes I+g\sigma_z\otimes h_s,
 \qquad h_0=X+8C,\quad \Omega=1,\quad g=1/8.
\]

Here \(X\) and \(C\) are the seam exchange and squared contrast. The three tested couplings are \(h_s=h_0,8C,0\). The tuned contrast control and the whole-seam probe give the same records on the vacuum and the symmetric one-load state. For an eigenvalue \(e\) of \(h_s\) in a joint seam eigenstate, the yes probability is

\[
 \tfrac12+\frac{ge}{\sqrt{\Omega^2+4g^2e^2}}
 \sin\!\big(t\sqrt{\Omega^2+4g^2e^2}\big).
\]

The [four prepared protected words](../../../../code/matched_einstein_limit/seam_overlap/README.md) have energy table \(((0,0),(33,32))\). The high ensemble weights them by \((1,3,3,1)/8\), the low ensemble by \((3,1,1,3)/8\). Each logical bit has the same marginal in the two arms. Their total energy differs by \(1/4\), all in exchange. The whole-energy prediction for the differential yes slope at release is \(g/4=1/32\). Direct evolution on all 198 seam and probe states gives the following high-minus-low probabilities. The contrast-only and idle slopes and probabilities vanish exactly by equal coupling-eigenvalue histograms.

| Exposure | Whole energy | Contrast only | Idle |
| ---: | ---: | ---: | ---: |
| 1/64 | 0.000484332 | 0 | 0 |
| 1/8 | 0.002051651 | 0 | 0 |
| 1/2 | −0.008974061 | 0 | 0 |

The sign change comes from coherent precession at finite exposure. Records here are Born probabilities for the specified detector, without shot noise. The four-word preparations and the probe Hamiltonian are supplied. A geometrical identification would require an independently recovered proper metric and a signed area response at the ordinary-matter coupling.

From the repository root, run

```sh
python -m model.experiments.interval_echo_einstein_station.native_rate_response.measure
python -m model.experiments.interval_echo_einstein_station.native_rate_response.energy_probe
python -m pytest -q model/experiments/interval_echo_einstein_station/native_rate_response
```
