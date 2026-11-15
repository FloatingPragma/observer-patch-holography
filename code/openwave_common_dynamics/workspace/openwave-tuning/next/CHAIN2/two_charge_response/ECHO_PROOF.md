# Connected phase under the unchanged full source

Choose the physical transport edge joining ports 0 and 7, edge index2.
It is outside the three-edge support of face0. At each endpoint take
the held event immediately preceding its moving event, and form an
equal superposition of these two modes. Both coherent event qubits
contain exactly one charge, so all four joint branches have total
charge2 and identical port loads. The complete initial field is
\(F=P_0+P_7\), which has value1 on this edge and zero elsewhere.
Every branch therefore has the same Gauss background \(2\delta_0\).

Apply the literal inverse local \(\Gamma(C_m)\) words at the two
endpoints as part of preparation, where
\(C_m=R\exp(-ihmR)\). The unchanged full source is
\(V_m=U_E(.17)U_B(.41)\Gamma(S)_F\Gamma(C_m)\), with all207
other face couplings zero. Its local coins cancel these preparation
inverses. All other local matter blocks are vacuum. Only the chosen
edge can receive particle current.

Let \(l,r\in\{0,1\}\) indicate the moving event at the left and
right endpoints. The edge current is \(l-r\), including zero for
the two-occupied-mode exchange. The shared edge electric phase is
\[
 z_{lr}=\exp[-i\alpha(1+l-r)^2/2].
\]
The face0 plaquette acts on different edge registers. The electric
energy is a sum over edges. Hence the remaining287-edge field evolves
as a common quantum factor, independent of \(l,r\). Starting from
zero on the face, its exact sixteen-orbit amplitudes are
\[
 \phi_q=\exp[-3i\alpha\bar q^2/2]
 [\exp(i\beta\cos\hat\theta)]_{q0},\qquad
 \bar q=(q+8)\bmod16-8.
\]
All sixteen amplitudes are retained, including those crossing the
centered modular branch. The construction makes no integer-field
cutoff substitution.

The terminal apparatus applies the literal inverse of the complete
native transport/current word. The chosen edge field returns to its
initial value and the matter occupations return to their original
ports. The canonical anticommutation sign from the forward exchange
is retained and cancelled by its inverse. All other transport factors
are identity on the zero-occupation modes. The rest-field factor is
unchanged. The resulting source state is
\[
 \frac12\sum_{l,r=0}^1 z_{lr}|l,r\rangle
 \otimes |F_{\rm chosen}=1\rangle\otimes|\phi\rangle.
\]
Its connected phase is
\[
 \frac{z_{00}z_{11}}{z_{10}z_{01}}=e^{i\alpha}.
\]
Any additive pair of one-charge phases has ratio1. Separate electric
registers also give ratio1. This cross phase follows from the shared
native current and electric source; no classical pair potential is
supplied.

Native occupation phases remove the one-body phases determined by
the two single-moving branches. Native local Givens rotations analyse
both event qubits. The complete event probabilities are
\[
 p_{00}=(10+6\cos\alpha)/16,\qquad
 p_{01}=p_{10}=p_{11}=(2-2\cos\alpha)/16.
\]
Their TV distance from the additive-phase control is
\(3(1-\cos .17)/8=0.005405712408914715\).
Each event probability is multiplied by every \(|\phi_q|^2\),
giving the complete64-outcome event/field law. Occupation, event,
field, and Gauss records are explicit in `FULL_V_ECHO.json`.

`edge_echo.py` executes the literal native current block, the full
fermionic swap with its occupied spectator parity, the frozen electric
word, the literal inverse transport, and the local preparation and
analysis. The native echo residual is about \(1.10\,10^{-13}\).
`full_v_echo.py` checks the outside-face support condition and the
common field factor through actual frozen native plaquette and
electric words. Its field residual is about \(6.93\,10^{-15}\).
The exterior-square source replay in `NATIVE_CONTRACT.json` separately
checks the full local coin and mass words, including vacuum, padding,
and retained global phases.

The full source step, inverse local coin preparation, inverse transport,
phase calibration, analysis, all3,203 source resets, all3,200 pointer
resets, and the complete3,200-bit native detector enter the rational
error bound. At either supplied mass its complete native TV bound is
about \(3.29647\,10^{-5}\). No ideal projector or free inverse is
inserted. The bound concerns the exact ideal law and arbitrary
original spectator references; displayed decimal laws use numerical
residual checks rather than formal interval enclosures.

This is a connected interaction record under one unchanged full
source step and an explicitly charged inverse/readout intervention.
The mass-dependent inverse-coin preparation makes the event law
independent of the supplied mass; it does not demonstrate mass
independence of arbitrary preparations. The echo is distinct from
a freely observed displacement or an unperturbed multistep force
trajectory. It gives no Coulomb law, physical mass, species identity,
reciprocal labelled acceleration, or continuum conclusion.
