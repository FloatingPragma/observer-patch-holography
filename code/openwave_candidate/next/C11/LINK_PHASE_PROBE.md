# External link-phase probe on the one interacting source

## Reuse and remaining adapter

| Needed | Existing result/API | Action and remaining difference |
| --- | --- | --- |
| Actual common H/reference/native instrument | C10 `JointNativeSource`, `ModularField.halfhops`, complete C/E/MIX | Reuse directly, full384-state interacting all-number source; no bare C9 gap. |
| Regional continuity current | `oph_openwave.native_currents.regional_current` | Reuse directly on each actual orbital-occupation projector, then sum for a region. |
| Current as derivative and gauge-paired preparation | RER `openwave_magnetic_response.response.current_and_gauge_controls` | Reuse its derivative/gauge-control method with a source-specific Peierls adapter; its spin-curvature/g≈2 result belongs to another source and is not transferred. |
| External native transport phase | C9 all-character Kraus phase and CAR transport convention | Reuse with external real phases in place of the character labels; this is an added apparatus intervention, not a dynamical internal field identification. |

For the actual oriented port→transit signed half-hop T_h, define `H(A)=H(0)+sum[(exp(−iA_h)−1)T_h+h.c.]`, keeping the same periodic electric and cyclic face terms. The intervention current is `J_h=−dH/dA_h=i(exp(−iA_h)T_h−exp(+iA_h)T_h†)`. With this sign convention, the regional charge derivative is `sum_h(1_R(port)−1_R(transit))J_h`. The actual matrix derivative, regional commutator and expectation derivative are checked independently in code.

The single development probe puts phase .17 on half-link4, all others zero. It begins in the actual C10 global reference, quenches to H(A), evolves for .37 native hop units using sparse exponential action, and applies the canonical probability .29 native event. MIX is selected; C and E are separate controls. All outcome probabilities, applied/proposed integer currents, half-record increments, conditional input and final regional charge weights, and final energies are retained. Their complete sums give mean native charge change and apparatus work. Conditional branch charge differences equal the applied current across that region's boundary. No initial-energy PVM is inserted or claimed in this work ledger.

For a local orbital phase G=exp(i sum lambda_z n_z), `A'_h=A_h+lambda_port−lambda_transit` gives `H(A')=G H(A) G†`. The state/reference must transform by G as well. The complete moving native branch is dressed by `exp(−ij(A_ea−A_eb))`; its transformed matrix is exactly GKG†. Every C/E/MIX branch probability, charge and energy, plus the full poststate, is gauge covariant. Leaving the preparation untransformed is a rejecting control. The dressing is a common scalar phase within each classical record outcome, so it cancels in that branch's density map; an undressed matrix can differ while representing the same classical-record instrument. No false observable channel difference is attributed to that phase convention.

The static perturbed ground is recomputed across **all** source numbers and remains in N3. Its derivative energy and persistent-current expectation obey the Hellmann–Feynman relation. This is a separately declared ideal static response, not relaxation of the quenched initial state. The initial phase-quench work, subsequent coherent current and native monitoring work are distinct quantities.

The local continuous gauge equivalence here concerns an **external** Peierls apparatus. It does not restore the invalid integer `[F,T]=T` identity at internal Z6 wrap or establish a physical U(1) field. The gauge-invariant loop phase is an imposed probe, distinct from the internal cyclic face operator W. No physical magnetic field, sourced dipole, spin, g=2, SI units or electron identification follows. No reservoir pulse is used in this experiment. The frozen C10 sources are imported unchanged; all new files stay under C11.
