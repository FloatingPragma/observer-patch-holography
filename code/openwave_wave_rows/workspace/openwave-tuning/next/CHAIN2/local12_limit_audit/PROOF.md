# Independent LOCAL12 clock and limit audit

This package accepts the stated bounded scope of the frozen clock and continuum results. It adds a finite-member native compiler with a precision parameter. The addition proves existence of a native word family with controlled approximation error. It does not identify the completed-word clock with physical seconds, implement every bandlimited preparation, or derive an interacting or massive continuum field.

The frozen source is `native-local12-current-v0-proposed`, with W = SR and the original twelve event labels. The two reviewed manifest hashes are `6f80b0fe899cdeffd59d3742c318fe627817c19685d8f164abed69721f0b0992` for clock and `3067f7bdea43ca5137d6b115f7e1bd92b4d753ba99684ed198ba2da170240adf` for continuum. `audit.py` verifies their transitive source files before and after its checks. No frozen producer file is changed. No reserved configuration is run.

## What the walk limit contains

The independent reconstruction starts from the five internal matchings and the external links of the actual cubic-routed graph. It builds all 144 internal/event coordinates, the reflection R and transport S. At zero momentum it constructs every spectral projector of C = (W + W*)/2 by rational polynomial interpolation. The ranks are 26 at C = 1, 96 at C = -1, 6 at C = 2/3, 4 at C = 1/2 and 12 at C = 1/3. Their sum is the identity. Thus the complementary phase-pi sector and all optical sectors are retained.

The compression A*SA = I - H_graph/12 gives the spectral mapping cos(theta_j) = 1 - lambda_j/12. Each interior eigenvalue of the twelve-dimensional compression contributes two walk eigenphases. At small nonzero momentum there are two acoustic modes, 22 other inherited optical modes, 24 exactly stationary walk modes and 96 phase-pi modes. At zero momentum the acoustic pair collapses into the rank-26 phase-zero sector. None of these counts is a count of particle species.

Let P project on that rank-26 sector. With W'_a(0) = i D_a R, the independent integer calculation verifies H_a = P D_a R P = g_a e* + e g_a*, where e has 144 entries 1/12. It also verifies e*g_a = 0 and g_a*g_b = delta_ab/48. There are 22 vectors annihilated by every H_a. For nonzero p, the two vector amplitudes transverse to p add two zero eigenvalues, leaving generator eigenvalues +/- i|p|/sqrt(48) and zero with multiplicity 24. The terminal native displacement has P J_a P = -H_a, with the minus sign fixed by arrival-port records.

For bandwidth |p|_infinity <= K, spacing hK <= 1/2 and t = n tau0 h, the common theorem bounds full-state trace distance and every complete terminal-instrument total variation by

    min(1, (5/2) h K + (77/32) n h^2 K^2).

The corrector is X = -[(W(0)-I)|_Q]^-1 Q W_1 P. The complementary spectral gap is sqrt(2/3), so ||X|| <= (5/4)K. The single-step Taylor estimates give 5/4 + 1/2 + 1/32 + 5/16 + (5/8)hK <= 77/32. Telescoping and the two endpoint correctors give the displayed bound. The supremum over momentum also bounds the direct-integral operator tensored with any spectator reference. Trace distance is bounded by the vector bound; it is not equal to it. The independent small matrix checks support this algebra and are not used as the certified upper bound.

The uniform internal preparation retains both acoustic signs. Its optical-weight bound is r^2/(4-r-r^2/2)^2 for r = hK. The inherited scalar return probe is a cosine of the two branches, not a selected one-way eigenbranch. The phase-pi and optical amplitudes in the actual output are not deleted by the full-state theorem. A continuum density integrated over detector cells needs the additional 3hK/4 binning error stated by the frozen theorem.

## Preparation and interventions

The frozen `initial_cell_state()` prepares one cell with uniform internal/event amplitudes. As a family on a fixed torus, this cell-localized state is not a fixed-bandwidth preparation. A general bandlimited state, possibly reference-entangled, remains an explicit theorem-domain input whose native preparation and error must be supplied. This audit does not turn a low-mode vector installed in a simulator into a native preparation procedure.

A terminal instrument can be arbitrary because trace distance contracts under a complete channel. The same statement does not admit its output as a fresh initial state for another continuum interval. For the zero-momentum uniform vector e, a complete internal/event label measurement yields I_144/144. Its weight outside P is exactly 1 - 26/144 = 59/72. `audit.py` also computes a strictly positive leakage for the binary native positive-axis displacement measurement. A sharp cell measurement, with internal coordinates unread, can preserve P while destroying the fixed bandwidth. These are different failures of the next interval's hypotheses.

The signed displacement identity concerns its mean, not identification of a vector amplitude with an electric field. A signed variable in [-1,1] changes in expectation by at most twice its total variation error. At t = n tau0 h, physical displacement h divided by tau0 h gives the corresponding probability-flux coefficient. This proves neither a Maxwell radiation sector nor a charge-current law of the augmented field source.

A common scalar sign also illustrates why norms and record distances must be separated. The vectors psi and -psi have norm distance two and trace distance zero. Controlled application makes that sign relative to the bypass branch, distinguishing plus and minus probe states with trace distance one. A claim about that experiment therefore requires the actual controlled source word and its apparatus certificate. An unimplemented controlled-W cannot be inferred from a global-sign operator diagnostic.

## Clock and retained records

The clock pulse is exp(-i h Z/16). Its relative Bloch phase is h/8 per completed W, giving frequency 1/(8 tau0). The exact clock time is t = n tau0 h. The frozen mean phase bias is a bound on systematic law error. It is not the resolution of one native bit or of a finite sample. With N independent complete trials in each quadrature, the inherited Hoeffding radius is sqrt(2 log(4/alpha)/N) per sign contrast, plus twice the systematic total variation. Failure receives sign zero and remains in the denominator. Zero contrast and all-failure counts are inconclusive. The executable checks these cases.

The thinning law is classical: a uniform global seam contacts the load with probability 6/M, while all nonincident seam, coin and law records remain. The capped geometric distribution includes its timeout and sums to one. At a fixed per-seam rate r, the local contact rate is 6r at each development size. Coherent repeated W is a different continuation: its two-step return is 25/72 versus the classical 7/24, a difference 1/18. No quantum-channel equality follows from the thinning identity.

The frozen source has 432,425 native pulses and 864,850 load/export services per base W. A serial duration requires a control-rate cap, service latency and attempt latency. The phase clock does not select those values. A retained-current source W_F needs its own entire word, register allocation, reset budget and controlled wrapper. Clocking or controlling bare W alone does not establish the corresponding W_F experiment. The family adapter's fresh clock wire is outside all its source, record and scratch wires.

## Constructive native family adapter

`refinement_adapter.py` supplies `reversible_program(L)`, `certificate(L, step_tolerance)`, `pulse_stream(L, precision_bits)`, `execute_full_response(event)`, `clock_tick_pulse(...)` and `reset_repetition_plan(...)`. The frozen compiler itself is size-two and casts pulse angles and integrated controls to binary64. It supplies neither an arbitrary-size generator nor an asymptotic precision guarantee. The additive adapter deliberately introduces a new exact-rational control serialization. It uses the same twelve native response generators and proper-code load/export services, not another Hamiltonian source.

For L >= 2, encode a vertex by (cell << 4) | port, with port below 12, and the event by 2 direction + sigma below 12. There are d = 8 + ceil(log2 L^3) data bits. Unused cell, port and event addresses remain explicit. The word has three work-source bits, three raw-record bits, one law bit and d-3 scratch bits, giving 2d+4 memory bits. The ideal reflection acts on all source addresses; event-padding entries receive its exact minus sign. Transport fixes every padded data address.

The adapter enumerates the actual graph's oriented seams. For every endpoint and direction, a reversible equality test toggles the corresponding one-load work-source bit. All zero-controls are implemented by paired flips. It then uses the existing native event and record-recoder permutations, with the original placement-coin conversion. The transport permutation swaps the two moving endpoint addresses on every seam. Gray-path transpositions compile that permutation to the existing complete Toffoli ladders. Reversing the extraction word clears the work source at the new position. The inverse record recoder is exposed separately for native raw-label interventions.

On admitted clean work inputs, the extracted source has code two at a left endpoint and code one at a right endpoint. The native event plus recoder exchanges that occupation precisely when sigma is one and clears the raw record. Routing changes the position by the same seam. The final extraction therefore clears the work register. On data padding extraction is zero and native vacuum is unchanged. This proves the complete ideal S on all data/law inputs. Each generated primitive is a reversible gate on every dirty scratch, raw and source address, so those inputs have a defined unitary extension too. No clean-work equivalence is claimed for initially dirty work.

The development run executes every clean data/law column at L = 2, selected data/law columns at L = 3 and 4, and dirty-input forward/inverse controls. It does not enumerate exponentially many dirty columns or any larger carrier. The all-input family statement follows from the equality tests, Gray-transposition and ladder proofs, not from the selected-column check.

The reflection is prepared by a two-bit three-component state and two Hadamards. Three ordered two-plane rotations have angles -atan(sqrt(2)), -pi/4 and -pi/2. Their exact Pauli expansions retain all scalar phases. Conjugating diag(1,-1,...,-1) by that preparation gives R on all sixteen event addresses. The preparation, reflection and reusable primitive words are evaluated on complete small proper codes using emitted twelve-port controls, with all complement Kraus branches and the first-failure Gram retained.

## Precision, complete errors and resources

The precision parameter selects rational outward intervals for sqrt(2), sqrt(5), pi and atan(sqrt(2)). Square roots use integer square root. Pi uses Machin's identity with alternating arctangent series. The nontrivial preparation angle uses pi/4 + atan((sqrt(2)-1)/(sqrt(2)+1)). Every series tail and each arithmetic rounding is bounded outward. Interval precision can grow without a fixed mantissa ceiling.

For each native pulse, the exact source control identities give twelve ideal integrated coefficients. The adapter multiplies their intervals by the desired angle interval and emits their exact rational midpoints. If D_p is the full six-level anti-Hermitian source response, its independently bounded row norm gives

    ||G_emitted - G_ideal|| <= sum_p radius_p ||D_p||.

Both generators are anti-Hermitian, so Duhamel bounds their exponential difference by the same number without an exponential growth factor. The existing complete load/export dilation preserves this norm bound. Every outside-code Kraus output is retained. Telescoping therefore bounds the complete word Stinespring error, half diamond distance and terminal total variation by the sum of pulse errors. The certificate is valid for clean-work inputs entangled with arbitrary references. Numerical six-level residuals are diagnostic only.

For a member's finite pulse count, `certificate` increases precision until the exact rational sum is below the requested tolerance. The intervals converge, so this search terminates. Choosing delta_W(h) <= h^2 gives n delta_W <= (t/tau0) h. The common clock tick has the same scalable rational control API. Preparation, final detector and timeout errors must also vanish; they are separate terms rather than silent deductions from the source-word bound.

The pinned native blank-reset certificate is uniform over all 78 input states. If one reset block has failure bound q < 1, r sequential blocks fail with probability at most q^r, even after conditioning on previous failed blocks. A failed block leaves a native state to which the same bound applies. For Q registers, the union bound is Q q^r. `reset_repetition_plan` chooses r for a supplied total budget. It stops each register after its first heralded success, retains every fine path and final timeout, and never normalizes away failures. Other apparatus resets, if present, need their corresponding uniform certificates.

There are 36L^3 disjoint moving pairs and 72L^3 extraction endpoints. An equality-test ladder uses at most 4d-7 reversible gates; a Gray transposition uses at most (2d-1)(4d-7). Two extraction words and transport therefore have an explicit O(L^3 log^2 L) gate bound, plus the fixed local native event/recoder words. Each reversible gate has a fixed finite native template. `resource_upper_bound` evaluates this symbolic bound without constructing a graph. This is a serial address-encoded implementation. Intermediate processor states can be far apart in graph coordinates; only completed words have the original graph light cone.

At fixed serial bandwidth and positive service latency, this growing pulse cost cannot fit into tau0 h as h decreases. The adapter closes a formal finite-word and approximation-error existence gap. Physical-second convergence, distributed graph-local services, hardware noise, general native bandlimited preparation and the refined detector remain explicit obligations. The new tape must not be passed through the old float-only API and then claimed to retain its precision certificate.

## Field and matter-wave boundary

The audited free result gives omega^2 = |p|^2/(48 tau0^2). After normalizing the shared propagation speed, it is the m = 0 instance of the literal row-30 expression “emergent ω²=k²+m² for matter waves.” It does not derive a nonzero matter mass, observed species masses or a common interacting field dispersion. Pending `common_dynamics/particles` mass-coin and coupled-word packages are separate proposed extensions. Their presence does not transfer their claims into these frozen manifests or into this audit.

The reviewed finite flux algebra gives [E,U] = U - m|0><m-1|. Modular Gauss preservation is not the integer U(1) commutator. A bounded integer-flux update with boundary holds is an explicit source amendment whose operator-norm distance from the cyclic shift can be two. A uniform finite affine-cycle state has no normalizable infinite uniform integer counterpart. Field preparation, energy and interacting particle transfers belong to their owning source-specific audits. This package changes no official OpenWave row.
