# Reusing one interacting radial field for complete native readout

## Source and controlled interval

Use the complete 27-site radial Hamiltonian, its normalized actual ground
projector \(G/2\), or either unique lowest charged state from
[`radial_interacting_core`](../radial_interacting_core/PROOF.md), with
spacing and velocity one, radial mass 64, and centered contact coupling
\(U=1/65536\). The entire 54-mode Fock space is retained. At the origin put
\(D=(n_\uparrow+n_\downarrow-1)^2\), \(P_1=D\), \(P_0=I-D\).
The same-field native instrument of
[`radial_geometry_join`](../radial_geometry_join/PROOF.md) is
\[
 K_0=(7P_0+P_1)/\sqrt{50},\qquad
 K_1=i(P_0+7P_1)/\sqrt{50}.                                  \tag{1}
\]
Each round has fresh blank probe registers and nine independent fair coins,
but it acts on the **same field**. During the controlled interval the field
undergoes only (1): Hamiltonian evolution is absent or compensated, and no
field repair or additional field operation is interleaved. The source
preparation, addressable controls and
drift compensation are supplied inputs. The statement is a composition
theorem for that declared instrument, not a claim that \([H,D]=0\).

The native work register is the full 78-state, total-two load carrier on
twelve ports. It starts independently of the field and its untouched
reference, with arbitrary classical law. Its initial mean is \(x_0\).
Every pointer result, scheduler word, placement coin, seam or idle label,
current and terminal native load is retained. Quantum probe outputs are
retained or committed as in the parent instrument; they are not silently
reset. The fixed detector schedule depends only on these declared records.

## Complete history instrument without initial dephasing

For an ordered binary history \(s=(b_1,\ldots,b_N)\), let \(k=\sum b_j\).
Since the two operators in (1) commute,
\[
 K_s=K_{b_N}\cdots K_{b_1}
 ={i^k\over50^{N/2}}(7^{N-k}P_0+7^kP_1).                    \tag{2}
\]
This equality acts on the complete field and extends with any untouched
reference. For an arbitrary joint input \(\rho\), its unnormalized
conditional state is
\[
 K_s\rho K_s^*=50^{-N}\bigl[
 49^{N-k}P_0\rho P_0+49^kP_1\rho P_1
 +7^N(P_0\rho P_1+P_1\rho P_0)\bigr].                       \tag{3}
\]
The cross terms are present. No parity measurement or dephasing is inserted
at the input. To include a complete native history, multiply (3) by its
independent coin probability \(512^{-N}\) and retain the native deterministic
trajectory from its initial load; histories incompatible with that trajectory
have zero weight. This describes the joint field, reference and record
instrument.

Put \(q=\operatorname{Tr}(D\rho)\). Taking the trace of (3) gives
\[
 \Pr(s)=(1-q){49^{N-k}\over50^N}
              +q{49^k\over50^N}.                            \tag{4}
\]
Thus record probabilities are a mixture of two independent Bernoulli laws,
with a single latent parity fixed through the block. This probabilistic
representation does not replace the coherent conditional state (3).
Summing all \(2^N\) ordered histories proves trace preservation and
\[
 \Phi^N(\rho)=P_0\rho P_0+P_1\rho P_1
 +(7/25)^N(P_0\rho P_1+P_1\rho P_0).                        \tag{5}
\]
The integer count representation in the certificate is an exact calculation
of (2)--(5); physical ordered records are not discarded or identified.

## Native mean and measurement work

For fixed latent parity \(d\), each round uses an independent pointer with
\(\Pr(b=1\mid d)=1/50+24d/25\), and fresh independent scheduler and
placement coins. The source's exact one-step conditional native mean is
\(M_d x\), where
\[
 M_d=I-\{L+(3/25)(1-2d)(D_{(0,1)}-D_{(0,5)})\}/64.
\]
The parent proves this identity on all 78 native inputs, including idle and
zero-current histories. Conditioning on latent parity and iterating that
linear identity gives the repeated-field mean
\[
 \mathbb E x_N=(1-q)M_0^N x_0+qM_1^N x_0.                  \tag{6}
\]
This is generally different from \(M_q^Nx_0\). For two rounds the
matrix difference is
\[
 q(1-q)(M_1-M_0)^2,                                        \tag{7}
\]
which is nonzero for this source. The checker independently enumerates all
79,872 one-round native histories, verifies the endpoint mean identities,
and propagates the complete 78-state transition kernels for two rounds.
There is no fresh field preparation hidden in (6).

Every local mass or centered contact term commutes with \(D\). Precisely
the hopping terms incident at the origin, denoted \(H_\star\), anticommute
with \(I-2D\). Equation (5), on the entire Fock space, therefore gives
\[
 \Phi^{N*}(H)=H-[1-(7/25)^N]H_\star,
 \quad W_N=-[1-(7/25)^N]\operatorname{Tr}(H_\star\rho).       \tag{8}
\]
This is the energy change caused by the controlled measurements. The work
of setup, compensation and reset is not assigned by (8). The postmeasurement
field need not lie in its initial spectral sector; (8) includes its virtual
and high-energy components. The work is not labeled reservoir heat.

## One-preparation ground versus charged classifier

The actual source states obey
\[
 q_G\le e,\qquad q_{C_-},q_{C_+}\ge1-e,\qquad
 e={188263256928\over480690783803197}<1/1000.                 \tag{9}
\]
These are the parent full-Fock interacting bounds, rather than projected
core probabilities. Use \(N=32\). Return `charged` when at least 16 pointer
bits are one, and `ground` otherwise; the tie is assigned to `charged`.
For latent parity zero, the error is \(\Pr\{\mathrm{Bin}(32,1/50)\ge16\}\).
For latent parity one it is \(\Pr\{\mathrm{Bin}(32,49/50)\le15\}\).
Hoeffding's one-sided bound, using a weak deviation event of size \(12/25\)
for both tails, gives
\[
 p_{\rm tail}\le\exp[-2\cdot32(12/25)^2]
 =\exp(-9216/625)<10^{-6}.                                 \tag{10}
\]
The final strict inequality follows by summing the first 21 nonnegative
Taylor terms of \(\exp(9216/625)\), whose rational sum exceeds one million.
The checker also computes both binomial tails exactly, independently of
this exponential estimate. Equations (4) and (9) imply error at most
\(e+p_{\rm tail}\) on each allowed actual spectral input.

The parent certifies one-half trace distance below \(10^{-9}\) between
each complete native round and its ideal instrument, for every untouched
reference. This includes all quantum-control and nine-coin failure outputs.
Composing the instruments with fresh apparatus and retaining the whole
archive gives error below \(32\cdot10^{-9}\), by channel contractivity and
the telescoping sum. It does not require conditional input states to be
independent or close to their initial spectral state.

For physical control failure, retain its round, pulse, complement label and
record prefix, halt the interaction, and assign zero to the failed and unused
pointer slots. Every nonlogical pointer outcome receives score zero, is
retained, and also halts the block. The ideal instrument has zero probability
for those labels; the
same complete-channel bound controls any classifier with this total scoring
rule. Hence the full classifier satisfies
\[
 p_{\rm error}<e+10^{-6}+32\cdot10^{-9}<1/2500.               \tag{11}
\]
No success branch is postselected. The result classifies the promised input
as ground or charged without resolving the sign of charge. It measures one
field's persistent parity. Repetition cannot consistently estimate the
ensemble parameter \(q\) from one field: as \(N\) grows the sample frequency
approaches \(1/50\) or \(49/50\), with weights \(1-q,q\).

## Finite resources and boundary

The block uses one field preparation, 32 rounds, 40 parity-instrument pulses
and nine four-pulse native fair-coin words per round. Its pulse budget is
2,432, with at most 32 integer native repair attempts. Without recycling
apparatus, each round has five fresh probe work
bits and nine fresh two-bit coin registers, totaling 736 fresh logical bits,
in addition to the 54 field modes and the 78-state classical native work
register. The archive has 32 finite frames with the complete fields specified
above. Native complement records have the finite alphabets of the 76-pulse
round, and the halted suffix has a fixed label. This gives a finite explicit
instrument without executing an exponentially large history table. The
receipt executes independent short-history and complete-source controls;
it does not claim a physical 32-round run or preparation of the spectral state.

This removes independent field re-preparation within the specified controlled
readout block. It supplies neither an autonomous readout under simultaneous
Hamiltonian/repair evolution, a reconstructed Gram matrix from one preparation,
a thermalization theorem, a source-selected radial texture, a continuum
interacting field theory, nor Einstein dynamics. Fresh apparatus, source
preparation and the stated control interval remain explicit inputs.
