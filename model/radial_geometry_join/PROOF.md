# One interacting defect source and its geometry readout

## Source and states

Use the full fermionic source on the 27-site cube `{-1,0,1}^3` from
`radial_interacting_core`: spacing and velocity one, bounded radial mass
64, centered contact coupling `U=1/65536`, and the complete `2^54` Fock
space. The texture, parameters and observation schedule are specified
inputs. The one-body gap lower bound is `d=62`. Its exactly two-dimensional
zero space has an origin-normalized frame with

\[
 Z_{00}=\alpha I_2,\qquad \alpha\geq a_0=7688/7691.
\]

The full interacting source has its actual half-filled Kramers ground
doublet and two charged singlets, with the gap bounds proved by that
package. Put `w=27U/4` and `t=w/(d-2w)`. Each charged projector is within
operator norm `t` of its quadratic reference. The normalized ground
mixture has half trace-norm distance at most `t` from the quadratic
single-zero-occupation mixture. These are statements on the full source,
including the virtual complement; no four-state truncation replaces it.

## Local energy separates these states

Let `n0=n_up+n_down` at the origin and define the even-occupation projector

\[
 D_0=(n_0-1)^2=\operatorname{diag}(1,0,0,1).
\]

The actual contact-energy term is exactly

\[
 W_0=U(n_\uparrow-\tfrac12)(n_\downarrow-\tfrac12)
     =\frac U2D_0-\frac U4 I.                         \tag{1}
\]

Inversion anticommutation at the origin gives the sea covariance
`(I-Z00)/2`. Wick's identity
`<D0>=1-tr(C00)+2det(C00)` gives `(1-alpha^2)/2` for a singly
occupied zero mode and `(1+alpha^2)/2` for zero occupation zero or two.
The single-occupation answer is independent of its spinor direction,
so the same formula holds for its ground mixture. Since `0<=D0<=I`,
the interacting projector bounds imply

\[
 0\leq n_G:=\operatorname{Tr}(\rho_GD_0)\leq e,
 \quad 1-e\leq n_\pm:=\operatorname{Tr}(\rho_\pm D_0)\leq1,
 \quad e=\frac{1-a_0^2}{2}+t<\frac1{1000}.             \tag{2}
\]

Here the numerical value of `e` is about `0.000391652`. In particular,
the local contact-energy contrast is at least `U(1-2e)/2>0`.
The exact particle-hole symmetry sends `n0` to `2-n0`, preserving `D0`.
It interchanges the two charged singlets, so `n_+=n_-` and their entire
geometry responses are equal. Equation (1) does not identify the local contrast with the total charged
excitation energy: the hopping and all other contact terms contribute
to that total. The preparation of these exact spectral states is an input,
not an executed ground-state preparation protocol.

## Complete controlled quantum instrument

Use two field occupation bits, four source probe bits and one pointer bit.
The field bits can be entangled with every other Fock mode and an arbitrary
untouched reference. A phase-exact X on the blank charge probe followed by
two CNOTs computes `D0=1 XOR n_up XOR n_down` into that probe. All these
operations are even functions of the field occupations, so the map needs
no omitted Jordan-Wigner string. The remaining charge and scalar spectator
bits are blank. Prepare the active scalar in `(|0>+i|1>)/sqrt(2)`.

Apply exactly the retained source interaction subword
`H_int=X_scalar Z_charge/4` for time `4 atan(3/4)`, copy the scalar bit
to the pointer and invert the parity-copy word. The charge probe returns
to its blank; both scalar and pointer outputs are retained. On the original
field the two complete logical Kraus operators are

\[
 K_0=(7I-6D_0)/\sqrt{50},\qquad
 K_1=i(I+6D_0)/\sqrt{50}.                            \tag{3}
\]

Thus `sum K_b^*K_b=I` and

\[
 p_1=\frac1{50}+\frac{24}{25}\langle D_0\rangle.       \tag{4}
\]

The identity holds on coherent and reference-correlated inputs. An
independent routing template verifies every pulse in the source subword.
Its exterior quarter-turns conjugate `Z_scalar` to `-X_scalar`, and its
central CNOT pair conjugates `Z_scalar` to `Z_charge Z_scalar`.
For its one nonrational coefficient, integer Taylor enclosures compare
the emitted sine and cosine to `3/5` and `4/5`. Their absolute errors
bound the operator error of the entire Pauli rotation. The fixed-angle
errors and the complete native-port bounds from `certified_reference`
and `centered_native_preparation` are added by isometry telescoping.

The forty emitted pulses execute on all 128 work rows with four retained
reference columns. The execution retains all 146 native first-failure
slots; neither survival conditioning nor renormalization is used. The
four reserved physical pointer outputs have zero ideal weight. Nine
fresh fair-coin words supply eight scheduler bits and one placement bit.
The sum of the quantum and coin isometry bounds is a rational number
less than `10^-9`. Consequently the complete single-attempt output,
including native repair, pointer, idle and failure records, has half
trace-norm distance less than `10^-9` from its ideal counterpart. Applying
the same exact integer-repair instrument to both outputs is contractive.
This assumes the admitted native-port specification, including its memory
and controller services; it does not bound unrelated hardware errors.

## Measurement work and continuing use

The channel of (3) fixes the diagonal parity blocks and multiplies their
off-diagonal blocks by `7/25`. Let `H_star` be the sum of the actual
Hamiltonian hopping terms incident to the origin. Each of these terms
changes origin parity. Every other hopping, onsite mass and contact term
commutes with `D0`. Hence, on the full Fock algebra,

\[
 \Phi^*(H)=H-\frac{18}{25}H_{\rm star},\qquad
 \Delta E=-\frac{18}{25}\langle H_{\rm star}\rangle.   \tag{5}
\]

The measurement preserves total particle number and the unconditional
local contact energy. It need not preserve the total Hamiltonian energy
or its ground projector. Equation (5) explicitly accounts for its field
energy injection. Native integer-load redistribution follows the recorded
controller; it is not asserted to compensate this energy. The controlled
pulses require their specified schedule and do not include simultaneous
uncompensated free-Hamiltonian drift.

Iteration with the same mean matrix uses independent copies of the
declared field input, initially independent of the integer-load work
system and its past records. Reusing a correlated field without this
condition does not give that iteration. Blank probe memory, coin services
and retention of the full transcript are explicit resources.

## Same-state geometry contrast

Feed (4) to the unchanged eight-bit scheduler and native placement law
of `matter_geometry_response`. Sum all 1024 complete histories for each
of its 78 integer-load inputs. With `n=<D0>`, their actual mean is

\[
 M_n=I-[L+\theta(D_{(0,1)}-D_{(0,5)})]/64,
 \qquad \theta=\frac3{25}(1-2n).                     \tag{6}
\]

The idle probability is `1/16` and every seam has positive probability.
The exact readout (4) supplies the parameter of this same source law;
it is not a fit to independently selected field and geometry energies.
Its rank-three spectral band and normalized Gram matrix `G(n)` are the
ones certified in the parent package. That proof gives
`G01(0)-G01(1)>3/40` and `|dGij/dtheta|<=80`, hence
`|dGij/dn|<=96/5`. Applying (2) to the actual ground and either actual
charged singlet gives

\[
 G_{01}(n_G)-G_{01}(n_\pm)
 >\frac3{40}-\frac{192}{5}e>\frac{59}{1000}.          \tag{7}
\]

The rational lower bound is approximately `0.05996058`. This proves a
finite energy/geometry response for the same charged source states.
The response Gram is a spectral observation geometry; no Lorentzian
metric, conserved relativistic stress or Einstein equation is identified.

For a concrete finite observation specification, take 48 million fresh
attempts for each of the ground and one charged preparation. Score
record one as one and every other record, including failure, as zero.
Hoeffding with `delta=1/4000` gives joint sampling failure at most
`4 exp(-6)<1/100`; the latter inequality follows from a finite rational
Taylor lower bound for `exp(6)`. Complete-control error contributes less
than `10^-9` to either probability. Clipping each estimate to the ideal
probability interval `[1/50,49/50]` cannot increase its error. Since
`|dG/dp1|<=20`, the inferred ideal Gram error is bounded by
`20(delta+10^-9)` per preparation. Equation (7) therefore retains a
contrast greater than `0.04899996` at joint confidence at least 99 percent.
These 96 million attempts are a costed specification, not an executed
experiment. This inference bound does not assert that failed or rounded
long-time dynamics has exactly the ideal matrix (6).
