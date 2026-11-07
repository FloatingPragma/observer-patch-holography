# State, clock and motion of a protected native charge

The binary native source has 66 occupation states. Its endpoint weights
are read from `native_program.interpreter.program`: two ports have weight
25 and ten have weight 8. Compose this declared source with the existing
ordered-cover memory in `defect_protection`. This specifies a 132-dimensional
process with one persistent orientation bit. The construction is an explicit
composition of two existing source laws; the unextended record-clock member
does not contain this bit.

This process has a complete classification of its stationary states, a
controlled preparation in either protected charge sector, and an exact
relation between its recorded motion and its modular energy changes.
Protection of the memory charge supplies no spatially localized particle or
nonzero rest energy.

## The complete biased ordered-cover instrument

For raw record `h=(law,seam,coin,raw)`, the existing ordered-cover map is

```math
L_h=\sum_{x\in D_h}|f_h(x)\rangle\langle x|
           \otimes X^{p_h(x)}.
```

Here `f_h` is the actual binary endpoint update and `p_h` is the parity of
the transported ordering of the two occupied ports. The map is a partial
injection. Each raw slot is retained, including zero and abort slots.
The biased complete Kraus map is `K_h=sqrt(w_h)L_h`, with

```math
w_h=\pi_{\rm seam}(e)\pi_{\rm law}(\ell)
       \frac{\theta_{a\text{ if coin}=1\text{ else }b}}{\theta_a+\theta_b}.
```

The named member uses uniform seam selection and equal tagged C/E mixing.
The source coin depends on the declared endpoint weights, independently
of the state. Replacing fair amplitudes with these positive amplitudes
changes no branch domain, ordered lift or generated star algebra.

`reference.py` independently checks each lifted column against the native
endpoint function and the raw-record encoding. It verifies complete
effects, all stationary populations and the signed moving columns. The
proof is for this exact mathematical instrument. The inherited fair-coin
source execution does not, by itself, certify a complete literal biased
ten-qubit implementation; no such execution is claimed here.

## All stationary states on the same carrier

Let

```math
\rho_*=Z^{-1}\sum_x\prod_i\theta_i^{x_i}|x\rangle\langle x|,
\qquad Z=7505,
\qquad Q=I_{66}\otimes X,
\qquad P_\pm=(I_2\pm X)/2.
```

The sum runs over binary occupations with total two. The existing endpoint
balance identity proves that `rho_*` is stationary. Every ordered branch
commutes with `Q`. On a diagonal source input, its ordering operation fixes
each `P_+` and `P_-`; hence

```math
\rho_a=\rho_*\otimes[aP_+ +(1-a)P_-],\qquad 0\leq a\leq1,
```

is stationary for the full biased channel. In particular
`rho_0half=rho_* tensor I2/2` is faithful.

There are no additional stationary density operators. The exact
`defect_protection.protection.algebra_certificate` constructs the branch
algebra `M66 tensor span{I,X}` from separating effects, a connected tree
of movements and an odd native exchange loop. Its commutant is `span{I,Q}`.
All branches used in that construction have positive weights here.

For completeness, let `T(A)=sum_h K_h^* A K_h` be the unital adjoint
channel. If a Hermitian `A` satisfies `T(A)=A`, then

```math
T(A^2)-A^2
 =\sum_h(AK_h-K_hA)^*(AK_h-K_hA)\geq0.
```

Its expectation in the faithful invariant `rho_0half` is zero. A positive
matrix with zero faithful-state expectation is zero; each summand therefore
vanishes. Thus `A` commutes with every `K_h` and its adjoint. Conversely each
element of the commutant is fixed. Real and imaginary Hermitian parts give
the full fixed space `span{I,Q}`.

In finite dimension a linear map and its adjoint have equal nullities.
The Schrödinger fixed space is therefore two-dimensional. The independent
fixed matrices `rho_* tensor P_+` and `rho_* tensor P_-` span it; positivity
and trace one give exactly `rho_a`. This argument covers all source and
memory coherences, rather than assuming a classical stationary ansatz.

For every initial state, `a=Tr[(I66 tensor P_+)rho]` is conserved. Ordinary
repair cannot select a different value. If any marginal limit exists, this
conserved value fixes its member of the stationary family. The theorem
does not assume the full growing classical archive becomes independent.

## Finite preparation in a sharp protected sector

Initialize the orientation bit in a specified eigenstate of `X`, using
the existing one-pulse preparation service. The joint input lies in a
sharp sector `q=+1` or `q=-1`; its source may be arbitrary and entangled
with an untouched reference. In that sector the branch is

```math
K_{h,q}=\sqrt{w_h}\sum_{x\in D_h}
             q^{p_h(x)}|f_h(x)\rangle\langle x|.
```

The signs have absolute value one. Every pair of distinct source addresses
that survives one common record stays distinct, because each branch is
injective. The existing complete-record off-diagonal contraction proof in
`local_repair_family.mixing` therefore applies with identical constants:
phases cannot increase the trace norm of a surviving reference block.
The diagonal population kernel is unchanged. The same Doeblin bound then
gives, after `n` attempts and after tracing the preparation archive,

```math
\frac12\|\rho_{SMR}^{(n)}
 -\rho_*\otimes P_q\otimes\rho_R^{(0)}\|_1
\leq \min\{1,\tfrac{65}{2}r^n+(1-\beta)^{\lfloor n/20\rfloor}\}.
```

The exact count certificate gives `beta>=33/200` and `r<=893/990` for
the declared weights. At 1600 attempts the bound is below `10^-6`.
Every preparation record exists in the complete process; only this
source/reference marginal comparison traces those records. The result
supplies neither free record erasure nor an independent full archive.

## One modular reference for the charge and the motion

For `0<a<1`, the full state is faithful and its modular generator is

```math
K_a=-\log\rho_a
 =\left[(\log Z)I-\sum_i(\log\theta_i)n_i\right]\otimes I_2
   -I_{66}\otimes[\log(a)P_+ +\log(1-a)P_-].
```

Tensor ordering in this formula is source then memory; the executable
matrix addresses use memory then source. At `a=0` or `a=1`, the logarithm
is defined only on the supported charge corner. On that corner it is
exactly the source occupation potential, up to its scalar normalization.
No finite energy of the absent charge sector follows from a corner log.

For seam `(u,v)`, let `j` be the actual signed current from port `u` to
port `v`. Every branch column satisfies

```math
n_k(f_h(x))-n_k(x)=j(\delta_{kv}-\delta_{ku}),
\qquad
K_a(f_h(x),q)-K_a(x,q)=j\log(\theta_u/\theta_v).
```

`reference.py` verifies the energy identity through the exact rational
ratio `rho_*(x)/rho_*(f_h(x))=(theta_u/theta_v)^j`, without numerical logarithms.
On the negative exchange sector, moving columns are precisely the existing
CAR bilinears `c_v^dagger c_u` or `c_u^dagger c_v`, restricted to two
occupations. Creation or annihilation separately leaves this fixed-number
carrier. This equality identifies the moving operator and its signs;
it does not identify a recorded stochastic jump with coherent hopping.

Since the raw record fixes `j`, every branch is a modular eigenoperator:

```math
[K_a,K_h]=\Delta_h K_h,
\qquad \Delta_h=j\log(\theta_u/\theta_v).
```

Consequently its complete record map is covariant under
`U_s=exp(-is K_a)`: the phase in
`U_s K_h U_s^*=exp(-is Delta_h)K_h` cancels between a Kraus operator and
its adjoint. This is a same-state, same-clock compatibility theorem for
every retained outcome. It composes under finite record-adaptive choices
of these branches. Stationary preparation and the fixed-space
classification, in contrast, use the specified fixed all-seam selector.

Along any retained history the energy changes telescope between its
initial and final configurations. The actual three-move odd exchange
loop returns to its source pair, so its modular energy change is zero,
although its ordering operation is `X`. The exchange character and the
modular energy change are therefore distinct observables.

The occupation potential commutes with every occupation projector.
Its unitary flow has no occupation propagation. Native repairs do move
occupations and dissipate coherences. Covariance proves their compatibility
with this reference; it does not equate repair-attempt count with modular
time or produce the scalar, gauge and hopping action of a Higgs field.

## Charge preparation, readout and the energy offset

The exact memory operation `F=I66 tensor Z` obeys `QF=-FQ`. It exchanges
the two sharp stationary states and leaves the source untouched. The
inherited coherent memory service realizes its target with a `Z(pi/2)`
pulse and the scalar correction `I(-pi/2)` in the compiler's
`exp(+i angle Pauli)` convention. Measurement of `X` is the inherited
`Y(pi/4)` rotation followed by both computational outcomes. These are
declared observer operations; unassisted native repair commutes with `Q`
and cannot perform a charge flip.

For the faithful reference `rho_a`, the mean modular energy difference
between the sharp minus and plus states is exactly

```math
\operatorname{Tr}K_a(\rho_- -\rho_+)=\log\frac{a}{1-a}.
```

The value is zero at `a=1/2`, negative at `a=1/4` and positive at `a=3/4`.
It is set by a conserved preparation weight, with no change in the native
repair weights or protected charge. This is a dimensionless reference
energy difference, not a physical mass or an energy extracted without
accounting for the charge-changing apparatus.

For `rho_+` and `rho_-`, every complete classical native history has the
same unnormalized source output. Source diagonality removes interference,
and the two branch signs have equal squared modulus. This remains true
for finite adaptive classical policies depending on earlier records.
The full source/memory states have trace distance one; source plus native
tape has trace distance zero. A readout that accesses memory can distinguish
the charges. This proves that the protected information resides in the
memory factor; it supplies no spatial local-core identification.

The complete state and energy connection to an interacting field requires
a map from this source/memory process to a field algebra, preparation,
propagating action and readout. The configured cube and Higgs source
experiments specify such field actions separately. Their literal
preparation and detector results do not select that map for this native
protected charge.

## Weighted reversal of the same complete records

`REVERSAL.md` proves that the Petz reverse of every complete native branch
is an existing reverse-labelled branch of this same biased instrument.
For the identical supported reference, the exact relation is
`rho_a^(1/2) K_h^dagger rho_a^(-1/2)=K_reverse(h)`. Both holding and zero
outcomes remain in the full 720-label involution. It gives quantum detailed
balance, KMS self-adjointness of the Heisenberg channel, opposite modular
transfers and equal probabilities for every complete reversed history
under the stationary fixed-selector experiment. The all-length history
identity follows by multiplying the operator identities; it is not an
extrapolation from a bounded enumeration.

The actual unread channel contracts the trace distance of an explicit pair
of orthogonal inputs from one to `949/990`. Its Petz self-duality therefore
provides no inverse on arbitrary states, no samplewise rewind and no
identification of the modular clock with physical time. This is a quantum
operator statement on the same protected carrier, including arbitrary
untouched references, rather than a substitute population dynamics.

## Finite thermodynamics on the same gauge corner

`GAUGE_JOIN.md` proves that the declared two-singlet, fixed-charge corner
has the complete unread quantum channel
`A -> Tr(A) diag(25/33,8/33)`. Its four retained outcomes give, on source
diagonals, the kernel `K(x,y)=p(y)` with `p=(25/33,8/33)`. These are the
native endpoint coin probabilities and also the native stationary weights
conditioned on the active pair. No energy fit determines them.

`Lean/Thermodynamics/ProtectedGaugeCorner.lean` counts 25 and 8 tickets
in a 33-ticket realization of that same coin, derives the diagonal kernel
by summing all four outcomes, and proves it equals the conditional heat
bath of `RepairLawData` with a constant visible datum. Its energy is the
supported modular energy `E=(log(33)-log(25),log(33)-log(8))`. The thermal
identity `-log p=E` supplies `CalibrationData` at multiplier one and
log partition zero in modular units. The theorem `corner_four_laws`
then instantiates the existing `FiniteFourLawCore` on this same kernel
and energy. Its clauses include the exact first-law split, relative
entropy contraction, detailed-balance fluctuation relations and the
finite low-temperature laws. The exact positive gap is `log(25/8)`;
the mathematical Gibbs family of this fixed two-state spectrum has entropy
tending to zero, and the nondegenerate energy makes its Gibbs multiplier
identifiable. The fixed native coin implements multiplier one, with no
implemented native variation of that multiplier. The theorem
`corner_same_kernel_four_laws` joins the finite conclusion to the equality
between its heat bath and the complete four-outcome diagonal kernel.
The visible datum is constant while this energy is nonconstant. The channel
therefore exchanges modular heat; the library's conservation theorem for
fibre-measurable quantities does not make this energy conserved. Its
zeroth-law conclusion identifies the Gibbs multiplier for a specified
spectrum, without constructing a thermal-contact network.

The executable certificate checks the native coin, conditional reference
and all five gauge-block corner channels against the Lean adapter's exact
probabilities. The quantum matrix-unit calculation binds the diagonal
kernel to the complete instrument, including arbitrary untouched
references; the Lean theorem concerns this diagonal marginal. The outcome
archive is retained in the full process and traced only for this marginal
comparison. Corner admission and the dual reference are declared inputs.
The modular calibration does not identify laboratory energy or temperature,
and the fixed finite entropy limit supplies no refinement-uniform capacity
family. No canonical physical premise is discharged by this adapter.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/protected_reference
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m model.protected_reference.reference --check
cd Lean && lake build ProtectedGaugeCorner
```
