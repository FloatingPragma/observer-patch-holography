# One native reference with gauge transport

The protected native state and the finite gauge-link transport admit one
explicit joint process. Its recorded successful moves carry the native
ordering flip, preserve the endpoint Gauss actions, and have energies
relative to the same native modular reference. Its forward and reverse
probabilities are derived from the endpoint coins. The matter and link
registers, the active edge and the flux-dependent held outcome specify a
declared extension. This is a finite common-model construction, with no
claim that the unextended observer selects this gauge law or a coherent
field Hamiltonian.

## The actual native edge and the declared fiber

Use the native seam `(0,11)` and a spectator occupation at port 5. The two
source configurations are `x={0,5}` and `y={5,11}`. A successful move
between them crosses the spectator once in the sorted-pair chart. The
existing ordered-cover map therefore applies `X` to its orientation bit
in either direction. Both C and E native columns are checked directly.

The native weights are `theta_0=25`, `theta_5=theta_11=8`. Hence

```math
\pi_x=200/7505,\quad \pi_y=64/7505,\qquad
\alpha=\Pr(\mathrm{coin}=1)=25/33,\quad
\beta=\Pr(\mathrm{coin}=0)=8/33,
\quad \alpha\pi_y=\beta\pi_x.
```

For one of the five retained matter blocks let `d` be its dimension.
The literal sixty-state link and its normalized coordinate `B` come from
`gauge_defect_transport.transport.literal_transport(charge, completed=False)`.
This uses the literal raising coordinate without the separate rotor
completion. The exact identities are `B†B=P`, `BB†=R`, with `P,R`
orthogonal projectors. The fiber `F=V_r tensor L` has dimension `60d`.

The full declared carrier is
`H66 tensor C2_orientation tensor F`. Attach the internal fiber to the
designated moving token on the active pair. The fixed spectator supplies
the native ordering parity and receives no additional matter fiber in
this construction. This is not the two-particle colored code with a
separate internal label per particle. Protection of `Q=I66 tensor X tensor I`
therefore proves a memory-charge statement, without identifying that
charge with an individual matter particle or deriving its statistics.

## Four complete outcomes on the full carrier

Define

```math
T=|x\rangle\langle y|\otimes X\otimes B,\qquad
P_y=T^\dagger T,
\qquad P_x=TT^\dagger.
```

The forward coin has Kraus maps `sqrt(alpha)T` and
`sqrt(alpha)(I-P_y)`. The reverse coin has `sqrt(beta)T†` and
`sqrt(beta)(I-P_x)`. These are four separate retained records. Since
`T` is a partial isometry, their effects sum to
`alpha(P_y+I-P_y)+beta(P_x+I-P_x)=I`.

Both moved outcomes vanish on the other 64 native configurations. Both
held outcomes act as identity there. All those configurations and every
link-flux state are retained. The executable stores the nontrivial
`4*60d` active block and states the identity complement of dimension
`128*60d`; this is an exact direct-sum factorization, not a projection or
success normalization. Cross-coherences between the active block and the
complement are acted on by the same full Kraus maps.

The record distinguishes a successful motion from a blocked attempt.
Forward success has native current `j=-1`, reverse success has `j=+1`,
and either held outcome has actual current zero. The coin alone cannot
be called a motion record. Link support can block a move that the bare
native endpoint rule would perform, and the complement uses an explicit
hold policy. This process therefore differs from the original C/E law.
No original raw label is reused to conceal its additional flux readout.

Every map commutes with `Q`: the moved map contains `X` and its support
projections are orientation independent. The existing endpoint Gauss
intertwiners for `B` lift by tensoring with `X`. At `y` and `x` use the
incoming and outgoing endpoint representations defined in the link
package. On the other 64 source configurations choose a fixed declared
representation, for example the trivial representation of these two
endpoint groups. The direct sum is a representation on the full carrier,
and identity holds commute with it. This constructs one edge operation;
it does not supply a common gauge network on all thirty native seams.

Thus every retained outcome is gauge covariant and charge preserving,
including for arbitrary untouched references and finite adaptive
compositions of these instruments. It is an exact finite algebraic
construction. Literal source pulses, apparatus accuracy and a universal
many-particle collision extension are separate inputs.

## The native reference is stationary on the same extended process

For any `0<=a<=1`, define the normalized full state

```math
\Omega_a=\rho_*\otimes\sigma_a\otimes I_F/(60d),
\qquad \sigma_a=aP_++(1-a)P_-.
```

The maximally mixed internal/link factor is a declared preparation; its
stationarity follows below and does not prove unique selection. Since
`X sigma_a X=sigma_a`, the changes caused by a complete forward or
reverse instrument are

```math
\Phi_f(\Omega_a)-\Omega_a
 =\frac{\pi_y}{60d}\left[|x\rangle\langle x|\otimes\sigma_a\otimes R
                         -|y\rangle\langle y|\otimes\sigma_a\otimes P\right],
```

```math
\Phi_r(\Omega_a)-\Omega_a
 =\frac{\pi_x}{60d}\left[|y\rangle\langle y|\otimes\sigma_a\otimes P
                         -|x\rangle\langle x|\otimes\sigma_a\otimes R\right].
```

Their native weighted sum is zero by `alpha pi_y=beta pi_x`. This proves
stationarity on the entire carrier, including the inactive configurations.
The calculation retains memory coherences in the computational basis;
it is not only a population check. No uniqueness or mixing statement
holds by this argument, and the inactive complement gives explicit
additional fixed states. The 1600-attempt preparation theorem for the
ungauged all-seam process is not transferred to this edge law.

This state commutes with the endpoint gauge action. That is
gauge invariance of a density, not support on the kernel of every Gauss
generator. A faithful density on a nontrivial representation cannot have
that kernel support. No physical Gauss-singlet preparation follows from
the tracial fiber alone.

For `0<a<1`, the full reference is faithful. Its modular generator is

```math
K_{\Omega_a}=K_a\otimes I_F+\log(60d)I.
```

The same protected-state generator therefore assigns

```math
[K_{\Omega_a},T]=\log(8/25)T,
\qquad [K_{\Omega_a},T^\dagger]=\log(25/8)T^\dagger.
```

The held maps commute with this generator. All four outcome maps are
covariant under the same modular clock, with all records retained. At
`a=0` or `a=1` these statements hold on the supported charge corner.
The common reference and these dimensionless energy transfers are derived
from the stated joint law. Modular time is distinct from the native
attempt counter and a laboratory clock.

The operation `T+T†` permits a separately parameterized coherent transport
control. Its coefficient is not fixed by completeness, Gauss covariance,
or the detailed-balance coin. Indeed it need not commute with the
nonuniform native reference. No coherent Hamiltonian, vacuum particle
pole, chiral continuum or measured mass is inferred from this stochastic
common-model connection.

## A stationary Gauss-singlet corner

A separate declared preparation can satisfy the stronger constraint.
Adjoin a reference in the dual of `V_r` at endpoint `y` and use the
existing normalized singlet

```math
|\psi_y\rangle=d^{-1/2}\sum_b|y,b,0_{\rm link}\rangle
                         \otimes|\bar b\rangle,
\qquad |\psi_x\rangle=T_0|\psi_y\rangle,
```

where `T_0=|x><y| tensor B` omits the orientation factor. The literal
endpoint equations prove that all endpoint Gauss generators annihilate
both vectors. They are orthonormal, and `T_0† psi_x=psi_y`. The four
complete maps preserve their two-dimensional span; their held maps and
the reverse annihilations are checked as well. This is a reducing corner
of the full instrument, with every orthogonal output retained.

Tensoring the corner with `sigma_a` on the orientation bit, the full move
applies `X` to memory and `T_0` to the singlet. Since `X sigma_a X=sigma_a`,
the state

```math
\Omega_{\rm singlet,a}
 =\frac{25|\psi_x\rangle\langle\psi_x|
          +8|\psi_y\rangle\langle\psi_y|}{33}\otimes\sigma_a
```

is exactly stationary under the same native-biased four-outcome
instrument. Its weights are the original `pi_x,pi_y` normalized within
the active pair; they are not fitted energies. It is supported on the
endpoint Gauss singlets, unlike the faithful full tracial-fiber reference.
Its source marginal is this conditional two-configuration marginal,
rather than the original full 66-state reference. Both the extra dual
reference and this preparation are explicit inputs; no native pulse
preparation of them is executed by this package.

There is an exact preparation statement within this admitted corner.
Fix a sharp memory character `q=+1` or `q=-1`. On the ordered logical
basis `(psi_x,psi_y)`, the four Kraus maps are

```math
\sqrt\alpha q|x\rangle\langle y|,\quad
\sqrt\alpha |x\rangle\langle x|,\quad
\sqrt\beta q|y\rangle\langle x|,\quad
\sqrt\beta |y\rangle\langle y|.
```

Their complete unread channel sends any operator `A` to
`Tr(A)[(25/33)|x><x|+(8/33)|y><y|]`. This is checked on all four
matrix units. Tensoring the identity on any untouched reference proves
exact replacement after one attempt on arbitrary reference-entangled
inputs in the two-singlet, fixed-charge corner. The initial corner and
dual reference are declared resources; no admission of the unused
orthogonal sectors into this preparation theorem is inferred. All four
records exist in the complete output. Only the stated marginal
replacement comparison traces them, and their correlations are not
erased or assigned an independent product law.

The supported modular generator has the same forward gap `log(8/25)`
and the same charge offset as the full reference. Its logarithm makes
no energy assignment on the unused orthogonal complement. Consequently
this construction supplies a common gauge-constrained state, complete
transport process, protected memory charge and modular clock on a stated
finite domain. It does not supply a vacuum or a physical particle spectrum.

```sh
OPENBLAS_NUM_THREADS=1 python -W error -m pytest -q model/protected_reference
```
