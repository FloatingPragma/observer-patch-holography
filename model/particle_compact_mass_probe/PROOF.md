# A four-outcome readout of the same finite filtered source response

## Fixed target operator and a different apparatus

Use the selected open Wilson Hamiltonian H, physical polynomial state rho,
full conjugate reference with idle Hamiltonian zero, and source clock of
the [primary contract](../particle_common_contract/PROOF.md). The
[finite survival theorem](../particle_primary_survival/PROOF.md) specifies
an algebraic Hermitian contraction

```
B_tilde=sum_j w_j C_R(t_j)* A C_R(t_j),
w_j>0,     sum_j w_j=1,     B_tilde*=B_tilde,     ||B_tilde||<=1.
```

Here A is its actual plaquette probe, the weights and nodes are its exact
positive triangular quadrature, and C_R are its actual original-H Cayley
gates. M and R retain their full factored values. This construction changes
neither B_tilde nor the source, state, reference, energy origin or clock.

The [filtered-probe theorem](../particle_filtered_probe/PROOF.md) realizes
B_tilde as one Kraus block of an M-state coherent selector with a neutral
probe qubit. The following construction applies that same theorem's
Hermitian-contraction dilation directly to B_tilde. It is a distinct
declared exact algebraic apparatus service, without an M-state selector.
It is not an autonomous repair implementation or a compression of the
selector's complete output channel.

## Direct dilation and four complete Ramsey records

Write B=B_tilde and C=sqrt(I-B^2). Reuse
`filter.hermitian_dilation(B)`:

```
U_B = [[B,C],[C,-B]],       U_B*=U_B,       U_B^2=I.
```

The root is positive and commutes with B by functional calculus. A neutral
qubit initially in 0 gives the complete two-outcome instrument

```
K_0=B,        K_1=C,        K_0*K_0+K_1*K_1=I.         (1)
```

Every source gauge commutant of B is preserved by C and U_B. The reference
is untouched. Exact algebraic matrix functions are finite computable
objects on each finite carrier. Their computational or physical service
cost is not bounded by the qubit count.

Use exactly the inherited Ramsey construction, with u either the ideal
source evolution or the specified Cayley echo, and phase phi=0 or -pi/2:

```
K'_r=u* K_r u,
L_(r,b)=[K_r+(-1)^b exp(i phi)K'_r]/2,
r,b in {0,1}.                                          (2)
```

The four effects sum to identity. Zero-probability labels are retained.
For the actual stationary source, rho commutes with u, and

```
F0=Tr rho B*B,
F_u=Tr rho B* u* B u,
P(0,b)=[F0+(-1)^b Re(exp(i phi)F_u)]/2.                (3)
```

Score (-1)^b on r=0 and zero on r=1. Its original mean is Re F_u or
Im F_u. The total blank probability is F0. No success division is used.
The finite Cayley echo has the same stated comparison error to the
unchanged source-time correlation as the inherited apparatus; it does not
define a replacement Hamiltonian.

The selector apparatus has the identical K_0=B. Its two blank Ramsey
Kraus operators in each phase therefore equal (2) exactly, on every
source input including every retained reference. Original successful
probabilities, source/reference blank outputs and scores agree. Its other
Kraus operators generally differ from C. Their summed quantum failure
channel and their Ramsey failure-bit distribution need not equal those
of (1)-(2), even after their many labels are grouped as failure. No
whole-channel equality or inherited whole-channel locality estimate is
asserted between these two different instruments.

## Same original spectral and excitation weights

The transfer measure of `Tr rho B* alpha_t(B)` depends only on H, rho and
B. These are identical to the finite survival construction, so its proved
positive original band `[delta_n,W_n]` and factored lower weight are
unchanged. After (1)'s blank outcome, the unnormalized source state is
B rho B*. The
[mass-readout transport theorem](../particle_primary_mass_readout/PROOF.md)
therefore gives the same final ground-relative excitation band, with no
additional weight loss on this full-width interval.

Summing both blank Ramsey bits gives a mixture of B rho B* and its
u-conjugate. Since rho commutes with u and u is a function of H, any
final H-energy spectral event has the same original probability in this
mixture as after (1). This does not identify their whole joint outputs.
The positive finite band supplies neither fixed-band noncollapse nor a
rest mass; its lower endpoint and weight can vanish and its upper endpoint
can diverge along refinement.

## Complete approximation errors on the compact instrument

Suppose two Hermitian contractions B,D satisfy `||B-D||<=delta`. Then
`||B^2-D^2||<=2delta`. The reused square-root continuity bound gives

```
||sqrt(I-B^2)-sqrt(I-D^2)|| <=sqrt(2delta).
```

Their two Kraus columns differ in operator norm by at most

```
s(delta)=sqrt(delta^2+2delta).                         (4)
```

For the same echo, summing the two Ramsey bit differences cancels the
cross terms. The four-column squared difference is the average of the
original two-column squared difference and its u-conjugate. Hence it
also has norm at most s(delta). Purifying an arbitrary input and using
trace-norm contraction gives complete output half-trace error at most
`min(1,s(delta))`, uniformly over all untouched references.

An echo-unitary operator error d contributes at most 2d by the unitary
product bound. A complete four-outcome compiler error sigma and initial
preparation error zeta give the sufficient bound

```
epsilon_complete <=min(1,s(delta)+2d+sigma+zeta).       (5)
```

Equation(5) compares compact instruments to compact instruments. It does
not compare their failure maps with the selector instrument. The exact
declared root service has delta=0; its inherited finite B definition and
Cayley readout allowances are kept separately.

For phi=0 the contrast is the single bounded record score
`2*1_(r=0,b=1)`. Consequently a whole compact-instrument half-diamond
allowance sigma changes that contrast by at most 2sigma. If the inherited
finite apparatus lower contrast is d_n/2 and sigma=d_n/16, the compiled
contrast is at least 3d_n/8, with exact preparation. Additional preparation
or sampling allowances must be included explicitly.

The optional compiler reuses the existing arbitrary-instrument theorem
under its coherent source-code load/export premise. With q padded source
bits, four fine outcomes and p=q+2, its sufficient pulse count is

```
20*2^(2p)+4*2^p+10q.
```

Writing `d_n=c_n 2^(-N_n)`, its precision bounds at sigma=d_n/16 are

```
port bits  =2p+N_n+ceil(log2(1972*16/c_n)),
angle bits = p+N_n+ceil(log2(8*(p+2)*16/c_n)).           (6)
```

These exact factored counts retain the density-floor exponent N_n.
It is not the carrier capacity ln D. Compiling the entire four-outcome
instrument avoids assuming that a close dephased instrument can be
coherently controlled. Code-complement and compiler-failure records stay
in the comparison. The source embedding, gatewise Gauss preservation and
a full spatial compiled execution are not supplied by this adapter.

## Quantum carrier size and unresolved implementation cost

The neutral dilation qubit and Ramsey qubit provide all four fine labels.
With the source and its declared full conjugate reference, one invocation
uses the quantum carrier

```
D_joint=4 D_source^2,
N_joint=2 ln D_source+2 ln2 nats.                      (7)
```

Using the [source-capacity bound](../particle_primary_capacity/PROOF.md),

```
N_joint <= (120V+2)ln2+6 ell V ln451,
qubit padding <=2V(60+27ell)+2.                       (8)
```

At degree one, replacing `ell ln451` by ln119 gives the exact ambient
carrier expression. At levels 1 and 2 the padded quantum counts are 116738
and 9240578 bits, respectively. The larger M of the quadrature remains
in the definition of B, but is not an allocated coherent selector register.

Equation(7) counts one invocation of the exact direct-service target.
Extra code-complement or export-failure flags of the conditional compiler
need their own allocation. This quantum-carrier result does not bound
the size of B's exact coefficient description, root computation, control
precision, controller and compiler workspace, repeated record archives,
preparation implementation, or empirical sample count. A compact symbolic
formula is not an execution of those services. The source operation can
be expensive or nonlocal even with two neutral qubits. No efficiency,
total self-contained capacity or necessary detector lower bound follows.

## Exact independent control

The executed control uses H=diag(0,4), rho=diag(3/4,1/4), A=sigma_x,
triangle half-width 1/2, two quadrature cells and one Cayley step. Its
exact finite target is `B=(3/5)sigma_x`, so `C=(4/5)I`. The same Cayley
readout at time 1/2 gives `F0=9/25`, `ReF=0`, `ImF=9/50`. The eight-label
selector/contraction Ramsey instrument and four-label direct instrument
have identical blank Kraus maps and these original scores.

On source input |0>, the old aggregated filter failure state is
`(16/25)|1><1|`, while the compact failure is `(16/25)|0><0|`.
This explicitly excludes complete failure-channel equality. Exact
matrix completeness, both Ramsey phases, zero/involution limits and
entangled-reference tests retain every output. These small matrices are
an independent control, not a numerical execution of the primary source.
