# A finite charged-energy detector on the same source clock

## Source and apparatus domain

The source is the unchanged finite Hamiltonian of
`model/affine_gauss_mobile`, with the complete spectral certificate in
`model/affine_gauss_full_spectrum`. Its exact neutral vacuum is the initial
state. All thirteen source CAR modes, fourteen static compensator charges,
twelve longitudinal bits and 5456 transverse states remain present. The
physical carrier has dimension 44695552. The conserved one-charge sector
has dimension 70928. The reference charge stays at its original vertex.

Write Hraw for the source Hamiltonian, Evac for its exact neutral energy,
and H=Hraw-Evac I when reporting excitation energies. This notation does
not change the autonomous clock exp(-it Hraw). The independently replayed
source certificate gives a unique one-charge ground with

```
E1 in [lambda-10^-12,lambda+10^-12],
lambda = 8438888914228825/4503599627370496,
E2 > 2.00219970691,          0 < H|_(N=1) < 12 I.
```

The charged energy includes one unit of static-reference energy. The
normalized zero-photon part z of the source's exact binary trial has
ground probability f in a certified interval whose lower endpoint exceeds
0.999927753447. This preparation is source-derived and uses no measured
mass. Every source coefficient, the full trial residual and the complete
parent dependency closure are replayed before apparatus bounds are admitted.

## Complete CAR rotation and number record

On the full physical carrier, the compensated creators are
`D_y=V c_y† V†`, where V is the exact Gauss embedding. For the normalized
spatial mode z, `D_z=sum_y z_y D_y` obeys

```
D_z²=0,        {D_z,D_z†}=Iphys.
```

All 8192 occupation words and every one of the 169 mixed CAR pairs are
checked with exact integer signs. Normalization of z is the exact square
root of a sum of thirteen squared binary rationals. Consequently
`A=D_z-D_z†` satisfies `A†=-A` and `A²=-Iphys`, so

```
Upump = sqrt(3)/2 Iphys + A/2
```

is exactly unitary. It sends the neutral vacuum to
`sqrt(3)/2 |vac> + |z,0ph>/2`. Measure total source number with a retained
fourteen-state pointer, using the complete unitary
`sum_(n=0)^13 P_n tensor S_14^n`. Pointer n=0 occurs with probability
3/4, n=1 with probability 1/4, and all other n have probability zero.
Reference number R equals n on every nonzero branch. This number record
also removes coherence between the preparation branches; no success-only
normalization is used for unconditional probabilities.

The literal boundary-zero extensions of D_y to arbitrary kinematic
registers need not obey CAR. The complete kinematic pump is therefore

```
Pphys Upump Pphys + (I-Pphys),
Pphys = product_v 1_{G_v=0}.
```

Here Pphys is a finite register-equality test, since physical states satisfy
`R=N` and `ell_y=n_y`. It is not an energy spectral service. The pump is
identity on every nonphysical complement and commutes with all Gauss
generators. The source and reference parities change together, so the
joint pump is even. Its external apparatus can supply work.

## Finite phase pointer and every outcome

Set `M=2^20`, `t0=2 pi/16` and prepare an M-state pointer in
`M^-1/2 sum_j |j>`. On the retained n=1 branch apply

```
C = sum_(j=0)^(M-1) |j><j| tensor
       [exp(-i j t0 Hraw) exp(+i j t0 Evac)].
```

The first factor uses the actual source clock. The second is an explicit
source-derived pointer phase, essential because the certified excitation
energies are measured above Evac. It changes no source generator or
reference energy. Binary control realizes C with twenty controlled powers.
Each extra control term is `Pphys P1 (Hraw-Evac) tensor |1><1|`, in
addition to autonomous Hraw. Over the corresponding controlled duration,
the exact total propagator is common free source evolution times the
displayed controlled gate. All such factors commute. Thus the physical
experiment retains `exp(-i T Hraw) C` for the total elapsed duration T.
Every control term commutes with Gauss and N. On other number branches
and the nonphysical complement the extra apparatus interaction is zero,
its pointer remains |0>, and autonomous source evolution continues.
Identity in the displayed controlled gate is an interaction-picture
statement. It does not silently pause, reset or cancel the source.

Apply the positive Fourier transform
`Fplus[k,j]=exp(+2 pi i k j/M)/sqrt(M)` and retain the outcome k. The sign
is fixed by the negative source evolution: a grid energy `E=16r/M` gives
k=r. The Kraus operators on the charged sector are

```
K_k = (1/M) sum_(j=0)^(M-1) exp(2 pi i j k/M) exp(-i j t0 H),
sum_(k=0)^(M-1) K_k† K_k = I.
```

The second identity follows from discrete Fourier orthogonality. It holds
on the entire charged carrier, with every photon occupation and every
charged energy. On all other retained sectors use `K_k=delta_(k,0) I`.
Together with the number measurement and pump this defines 14M complete
Kraus records. The actual maps include the common left factor
`exp(-i T Hraw)` after these operations. It preserves every norm and
ground fidelity used below, so it is suppressed only in the probability
formulas. For the specified initial state their exact law is

```
p(0,0)=3/4,
p(1,k)=||K_k z||²/4,          0<=k<M,
p(n,k)=0                    for every other record.
```

All M charged outcomes and every zero-weight number outcome are part of
the apparatus. No quantum source state, compensator or complement is
erased. The exact formula specifies every outcome without enumerating a
70928M-dimensional source-pointer matrix.

## Uniform circle-grid bound

For an exact energy E, put `theta=E/16` modulo one. The pointer probability
at k is the squared geometric sum

```
q_k(theta) = |M^-1 sum_j exp[2 pi i j(k/M-theta)]|².
```

Let d be the circular distance from theta to k/M. For `0<d<=1/2`,
`|sin(pi d)|>=2d`, while the numerator of the Dirichlet quotient has
absolute value at most one. Thus `q_k<=1/(4M²d²)`. Consider outcomes with
d at least delta, where `0<delta<=1/2`. Split the circle at theta into its
two oriented semicircles. On each side distances are separated by 1/M,
so their inverse-square sum is bounded by

```
sum_(r=0)^infinity (delta+r/M)^-2
 <= delta^-2 + integral_0^infinity (delta+x/M)^-2 dx
 = delta^-2 + M/delta.
```

The antipodal point can be counted twice, which only enlarges the bound.
Consequently the complete tail satisfies

```
sum_(distance>=delta) q_k <= B_M(delta),
B_M(delta) = min(1, 1/(2M²delta²)+1/(2Mdelta)).
```

This proof includes off-grid phases and arbitrary alignment with the
first retained grid point. It uses no continuous pointer or limiting
Fourier transform.

Accept the integer labels satisfying `|16k/M-lambda|<=.03`. Exactly
k=120836 through 124768 are accepted, a set of 3933 outcomes. For every
ground energy in the certified radius, rejected labels have circular
distance at least `delta_g=(.03-10^-12)/16`. For any other actual charged
energy, the whole-source bounds `E>=L2` and `E<12` imply distance from the
accepted labels at least

```
delta_e = min((L2-lambda-.03)/16, 1-(12-lambda+.03)/16)>0.
```

The second term controls phase wrapping. The apparatus runs only on the
n=1 record; a norm bound from that sector is never applied to other
charge sectors. Define `alpha=B_M(delta_g)` and `beta=B_M(delta_e)`.
Then ground rejection is at most alpha and excited acceptance at most
beta. For this source alpha is below 0.000254442502 and beta below
0.000077554442.

## Finite rational source gates with a systematic error

The exact controlled exponentials above are explicit finite apparatus
actions. An additional exact finite-matrix realization replaces them by
rational operations on the same source. The general finite-clock identity
in `model/coherent_defect_wave/chiral_joint/finite_clock/PROOF.md` is
independently replayed here; none of that package's different source
probabilities is imported. For self-adjoint X put

```
a=I-X²/10,       b=X/2-X³/120,
R33(X)=(a+i b)^-1(a-i b),
a²+b²=I+X²/20+X⁴/600+X⁶/14400 >= I.
```

The inverse exists and R33 is exactly unitary. With
`p(z)=1+z/2+z²/10+z³/120`, exact integration by parts gives

```
exp(z)p(-z)-p(z)
 = -z^7/720 integral_0^1 exp(uz) u³(1-u)³ du.
```

The weight integral is 1/140 and `|p(-itE)|>=1`. The finite spectral
theorem and unitary product telescoping therefore give

```
||R33(t A/n)^n-exp(-it A)|| <= |t ||A|| |^7/(100800 n^6).
```

The actual neutral energy is `Evac=Tr(Omega)/2<45`, because all thirty
frequencies are below three. On the full charged sector
`||Hraw||<Evac+12<57`. For n=4096 implement the base gate as

```
Uhat=R33(t0 Hraw/n)^n R33(-t0 Evac/n)^n.
```

The scalar second factor is the explicit finite pointer phase correction.
The gate is extended by identity off `Pphys P1`; common autonomous source
evolution on all sectors is retained as above. Bounds 12 and 57 are used
only on the charged physical sector. They are not bounds on higher
charge or nonphysical sectors. The ideal and approximate extra apparatus
gates coincide there, so their full difference has the stated sector norm.
All coefficients, inverses, the finite Fourier matrix and the pump's
normalization are declared exact computable finite-matrix services. An
efficient local gate synthesis and finite-precision service implementation
are not asserted. The rational approximation error itself is proved,
rather than assigned to an unspecified simulator tolerance.

Using `pi<22/7`, base-gate error is at most
`[(627/28)^7+(495/28)^7]/(100800 n^6)`. Each controlled branch has j
applications. Since j<=M-1, the full controlled-unitary error is at most

```
epsilon_c=(M-1)[(627/28)^7+(495/28)^7]/(100800 n^6)
         < 7.409*10^-12.
```

Fourier operations and complete pointer measurements cannot increase
trace distance. Thus the total variation difference of the complete
phase record law is at most epsilon_c on every state, including a
retained reference. All factors are functions of actual Hraw, so its
eigenvectors, every Gauss sector and N are preserved exactly. Applying
the same bound to each eigenvector separately gives ground rejection at
most `alpha+epsilon_c` and excited acceptance at most
`beta+epsilon_c`. This is stronger than subtracting a global state error
from a rare false-acceptance probability. Rational eigenphases are not
identified with exact source energies; their deviation is the explicit
apparatus error above.

Twenty binary controlled powers require at most
`2n(M-1)=8589926400` rational factors in the straightforward serial
realization. The largest coherent source duration is `(M-1)pi/8`.
These are declared finite resource bounds. They describe the controlled
extra interaction; the common source evolution is not removed. The large source-apparatus
matrix is not numerically executed, and these counts are not derived
from an observer memory law.

## Acceptance, false positives and retained conditional state

Write `f>=f_lower` for the actual ground weight of the prepared z state.
Every Kraus operator is a function of H, for both exact and rational
gates. The accepted ground contribution is at least

```
g=f_lower(1-alpha-epsilon_c),
```

and the accepted excited contribution at most

```
b=(1-f_lower)(beta+epsilon_c).
```

Ground and excited spectral subspaces are orthogonal, so no interference
term changes their probabilities. The conditional detector acceptance is
at least g. Its upper bound is
`f_upper+(1-f_upper)(beta+epsilon_c)`. Including the pump, the complete
three aggregate records are `3/4`, `p_accept/4` and
`(1-p_accept)/4`. The receipt retains outward intervals for both detector
aggregates and the exact full 14M-label law.

The actual unconditional false-positive probability is at most b/4.
Conditional on acceptance, with the phase outcome still retained, the
source ground fidelity is at least `g/(g+b)`. This bounds the source
marginal of the accepted joint state; it does not assert an exactly pure
ground or an erased phase record. The bounds give

```
Pr(accept | pump success) > 0.99967332932,
Pr(pump success and excited acceptance) < 1.400761*10^-9,
accepted source ground fidelity > 0.99999999439512.
```

The first probability is conditional on the explicitly retained pump
success. Unconditional acceptance is greater than 0.24991833233. Every
failure, rejected phase and zero-weight register outcome is kept.

## Capacity and physical interpretation

The phase pointer has M states and the number pointer fourteen. The
unchanged full kinematic source has `2^13*14*2^12*5456` states, including
its charge reference. Their combined log dimension is

```
ln(2^13*14*2^12*5456)+ln(14*2^20).
```

The receipt gives an outward interval for this logarithm. This is a
logical capacity, not the entropy of the prepared pure state.
Gauss reduction lowers the source dimension to 44695552. The twenty
control qubits are the same phase pointer, not an additional register.
Extra gate-synthesis workspace and autonomous cooling are outside this
accounting; blank-register preparation and exact finite coefficient
services are declared architecture choices.

The result supplies a finite complete energy measurement with quantified
ground preparation fidelity on the fixed interacting source. It supplies
no exact spectral projector, translation-invariant mass shell, Lorentz
spin, Standard Model species or GeV conversion.
