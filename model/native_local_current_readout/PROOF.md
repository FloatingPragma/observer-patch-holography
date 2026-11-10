# A local current pointer with native evolution retained

## Source and apparatus

The source, actual depth-three preparation, elementary/Higgs pumps, spectral
bin, spin environment and old/dual/port references are those of
[native Higgs current response](../native_higgs_current_response/PROOF.md).
The first three source-clock intervals are unchanged. The current source
operator is the original W, not its charge-conjugate Wbar. This package
replaces only the final measurement service by a specified local interaction.
It neither subtracts Hext from the drive nor changes the native h.

On the complete five-scalar/four-weak-mode carrier W has dimension 80 and
spectrum {-1,-1/2,0,3/2,2}, with projector ranks {4,12,58,4,2}. The other 26
native CAR modes are spectators in this even operator's graded tensor
representation. The provider checks the full 30-mode polynomial and all 12
native Gauss commutators. These statements hold on all occupation sectors.

For each eigenvalue w, take the exact rational Lagrange projector

```
Pi_w = product_(v != w) (W-vI)/(w-v).
```

The executable checks their orthogonality, Hermiticity, completeness, ranks
and W Pi_w=w Pi_w. Give the five outcomes the centered integer labels
{-2,-1,0,1,2}. The same physical W values are recovered by the retained decoder.
Their polynomial is

```
f(W)=(4W^4-4W^3-13W^2+25W)/15,       ||f(W)||=2.
```

The pointer has five neutral states, Hamiltonian zero and cyclic shift
S|r>=|r+1 mod5>. Its Hermitian principal generator is

```
P_ab = (4 pi i/25) sum_(k=1,2) k sin[2 pi k(a-b)/5],
exp(-iP)=S,                         ||P||=4pi/5.
```

Thus exp[-i f(W) tensor P]=sum_w Pi_w tensor S^j(w). On pointer blank 0,
all five pointer outcomes give exactly the five current projectors. Every
column of this 400-state unitary is included, not only the declared input.

## Exact rational control

The two independent real coefficients of the imaginary antisymmetric P
are rounded to a common 40-bit dyadic grid. Opposite entries receive opposite
values, so P40 is exactly Hermitian. Integer interval arithmetic encloses pi
and each sine; the outward Frobenius error epsilon_P bounds operator norm.
It is below 9.078e-13. Consequently

```
||P40|| <= 4pi/5+epsilon_P,
||exp[-i f(W) tensor P40]-exp[-i f(W) tensor P]|| <= 2epsilon_P.
```

The represented coupling G=f(W) tensor P40 is an explicit rational
400-dimensional matrix. Its exponential defines the mathematical service.
An independent rational degree 64 Taylor evaluator uses the exact five
spectral projectors and five 5-dimensional pointer exponentials. Since
||jP40||<6, its uniform operator remainder is at most

```
6^65 / 65! / (1-6/66).
```

The code also exponentiates the complete 400 matrix numerically and checks all
five rectangular Kraus maps. Those residuals are backend diagnostics; they
are separate from the rational control and Taylor bounds. No state is
renormalized and no pointer label is discarded.

The added apparatus primitive is this bounded-support, Gauss-invariant
polynomial-current/pointer Hamiltonian, with switching controlled by retained
success flags. It is not a compiler from the integer C/E rule. Polynomial
powers of W are local many-body terms; an efficient synthesis of that
interaction is not supplied. Because every polynomial in W commutes with
the native Gauss generators and the pointer and flags are neutral, the full
driven Hamiltonian preserves all Gauss constraints throughout the pulse.

## Uninterrupted source clock

Let tau=2^-20 source units and apply V=G/tau on the admitted pump/bin branch.
The success projector acts on retained zero-energy apparatus flags and
commutes with Hext. All other branches have V=0 and continue freely for the
same duration. The complete elapsed clock is 3+tau; the parent's exact-time 4
instrument and receipt remain separate.

The actual unitary is exp[-i tau(Hext+V)]. For any finite Hermitian H and G,
a double Duhamel integral and unitary norm one give

```
||exp[-i(tau H+G)]-exp(-i tau H) exp(-iG)||
    <= integral_0^1 integral_0^s tau ||[H,G]|| du ds
     = tau ||[H,G]||/2.
```

This bound includes arbitrary untouched references. It does not assume that
the current or admitted state is stationary.

Only terms touching the measured cell contribute. The scalar penalty and
all link-energy terms commute with W. The exact native onsite Yukawa word,
including up, down and lepton terms sharing the scalar, has norm at most 56.
The two quark commutators with W are explicitly nonzero, so omitting these
terms would be invalid. At side two the actual kinetic schedule has six
incident directed terms with total coefficient modulus 3/2. Its full-CAR
row/column and link-contraction bound gives touched kinetic norm at most 90.
All other source cells, all other preparation copies and every reference
commute with the detector. Therefore

```
||[Hext,f(W)]|| <= 2*2*(56+90*2^-34) =224+360*2^-34,
||[Hext,G]|| <= (224+360*2^-34)(4pi/5+epsilon_P) = C.
```

No lepton-only approximation replaces the full onsite Hamiltonian in C.
Combining source drift and rational-pointer approximation gives the exact
uniform output-operator bound

```
delta = tau C/2 +2epsilon_P < 0.000268446640.
```

Half diamond distance, and hence total variation of the complete retained
record law, is at most delta. Coherent flag control makes the failure block
identical in both comparisons; it does not multiply the bound by the number
of branches. The ideal trailing exp(-i tau Hext) acts trivially on the pointer,
so its pointer probabilities refer to the current at source time 3.

## Retained response and conditional comparisons

The complete setting/pump/bin/current alphabet has 80 labels, in addition to
all fine preparation flags and partners. Disabled-pointer zero denotes
absence of that measurement, not a zero source-current reading. The actual
instrument retains its full source output after every outcome.

On either actual admitted preparation, each ideal conditional pointer
probability interval widens by delta and is clipped to [0,1]. A lower bound
on the original-attempt joint probability is its original bin-success lower
bound times the widened conditional lower bound. The fair setting adds its
factor one-half. Both W=-1/2 original-attempt probabilities remain positive.
There is no conditioning on a successful final readout.

Since ||W||=2, the decoded mean changes by at most 4delta. Both prepared means
remain strictly negative. These are pointer-decoded mean intervals; no claim
is made that the disturbed final source has that instantaneous W expectation.
Applying this same complete apparatus to the two normalized admitted states
also contracts their parent's trace-distance bound. That comparison conditions
on each chosen setting and retained bin success and compares corresponding
histories. It does not equate the unequal unconditional A and B rates or
compare their distinct classical setting labels as equal quantum states.

## Resources

The five-state pointer, retained registers and flag records are the same as
the parent current construction. Its quantum capacity formula is retained.
The two signed pointer numerators and common denominator are explicit;
coefficient synthesis and controller workspace are not counted as part of
that quantum-register formula.

The final local drive has norm at most 2(4pi/5+epsilon_P)/tau, below 5,270,718.
The two switches have gross work bound twice that norm, below 10,541,436.
For the complete unconditioned source state,

```
||U_actual^dagger Hext U_actual-Hext|| <= C,
```

by integrating i[V,Hext]. This bounds the mean bare-source energy change for
any initial state. It is not a separate measurement of individual work values.
Using the parent's full raw bound B_H, the complete four-interval net energy
change is also bounded by 2B_H and gross switching by
6(B_H+pi)+2||V||. Source partners, references and failure branches enter these
bounds. The short local pulse has a substantial control-cost tradeoff.
Hardware accuracy and efficient realization are distinct from this finite
mathematical apparatus.

The preceding full-matrix preparation and exact spectral-bin services remain
assumptions. Relativistic propagation, weak scattering, a W boson, Standard
Model particle identity and a GeV mass are not consequences of this readout.
