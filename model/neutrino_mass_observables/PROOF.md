# Conditional neutrino mass observables and complete phase records

The input is the optional calibrated Weinberg kernel in
`model/common_source/calibrated_matter.py`. Its normal ordering,
zero lightest mass and two printed Majorana phase settings are declared
inputs; its mass splittings and PMNS angles come from its retained
oscillation packet. With rank two there is one independent relative
Majorana phase. The finite result retains this ancestry and the unproved
running/matching step to the common Higgs background.

## Same-source masses and current weights

Let `U` be the actual unitary PMNS matrix and
`D=diag(0,m2,m3)`, with masses in GeV. The symmetric all-left kernel is

```
M=conjugate(U) D U*,       U^T M U=D.                 (1)
```

The executable extracts this block from all 45 native internal/family
coordinates and verifies that it has no matrix elements to the 42 other
coordinates. The actual weak raising current in the mass frame has
neutral-to-charged-lepton block U*. Thus its electron weights are
`|U_ei|²`, not a separately supplied family assignment. It follows that

```
m_beta² = (M* M)_ee = sum_i |U_ei|² m_i²,
m_ee = |M_ee|.                                        (2)
```

These are finite current-weighted mass norms. Naming their conventional
comparison observables does not calculate an endpoint beta spectrum or
a nuclear double-beta decay rate.

Put `a=m2 sin²(theta12) cos²(theta13)` and
`b=m3 sin²(theta13)`. For this packet a>b>0. The complete phase family
satisfies

```
m_ee²=a²+b²+2ab cos(alpha31-alpha21-2delta),
a-b <= m_ee <= a+b.                                   (3)
```

The code constructs the two endpoint kernels by multiplying the third
PMNS column by `exp(i delta)` or `i exp(i delta)`, starting from the
packet's zero Majorana phases. It verifies their electron matrix elements
against both exact endpoints. Column phases cancel from
`U D² U*`, so the oscillation generator and m_beta are unchanged.
The massless column and the common phase of the two massive columns do
not supply additional independent mass observables.

## A complete mass-contraction instrument

Use the actual spin-resolved creation coefficient
`Delta90=epsilon_spin tensor conjugate(K45)` and its zero-momentum
180-dimensional Nambu generator. Its invariant six-coordinate subspace
consists of the upper spin-zero neutral coordinates and the lower
spin-one neutral coordinates. Direct extraction gives

```
Hnu=[[0,M*],[M,0]],       C=Hnu/m3.                    (4)
```

Here M is symmetric, so `M*=conjugate(M)`. The exact Takagi frame
`W=diag(U,conjugate(U))` transforms Hnu to
`[[0,D],[D,0]]`. Since 0=m1<m2<m3, C is a Hermitian contraction of norm
one. Define the positive defect `F=sqrt(I-C²)`. In the same frame it
is diagonal with entries `sqrt(1-(mi/m3)²)`, repeated twice. Consequently
`C²+F²=I` and CF=FC exactly, and

```
J=[[C,F],[F,-C]]                                       (5)
```

is a unitary dilation on a binary flag. The six success operators are
`K_z=|z><z| C`; the failure operator is F. An eighth operator is the
identity on the complete 174-dimensional spectator complement. Their
squared effects sum to the full identity. No mode or unsuccessful trial
is removed. The resulting isometry preserves arbitrary input coherences
and every idle quantum reference in its complete record isometry.
Forgetting the records need not give the identity channel on the input.
This is a declared logical apparatus
on the Nambu coefficient carrier, not an asserted native pulse word or
a derivation of its many-body laboratory implementation.

On the electron input in the upper block the original probabilities are

```
P(success, lower alpha)=|M_alpha,e|²/m3²,
P(success, upper alpha)=0,
P(failure)=1-m_beta²/m3²,       P(spectator)=0.          (6)
```

Thus m_beta² equals m3² times the total success probability, and m_ee²
equals m3² times the lower-electron success probability. Equations
(3) and (6) give a phase-sensitive finite readout while retaining the
same failed trials. The normalization m3 is a conditional mass input;
it is not an independently measured rate or source-selected coupling.
No physical absolute neutrino mass, Majorana phase, sterile neutrino,
cosmological mass bound, nuclear matrix element or decay rate is derived.
