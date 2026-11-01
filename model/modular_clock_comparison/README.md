# Comparing the actual interacting regional clocks

The same two-site Weyl/contact field now has source-executed comparisons
between its global modular clock and both one-site modular clocks. A
comparison between the two regional clocks costs **20 supplied native
pulses**. Each global-to-regional comparison at parameter `s=1/2` uses
**2,826 pulses**, with ideal splitting bound `11/32768`. Every input column,
number sector, untouched quantum reference and first export failure is
retained. No re-preparation or postselection is inserted.

This closes a finite operational comparison premise of CP-02 / [#1138](https://github.com/FloatingPragma/observer-patch-holography/issues/1138).
It uses the accepted coherent service class and the existing supplied Gibbs
reference. It does not calibrate physical duration or imply that a regional
clock is the restriction of the global clock. The latter assertion is
actually impossible on this fixed reference if it must transmit a moving
regional clock hand: the exact certificate below proves that every
common-parameter linear intertwiner has stationary regional range.

## One source and its actual references

The source is the [complete-Fock Weyl/contact dimer](../full_fock_field/README.md),
with four CAR modes, all sixteen occupation states, speed `c=1/8`, contact
`U=1`, and the faithful reference `rho=exp(-H)/Z` at `beta=1, mu=0`.
The existing [thermal preparation](../field_thermal/README.md) and
[modular history](../modular_flow/README.md) supply its preparation and
global source-word interface. The two reduced densities are obtained by
the existing signed CAR restrictions, not separately chosen temperatures
or regional bare Hamiltonians.

Their common spectrum in occupation order is approximately

```
rho_A = diag(0.29691079, 0.29687342, 0.29687342, 0.10934237).
```

The small first splitting is retained. Treating the first three entries as
one degenerate eigenvalue would erase a real finite modular frequency.
The regional generator is `K_A=-log(rho_A)`, up to a scalar identity gauge.
Its source actuator uses all four occupation states, including coherences
between distinct local number sectors. Physical superselection or detector
access may impose further restrictions; they are not silently added here.

## Exact relative-unitary construction

Use the inherited observable convention
`sigma_s(X)=rho^(is) X rho^(-is)`. Its state unitary is
`W_S(s)=exp(+is H)`, after dropping the scalar partition phase. For either
region, embed `W_A(s)=rho_A^(-is)` into the full Fock carrier, again in a
fixed scalar gauge. Define

```
C_AB(s) = W_A(s) W_B(s)*.
```

This is an explicit active change between the two implemented actions:

```
C_AB(s) W_B(s) = W_A(s),
C_AB(s) C_BD(s) = C_AD(s),
C_AS(s+t) = C_AS(s) W_S(s) C_AS(t) W_S(s)*.
```

The last identity is a time cocycle, not a group law for `C_AS` alone.
All three follow by unitary cancellation and the group laws of the two
flows. They hold on the whole carrier, hence with arbitrary untouched
references. [The Lean module](../../Lean/QFT/ModularClockComparison.lean)
kernel-checks these group identities and the finite spectral theorem below.
It does not encode the concrete Gibbs matrix or transcendental intervals.

The actuator executes the inverse origin word followed by the target word.
Regional words are the existing two-bit diagonal compiler applied to the
actual reduced logarithm. Global words use the existing full-field source
split. There is no arbitrary-unitary primitive. Counts include each
proper-code load and export; initialization and physical apparatus costs
retain their existing service contracts.

For segment ideal operator errors `epsilon_A, epsilon_B`, unitary
telescoping gives comparison error at most `epsilon_A+epsilon_B`. The
complete source executor also retains every first-failure Kraus map `F_j`
and successful matrix `S`, with

```
S* S + sum_j F_j* F_j = I,
complete dilation error <= sqrt(||S-C_AB||^2 + ||sum_j F_j*F_j||).
```

The ideal splitting bound is rational. Executed primitive/logarithm
residuals are separate floating diagnostics, not certified hardware bounds.
The baseline global-to-regional error is approximately `2.03e-5`, below
the ideal bound `11/32768`; local-to-local splitting error is zero.

These corrections have a real state cost. At `s=1/2`, the ideal
global-to-regional correction changes the common Gibbs reference by half
trace distance about `0.00444011`, and the regional-to-regional correction
by `0.00733040`. Each regional flow preserves its own reduced density but
generally changes the global correlations. Thus this is a working
operational comparison, with explicit disturbance, rather than a common
stationary regional automorphism or a geometric modular net.

## A complete finite obstruction to uncorrected clock transmission

For a faithful finite reference, diagonalizing its modular generator
decomposes the matrix algebra into eigenoperators. If `L` intertwines the
source and target modular groups with the **same parameter**, then an
input eigenoperator with frequency `omega` can only have output components
of that same frequency. Indeed, in the two eigenoperator bases, each
coefficient obeys

```
(nu_j - omega_i) L_ji = 0.
```

Therefore, if every nonzero target frequency is absent from the source
spectrum, `L` has range in the target modular centralizer. This applies to
every linear map, so imposing positivity, complete positivity or unitality
cannot restore a missing frequency. The Lean theorem
`nonresonant_output_stationary` proves this finite coefficient implication.

The [certificate](comparison.py) verifies the hypothesis for this actual
interacting dimer, rather than assuming a fitted numerical spectrum.
Put `t=1/16` and `D=sqrt(17)/4`. Its nine distinct energies are

```
-t, t, 0, 1-t, 1+t, 1, 2, (1-D)/2, (1+D)/2.
```

The exact characteristic polynomial of the actual sixteen-state matrix is

```
x^4 (x^2-t^2)^2 ((x-1)^2-t^2)^2
    (x-1)(x-2)(x^2-x-4t^2).
```

For each energy, the partial trace of its spectral projector is diagonal.
In the displayed energy order, its four diagonal weights are

| Energy | Regional spectral weights |
| --- | --- |
| `-t`, `t` | `(1, 1/2, 1/2, 0)` each |
| `0` | `(1, 3/2, 3/2, 0)` |
| `1-t`, `1+t` | `(0, 1/2, 1/2, 1)` each |
| `1` | `(1/2, 0, 0, 1/2)` |
| `2` | `(0, 0, 0, 1)` |
| `(1-D)/2` | `((1-1/D)/4, (1+1/D)/4, (1+1/D)/4, (1-1/D)/4)` |
| `(1+D)/2` | `((1+1/D)/4, (1-1/D)/4, (1-1/D)/4, (1+1/D)/4)` |

Nine exact partial-trace identities for `H^k`, `k=0,...,8`, verify these
weights against the source matrix. Distinct eigenvalues make the associated
Vandermonde system invertible, so the moment checks determine every
spectral weight. Multiplying by `exp(-E)` and normalizing gives the actual
reduced Gibbs density.

The positive regional modular frequencies are approximately

```
0.0001258905675, 0.9988218345429, 0.9989477251104.
```

The rational interval kernel bounds `sqrt(17)`, exponentials and logarithms
with explicit series remainders and outward rational rounding. All **243**
comparisons of these three frequencies with the `9 x 9` global energy
differences have separation **strictly greater than `1/10000`**. Negative
frequencies follow by sign symmetry. Hence every same-parameter
intertwiner has range in the six-dimensional complex centralizer
`C direct_sum M2(C) direct_sum C`; none transmits nontrivial local modular
phase on this full regional algebra. No floating equality threshold enters
that conclusion.

This is an exclusion of an optional strong clock identification on one
fixed finite reference, not of CP-02, the accepted service class, other
references, different clock rescalings, approximate maps or a field limit.
Deleting the boundary hopping gives a product reference and a positive
control with nontrivial covariant Petz observation. The active comparisons
above supply a permissible alternative on the unchanged interacting source.

## Reproduction and remaining premises

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.modular_clock_comparison --check model/modular_clock_comparison/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/modular_clock_comparison
cd Lean
lake build QFT.ModularClockComparison QFT
```

Tests independently check exact moments, numerical spectral/logarithm
comparisons, arbitrary complex reference columns, both modular signs,
all complete source outcomes, deliberately injected source leakage,
the non-group cocycle, a product-reference positive control and malformed
contracts. The receipt retains the exact frequency intervals and all
source resources.

The source action, Gibbs preparation, inverse temperature, access to
initialized coherent work bits and interpretation of the parameter remain
declared specification inputs. These are permitted readout/service data in
the executable hypothesis, but they are not physical clock calibration.
Common native reference selection, scalable preparation and interaction,
physical energy/time normalization and geometric modular/Einstein premises
remain open. No canonical physical premise or integrated gate is promoted.
