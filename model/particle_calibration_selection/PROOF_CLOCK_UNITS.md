# Native record rates without a fitted time unit

The declared source is the binary total-two carrier with uniform selection
of thirty seams and a fair placement coin. All sixty branch outcomes,
including no movement, count as one attempt. The calculation imports its
integer transition counts from
[`native_quantum_reconstruction`](../native_quantum_reconstruction/PROOF.md).
It does not supply the adaptive observer's schedule, a physical clock, or
the time dependence of the calibrated charged action.

## Exact rates from complete records

Let B be the 66 by 12 occupation matrix: its row is the twelve binary
occupations of one configuration. Let L be the actual icosahedral graph
Laplacian and P the symmetric complete-attempt transition matrix. Direct
enumeration of every seam and coin gives

    P B = B (I-L/60),        B^T B = 10 I + 11^T.

The first identity also follows by summing the local occupation changes:
an edge incident on a port replaces its expected occupation by the mean
of its two endpoints. Hence, for every integer n>=0,

    P^n B = B (I-L/60)^n.

The exact characteristic polynomial of L is
`z (z-6)^5 ((z-5)^2-5)^3`. Let Pi_lambda be one of its nonconstant
orthogonal spectral projections, of rank d. For a uniformly prepared
configuration s, score its endpoint by the twelve-vector
`B(s,:) Pi_lambda`. This is a specified function of the complete occupation
record. The stationary correlation, with every history retained, is

    C_lambda(n) = (1/66) Tr(Pi_lambda B^T P^n B Pi_lambda)
                = (5d/33) (1-lambda/60)^n.

The squared-score weights are 5/11, 25/33 and 5/11 for the slow, middle
and fast modes. They are nonzero. They need not sum to one, since this is
a weighted correlation rather than a probability partition. Current-only
reconstruction of the endpoint configurations has the preparation and
census bounds of
[`native_thermal_return`](../native_thermal_return/PROOF.md); the formula
above is the ideal complete-law value. It does not discard failed census
records or promise an error-free physical detector.

The three logarithmic rates per attempt are therefore

    e_s = -log((55+sqrt(5))/60),
    e_m = -log(9/10),
    e_f = -log((55-sqrt(5))/60).

Their ratios use no measured mass, phase or energy unit. These are native
record-relaxation rates. They are also gaps above the constant eigenvector
of the reconstructed H=-log P. This functional-calculus identification
does not turn the stochastic P into an executed unitary channel.

## The common scale cancels

The closed-history reference `rho_L=P^L/Tr(P^L)` gives
`-log rho_L=L H+log Tr(P^L) I`. Its spectral gaps are L times those
of H; the additive scalar cancels. The modular action uses the signed
Heisenberg conversion `tau=-Ls`, as proved in
[`native_modular_dynamics`](../native_modular_dynamics/PROOF.md).
Grouping b complete attempts gives `-log(P^b)=bH`. For any common
positive conversion a, the ratio `(ab e_i)/(ab e_j)=e_i/e_j`.
This removes the arbitrary common unit and block length from a ratio on
this source. Different references generally have different expectations;
their common H alone is not a state-identification theorem.

For positive g_b and `epsilon_b<g_b`, bounded gap errors obey

    |(g_a+delta_a)/(g_b+delta_b)-g_a/g_b|
      <= (epsilon_a g_b+g_a epsilon_b)/(g_b(g_b-epsilon_b)).

The program computes rational interval bounds on every displayed rate
and ratio. For rational `0<p<=1`, set `y=(1-p)/(1+p)`. Summing the first
m terms of `2 sum y^(2k+1)/(2k+1)` underestimates `-log p`; the remainder
is at most `2 y^(2m+1)/((2m+1)(1-y^2))`. Rational integer-square bounds
enclose sqrt(5), so no numerical logarithm or eigensolver supplies the
certificate. The all-horizon result follows by induction, independently
of the finite implementation controls.

## Physical scope

A common energy conversion is unnecessary for dimensionless source-gap
ratios. It is required for a value expressed in GeV. Physical mass ratios
also require a common interacting state, identified spectral observables,
the appropriate dispersion and the particle dictionary. None is supplied
by unit cancellation. A scalar record transfer and a working-state modular
clock can produce different complete readouts, as demonstrated by
[`native_clock_comparison`](../native_clock_comparison/PROOF.md).
Independent sector clock rescalings do not cancel, and dropping no-movement
attempts changes the tested kernel. The charged-lepton controller has no
proved attachment to these particular carrier rates.
