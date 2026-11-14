# Independent review of the finite modular-correlation experiment

**Outcome: mathematical and numerical checks pass, with scope notes below.**
The independent audit script is `review_modular_time.py`; its receipt is
`INDEPENDENT_REVIEW.json`. It uses the closed-form projectors of a two-by-two
Hermitian matrix, without importing the producer's spectral-flow or matrix-log
implementation. No producer source was edited.

## Conditioning, algebra and faithfulness

The experiment declares a two-port corner algebra with identity P and conditions
on its inside outcome. Its probability plus the retained complementary probability
agrees with the full native trace. This is postselection, explicitly named in the
method; it is neither a partial trace nor faithfulness of the full 87-position
state. The displayed two conditional eigenvalues are positive. No floor, Gibbs
preparation or hidden pseudoinverse appears.

Strict positivity is numerical: the matrix-free native propagation has no
certified enclosure. Its minimum conditional eigenvalue is well above the stated
threshold, but this audit supplies no independent interval theorem for the native
output. The source-bound `native_state.npz` contains the full 87-position density,
P and Q, both unnormalized terminal branches, their nonselective sum and the
conditional M_2 density. The independent review verifies positivity, trace,
projector/complement identities and every branch against P rho P and Q rho Q.
This resolves the previous missing complementary-state custody. The terminal
P/Q instrument is computed separately from the native evolution: its disturbance
must not be presented as a sequence of measurements during that evolution.

## Modular conventions and KMS sign

Column vectorization gives L_A=I tensor A and
Delta=(rho^-1)^T tensor rho. Direct action on a generic complex matrix agrees
with X -> rho X rho^-1, and vec(sqrt(rho)) is fixed. These checks independently
confirm the finite standard-form/GNS convention. For
sigma_s(A)=rho^(is) A rho^(-is), K=-log(rho) gives
sigma_s(A)=exp(-isK) A exp(isK).

The minus-i KMS shift is correct:

```text
Tr(rho A sigma_(s-i)(B)) = Tr(rho sigma_s(B) A).
```

The review tests this with distinct generic complex A and B at three real
parameters. The plus-i shift fails. This is stronger than checking only the
selected Pauli observable squared. The producer's real-flow sign and closed
correlation formula agree with this convention. A separate Gibbs interpretation
would give physical Heisenberg time t=-beta*hbar*s; none is established here.

## Observable, correlations and interpretation

The state-dependent eigenframe sigma_x is a declared apparatus choice in the
full corner matrix algebra. It is not one of the raw classical C/E letters.
The complex function Tr(rho A sigma_s(A)) is an algebraic quantum correlation,
not a classical joint probability for two projective measurements. A direct
measurement would require a specified phase-sensitive correlation instrument
and its backaction/reference records. The method does not execute such an
instrument or identify an autonomous subjective clock.

The state is invariant under its modular automorphisms while noncommuting
observables evolve. A commuting state observable and a maximally mixed state
have trivial flow as stated. Restricting a singular pure state to M_1 changes
the algebra and removes off-diagonal observables; it cannot preserve the M_2
experiment. These distinctions are correctly represented in source/tests.

## Plot and source audit

All 161 CSV values match the JSON and the independently evaluated analytical
flow/correlation within the reported numerical errors. Current plot-source pins
and native source pins match. The PNG has legible axes and explicitly labels
s as dimensionless, with no physical-time conversion. The plot is an illustration
of this conditional algebraic flow, not an image of experienced time or physical
motion.

Two small presentation/provenance recommendations are communicated separately:
the producer corrected METHOD's clock wording to "declared native evolution
parameter"; store numerical library/environment metadata plus a semantic
CSV/JSON consistency guard in any maintained plot regeneration path. The latter
guard passes for the present artifacts; the recommendation concerns future drift.
No source changes are made by this review.
