# Same-action full-Fock reference and modular observations

This attachment prepares a faithful reference of the **same four-mode
Hamiltonian** used by `defect_field_energy` and the selected source factors
of `defect_field_bridge.joint`. It then applies the actual complete native
repair and field histories to that reference. It supplies finite quantum
entropy and regional-observation inputs to the Einstein premise chain;
it does not construct `EinsteinRegisterData`.

This is a new reference-state experiment using the same complete field
instrument. Its Gibbs preparation differs from the original joint witness's
localized input; that earlier witness is unchanged. The entropy identities
refer to the field and its raw records. The clock factor may be traced
consistently or retained as reported below. No stationarity of the full
field/clock/reference state is asserted.

The declared inverse temperature is beta=1 GeV^-1. On the complete16-state
Fock algebra, rho=exp(-beta H)/Z is faithful and
K=-log(rho)=beta H+log(Z) I. No total-number restriction or low-weight
spectral deletion is made. Adding a common scalar to H changes Z and
leaves rho and K unchanged. The modular convention is
sigma_s(A)=rho^(is) A rho^(-is), corresponding to Hamiltonian time
**t=-beta s**. This Gibbs parameter is not an extra physical clock. Field
intervals still use the common Omega_ref=1/4 GeV source clock and its
lambda_j multiplier.

The existing complex preparation compiler emits an eight-qubit word for
|sqrt(rho)>. All sixteen purification columns remain. Giving the purifier
the dual charge -Q^T makes the ideal complete state neutral because
[Q,rho]=0. Preparation costs and every first-export complement output
remain accounted. This is a declared preparation service, not a claim that
passive repair thermalizes to rho.

The signed full-Fock CAR trace gives the first-site reference rho_A and a
nested one-mode reference rho_a. All local occupation sectors remain and
both references are faithful. The strict bound
Tr_B(rho)>=dim(B)lambda I whenever rho>=lambda I proves this without a
support projection. Existing `ReferenceObservation` implements

    E_A(X)=rho_A^(-1/2) T_A(sqrt(rho) X sqrt(rho)) rho_A^(-1/2).

It is completely positive, unital and reference preserving. Its adjoint
recovers rho from rho_A. The nested maps compose exactly because the
intermediate square roots cancel; the code checks every matrix unit and
an untouched ancillary operator. This is not recovery of every state.

The existing bounded `RecoveryService` compiler also emits a native
four-input,64-output isometry for this same Petz recovery. It retains the
16-state field and four-state environment, all input/reference columns,
and every first-export failure. The preparation and recovery words are
frozen in `source_controls.npz`: they are explicit emitted programs, not
new source primitives. Each replay independently checks their full targets
against the current H and rho. Freezing the words keeps actual costs and
otherwise arbitrary QR phases reproducible across platforms; controls are
regenerated only with `--freeze-controls`.

There is an exact obstruction to calling this observation a local
reference-preserving conditional expectation. In the bare chart the
cross-boundary hopping is (i/2) times an integer CAR matrix. Both partial
traces vanish and its squared Hilbert--Schmidt norm is exactly4. The
actual chart change is a product of local unitaries, so this nonadditive
component remains. For beta>0 a product Gibbs reference would imply an
additive Hamiltonian by taking its logarithm. Therefore the reference is
not a product, and no UCP retraction onto the full first-site factor can
preserve it. The compatible non-idempotent observation remains valid.
Numerical controls retain its retraction and modular-covariance defects,
and compare the actual reduced K_A with beta times the bare local energy.

For each C/E law and each coherent/flagged collision policy, the two-event
history retains all1024 native records and the original purifier. These
entropy values concern the unnormalized **continuing source component**.
First-export abort weights, Grams and reconstructible prefix outputs are
retained separately; their entropy is not included in these numbers.
They have zero weight on the exact proper code and remain subject to the
declared numerical admission. Neither continuing outputs nor tiny branches
are normalized or removed. For a normalized ideal unread output sigma,
full matrix logarithms give

    S(sigma)-S(rho)=beta Tr[H(sigma-rho)]-D(sigma||rho).

For the retained classical-quantum record, the corresponding identity is

    S(CQ)-S(rho)=H(p)+beta DeltaE-sum_h p_h D(sigma_h/p_h || rho).

Conditional normalization appears only inside this mathematical identity;
stored quantum outputs remain subnormalized. Numerical trace drift is
retained explicitly: both identities gain (Tr(sigma)-1)log(Z). The
first-order entropy tangent has the corresponding (log(Z)-1)Tr(delta)
term. The second-order Bregman remainder is
D(rho+epsilon delta||rho)-epsilon Tr(delta), not raw relative entropy
when the computed trace drifts. Its Hessian bound uses the positive lower
bound (1-epsilon)lambda_min(rho). All numerical admissions are distinct
from outward-rounded physical or hardware certificates.

The same final clock outputs are tensor-composed with every native
history, retaining all4096 continuing native/clock records per policy/law.
Within that component, work entropy and entropy including the purifiers
are reported separately. Neither includes the separately stored abort family.
This prevents entropy of an unread field, a retained record and a complete
purification from being silently identified. Reference stationarity under
the driven repaired history is a measured diagnostic, not an assumption.

The new construction closes a finite same-action/reference/preparation
and compatible-observation join. It does not identify a geometric cap
kernel, null stress, a bulk/edge split, fixed-volume area variation,
shrinking balls or a positive measured Newton normalization. The two-site
source factors are not an invariant restriction of the full ring.

```sh
python3 -m model.defect_field_modular --check
python3 -m pytest -q -W error model/defect_field_modular
```

Intentional regeneration uses `--freeze-controls`, then `--write`. Retained
controls and outputs are inspected independently by `--check`.
