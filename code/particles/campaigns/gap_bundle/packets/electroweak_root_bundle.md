# Electroweak Root Bundle Packet

Close or obstruct the coupled package: compressed `P` trunk, Ward-projected
Thomson endpoint, RG/matching/scheme clauses, interval certificate, and
certified code adoption.

Return `closure_candidate` only with a source-emitted `Delta_Th(P)` object,
matching clauses, interval certificate plan, and exact follow-up surfaces.
Otherwise return `blocked` with the minimal obstruction list.

Measured `alpha(0)` must not enter the solver.

## Worker Result

`constructive_endpoint_contract_emitted`.

The endpoint package closes #223 as blocker isolation. The remaining root
package is #235/#32/#224, so the result is a constructive local contract rather
than a no-go.
The compressed `P` trunk can remain a candidate/audit artifact, and the next
implementation target is now explicit:

```text
code/P_derivation/runtime/thomson_endpoint_contract_current.json
```

That artifact requires the same D10 source branch to emit a theorem-grade
Thomson transport object:

```text
T_Th(P) =
  (
    Delta_Th(P),
    Sigma_Q(P),
    C_match(P),
    C_quad(P),
    C_I
  )
```

with

```text
alpha_Th^-1(P) = a0(P) + Delta_Th(P)
Delta_Th(P) = Delta_lep^src(P) + Delta_had^src(P) + Delta_EW^src(P).
```

The admissible closure must include:

- a Ward-projected `U(1)_Q` scheme `Sigma_Q(P)` shared by `a0(P)` and the
  zero-momentum endpoint;
- source-emitted charged-lepton transport, not charged physical values imported
  as endpoint data;
- a source-emitted hadronic spectral density `rho_had(s;P)`;
- either a proof that `Delta_EW^src(P)` vanishes in the source scheme or a
  certified bound for it;
- OPH-derived RG matching, threshold placement, decoupling, and scheme
  conversion;
- quadrature/remainder certificates for every transport integral;
- an interval-level fixed-point certificate for the final closed map.

The current implemented continuation,

```text
Delta_impl(P) =
  Delta_lep^1loop(P)
  + screening(P) * Delta_quark_free^1loop(P),
```

is useful as a reproducible diagnostic, but it is not an admissible theorem-path
replacement for `Delta_Th(P)`.

## Constructive Work Items

- Populate `sigma_q_scheme_lock`.
- Populate `delta_lep_source_transport`.
- Populate `rho_had_spectral_measure` from the hadron spectral-measure schema.
- Populate `delta_ew_remainder`.
- Populate `full_endpoint_interval_certificate`.

Any future worker response on this packet must return one of those populated
objects, a code/schema patch that moves one object forward, or a runnable local
target.  Obstruction-only output is not an accepted result.
