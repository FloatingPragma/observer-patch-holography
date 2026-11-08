# Recorded charge transport on a local window

Fix an admitted fair binary source and one connected retained port set `B`.
The local readout's source attachment supplies every incident E record,
including crossing seams, in the original native attempt order. There is no
extra state inspection. The public global total may resolve an empty/full
ambiguity exactly under the owning readout's feasibility test.

For an oriented seam `(a,b)`, the actual recorded current is `j` and the
endpoint update is `n_a'=n_a-j`, `n_b'=n_b+j`. Extend a rational profile
`phi` on B by zero outside B. Every live branch therefore satisfies

```
A_phi(n') - A_phi(n) = (phi_b-phi_a) j,
A_phi(n) = sum_(v in B) phi_v n_v.
```

Summing over any prefix gives its exact charge balance. For an indicator of
a cut, the change of enclosed occupation is minus the outward boundary
flux. Internal currents cancel. An entire source has zero boundary and
conserves its total occupation. This source charge has no assigned electric,
weak or color normalization.

The statement is an operator identity, including arbitrary untouched
references. Each fine native history is a scaled partial permutation `K_h`.
Write `Delta_h` for its summed profile increment. On every supported input
configuration, `A_phi K_h = K_h A_phi + Delta_h K_h`. Consequently

```
K_h* A_phi K_h = K_h* K_h A_phi + Delta_h K_h* K_h.
```

The effects `K_h* K_h` commute with the diagonal `A_phi`. Summing all branches
gives the finite expectation identity. Zero branches contribute zero, and
erasure branches are included. An arbitrary source/reference density matrix
is allowed; no diagonal input assumption enters the operator proof. This is
the discrete continuity identity for occupation, not a stress-energy Ward
identity or a relativistic field equation.

## Complete reconstruction and its error

The local detector infers an initial mask on a decoded branch. Its incoming
contrasts and recorded currents then determine every local snapshot through
the final attempt. Remote events leave that snapshot unchanged. Prefix masks
and the terminal mask use exactly the same clock as the detector.

For decoded label `x`, the readout theorem gives `K_h=K_h P_x`. Attaching the
whole reconstructed local history therefore does not add measurement
back-action. The original terminal source/reference block, FIFO and archive
are retained. Coverage and constant erasures remain explicit labels; no
trajectory is invented for them. Currents still give occupation differences
on live erasure branches, without fixing the initial occupation.

Let D be this actual complete history instrument. Let O first measure initial
local occupation and then run the same native instrument, attaching the
corresponding history. The original readout theorem and the orthogonal
erasure labels give

```
d(D(rho), O(rho)) = p_coverage + (1-p_coverage) eta_amb
```

on its stated independent-selector domain; the global empty/full singleton
exception has distance zero. All references and terminal disturbance are
included. The label retains the initial mask, so the equality survives this
history attachment. A further lossy postprocessing gives an upper bound.
There is no factor equal to the number of snapshots: one complete
instrument comparison controls their joint distribution. Composing with
the local preparation gives the owning `2^-b + gamma` bound on its stated
comparison output. That comparison traces the preparation mixing-record
values; those values cannot be used as extra feedback under this bound.

For public total two on the growing fair source, a twelve-port cell's full
configuration is infeasible. Eleven fixed E reads resolve every feasible
local number, including zero. After 64 preparation blocks, 3862 native
attempts supply the complete local charge history during the final eleven
detector attempts within `2^-64` of its conditional mixed reference. The
comparison does not reconstruct the marginalized preparation mixing history.
No global canonical preparation is inferred.

The exact current law does not label individual identical units through a
collision, produce a stable localized particle, or establish spin,
statistics, mass or a physical energy scale. The actual coherent field
operations and a same-state relativistic dictionary need separate maps.
The evaluator's fixed-size raw-word query and its expanded resource stops
retain their owning meanings.
