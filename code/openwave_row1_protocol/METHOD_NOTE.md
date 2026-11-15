# Current-derived spatial charge and native readout

For an oriented seam e=(a,b), let B_e have entries +1 at a and −1 at b. The literal native endpoint event returns current j and updates

\[
n'=n-jB_e,\qquad F'=F+je,\qquad BF=G-n.
\]

Here G=n(0) is retained and F(0)=0. Thus regional spatial degree and relative charge satisfy

\[
w_A=\sum_{v\in A}(BF)_v=G(A)-n(A),\qquad q_A=-w_A.
\]

On the actual spherical carrier, its dual graph has twenty vertices, thirty edges and twelve pentagonal boundaries. The integer quotient of edge labels by dual vertex coboundaries is

\[
\mathbb Z^{30}/\operatorname{im}A
\cong\{q\in\mathbb Z^{12}:\sum_vq_v=0\},
\]

where A is the oriented primal face-boundary matrix and the map is induced by B. Nineteen unit Smith factors establish the full integral image, and a unimodular spanning-tree minor certifies the primitive neutral lattice. The proof does not infer integral equality from real ranks alone. A continuous circle-valued representative winds F_e times along each oriented dual edge. Its contour degree is BF. A positive primal current runs from the lower to the higher port; `dual_boundaries` fixes a consistent surface orientation and checks that signed dual pentagons reproduce each row of B. Shared vertex phases and nonlinear interior homotopies preserve the degree. The equivalence is free homotopy: shared vertex phases may vary during it. The domain is fixed, while the faces of nonzero degree change when the load moves. This quotient is for one carrier. Glued carriers require their specified additional cells and can retain independent cycle periods.

The class is exactly the class of the native current record on this dual graph. The phase map is its chosen observer encoding. The contour diagnostic re-reads that encoding; it is not a separate field sensor. The additive law controls fix 32 subdivisions per edge before execution on |F_e|≤3, instead of selecting density from an observed tape. The native pointer reads the counters described below and does not use phase samples.

Half a turn fails closed endpoint gluing; a unit class cannot split into two closed half-unit classes. This integrality does not by itself establish that the source preserves one physical unit. That separate finite statement follows from the nonnegative integer endpoint rule, which maps a total-one state to another total-one state. The executed half-current law retains BF=G−n while its phase curves fail to close. An executed random-integer law retains integral degree while violating whole-unit occupancy and regional bounds. The separate native charge-two example can split into two unit classes. These controls distinguish the identity, the integer input alphabet and the nonfragmentation property of the declared law.

For the one-load source, every regional degree belongs to {−1,0,1}. Current-driven counters modulo six decode these values injectively. The meter discovers an unknown one-load pointer using emitted native currents, routes it to a fixed cycle, and completes six signed native hops per nonzero degree. Zero degree requires discovery/reset too. No coherent phase-control service or ideal terminal projector is used.

With R retries per requested hop, write s=1−2^(−R), and let p_x be the exact discovery/reset success probability for initial pointer port x. For source history h with degree w(h),

\[
\Pr(h,\mathrm{success}\mid x)=2^{-|h|}p_xs^{6|w(h)|},\qquad
\Pr(h,\mathrm{timeout}\mid x)=2^{-|h|}(1-p_xs^{6|w(h)|}).
\]

Every complete fine history has weight 2^(−number of actual attempts). Configured source, region and pointer inputs carry no invented prior. At R=16, the tested pointer cap is 624 attempts and the longest tested joined program has 627. The largest unconditional meter timeout over the primary input domain is 0.00013732048715552606. Success probability depends on degree; source statistics must retain timeouts rather than condition them away.

An omitted moving record changes the decoded class. The omission is detected by comparison with the source loads, which the controller does not read. The readout has no internal check for a lost record: a damaged integer tape still defines closed phase contours with integral degree. `law_controls.py` executes this distinction on every affected development history. Original source, receipt and review bytes remain frozen; this method note supplies the corrected interpretation.

## Equation-to-code map

| Object | Commit-pinned executable |
| --- | --- |
| Literal endpoint rule | [native](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/native_kernel.py#L2-L9) |
| Native source and current counter | [source_prefix](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/joined_readout.py#L18-L42) |
| Dual cell boundaries | [dual_boundaries](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/native_charge_winding.py#L34-L96) |
| Continuous representative | [curve](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/native_charge_winding.py#L99-L112) |
| Contour degree | [winding](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/native_charge_winding.py#L115-L120) |
| Spatial evaluation | [evaluate](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/native_charge_winding.py#L123-L232) |
| Native discovery/reset/cycle | [execute](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/native_cycle_meter.py#L45-L120) |
| Joined source and pointer | [execute_joined](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/joined_readout.py#L45-L76) |
| Complete primary execution | [execute](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/run.py#L34-L127) |
| Integral lattice equality | [Exact proof](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/provenance/SPATIAL_LIFT_PROOF.md), [arithmetic audit](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/provenance/SPATIAL_LIFT_AUDIT.json) |
| Independent full-join reviews | [Phase/controller review](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/provenance/JOINED_READOUT_INDEPENDENT_REVIEW.json), [lattice/implementation review](https://github.com/FloatingPragma/observer-patch-holography/blob/19288e2b53b680e244e8ec5f4c5d45e000b2ef51/code/openwave_row1_native/provenance/JOINED_READOUT_REVIEW.json) |

**Scope**: This finite known-reference one-load construction supplies current-derived spatial integer charge and complete native readout, with carrier, observer memory, addressed schedule and auxiliary pointer declared. The circle-valued representative and pointer cycle carry the specified diagnostic meaning; no independently dynamical electromagnetic field, coherent-source transfer, physical mass, spin, force, electric normalization or official platform award is inferred.
