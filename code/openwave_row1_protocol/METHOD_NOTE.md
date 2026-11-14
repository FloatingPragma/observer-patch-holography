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

where A is the oriented primal face-boundary matrix and the map is induced by B. Nineteen unit Smith factors establish the full integral image, and a unimodular spanning-tree minor certifies the primitive neutral lattice. The proof does not infer integral equality from real ranks alone. A continuous circle-valued representative winds F_e times along each oriented dual edge. Its contour degree is BF. Shared vertex phases and nonlinear interior homotopies preserve the degree. Half a turn fails closed endpoint gluing; a unit class cannot split into two closed half-unit classes. The separate native charge-two example can split into two unit classes.

For the one-load source, every regional degree belongs to {−1,0,1}. Current-driven counters modulo six decode these values injectively. The meter discovers an unknown one-load pointer using emitted native currents, routes it to a fixed cycle, and completes six signed native hops per nonzero degree. Zero degree requires discovery/reset too. No coherent phase-control service or ideal terminal projector is used.

With R retries per requested hop, write s=1−2^(−R), and let p_x be the exact discovery/reset success probability for initial pointer port x. For source history h with degree w(h),

\[
\Pr(h,\mathrm{success}\mid x)=2^{-|h|}p_xs^{6|w(h)|},\qquad
\Pr(h,\mathrm{timeout}\mid x)=2^{-|h|}(1-p_xs^{6|w(h)|}).
\]

Every complete fine history has weight 2^(−number of actual attempts). Configured source, region and pointer inputs carry no invented prior. At R=16, the tested pointer cap is 624 attempts and the longest tested joined program has 627. The largest unconditional meter timeout over the primary input domain is 0.00013732048715552606. Success probability depends on degree; source statistics must retain timeouts rather than condition them away.

## Equation-to-code map

| Object | Executable or proof in `../openwave_row1_native/` |
| --- | --- |
| Literal endpoint update | `native_kernel.py::native`, exact extracted source ancestry in `source_packet.json` |
| Actual source/current counter | `joined_readout.py::source_prefix` |
| Dual contours and deformed circle representatives | `native_charge_winding.py::dual_boundaries`, `curve`, `winding`, `evaluate` |
| Integral lattice equality | `provenance/SPATIAL_LIFT_PROOF.md`, `SPATIAL_LIFT_AUDIT.json` |
| Discovery/reset/native cycle | `native_cycle_meter.py::execute`, with original body retained in provenance |
| Whole source-to-pointer process | `joined_readout.py::execute_joined` |
| Independent composition and complete grouped laws | `run.py`, `PRIMARY_RESULT.json` |
| Independent joined-runtime reviews | `provenance/JOINED_READOUT_INDEPENDENT_REVIEW.json`, `JOINED_READOUT_REVIEW.json` |

**Scope**: The result is a finite static spatial integer-charge mechanism and a complete native readout on the known-reference one-load source. Carrier geometry, native integer state alphabet, reference preparation, finite observer memory, addressed schedule and auxiliary pointer are declared architecture. The circle-valued representative interprets actual current cohomology; no independently dynamical electromagnetic field is asserted. The pointer cycle is a readout path, not a claim that it is noncontractible on the unpunctured carrier sphere. No coherent-source, physical mass, spin, force, electric normalization or official platform award is inferred.
