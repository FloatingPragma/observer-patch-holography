# M10.5 · Assess free stability of the one-load Interval Echo excitation

Requested decisions are adoption of the one-load state and its home-port detector at proposed OPH pin `69ad9f5eea6fa9b2553f08c4e49e5303b83114c9`, and full row 4 credit for that candidate. The current M10.4 pair remains partial. Its eleven written seams were selected with the stability target known. The M10.5 release has one material load, blank seam memories and all 630 contacts open on the scored board. The port filter was designed with the localization target known. [M10.4 decision](https://github.com/openwave-labs/openwave/pull/749#pullrequestreview-6060088147) · [Public candidate physics and executable view](https://github.com/FloatingPragma/observer-patch-holography/blob/09f062f77f7f540f2783f128ae6314c5d8b818ea/model/candidate_laws/interval_echo/openwave/row_04_compact/gateway_view.ipynb)

A cell starts with one sharp load. The observer prepares a coherent excitation by querying its three gateways and stops on acceptance or exhaustion. During the subsequent free flight the two routes from opposite-sign amplitudes can cancel at exits. The native contrast contributes a common phase at every coupling. The charge, mass and species of the excitation are unassigned.

The unconditional all-time home probability is at least 0.999940596796 on all 420 cells and for every complete sweep order. The 20-cell development record covers all 11,520 possible confirmation settings. At the historical cell-13 exposure of 37, home recovery is 0.999941, falls to 0.570144 after dephasing, and is 0.151119 under random swaps. The platform script reproduces the four complete historical distributions from the public board. The task was designed with these outcomes known. The first official commit containing its M10.5 text selects a separate future comparison of all 56 bins. Its result is pending. The proof and development measurements are presented for the scope they establish.

**Proposed cell wording**

A gateway-prepared one-load excitation on blank memories persists in free flight. Across 420 cells its unconditional home-detection probability exceeds 99.994% at all exposures and couplings. The bounded filter counts exhausted attempts. On the 20-cell pilot at exposure 37, dephasing gives 57.01% home detection and random swaps give 15.11%. Species is unassigned.

**Reproduction.** Use the task, method note, raw-board platform script and public board linked in the patch. The exact source and tagged candidate accompany the local owner packet. The proposed change requests row 4 only.
