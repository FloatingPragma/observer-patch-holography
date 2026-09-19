# Post-run audit: record-metric control families

This addition addresses the issue #776 Exit requirement that the paired
receipt contain both readouts **and their controls**. The original receipt
carried matched W3/isolated controls, but omitted the record-metric family's
existing lower-dimensional scientific controls. Mutation tests establish
verifier behavior; they do not replace those scientific comparisons.

This is a post-run audit addition, not a claim of blind preregistration. The
historical family's results are already public. The original specification,
W12 run, selected intervals and q13 three-dimensional trace remain frozen.
Declare the additional procedure here before generating its new read tapes:

* Run both existing q=13 control families, with one and two spatial axes.
  Use the same golden orbit, exact metric threshold, site ordering, initial
  values `site+1`, self reads, four rounds and `1+sum(read values)` law as
  the three-dimensional comparison. Unused readback axes are zero.
* Retain every actual writer version, integer value and read offset. Choose
  the existing central site and all four lags, with inclusive endpoints.
  Report counts, ordering fractions, MM inversions, longest-chain heights,
  count ratios, spatial summaries and exact clipping flags at every lag.
* Compare exact interval counts with the corresponding archived family.
  The flat-space reference fractions are 1/2 and 8/35 (MM dimensions 2 and 3).
  They are comparison values, not acceptance thresholds. Do not retune a
  horizon, seed or interval to make the measured dimensions match them.
* Independently reconstruct each metric graph and replay each complete tape.
  Require both families in the paired receipt and its causal-poset mirror,
  and include all additional arrays and readouts in fresh CI reproduction.

The q13 lag-four geometric diamond is clipped in all three spatial families;
the flag must be computed and shown, rather than suppressed for a control.
No control establishes a continuum limit or selects the supplied read law.
