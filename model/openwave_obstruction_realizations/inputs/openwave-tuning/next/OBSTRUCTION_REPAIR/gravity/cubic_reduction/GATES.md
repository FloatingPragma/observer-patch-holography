# Exact cubic reduction gates, before large-ball execution

The small-ball full-versus-reduced test must pass before the large run.
All copies of G1 and its two-dimensional multiplicity are retained, including
both parities. Excitation support outside G1 must vanish to numerical error;
only then do untouched irreps cancel from the entropy difference.

The large run compares ell/R4 and6 on the same improvedR1.5,L96 source,
with the same complete neutral preparations at widths27 and36. Actual
energies, box ratios and profile widths are retained. Desired evidence is
decreasing three-dimensional modular error as ell/R increases, checked on
both preparations and under eigenvalue-floor variation. No full modular
operator identity or joint gravity stationarity follows from one ratio.

No matrix exceeding5000 by5000 is allocated. The unreduced large ball is
never assembled. Covariance blocks are built orbit by orbit; its reduced
matrix is saved before preparation and entropy work. Each run remains below
five minutes CPU; if needed, later stages resume the saved matrix.

The G1-only reduction applies to this specific matter change. Geometric
entropy variations generally affect all irreps and must not use it without
a separate proof. The positive-energy periodic mean mode remains open.
