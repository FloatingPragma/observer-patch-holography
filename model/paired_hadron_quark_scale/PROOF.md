# Conditional hadron scales and bare Wilson coordinates

The source consists of the retained 32 paired stream-deletion tables for
the quenched 12³×24 lattice, the fixed full-covariance hadron fits, and the
original Wilson-potential fits. No gauge configuration is reconstructed.
Every retained fit window, sink and admitted chiral variant is reused
without an additional target-dependent rule. The lightest unresolved nucleon row
remains in the original fit table; its inclusion is determined by the
parent's pion admission rule, not its agreement with a desired mass.

For every deletion j, the same label j is used in all hadron channels,
the critical and ratio-tuned hopping parameters, and the potential.
The entire chiral chain is recalculated from the frozen fitted masses.
It agrees with every central value of the parent receipt. Writing R for
the fitted ratio-tuned vector mass, S for sqrt(a² sigma) and L for r0/a,
the three conditional inverse spacings are

    A_rho = m_rho_input/R,
    A_string = sqrt_sigma_convention/S,
    A_r0 = L hbar_c/r0_benchmark.

The first fixes an empirical mass input. The second uses the archived
0.440 GeV string convention. The third is an explicitly supplied 0.5 fm
benchmark, with results proportional to 0.5 fm/r0_physical for a different
choice. Neither potential convention is asserted to be an exact measured
physical truth. The retained hbar_c decimal only converts units. External
anchor errors, their covariance and conversion rounding are not propagated;
all statistical errors are conditional on these declared inputs.

Scale-free diagnostics R/S, RL and LS and their paired covariance are
reported before these assignments. Comparing A_rho with A_string or A_r0
therefore tests the compatibility of these assignments on this retained
quenched source. Their differences are not a physical discovery significance.
The finite-volume, cutoff, sea-quark and possible within-stream correlation
uncertainties are not bounded by these tables.

For a Wilson hopping coordinate kappa define

    a m_q(kappa) = (1/kappa - 1/kappa_c)/2,
    m_q^(scale) = A_scale a m_q.

The receipt gives this coordinate at every original hopping value, the
ratio-tuned value and kappa_c. Its zero at kappa_c is definitional. It is a
bare action coordinate, not a renormalized MSbar mass, a quark pole mass or
an identified up, down or strange particle. A mass renormalization factor,
scale and scheme matching and a physical source identification are absent.
The ratio-tuned parameter also obeys exactly within the declared regression

    a m_q(kappa_phys) = (a m_pi_phys)²/(2 slope_mpi²).

This cancellation makes it necessary to refit kappa_c and kappa_phys
together in each deletion. Fixing kappa_c while varying kappa_phys is an
explicitly executed invalid covariance control.

For a vector f_j of 32 derived deletion values the covariance is

    C = (31/32) sum_j (f_j - mean(f))(f_j - mean(f))^T.

Thus Var(A-B)=C_AA+C_BB-2C_AB. The code checks this identity directly and
compares it with ignoring correlations or cyclically mispairing one scale.
Every replica is retained, with no repair of inconvenient signs or values.
All crossed hadron/potential central choices supply a separate maximum
absolute displacement from the primary value. This fit-sensitivity allowance
is not a confidence interval and is not silently added in quadrature to
statistical uncertainty. The full vectors and covariance remain available
for consumers rather than only rounded individual errors.
