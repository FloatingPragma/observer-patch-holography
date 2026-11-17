# Independent canonical heat-trace calculation

Root will import no producer code. Reconstruct the complete zero-mass
L=4, S=2 ring spectrum from its analytic momentum and wall singular values.
At every one of 1024 holonomies compute the full eight-particle canonical
trace as the coefficient of z^8 in the sixteen-factor generating product,
using interval exponentials. This differs from the producer's enumeration
of 12,870 occupations. Retain every one-particle band and Fock contribution.

Consume the frozen Coulomb pointwise ground lower bounds, complete source
number-sector bounds, preparation residuals and earlier Green residual
bounds explicitly. Replace only the free ground contribution in the trace.
Independently bound the shifted-rotor heat trace and the holonomy quadrature,
then check the eigen-count and vacuum-fidelity consumers. The witness is
beta=5/4 with cut offset1/8, fixed before this independent run.

Also test the literal boundary-current derivative: its one-body spectrum
is expected to be +8 twice, -8 twice and zero twelve times. If verified
exactly, the complete N=8 matter derivative has norm16. Report this
stronger quadrature bound separately from the producer's16 sqrt(2) bound.
Neither is a source parameter or a physical mass choice.
