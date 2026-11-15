# Same-group native Coulomb spectral-limit review

Root independently read the exact native controller proof and the full shell
limit, checked all frozen leaves, and requested the three corrections now
present: positive-band finite time error below the binding scale, the dual
Birman-Schwinger spectral identity, and a precise eigenvalue/multiplicity
scope. This is a written mathematical review with no new native simulation
or large-shell eigenvalue run. It supplies no G6 pass.

The source conditional spatial jump probabilities are exactly the six-neighbor
walk. Its Green function has the quoted normalization1/6 relative to the
discrete-time walk, so the leading coefficient is1/(4pi). I checked the
primary Lawler-Limic Theorem4.3.1 and the hit identity in the linked author
book (https://www.math.uchicago.edu/~lawler/srwbook.pdf). The maximum-principle
boundary error divided by delta is O(1/m) for delta=m^-2,B=m³. The outer
Bohr box then expands rather than staying fixed.

The discrete periodic Sobolev inequality has the stated D^-2 remainder
under trilinear interpolation and finite-dimensional norm equivalence. The
local L^(3/2) norm of C/(r+delta) is O(rho), including the origin cell's
O(delta) contribution. This gives uniform infinitesimal form bounds and
excludes an uncontrolled deep defect at the native target. On fixed annuli,
free-resolvent consistency, local H1 compactness and self-adjointness give
norm convergence of localized resolvents. The disjoint inner/outer pieces
are bounded in the resolvent form norm, with cross terms controlled by
Cauchy-Schwarz. Thus the full compact Birman-Schwinger limit holds. The
positive nonzero spectra of A*A and AA* identify eigenvalue counts at each
fixed negative threshold. This proves every fixed shell and its spatial
n² multiplicity while allowing finite cubic splitting and levels near zero.

The radial Coulomb equation and its regular square-integrable solutions give
-beta²/(4n²), beta=lambda/(4pi G0). No measured coupling is used. The
same-source standing energy is E_a(s), and binding is correctly measured
from E_a(lambda). Replacing that reference by sqrt(lambda) would add a
same-order lattice shift, which the current proof explicitly retains.
The positive native band has S_a>=lambda/2 eventually, so H=2m³+1 makes
both temporal image and first-hit energy errors exponentially small,
stronger than o(a²). The massless uniform clock estimate alone would not do so.

The capped history construction has finite charged costs. Incident probes
have squared-distance drift3/11, expected absorption at most11B², and an
iterated Markov bound2^-j. Native relocation contributes at most6B route
edges per probe. The stated caps and union bound couple finite histories
to the exact-preparation marginal with arbitrary prescribed error. Choosing
this error after each finite moment bound avoids an invalid infinite-word
exchange of limits.

No blocker was found in the stated spectral theorem. It is a nontrivial
source interaction and standing-level derivation, with an explicit resource
family. It does not identify electron/proton charge or a laboratory Rydberg
constant, and the calibration is not an executed large-shell production
experiment. Those scope limits are preserved in the packet.
