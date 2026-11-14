# The discrete macro and autonomous Dirac flow have different evolutions

The calibrated L3 source8a64e359 applies ordered full spatial shifts,
followed by its field, Coulomb and clock factors. The separately declared
autonomous Dirac source44c00889 uses a continuous generator with signed
nearest-neighbour differences. These constructions do not have the same
evolution, even when their clock and mass coefficients are matched.

Consider the admitted neutral charge sector and momentum
theta=(2pi/3,0,0). The field evolution and independent clock factor out.
Let alpha and beta be anticommuting Hermitian involutions on the four spin
coordinates, and let h=1/3. The discrete spin operator is a product of
exp(-i theta alpha) and exp(-i h m beta), in the source's declared order.
Its normalized trace is

    tr(U_macro)/4 = cos(theta) cos(h m) = -cos(h m)/2.

The corresponding continuous signed-difference generator gives

    U_aut = exp[-i (h m beta + sin(theta) alpha)],
    tr(U_aut)/4 = cos(sqrt((h m)^2 + 3/4)).

The trace formulas are independent of the sign convention for the shift
and of the ordering of the two discrete factors. The script encloses both
traces by rational alternating Taylor bounds for each of the three mass
coefficients in the frozen L3 source. Each normalized trace difference
exceeds 57/50. Since |tr(A)|/4 <= ||A||, the operator-norm difference exceeds
1.14 in every flavor. No eigensolver, numerical quadrature or choice of
logarithm branch enters this check. Replacing the macro by the autonomous
formula is an explicit rejected mutation.

This counterexample belongs to the exact neutral reduction of the declared
family. It excludes an identity of those source families. It does not
calculate the charged interacting spectrum or exclude a separate
approximation theorem on a smaller preparation and observation domain.
A Floquet logarithm requires a specified branch and cannot be identified
with the displayed continuous sum by its name.

The canonical harmonic field preparation has Var(p)=1/2 in each of three
independent coordinates. Its harmonic energy H_h=sum(p_j^2)/4 has mean 3/8
and variance 3/32. In physical time at h=1/3 the variance is 27/32. The
supplied normalized harmonic packet is therefore not a free energy
eigenstate. A stationary perturbation calculation requires a suitable
reference or momentum-fiber contract. This variance calculation is not a
nonexistence theorem for interacting stationary states.

PARAMETERS.json retains the exact source mass coefficients, update factors,
clock and source hashes. Run `python3 model/openwave_macro_generator_identity/check.py`
from the repository root. RESULT.json contains exact rational enclosures
and the source-substitution control. Independent matrix checks and review
status are recorded in REVIEW.json.
