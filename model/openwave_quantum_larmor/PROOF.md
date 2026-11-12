# Quantum Larmor dynamics on the same complete source

## Source, current and state

Keep the cap-two complete matter/Maxwell source of
`openwave_quantum_rotational`, including its positive absolute kernel, trap,
reference and clock. A real uniform external field b along z changes

    H(b)=sqrt(16^2+(B+b E_z)^2)+N_matter+N_photon.

E_z is the native derivative of the symmetric-gauge vector potential.
Neither the Gaussian quantum field nor its energy is removed. Each source
has both signed Dirac branches, every Clifford copy and all fermion sectors.
The neutral reference and charge attachment to one complete native block
are specified by `openwave_native_quantum_charge`.

The isolated unperturbed reduced ground band P has energy E0 and exact spin
one-half. The exact derivative X satisfies

    K X+X K=B E_z+E_z B.

Rotation covariance and a pi rotation give PXP=mu sigma_z in its actual
Jz eigenbasis. The inherited static coefficient is g_static=4 E0 |mu|.
The same unit, field derivative, state and source time define the dynamical
coefficient below. No fitted Larmor clock is introduced.

## Uniform finite-time comparison

At cap two, ||B||<=a=5/2, ||E_z||<=e=4/5 and K>=m=16. The Sylvester inverse
has norm at most1/(2m). For field magnitude b this gives

    v=a e |b|/m+e^2 b^2/(2m) >= ||H(b)-H(0)||,
    L2=[e^2+((a+|b|e)e/m)^2]/m >= sup ||H''||,
    r=L2 b^2/2 >= ||H(b)-H(0)-bX||.

Let gamma=17-E0_upper and require gamma>2v. In either exact Jz sector the
unperturbed ground is simple and the next level is at least17. Temple's
bound on the unperturbed ground vector and the perturbation bound imply

    |E_s(b)-E0-s b mu| <= D=r+v^2/(gamma-2v), s=+1,-1.

The perturbed eigenvector can be phase aligned within
sqrt(2) v/(gamma-2v) of its unperturbed line. The two Jz sectors are
orthogonal, so the full ground-space propagator obeys

    ||exp(-itH(b)) P-exp(-it(E0+b mu sigma_z))P||
      <= eta(t)=2 sqrt(2)v/(gamma-2v)+|t|D.

This gives trace distance at most min(2,2 eta) for every ground-space state
and every spectator reference, including entangled input. Applying any
complete instrument cannot increase this distance. There is no photon trace
or erased field outcome in this assertion. The finite-field precession
coefficient and its static comparison are

    g_pre(b)=2 E0 |E_+(b)-E_-(b)|/|b|,
    |g_pre(b)-g_static| <= 4 E0 D/|b|.

The analytic bound vanishes for fixed scaled time |bt| as b tends to zero.

## Tight finite-field branch certificates

For rational b=1/100 and1/1000, numerical diagonalization supplies a trial
only. Its integer dyadic coordinates are projected using the exact spectral
polynomial of Jz with all half-integers from−9/2 through9/2. Dividing by
its enclosed exact norm gives an exact normalized vector in Jz=+1/2 or−1/2.
The finite matrices are evaluated in outward50-digit interval arithmetic.
A sixth-degree polynomial for sqrt(16^2+B_b^2), with positive scalar
remainder at most(33/2048)(a+|b|e)^14/16^13, encloses both its Rayleigh
energy and residual. In either projected sector, H(b)>=16+N_matter+N_photon
leaves only one vacuum line, so the next energy is at least17. Temple's
bound therefore encloses the actual branch energy. Subtracting the two
intervals and using the inherited E0 interval encloses g_pre directly.

## Complete executable records

The fixed compensated source seed undergoes its complete ground-line
filter. Its success weight is retained. On success a total angular rotation
prepares (up+down)/sqrt2. The original filter failure is evolved by the same
H(b), with its norm retained. The output jointly records photon number and
total Jx, which commute. Every joint probability is included, together with
all complex coefficients of the unmeasured joint matter/field success and
failure states. The opposite-charge apparatus reference is idle and retained.
Other native Fock sectors have zero probability for this invariant source
preparation and are not removed from the declared source.

These numerical records illustrate the proven uniform bound; their
floating-point traces are not substituted for the interval energy theorem.
The source is confined and finite in its matter and photon bases. The field
coefficient, mass law and calibration retain their declared ancestry. This
same-source precession result supplies neither a physical Tesla conversion,
a physical electron magnetic-moment error budget nor an infinite-photon
magnetic coefficient.
