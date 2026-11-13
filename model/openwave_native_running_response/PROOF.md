# Native current response and its scale dependence

The source is the complete 90-CAR Higgs-Wilson construction of
`openwave_single_valley_leptons`, with all seven native charged pairs in
all three families, the retained neutral modes, and its filled-negative-band
reference. It is distinct from the interacting 192-CAR transverse member.
The electron fixes the energy unit. Native pair charges are taken from the
actual invariant tensors, not supplied species labels. Their squared sum
per family is exactly 8/3. The neutral complement has zero electromagnetic
current. The complete Nambu trace includes its factor one half.

For each exact native Dirac block write

    h(k)=alpha.sin(ak)/a+beta[m+sum(1-cos(ak_i))/a].

A native electromagnetic Peierls probe couples through the conserved charge.
The uniform y-current is q[alpha_y cos(ak_y)+beta sin(ak_y)]. Its operator
norm is at most |q|. The full current on V sites has the conservative Fock
norm bound 90V. The actual 180-dimensional matrices independently reproduce
all occupied-to-empty transition sums; no pair or charge-conjugate branch
is discarded. The full symbolic charge and tensor identities supply the
reduction for every momentum.

## Gauge contact term

For a pure-gauge spatial profile theta, the entire finite source transforms
by H(epsilon)=exp(i epsilon Q theta)H exp(-i epsilon Q theta). Thus
J=i[G,H] and T=-[G,[G,H]]. The exact static energy derivative is

    Tr(P_negative T)-2 sum_(a occupied,b empty)|J_ab|^2/(E_b-E_a)=0.

This follows either by differentiating the unitary identity twice or by
inserting its complete eigenbasis. The executed 256-dimensional real-space
Dirac block checks both nonzero terms independently. Dropping the contact
term fails. The native block intertwiners transfer this identity to every
charged channel and its conjugate; neutral terms vanish. The subtraction
below cancels frequency-independent contact terms but does not pretend
that they were absent from the source.

## An operationally specified coupling

The Gaussian reference has exact retarded current spectral measure. At
imaginary frequency nu, its even response is
2 sum |J_ab|^2 Delta/(Delta^2+nu^2), Delta=E_b-E_a>0. Two subtractions
remove the static response and specify the zero-frequency electric
normalization. The resulting exact quadratic-response coefficient is

    Pi_R(nu)=Vphys^-1 sum 2|J_ab|^2 nu^2/[Delta^3(Delta^2+nu^2)].

The full Nambu expression carries one half. Every summand is positive and
strictly increasing for nu>0 when its transition current is nonzero. The
source has minimum charged pair energy 2m_e. Its small-frequency response
is quadratic; the transition denominators expose the onset scale rather
than inserting a beta function. No localized spatial core is inferred from
this pair-production threshold.

Define the photon-probe quadratic inverse kernel by

    alpha_response(nu)^(-1)=alpha_IR^(-1)-4pi Pi_R(nu).

The supplied alpha_IR fixes the zero-frequency normalization. This is an
explicit finite response scheme, realized by the second derivative of the
Gaussian fermion source with respect to its native probe. Positivity of the
reported inverse kernel is checked. It is not the all-orders interacting
scattering coupling of the separate transverse construction. In particular,
resumming this quadratic coefficient is an interpretation of a declared
probe kernel, not a proof that higher current cumulants vanish.

For one charge-one Dirac block, P_plus/minus=(I plus/minus h/E)/2 gives

    S=Tr(P_minus j P_plus j)
     =2[1-(p_y cos(ak_y)+M sin(ak_y))^2/E^2],
    Pi_R=Vphys^-1 sum_k S nu^2/[4 E^3(4E^2+nu^2)].

The code evaluates the identical positive expression
2[p_x^2+p_z^2+(p_y sin(ak_y)-M cos(ak_y))^2]/E^2, which avoids subtraction
near a current eigenvector. Complete sums retain every Fourier address.
Fixed physical box and varying cutoff are separate from fixed-spacing
volume comparisons. These are deterministic finite sources; their numerical
variation is reported without a certified continuum or volume extrapolation.

A direct continuum radial integration of the same spectral formula yields

    Pi_R(nu)=1/(2pi^2) integral_0^1 dx x(1-x)
                   log[1+nu^2 x(1-x)/m^2].

The radial and parameter integrals are independently executed. The
small-frequency coefficient is 1/(60pi^2 m^2) per unit-charge Dirac field.
This derived comparator is not inserted into the finite lattice sum. The
full native curve contains colour partners with the same family masses as
the leptons; those masses are not observed quark masses. Consequently it is
not a physical QED or QCD running-coupling prediction.

## Complete current records

The literal instrument J/(90V), with its positive square-root failure
complement, acts on the supplied Gaussian reference. Reflection of k_y
makes its vacuum current mean exactly zero. The success vector is the
complete native one-pair superposition. Resolving momentum and native pair
channel, while retaining its quantum spin subspace, gives the recorded
S times q^2 weights divided by (90V)^2. Every native pair, family, momentum,
neutral zero output and failure is retained. The spectral response is the
weighted sum of those transition records, not a fitted scale curve.
No physical photon emission, continuum scattering or Yang-Mills mass gap
is supplied by this quadratic-response experiment.
