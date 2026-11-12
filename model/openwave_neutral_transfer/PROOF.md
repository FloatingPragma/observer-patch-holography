# Neutral spin-half transport on the native four-cell torus

## Source and observation domain

The source is the literal `periodic_family.Configuration(side=4)` and
`peter_weyl_links` Hamiltonian, with spacing π and all five coefficients one.
It has 64 sites, 192 distinct oriented links of dimension119, 1920 matter
modes obeying the canonical anticommutation relations, and five scalar states
at each site. The antiperiodic seam signs in `kinetic_terms` are retained.
The Hamiltonian contains every native Yukawa word, every scalar and link
penalty, and every axial kinetic path. No scalar or link expectation replaces
an operator. The clock is exp(−itH) on this complete carrier.

The preparation is a local three-fermion gauge singlet, blank scalars and
links, and empty matter elsewhere. An idle neutral odd reference carries
opposite raw lepton number. A declared parity-even unitary prepares the
source/reference pair from its vacuum. It is a preparation service, not a
claim about the interacting ground. Measurement has all128 site/spin packet
projectors and their full-carrier complement. Reference outputs are retained.

In the native local ordering let

    a = |10,11,24>, b = |9,10,26>, c = |9,11,25>,
    χ_up = (2a − b + c)/sqrt(6), χ_down = S_minus χ_up.

All12 native gauge currents annihilate both vectors. Their Gram matrix is
I₂, their total spin Casimir is3/4, their electric charge is0, and their raw
lepton number is1. These are actual `e^c LL` states, rather than extra neutral
coordinates. The exact column computation checks these statements.

## Third-order transfer

For one complete link with unit Hermitian hopping coefficient, the exact
matrix from a packet at its initial endpoint to a packet at its final
endpoint is zero for powers0,1,2 and −6σ_a for power3. The calculation uses
all119 intermediate link labels and complete source columns. Each initial
spin vector has three columns; the complete successive supports have sizes
3,15,18,18. No projected intermediate two-state clock is used.

A useful equivalent proof starts with the uncompressed group-coordinate
multiplication operator. Three transfers are necessary to move the packet.
Acting once on a blank link produces only native matrix coefficients and
therefore lies in the retained link carrier. This also holds for the
receiver acted on once. Consequently

    <receiver|(PHP)^3|initial> = <receiver|H^3|initial>

for these minimum-length transport contributions: both outer actions already
lie in P. The fully transferred gauge-singlet polynomial is independent of
the group element. Haar integration preserves the coefficient. The argument
also applies to a path of distinct links because its first action lies in
the tensor product of the retained native coefficients.

For x=(0,0,0), y=(0,0,1), the literal source has two paths joining the
endpoints: the direct link and the three-link route around the torus. With
the receiver orientation, each has coefficient −i sqrt(2)/8. Mixed choices
among the three matter transfers integrate to zero: a proper nonempty subset
of `e^c LL` has no gauge-singlet component, and a link belonging to only one
path sees that nontrivial coefficient. Thus the third moment is the sum of
the two cubes, rather than the cube of their sum:

    <y,up|H^3|x,up> = −6 · 2 · (−i sqrt(2)/8)^3
                       = −3i sqrt(2)/64.

No other full-source term contributes at this order. Three hops must each
connect x to y; an intermediate matter site would require an additional hop.
A single onsite pair term cannot conserve the final matter number. Two
onsite pair terms and one hop would leave nonblank scalar states at the two
endpoints. Penalties cannot change the required matter occupations. These
observations also prove that the complete H⁰,H¹,H² receiver moments vanish.

Removing the antiperiodic seam changes the two oriented coefficients to
opposite signs and cancels this third moment. The boundary condition is an
explicit source choice. No ground, low-energy band, or periodic side-two
flatness theorem is transferred to this member.

## Complete probability bound

Every onsite word is a contraction. On the uncompressed link space, each
path coordinate is a unitary multiplication matrix. The Hermitian hopping
on its two endpoints has Fock norm at most30 times its scalar coefficient;
compression to the native link carrier is contractive. The looser factor60
therefore supplies the conservative bound

    ||H|| <= 64(sum_Yukawa |coefficient| + 1)
             +192+60 sum_paths |coefficient| < 12937 = B.

The spectral theorem gives the full-vector Taylor remainder at degree3:

    ||exp(−itH)ψ − sum_(k=0)^3 (−itH)^k ψ/k!|| <= B^4 t^4/24.

For t=10^−20 the leading receiver modulus is
L=sqrt(2)t³/128 and R=B⁴t⁴/24<L. Its original-input probability lies in
[(L−R)²,(L+R)²], numerically
[1.2181254471,1.2232835274]×10^−124. This intentionally small time certifies
strict positivity using a crude complete-carrier bound. It estimates neither
a transport speed nor a phenomenological cross section. The other128
outcomes contain the entire complement and restore probability exactly.

The finite theorem is a nonzero native neutral transport result. It does
not identify an observed neutrino, derive a relativistic dispersion law,
or establish a physical mass gap.
