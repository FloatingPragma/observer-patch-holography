# Source correlations bounded by the interactions across a cut

Use the unchanged [primary construction](../particle_common_contract/PROOF.md),
its physical Gauss carrier, selected polynomial density rho, retained
canonical reference and Hamiltonian clock. The region contains its matter
and every link with both endpoints inside it. Crossing links belong to
the complement. This is the [induced regional cut](../particle_primary_geometric_cut/PROOF.md),
including its boundary representations and supported partial traces.

The argument adapts the finite thermal mutual-information comparison of
[Wolf, Verstraete, Hastings and Cirac](https://arxiv.org/abs/0704.3906).
The primary state is a specified polynomial, not an exact Gibbs density.
The actual polynomial error and the physical boundary pairing are retained.
No regional modular generator is identified with a geometric energy.

## 1. A physical comparison with the same regional marginals

On the actual physical carrier write the Schur decomposition and source
blocks as in the [boundary information identity](../particle_primary_boundary_information/PROOF.md):

    H_phys = direct_sum_lambda M_R,lambda tensor M_out,lambda
                                 tensor span(Omega_lambda),
    Pi_lambda rho Pi_lambda = p_lambda tau_lambda tensor
                                      |Omega_lambda><Omega_lambda|.

Omega_lambda is the normalized invariant representation pair of dimension
d_lambda. The actual rho may have coherence between different lambda.
Its positive physical support makes every populated multiplicity block
faithful. Let tau_R,lambda and tau_out,lambda be the two marginal densities
of tau_lambda, and define the comparison density

    sigma = direct_sum_lambda p_lambda tau_R,lambda tensor tau_out,lambda
                                      tensor |Omega_lambda><Omega_lambda|. (1)

This density is normalized and physical. It retains exactly the original
center weights, inside marginal and outside marginal. It removes conditional
multiplicity correlations and global center coherence, while preserving
all mandatory representation-pair entanglement. Replacing the singlet by
two independent maximally mixed representation factors would leave the
physical Gauss space and is not used.

Equation(1) is a variational comparison state. The map from an unknown rho
to its product marginals is nonlinear; no universal channel implementing
that map is asserted. The actual experiment keeps its original preparation
and reference. All entropies here concern specified source marginals, not
the complete pure source/reference output.

Put J=sum_lambda p_lambda I(tau_lambda) and
C=D(rho||sum_lambda Pi_lambda rho Pi_lambda). The inherited entropy split
immediately gives

    X := S(sigma)-S(rho) = J+C >=0,
    I_ambient = H(p)+2 S_edge+X.                               (2)

There is no assumption of zero conditional correlations or zero center
coherence. The source already supplies a strictly positive finite C for
the two-site cut.

## 2. Exact polynomial free-energy comparison

The [selected polynomial slope](../particle_primary_gravity_premises/PROOF.md)
gives, on the whole physical carrier,

    K=-log rho = beta H+k I+R(H),
    0<=R<=r I,    r=beta delta(U-L),                          (3)

where every parameter is reconstructed from that same member. In
particular the original density, not an approximate thermal state,
appears on both sides. For the comparison(1), relative entropy gives

    D(sigma||rho)
       = beta Tr[(sigma-rho)H]+Tr[(sigma-rho)R]-X >=0.         (4)

Because both densities are normalized and R has range[0,r], the absolute
value of its expectation difference is at most r. Decompose the actual
Hamiltonian as H=H_R,int+H_out,int+V_cut. The first two terms belong to
the respective regional gauge-invariant algebras. Their expectations in
sigma and rho coincide by(1), hence

    X <= beta Tr[(sigma-rho)V_cut]+r
      <= 2 beta ||V_cut||+r.                                (5)

The energy difference in(4) is not a unitary cyclic-work experiment and
is not asserted positive. Equations(4)-(5) require neither a Gibbs
replacement nor an identification of K_R with beta H_R,int. The positive
work theorem for recording the center concerns a different, complete
pinching intervention.

## 3. Actual source boundary terms

For the induced factorization, every electric operator acts on one link
and belongs entirely to one regional algebra, including electric terms on
crossing links assigned outside. Every onsite term is also wholly inside
or outside. A Wilson hopping term belongs to V_cut exactly when its two
matter endpoints lie on opposite sides. A magnetic plaquette belongs to
V_cut exactly when it contains both an internal-region link and a link
assigned outside. A plaquette with all its links assigned outside is an
outside operator even if a vertex belongs to the regional cell set.
These statements concern the operator supports, not their energy-density
shares. Every term remains compressed to the actual cutoff and Gauss space.

Let b count crossing bonds and p_mix count mixed plaquettes. Full-CAR
hopping norm is at most30/a per bond; the signed magnetic term has norm
at most1/(7a). Thus

    ||V_cut|| <= (30 b+p_mix/7)/a,
    p_mix <=2b,
    ||V_cut|| <= (30+2/7)b/a.                                (6)

Every mixed plaquette has at least two crossing geometric edges. Every
crossing edge belongs to at most four elementary plaquettes; counting
incidences proves the middle inequality. The finite code enumerates all
mixed plaquettes adjacent to the actual inside links. Independent full
box enumeration tests include disconnected and nonconvex cell sets.
The two-site cut R={0,e1} has p_mix=2 at every admitted side, with b=4
at member zero and b=5 at larger members. Electric coefficients and their
cutoff dependence have not been hidden in a crossing norm.

## 4. A finite boundary-information bound

The [actual link capacity bound](../particle_primary_capacity/PROOF.md)
is dim H_link<=451^ell. The tensor product of the b crossing-link spaces
carries the outside action of all regional boundary gauge groups. Each
allowed tuple lambda must occur in this representation, by the physical
Gauss pairing. Its number is therefore at most451^(ell b). The inherited
representation dimension bound is d_lambda<=6^(ell b). Consequently

    H(p)<=ell b log451,    S_edge<=ell b log6,
    I_ambient <=ell b(log451+2log6)
                     +2 beta(30b+p_mix/7)/a+r.              (7)

The bound is uniform in the region's volume at fixed source parameters;
it has the explicit additive polynomial correction r. On the lattice
area b a², its boundary coefficient is bounded by

    ell(log451+2log6)/a² +2 beta(30+2/7)/a³.                 (8)

The upper coefficient(8) diverges along the chosen refinement; this
does not establish divergence of the actual source information per area.
Neither(7) nor the positive lower information floor proves a continuum area law,
a regional diamond kernel, generalized-entropy stationarity, a Newton
constant or an Einstein equation. Boundary cardinality is a source
incidence count; it has not been identified with proper area.

## 5. Why the polynomial correction cannot be omitted

Take a separate two-qubit control with no interaction,
H=diag(0,1,1,2)=H_R+H_out. Use the same positive polynomial form,
with beta=1, Taylor order1 and dyadic scale2, so y=(2-H)/4 lies in[0,1/2].
Its exact normalized density is

    rho_control=diag(1296,625,625,256)/2802.                  (9)

The product determinant of this two-by-two probability table is
(1296*256-625²)/2802²=-58849/2802², which is nonzero. The marginals
are correlated despite V_cut=0, so X=I(rho_control)>0. Equations(3)-(5)
include this correlation in r. An exact zero-boundary Gibbs inequality
without r would be false for the polynomial preparation. This control
uses a different carrier and is not substituted for the primary source.
The tests also evaluate(1)-(5) with nontrivial representation dimensions,
conditional quantum correlations and off-diagonal center coherence.
