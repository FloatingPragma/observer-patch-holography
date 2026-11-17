# Exact finite bridges and physical domains

The [specification](../../SPEC.md) defines H_S=sum_v h_v,
k_v=g_v+N_v(h_v−e_v)+e_v and K_v=A_v+M_v(k_v−F)+F.
The complete local-capacity Hamiltonian is H=sum_v K_v. Each vertex owns
its Phi and capacity factors. F is one common number. All normal-order
scalars and finite field-wall projectors remain present.

## Sectors, current and unitary dynamics

For a diagonal native charge D=sum_i d_i c_i†c_i and normal CAR word T,
[D,T]=(sum_created d_i−sum_annihilated d_i)T. The complete polynomial
certificate applies this identity to electric charge, baryon charge, Dirac
lepton charge and parity, on every word and every shard. Auxiliary and field
terms commute with these charges. The sector support identity therefore
passes to each K_v and its sum. A nonempty finite charge sector has a lowest
eigenvector; neither a positive mass gap nor particle localization follows.

For a dressed word with charge change s and shift exponent ell, B ell=s.
The current of the complete source is the phase derivative
J_e=sum_v M_v N_v partial_theta_e h_v,matter. Each word receives the mean
M_v N_v over its allocated support. Since these factors commute with charge,
the same exact local continuity identity holds after the join. Actual finite
shifts satisfy [E,U^ell]=ell U^ell, including at the wall. The proof never
uses U†U=I. Second derivatives supply the declared contact response.

Hermiticity and finite dimension yield unitary dynamics and finite Gibbs
identities. The abstract Lean interfaces expose extraction of the actual
finite matrix as an assumption. A stationary density commuting with H has
constant statistics under its own evolution; its global phase cannot supply
an observable clock.

## Complete coherence and stationarity countervariation

The specified triangle selection rule gives H_S matrix element 1/2 between
the native empty state and a charged pair/flux state. All capacity labels
zero and Phi index one give endpoint factors M_v=N_v=3/4. The complete
local-capacity matrix entry is z=9/32. Every output of its initial column
is retained and checked for Gauss compatibility.

The projectors on (|i>±i|f>)/sqrt(2), together with I−P_+−P_−, sum to identity
on the complete carrier and on every reference extension. The initial
probabilities are (1/2,1/2,0) and derivatives (−z,+z,0). Terminal dephasing
and deletion of the actual Dirac source term give zero for this response.
No two-state evolution approximation or finite-time native propagation is
used; the orthogonal complement is the entire complement.

At beta=1/1024 the faithful Gibbs state admits sufficiently small variations
delta rho=epsilon(|i><f|+|f><i|) of either sign. Trace and every diagonal
area/volume expectation have zero variation, whereas

    d(S+eta A)/d epsilon = 2 beta z = 9/16384.

Thus unrestricted full-Gibbs stationarity with diagonal area and fixed mean
diagonal volume fails for this source. This counterexample does not exclude
a justified regional area/entropy law with physical joint variations.

## Faithful reference and regional first law

At cap6 and zero reference charge, count the multiplicity n(q) from
product_i(1+z^(charge_i)) for all 96 local modes. Summing product_v n((BE)_v)
over all 13^|E| flux assignments gives the complete matter/field Gauss rank.
Multiply by 48^|V| to retain all independent Phi/capacity factors.

The full coefficient norm bound for H_S, local geometric bounds and the
actual capacity norms give ||H||<=M. With d=rank(P_G)>0,

    rho=P_G exp(−beta H)/Z >= exp(−2 beta M) P_G/d.

On a declared ordered region, D_R=Tr_out P_G is diagonal with integer
multiplicities. For P_R=supp(D_R) and lambda=exp(−2 beta M)/d,
rho_R>=lambda D_R>=lambda P_R, so the supported logarithm exists. Nested
partial traces compose. The exact entropy identity is

    S(tau)−S(rho_R)=Tr[(tau−rho_R)(−log rho_R)]−D(tau||rho_R).

For ||tau−rho_R||op<=lambda/2 the remainder is bounded by
||tau−rho_R||HS²/lambda. These finite floors need not be uniform under
refinement. No inverse recovery map, efficient native preparation, physical
vacuum identification or small-ball stress kernel is supplied by positivity.

## Local energy Ward identity

Define J_vw=i[K_v,K_w]. Hermiticity, gauge invariance and antisymmetry follow
from the same properties of K_v. Under H=sum_v K_v,

    dK_v/dt + sum_w J_vw = 0.

Taking expectations or time-integrated expectations gives the exact endpoint
energy identity. Summing over a region cancels every internal pair. Pushing
energies and currents along a vertex block map preserves continuity and
antisymmetry; two block maps equal their composition. These instantiate
[CoarseEnergyWard.lean](../../../Lean/QFT/CoarseEnergyWard.lean)'s
`coarse_continuity`, `coarse_antisymmetric`, `internal_block_flux_zero`,
`push_compose` and `coarse_compose`. The independently assembled complete
256-state local-capacity diagnostic checks the partition, every column,
nonzero current, cancellation and an omitted-current control. The algebraic
argument applies to the full native operator without dense Fock enumeration.

Disconnected components have disjoint even native supports and independent
auxiliaries, so their K_v commute. At bounded graph degree/face size and fixed
finite local dimensions the interactions have finite range. Physical
relativistic stress, a uniform refinement cone and covariant conservation
remain separate hypotheses.

## Geometry, coframe and theorem inheritance

The shell family has n=2^level, (2n²+2)(n+1) vertices, exact B1 B2=B2 B3=0,
compatible retained vertices and edge paths, and diameter at most
(8 sqrt(6)+1)/n. These facts do not establish native source/reference or
instrument naturality. Hard-wall plaquette cancellations require their own
cutoff analysis.

For constant det(V)>0, p=V^(−T)d changes the directional Dirac coefficients,
p.p changes the Wilson coefficient, b.p gives shift response, and the field
weight is d^T V^T V d/det(V), using the face normal for magnetic energy. The
origin is exactly H_S. Six symmetric strains and three shifts provide
mechanical source derivatives, including Wilson terms and mixed contacts.
They are not a spatial spin connection or an identified Lorentz stress.

The retained Einstein theorem requires the same regional source/reference,
proper metric, entropy/area data and o(ell^4) errors. Native charge-sector
existence does not instantiate its dark-sector hypotheses either. Universal
coupling, collar recovery errors, proper cut density and independent variance
composition remain required for the existing DarkSector, DeepProfileClosure
and PerCutCollarComposition implications. A stationary neutral finite state
is not an independently derived dark particle or cosmology.
