# Joint native charge and Gauss-flux records

## One declared interacting source

Use the complete eight-cell Peter–Weyl source of
[`interacting_charged_transport`](../coherent_defect_wave/chiral_joint/interacting_charged_transport/PROOF.md):
30 CAR modes and five scalar states per cell, 119 states on each of 24
oriented periodic links, and
\(H_\lambda=H_0+2^{-34}T\). The onsite coefficients are one. The static
30-mode CAR reference sits at cell zero and has Hamiltonian \(4N_{\rm ref}\).
The source and reference clocks, the physical interacting-ground preparation,
its retained environment, and every preparation failure are unchanged.

Measure the complete atomic-ground/blank-link intervention \(\{P_0,I-P_0\}\),
then use the actual \(e^c\) removal pump on success. The reference has the
same native representation as a source creation operator. On pump success
its one particle has \(Y_{\rm ref}=Q_{\rm ref}=+1\); the source has lost one
unit of this charge. The projector effect of the pump on \(P_0\) is exactly
\(P_0/2\). Its failure is retained. Successful packets evolve for
\(t=2^{-24}\) under the complete interacting clock.

For each of the six actual links incident to the origin, the arrival
projector \(P_e\) requires that link nonblank, every other link blank,
and the true local matter-plus-scalar ground at every site except its
other endpoint \(y\). The entire local carrier at \(y\) is retained.
The six effects and their complement form the actual arrival instrument.
Their two periodic routes to each neighboring site have different link
memories. They are distinguished in the output.

## Central current and oriented cut flux

Let \(Y\) be the actual native hypercharge generator and write
\(Q=T_3+Y\). The local source density includes both matter and scalar
currents,

\[
\rho_x^Y=d\Gamma(I_2\otimes Y)_x+Y_{H,x}.
\]

The native endpoint generators on a forward link block are
\((-Y^T,Y)\), with the reversed convention on its dual. Every retained
irreducible block has a constant hypercharge. Consequently the full
119-state identity is
\(L_{e,0}^Y+L_{e,1}^Y=0\). This identity includes all excited-to-excited
link channels. Hypercharge commutes with all twelve currents on matter,
scalar and both link endpoint factors.

For any subset \(S\) of the eight cells, define the outward flux by

\[
F_S=-\sum_{e\text{ crossing }S}L_{e,\,\text{endpoint in }S}^Y.
\]

Summing the actual Gauss generators gives the operator identity

\[
\sum_{x\in S}G_x^Y
 =\sum_{x\in S}\rho_x^Y+\mathbf1_{0\in S}Y_{\rm ref}-F_S.
\tag{1}
\]

Every interior endpoint contribution cancels. Thus on the joint physical
source/reference carrier,

\[
F_S=\rho_S^Y+\mathbf1_{0\in S}Y_{\rm ref}.
\tag{2}
\]

The calculation compares the actual Gauss compiler with every additive
factor basis: the vacuum, 240 single matter modes, 32 nonvacuum scalar
states, and all 2,832 single nonblank link coordinates. There are 24,840
cell-current comparisons. Additivity of the actual diagonal current
extends this check to every full occupation and simultaneous link/scalar
configuration. The complete endpoint Ward identities, the scalar Yukawa
identities and the reference-pump current identities preserve the physical
constraints under the source clock and every retained preparation branch.
No large Fock-space matrix is substituted for this finite operator proof.

Equation (1) is a central-current identity. The complete link operator
\(L_0^Q+L_1^Q\) is nonzero in weak-charged sectors. Summing that current
requires the interior gauge-field charge as well. A complete-basis
counterexample rejects its omission. The measured cut flux in this
construction is the gauge-invariant central \(Y\) flux; a general
nonarrival measurement of a gauge-dependent weak-current component is
not assumed.

## Exact arrival charge and flux

Let \(R_{e^c}\) project onto the one-particle \(e^c\) sector of the
complete static reference Fock space. The specified pump on reference
vacuum has range in this sector. The source clock with \(4N_{\rm ref}\),
its spectral projectors and the source-only arrival effects commute with
\(R_{e^c}\). The fixed charge and flux values below apply on the joint
physical subspace in \(R_{e^c}\); the general Gauss identity (2) applies
to every physical reference sector. The detector \(P_e\) itself is
identity on the full reference and imposes no reference-sector filter.

On \(P_e\), the origin source is a true local gauge singlet. All links
except \(e\) are blank. The joint origin Gauss equation therefore gives

\[
(L_{e,o}^A+T_{\rm ref}^A)P_e\Psi=0,
\qquad A=1,\ldots,12.
\tag{3}
\]

The full 119-state origin endpoint has exactly one basis vector with
\(Y=-1\). The other ten retained irreducible charges are distinct from
that one. This vector satisfies every equation (3), not merely its
hypercharge component. It is the one-dimensional flux sector conjugate
to the reference at the origin. Its arrival endpoint has \(Y=+1\), and
all colour and weak generators vanish at both ends. The source's actual
hopping compiler produces exactly this oriented link vector in each of
the six refilling branches; its coefficient norm squared is \(1/16\).

The arrival Gauss equation fixes the entire local matter-plus-scalar
state's current, even though that local state is unrestricted in the
detector definition:

\[
\rho_y^Y=\rho_y^Q=-1,\qquad
\rho_x^Y=\rho_x^Q=0\quad(x\ne y).
\tag{4}
\]

In particular the total weak current vanishes at \(y\), so \(Q=Y\)
on this attached sector. This is an operator consequence of the actual
joint constraint and complete detector, without choosing a particle label.
For every cut and every joint-physical vector in \(P_eR_{e^c}\),

\[
F_S=\mathbf1_{0\in S}-\mathbf1_{y\in S}.
\tag{5}
\]

The raw link endpoint contribution at the arrival site is \(+1\).
The outward flux enclosing that site and excluding the reference is
\(-1\), equal to its charge under the sign convention in (1).
The executor checks all 1,536 route/cut identities, including the empty
and full cuts. For the three half cubes with coordinate one, the flux
triples are \((-1,0,0)\), \((0,-1,0)\) or \((0,0,-1)\), according to
the arrival axis.

## Complete simultaneous readout

Every \(6F_S\) has an integer spectrum. Its spectral projectors are
literal diagonal projectors on the complete boundary-link tensor product.
The executable computes their exact multiplicities by integer convolution,
with total dimension \(119^{|\partial S|}\); no link state is omitted.
Projectors for different cuts commute. They also commute with the six
arrival effects and all local gauge currents. The corresponding total
central-charge projectors act on the matter, scalar and reference factors.
They commute with the flux projectors. On the physical carrier, (2) makes
every charge/flux mismatch outcome have probability exactly zero.

Let \(K_r\) denote the actual complete Kraus operator for a coarse record,
with its preparation environment retained. If \(\Pi_{S,k}\) projects
onto \(6F_S=k\), and \(R_{S,l}\) projects onto six times the right-hand
side of (2) equaling \(l\), the joint probability is

\[
p(r,k,l)=\operatorname{Tr}\!\left[
 R_{S,l}\Pi_{S,k}K_r\rho_{\rm in}K_r^\dagger
 \Pi_{S,k}R_{S,l}\right],
\quad p(r,k,l)=0\ (k\ne l).
\tag{6}
\]

The controller record is read first. This formula refines its output
state and requires no commutation of a prior preparation Kraus operator
with the flux. All spectral outcomes are retained. Summing (6) returns
the original coarse probability, and the sum over all records is one.

Write \(q\) for physical preparation success, \(a\) for conditional
\(P_0\) success, and \(p_e\) for the original unnormalized route
probability conditional on preparation. The ten coarse probabilities are

\[
1-q,\quad q(1-a),\quad qa/2,\quad
\{qp_e\}_{e=1}^{6},\quad q(a/2-\sum_e p_e).
\tag{7}
\]

Preparation failures, localization failures, pump failures and the complete
nonarrival output receive their own flux refinements. No failed output is
removed. On every route, (4) and (5) fix the charge/flux atom uniquely,
so its probability is exactly \(qp_e\). The inherited full-Hamiltonian
bounds, freshly calculated with the physical preparation error, give

\[
6.430004784\,10^{-38}<qp_e<1.274554272\,10^{-37},
\]

and a six-route union between \(3.858002870\,10^{-37}\) and
\(7.647325626\,10^{-37}\). Exact outward rational bounds occur in the
receipt. Individual flux probabilities inside the unselected outputs
are specified by (6) and constrained by (2); they are not numerically
evaluated here. The positive route signal is small because the declared
hopping and observation time are small.

A Gauss-preserving energy filter on the same source preserves (2)–(5)
whenever an arrival clicks in the one-\(e^c\)-reference sector, provided
the filter also preserves that sector, as the actual full-H spectral
projector does. This statement supplies no positive probability
for an energy-filtered arrival and assumes no commutation of energy and
arrival projectors.

## Physical interpretation

The result attaches a unit of the specified native current to its actual
finite source flux on the same transported packet, clock, reference and
complete record law. It requires no equality between current and texture
winding. A physical electromagnetic field normalization, a measured
fine-structure constant, an inverse-square field profile, and an autonomous
reference-free electron identity do not follow from this finite attachment.
