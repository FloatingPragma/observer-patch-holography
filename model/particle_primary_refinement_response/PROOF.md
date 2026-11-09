# Populated original responses at every finite capacity

The source is the selected open Wilson member, with all sixty CAR modes
per cell, the native fifteen-dimensional gauge representation, the
tensor product link-word cutoff, and the source Gauss constraints at
every endpoint. The coefficients are \(e=105\), \(\kappa=1/7\).
The source density is its specified positive polynomial \(Q(H)^2/Z\),
the full conjugate reference is retained with idle Hamiltonian zero,
and the clock is \(\exp(-itH)\). No source parameter, preparation,
energy zero or empirical conversion is changed.

At level \(n\), write \(L=2^{2n+1}\), \(V=L^3\), \(a=2^{-n}\),
and \(\widehat H=aH\). The selected mean ground excess is at most
\(\epsilon=2^{-n-4}\). The source has

\[
 B=3L^2(L-1),\qquad P=3L(L-1)^2                         \tag{1}
\]

interior bonds and plaquettes. This construction proves a positive
original response weight for every such finite member. Its energy
interval depends on \(a\) and the selected inverse-temperature
parameter. It supplies no fixed finite-energy continuum spectral weight.

## 1. An extensive physical trial bound

The [eleven-packet trial](../particle_primary_variational_readout/PROOF.md)
on the two-cell cube is a physical normalized vector. Its exact rational
dimensionless Rayleigh quotient \(u_0\) obeys

\[
 -723.767820<u_0<-723.767819,\qquad
       g=(u_0+720)/8\simeq-0.4709774928.                 \tag{2}
\]

Partition the larger cube into disjoint two-cell cubes. Put this vector
in each cube and put every link between cubes in its constant Haar
vector. The cube vector has constant outgoing ghost links. Its Gauss
invariance at a boundary vertex therefore persists when constant incoming
and outgoing link factors are attached. Its total fermion parity is
even, so the tensor product of the tile vectors defines an unambiguous
physical state on the finite graded tensor product. This is an auxiliary
trial vector; it does not replace the selected density.

The internal tile terms have exactly their two-cell expectations.
An intertile hopping term contains a nontrivial native matrix entry
on an independent constant link. Its expectation is zero by the Haar
integral of that entry. An elementary plaquette that is not internal
to one tile contains an intertile link. That link occurs once in each
character term, so its independent Haar integral also vanishes. This
argument uses a zero expectation, without discarding those terms from
the source Hamiltonian. Intertile electric terms vanish on their
constant vectors. Every trial factor has degree at most one per link,
so it belongs to all selected cutoffs. Compression preserves its
Hamiltonian expectation, including when a Hamiltonian term maps it
outside the trial span.

The full physical ground consequently satisfies

\[
 E_0(H)\le {V u_0\over8a},\qquad
 a\operatorname{Tr}(\rho H)\le -90V+gV+a\epsilon.       \tag{3}
\]

This bound includes all omitted occupations and link states. No
invariance of the tile trial span or equality with the physical vacuum
is required.

## 2. Native-dependent onsite retention and an exact local block

Use the all-occupation [electric Schur argument](../particle_primary_occupation_bound/PROOF.md).
The native block dimensions and Casimirs are

\[
 (d_r,C_r)=(3,4/5),(2,9/20),(6,21/20),(3,7/10),(1,3/10).
\]

Let \(D_x^r\) count beta-positive occupations and beta-negative holes
at a site in block \(r\). It is positive and has eigenvalues between zero
and \(4d_r\). Retain an onsite fraction \(r_r\) separately in each
block, and assign \(m_r=(1-r_r)/2\) to each endpoint of an interior bond.
The degree is at most six, so all unassigned onsite terms are positive.
The pointwise full-Fock negative shift in the link-complement block is

\[
 \delta=4\sum_r d_r\left(\sqrt{m_r^2+1/4}-m_r\right).
                                                               \tag{4}
\]

The exact rational square-root upper bounds in the code satisfy
\(\delta<63/4\). The full nonconstant electric block is therefore
positive. Its Schur complement on the constant-link block bounds each
native contribution by

\[
  m_r(D_x^r+D_y^r)-w_r Q_r,\qquad
            w_r={1\over105 C_r/2-\delta}.                \tag{5}
\]

For \(d_r>1\), the inherited four-defect inverse-order bound supplies
a rational \(c_r\) for which (5) is at least \(-c_r I\).
For \(d_r=1\), execute all eight local CAR modes. In a bond-adapted
Dirac basis there are two active and two inactive rows on each endpoint.
Let \(E_{ab}^x=c_{xa}^\dagger c_{xb}\), \(a,b=0,1\), and similarly
at \(y\). Then

\[
\begin{split}
 D&=4I+\sum_{z=x,y}\sum_{j=0}^1
       (c_{zj}^\dagger c_{z,j+2}+c_{z,j+2}^\dagger c_{zj}),\\
 Q&=N_{x,\mathrm{active}}+N_{y,\mathrm{active}}
                       -2\sum_{a,b=0}^1 E_{ab}^x E_{ba}^y.
\end{split}                                                     \tag{6}
\]

Both matrices preserve the two total endpoint occupations. Their
twenty-five sectors have total dimension256 and maximum dimension36.
The code constructs their integer entries directly from the CAR signs.
Rational LDL elimination proves every pivot of \(mD-wQ+cI\) positive
in every sector. Dyadic bisection chooses \(c\); no floating-point
eigenvalue is used as a certificate. Empty, filled, inactive and
unbalanced-number sectors are included.

Summing (5), keeping the retained onsite terms, and using the magnetic
norm bound before the Gauss restriction gives the full operator inequality

\[
 \widehat H\ge -90V-B\sum_r c_r-P/7
                          +3\sum_{x,r}r_rD_x^r.          \tag{7}
\]

Taking the actual selected-state expectation and using (3) is valid
without knowing the ground occupation sector.

For \(r_r=3/20\) in every block, the exact local computation gives
\(\sum c_r<1.461\). Equations (3) and (7), together with

\[
 B/V\le3,\quad P/V\le3,\quad a\epsilon/V\le1/128,
\]

give the uniform result

\[
 {\operatorname{Tr}(\rho D_{\mathrm{total}})\over V}
 \le {g+3\sum c_r+3/7+1/128\over9/20}
 <{97\over10}.                                          \tag{8}
\]

This is an average over the entire box, not a per-site estimate.
A separate allocation \(r_6=16/25\), with all other retained fractions
zero, yields a bound on \(D_{q=6}\) itself. At level one the code certifies

\[
 {\operatorname{Tr}(\rho D_{q=6})\over V}<1.821,
 \qquad {1\over2V}\sum_{x,\alpha}
 \operatorname{Tr}\rho\,
     n_{-,x\alpha}(1-n_{+,x\alpha})>{89\over1000}.        \tag{9}
\]

The last implication uses the all-occupation pair inequality from the
local response proof. The same q=6 allocation does not prove a positive
pair weight at every level. A failed lower estimate does not exclude
a physical response. The code obtains the stronger level-one bound
\(\operatorname{Tr}(\rho D_{q=6})/V<1.820\), hence a matched-pair
average greater than \(9/100\).

## 3. A uniformly populated original collective response

Use the native block of dimension six. At each site define the
gauge-invariant even operator

\[
 S_x^+=\sum_{\alpha=1}^2\sum_{i=1}^6
               c_{+,x\alpha i}^\dagger c_{-,x\alpha i},
 \qquad C_x=S_x^+/\sqrt{42}.                             \tag{10}
\]

The twelve beta pairs form the usual finite CAR spin algebra. On its
full Fock space, the maximum raising-operator norm squared is
\(6(6+1)=42\), including singly occupied, empty and doubly occupied
pairs. Thus \(\|C_x\|\le1\), and the exact commutator is

\[
 [C_x^\dagger,C_x]={12-D_x^{(q=1)}\over42}.
\]

Since \(D^{(q=1)}\le D_{\mathrm{total}}\), (8) implies

\[
 \chi={1\over V}\sum_x\operatorname{Tr}\rho
                    [C_x^\dagger,C_x]\ge {1\over20}.      \tag{11}
\]

All constituent native indices are contracted; the operator has gauge
charge zero. The constituent block's representation label is not a
particle assignment. A complete apparatus has \(V\) success Kraus
operators \(K_x=C_x/\sqrt V\) and one failure operator

\[
 K_f=\sqrt{I-V^{-1}\sum_x C_x^\dagger C_x}.               \tag{12}
\]

Equation (12) is a contraction instrument on the physical carrier and
retains the reference and every failure. It uses the declared exact
positive-root service. It does not preserve the input density in
general. Its unnormalized original success response is

\[
 \mu(d\omega)={1\over V}\sum_{x,i,j}
       p_i |(C_x)_{ji}|^2\delta_{E_j-E_i}(d\omega),
 \qquad \mu(\mathbb R)\le1.                              \tag{13}
\]

The selected polynomial eigenweights have the precise form

\[
 p(E)=Z^{-1}T_m\!\left({\beta(U-E)\over2s}\right)^{2s},
 \qquad 0\le -{d\over dE}\log p(E)
       =\beta {T_{m-1}\over T_m}\le\beta.                \tag{14}
\]

The argument is nonnegative on the entire source spectrum. Equation
(14) uses the actual positive polynomial, with no exact Gibbs
substitution. Reuse the [polynomial passivity window](../particle_primary_pair_channels/PROOF.md):
the negative-frequency part of
\(\chi=V^{-1}\sum_{x,i,j}(p_i-p_j)|(C_x)_{ji}|^2\)
is nonpositive, and its positive-frequency part below \(\delta_0\)
is at most \(\beta\delta_0\). If
\(\int(\omega-c)^2d\mu\le M_2\) and \(R>c\), then

\[
 \mu([\delta_0,R])\ge
                   \chi-\beta\delta_0-{M_2\over(R-c)^2}. \tag{15}
\]

The onsite commutator is \(6C_x/a\). Electric and magnetic terms
commute with \(C_x\). A native dimension-six bond hopping term has
norm at most12; each of the at most six incident bonds consequently
contributes at most24 to the commutator norm. Compression commutes
with the local matter operator. Hence, on the actual cutoff carrier,

\[
 \int(\omega-6/a)^2d\mu\le(144/a)^2.                    \tag{16}
\]

This conservative estimate requires no one-particle truncation. Choose
\(\delta_0=1/(100\beta)\) and \(R=1446/a\). Equations (11), (15)
and (16) give a bound uniform in capacity on the original record weight:

\[
       \boxed{\mu([1/(100\beta),1446/a])\ge3/100.}        \tag{17}
\]

For the dimension-one block use the scalar pair
\(C_x^{q=6}=\sum_{\alpha=1}^2c_{+,x\alpha}^\dagger c_{-,x\alpha}/\sqrt2\).
Its averaged commutator is \(1-\operatorname{Tr}(\rho D_{q=6})/(2V)\).
The native bond norm is at most two, so its centered second moment is
at most \((24/a)^2\). Applying (15) with the same lower endpoint and
upper endpoint \(246/a\) gives, at level one,

\[
 \mu_{q=6}([1/(100\beta),246/a])>7/100.                    \tag{17a}
\]

This is the original scalar-pair response, without a fermion-species,
spin or rest-momentum identification. The receipt supplies this window
only when its independently proved commutator lower bound is positive
enough to cover both frequency-tail allowances.

## 4. Final excitation energy and complete finite clocks

The [transfer-to-excitation coupling](../particle_primary_mass_readout/PROOF.md)
retains the original weights. The initial ground-energy excess has
mean at most \(\epsilon\). Discarding only its original tail above
\(100\epsilon\) costs at most \(1/100\), so the final excitation
measure satisfies

\[
 \nu([1/(100\beta),1446/a+100\epsilon])\ge1/50.           \tag{18}
\]

The existing [complete two-clock detector](../particle_primary_joint_spectroscopy/PROOF.md)
applies to (12), with all \((V+1)J^2\) fine records. Choose energy
resolution \(\eta=1/(400\beta)\), one-clock tail allowance \(1/400\),
and one-clock implementation allowance \(1/400\). Its exact rational
plan uses the same source clock and alias-free signed chart. Two tail
allowances and two complete implementation errors cost at most \(1/100\).
The reported energy-difference event therefore has original probability
at least \(1/50\) in

\[
 [1/(200\beta),1446/a+1/(200\beta)].                      \tag{19}
\]

The initial energy apparatus preserves only the stationary source
marginal; it changes the joint purification and retains that record.
Neither the full density matrix nor the enormous clock histogram is
numerically allocated. The positive-root and controlled-clock services
have their declared finite architecture and resource costs. A further
compiler/preparation error must be subtracted from (19).

## 5. Precise continuum and particle premises

The [primary field algebra](../particle_primary_field_algebra/PROOF.md)
supplies a common subsequential state, positive time generator, and
complete finite detector histories for its admitted uniformly bounded
time-regular probes. The [primary spatial algebra](../particle_primary_spatial_algebra/PROOF.md)
supplies regular spatial automorphisms and comparison of bulk moments.
Neither theorem makes the selected state spatially invariant or prevents
time-smearing from removing a capacity-dependent high-frequency response.
The generic consumers
`Lean/QFT/ColimitSelectedStateGNS.lean`,
`Lean/QFT/ColimitGNSRepresentation.lean`, and
`Lean/QFT/GroundRecordTail.lean` retain these state, domain and convergence
hypotheses. The CAR/contact programs in `continuum_field_net` and
`continuum_free_limit` have a different source and preparation; their
field dictionaries do not transfer to (10).

Equation (17) does not supply a fixed energy window: its lower endpoint
can tend to zero and its upper endpoint tends to infinity. An exact
two-level positive-polynomial control in the tests has vanishing mean
ground excess and success weight approaching one, while its excitation
energy tends to infinity. Thus the estimates used here do not imply
tight finite-energy spectral weight. A fixed positive band with a
uniform original lower weight, or an equivalent noncollapse and
tightness argument for a specified field normalization, is a separate
premise. Spatial implementation on the selected representation and the
momentum-origin control are also separate. An isolated spectral branch
and its dispersion must be attached to the declared physical clock
conversion before giving a mass interval. The existing source-unit
service supplies an SI/GeV conversion after the clock scale is declared;
unit declaration alone does not identify a physical rest energy. No interacting continuum theorem, named particle,
physical rest mass or GeV spectrum is asserted by (8), (17), or (19).
