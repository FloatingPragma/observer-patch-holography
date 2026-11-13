# Exact strong-source reduction and constrained potential

The same seven-input current action and its own conditional length and electric
coefficient define every finite chain. The exact native colour-singlet
isometry gives dimensions120,364,816,1540 for two through five cells. Every
charge/spin-projection block is retained. The complete underlying native
three-quark occupation labels and their colour-filter rank counts are
retained as a separate kinematic preparation stage. Physical evolution and
measurement act on the full Gauss-dressed singlet carrier. The pair-triplet interaction vanishes on all these pure
quark sectors. The two-cell intertwiner checks the complete common source.
Changing graph size changes the architecture; equality of physical predictions
across that change is tested rather than assumed.

## Complement memory and current vertices

Within a complete conserved-charge/spin block let P contain all states of
occupied span at most one cell, and Q its full orthogonal complement. Write
H=[[A,V],[V*,D]]. For every real E outside spec(D),

    T(E)=(E-D)^(-1)V*, R(E)=[I;T(E)],
    H_eff(E)=A+V(E-D)^(-1)V*, M(E)=R(E)*R(E).

The full eigen-equation is equivalent to H_eff(E)x=Ex together with
psi=R(E)x. No Q eigenmode is removed from the resolvent. The normalized
embedding is Z=R M^(-1/2). For the actual native current J between initial
and final charge sectors, the full matrix element is exactly

    y*R_f*J R_i x / sqrt[(y*M_f*y)(x*M_i*x)].

Replacing it by P_f J P_i drops a calculable vertex and normalization
correction. All position effects pull back through Z to a complete POVM;
the full-state complement is retained. The code reconstructs both charge
states and verifies the actual current against the full native matrices.
Equal weak-current labels alone are not used as evidence.

For a declared exact rational energy window below D, the source's exact
unnormalized integer Gram form gives a rational Gershgorin lower bound d0.
The exact sum of squared V entries bounds ||V||^2. With gap delta=d0-Emax,

    ||T|| <= ||V||_F/delta = tau,
    ||Sigma(E)-Sigma(F)|| <= ||V||_F^2 |E-F|/delta^2.

Since ||M^(-1/2)-I||<=tau^2/2 and ||M^(-1/2)||<=1,
||Z-P||<=tau+tau^2/2. Both native three-quark raising currents have norm at
most3, giving the displayed6(tau+tau^2/2) isometric-current error bound.
These bounds apply to the whole declared window. Floating lowest eigenvalues
and their residuals are explicitly diagnostic; a physical pole certificate
is not inferred from them.

Eliminating Q in real time gives the exact Volterra equation

    i dc/dt = A c - i integral_0^t V exp[-iD(t-s)] V* c(s) ds
              + V exp[-iDt] q(0).

The complete spectrum of D determines the memory kernel. The producer checks
its zeroth and first moments and evolves the entire finite H from a fixed
coarse basis preparation. The physical input is I/d on the complete singlet carrier, with a
charge-balancing reference correlated by sector. Use compressed effects
E_i=B*|i><i|B and their positive square-root Kraus operators. Each row of B
has at most one orbital column, so these effects sum to I and define a
unital physical instrument. After a rank-one filter and full source
evolution, the original-input records are |BU(t)psi|_i^2/d and
[(BB*)_ii-|BU(t)psi|_i^2]/d. Every Q, failure, native label and outside-carrier
identity effect is retained. All maps preserve charge and are unital; the
declared physical maximally mixed/reference input is unchanged
nonselectively, so its total mean source/reference work is zero. Raw
coloured occupation poststates are not used as physical Gauss outputs. No irreversible
relaxation or continuum QCD running is asserted.

The occupation record domain is determined by the exact epsilon isometry,
not by testing a floating amplitude for nonzero. A native occupation is
compatible precisely when it has one quark of each colour and its sorted
orbital label lies in the conserved charge/spin block. For Nu=1 and
Nspin-up=2 this gives exactly `9 L^3` native labels: three choices for the
up-flavour colour, three for the down-spin colour, and L positions for
each colour. Each associated epsilon coefficient is nonzero. Every row
outside this domain is identically zero for every vector in the block,
and therefore for all source times. The records enumerate this entire
fixed domain, including zero and arbitrarily small computed probabilities.
The source evolution, amplitudes and success/failure formulas are unchanged;
there is no numerical support threshold or removed outcome.

## A full-source linear inter-charge potential

On an eight-cell chain pin a quark at cell2 and a diquark at cell2+r,
r=2,3,4. Both occupied sites have degree2 and are nonadjacent. The complete
internal colour-singlet carrier has dimension4 times10=40. It is the same
native CAR carrier under a cell relabelling, not a selected spin state.
Every source term compresses as follows:

    H_pinned(r) = [3mbar+5A/9+(8/3)(U-kappa/2)+6t+(4kappa/3)r] I
                + delta(3-2Nu)/2
                + (Js+2lambda)(Sq^2+Sd^2)
                + (Jw+2eta)(Iq^2+Id^2).

Every off-diagonal hop leaves the pinned sector; its retained positive
diagonal cost is6t. Each local spin/isospin square occurs on two adjacent
edges, which explains the factors2. The global colour square annihilates
all40 singlets and the pair-triplet force annihilates all pure-quark states.
The exact epsilon colour identity gives the electric term4kappa*r/3. Thus
H_pinned(r+1)-H_pinned(r)=(4kappa/3)I on all40 internal states, verified as
an exact matrix identity. Every complete constrained eigenvalue has the same
linear shift. With the source's conditional a, the slope is exactly1GeV/fm.

This is a total-source constrained potential with explicit external position
constraints. It is not the unpinned ground energy, a complex resonance pole,
or a four-dimensional continuum scattering potential. Its1GeV/fm value is a
disclosed attachment target. The full pinned-occupation colour-filter rank counts describe a
kinematic preparation stage. A physical spectral instrument acts within
all40 Gauss-dressed singlets, retains every internal spectral rank and its
reference, and adds the complementary identity effect. No raw coloured
poststate is identified with physical Gauss evolution. No pure-gauge Yang-Mills mass-gap statement follows from this matter
string construction.
