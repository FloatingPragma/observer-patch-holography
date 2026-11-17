# A continuum proper-regulator source with finite actual regional entropy

This is a new declared free-field architecture. It establishes existence,
finite regional entropy and exact spatial coordinate covariance. It does
not establish an area asymptotic, a shape derivative, entanglement
equilibrium, a lapse equation, interacting chiral protection or physical G.

## 1. Source, spin doubling and domains

Let (M,g) be a closed smooth oriented three-dimensional spin manifold with a
specified spin structure. Let d_g be its self-adjoint intrinsic two-component
Dirac operator on L²(S_g,dvol_g). Use the doubled bundle E=S_g tensor C² and

    D_g=d_g tensor X,    beta=I tensor Z.

Then beta is a bounded parallel self-adjoint involution and beta D=−D beta.
The doubling is essential: an irreducible odd-dimensional spinor bundle has
no additional chirality matrix with this property. The Hilbert-space inner
product is the proper-volume spinor inner product; equivalent half-density
coordinates may be used.

For one independently declared proper length R>0, define

    H_g=D_g+beta R³D_g⁴,    Dom(H_g)=Dom(D_g⁴).

The leading term beta R³D⁴ is self-adjoint. D is relatively bounded with
respect to D⁴ with arbitrarily small bound, by the scalar inequality
|lambda|≤epsilon|lambda|⁴+C_epsilon. Thus H is self-adjoint on this domain.
Equivalently, the orthogonal eigenspaces of D² reduce H to finite-dimensional
Clifford blocks, which define the same self-adjoint direct sum. On the common
core of smooth sections, and hence for the associated self-adjoint closure,

    H²=D²+R⁶D⁸.

The cross terms vanish because beta anticommutes with D and commutes with D⁴.
No negative or complex energy-square branch is introduced.

Initially assume ker D=0. This assumption has concrete sources: the flat
three-torus with an antiperiodic spin structure in at least one cycle has
no zero momentum, hence no zero Dirac eigenvalue. Small smooth deformations
within an open spectral-gap neighbourhood preserve invertibility.

## 2. Exact projector and Schatten estimate

Write A=sign(D), T=|D| and

    f(T)=(1+R⁶T⁶)^(-1/2),
    b(T)=R³T³(1+R⁶T⁶)^(-1/2).

A and beta anticommute and both commute with f,b. The negative-energy
projector and the local product-reference projector are

    P=(I−A f(T)−beta b(T))/2,
    P_beta=(I−beta)/2.

Their difference Delta=P−P_beta has the EXACT modulus

    Delta²=(1−b(T))/2,
    |Delta|=[(1−b(T))/2]^(1/2).

For t>0, with x=R³t³,

    1−x/sqrt(1+x²)
       =1/[sqrt(1+x²)(sqrt(1+x²)+x)]≤1/(2x²).

Consequently |Delta|(t)≤1/(2R³t³), and asymptotically
|Delta|(t)~1/(2R³t³). Weyl spectral growth for a first-order elliptic operator
in dimension3 gives |lambda_n(D)| comparable to n^(1/3), with multiplicities
included. Therefore

    Delta belongs to S_q for every q>1.

Generically this argument does not give trace class S_1. It DOES give
Hilbert-Schmidt S_2, so the pure quasifree state P is implementable relative
to the local P_beta polarization on a compact M.

An explicit rotation removes any ambiguity about that implementation. Set
theta(T)=arctan[1/(R³T³)] and

    V=cos(theta/2)+A beta sin(theta/2).

A beta is anti-Hermitian with square−I. Hence V is unitary and
V beta V^dagger=beta b+A f, so P=V P_beta V^dagger. The off-diagonal part of
V with respect to P_beta has singular values sin(theta/2)=|Delta| and is
Hilbert-Schmidt. This meets the stated implementation criterion directly.

The spectral-growth input is the standard elliptic Weyl law; an original
primary source is [Hörmander, The spectral function of an elliptic operator](https://doi.org/10.1007/BF02391913).
The implementability criterion is discussed, with its operator-ideal
hypotheses, in [Avron, Bachmann, Graf and Klich, Fredholm determinants and the statistics of charge transport](https://arxiv.org/abs/0705.0099).
Neither reference supplies the regulator construction or the estimate above;
those are the explicit calculation of this packet.

## 3. Finite actual entropy of every measurable region

Let chi be multiplication by the indicator of a measurable region B. On
chi H define the actual compressed covariance C_B=chi P chi. P_beta is local,
so chi P_beta(1−chi)=0. Thus

    T_B=chi P(1−chi)=chi Delta(1−chi) belongs to S_q, q>1,
    C_B(1−C_B)=T_B T_B^dagger.

Choose any p with1/2<p<1 and q=2p. Singular-value calculus gives

    Tr [C_B(1−C_B)]^p = ||T_B||_(2p)^(2p)<infinity.

For h(t)=−t log t−(1−t)log(1−t), the elementary bound

    h(t)≤C_p [t(1−t)]^p,   0≤t≤1,
    C_p=2^p[1+1/(e(1−p))]

is sufficient. To verify it, use symmetry to take t≤1/2,
−(1−t)log(1−t)≤t, and
t^(1−p)(−log t)≤1/[e(1−p)], then t(1−t)≥t/2.
Functional calculus and the trace bound imply

    S(B)=Tr h(C_B)<infinity.

This is actual quasifree regional entropy, not a chosen area functional or
a subtraction of infinities. For clarity about the infinite-dimensional
representation, let Q_B=1_[1/2,1](C_B). The bound implies
Tr|C_B−Q_B|<infinity. Also C_B−(P_beta)|_B is Hilbert-Schmidt, so
Q_B−(P_beta)|_B is Hilbert-Schmidt. The Q_B Fock polarization is therefore
equivalent to the local product-reference polarization. Diagonalizing C_B
then gives a normal independent-mode density matrix, with summable mode
entropies exactly Tr h(C_B). This can also be obtained by increasing finite
spectral-mode truncations and monotone sums of the nonnegative h(c_j).

The finite-dimensional quasifree entropy identity and its mode-product
derivation are given in [Dierckx, Fannes and Pogorzelska, Proposition4](https://arxiv.org/pdf/0709.1061).
The trace and representation conditions needed here are established above;
the finite-dimensional result is not silently treated as an infinite one.

This bound does NOT prove proportionality to boundary area. It even applies
to arbitrary measurable B on the compact source; its upper bound is not a
sharp geometric estimate.

## 4. Flat infinite volume and a genuine zero-mode identity

On R³, with its usual spin structure, take D=alpha·(−i gradient) and the same
H=D+beta R³D⁴. The point p=0 has measure zero, so the spectral projector is
unambiguous almost everywhere. Delta(p) is bounded at low momentum and
O(|p|^−3) at high momentum. Hence Delta maps L²(R³) boundedly into H³(R³).

For a bounded region B, choose a smooth compactly supported function a=1
on B, with support inside a bounded smooth domain Omega. Then a Delta maps
L² boundedly into H³_0(Omega). The compact Sobolev embedding H³_0(Omega) to
L²(Omega) has singular values O(n^−1), so a Delta and chi_B Delta belong to
S_q for every q>1. Multiplication by chi_B on the output is bounded; its
boundary need not be differentiated. The argument of section3 gives finite
S(B) on flat infinite space as well. Global Delta need not be Hilbert-Schmidt
in infinite volume; only its localized version was used.

Now g_sigma=e^(2sigma)delta with constant sigma. The dilation F_sigma(x)=
e^sigma x is an isometry from (R³,g_sigma) to (R³,delta). It maps the
coordinate ball B_(e^−sigma ell) to B_ell. Naturality of the spin Dirac
operator, its functional calculus and the local region projection therefore
give EXACTLY

    S_(g_sigma)(B_(e^−sigma ell))=S_delta(B_ell).

In flat half-density coordinates one can also see the identity directly:
H_sigma(k)=H_0(e^−sigma k), so P_sigma(k)=P_0(e^−sigma k) and the covariance
kernel satisfies C_sigma(x,y)=e^(3sigma)C_0(e^sigma x,e^sigma y). The unitary
map from the moving coordinate ball to the fixed proper ball is
(U_sigma f)(y)=e^(−3sigma/2)f(e^−sigma y). It conjugates the compressed
covariance exactly to the reference one, proving equality of its spectrum.

This is uniform fixed-PROPER-volume cancellation of the actual continuum
source and actual local CAR region. No fitted K(0) subtraction is involved.
It is an invariance statement, so it does not require establishing a shape
derivative first.

## 5. General spatial diffeomorphisms and spin structures

For an orientation-preserving diffeomorphism F, pull back BOTH metric and
spin structure. The spin lift is unitary (two lifts differ only by sign),
and intertwines D, beta, H, P and the mapped local region projection. Entropy
is consequently unchanged under this simultaneous coordinate change.
If one insists on a fixed spin structure, only diffeomorphisms preserving it
act within that chosen source; diffeomorphisms isotopic to the identity do.
These spin-lift and Dirac-equivariance statements are established in
[Dabrowski and Dossena, Dirac operator on spinors and diffeomorphisms](https://arxiv.org/abs/1209.2021).

A uniform rescaling of the metric on a fixed compact torus changes its total
proper volume and is NOT such a diffeomorphism. One must also scale the
coordinate period to describe the same proper torus. Therefore the compact
covariance theorem does not eliminate finite-volume/Casimir effects in the
existing fixed-period lattice diagnostics. The infinite-space identity above
is a separate statement with its own local entropy proof.

Spatial diffeomorphism covariance is not Weyl/conformal symmetry. In
particular D_g⁴ is simply the fourth power of the geometrically defined
Dirac operator; no conformally covariant fourth-order operator is asserted.
Nor is this a spacetime Lorentz-covariant ultraviolet regulator.

## 6. Zero modes, stress and remaining obligations

If a compact source has ker D nonzero, it is finite-dimensional. Choose a
zero-energy filling explicitly; the covariant choice P_beta restricted to
ker D changes the above formula by a finite-rank term, so finite regional
entropy still holds. A spectral crossing may make a chosen vacuum family
nondifferentiable. Response calculations therefore need a gap neighbourhood
or a separately specified treatment of zero modes; no such crossing is
silently discarded.

Finite entropy is not finite vacuum stress. Per high-energy two-level block,
the energy difference of the P vacuum from the P_beta reference is

    −sqrt(lambda²+R⁶lambda⁸)+R³lambda⁴
      ~−1/(2R³lambda²).

In dimension3 its spectral sum has a leading volume ultraviolet divergence.
This packet has not supplied a common stress renormalization/counterterm
prescription across geometries. Finite-energy neutral excitations are a
separate, well-defined preparation, as in ../pure_neutral_slater/.

Still unproved: differentiability under physical curvature/shape changes;
the area-law coefficient and a controlled R/ell expansion; the full ball
modular-stress form; entropy equilibrium selecting one common geometry;
noncompact curved-source and boundary error bounds; lapse and all-observer
response; interacting chiral protection; a native architecture selection;
and transfer of the independently declared R into the physical clock map.
The periodic positive-energy static scalar mean-mode problem is not solved
by the coordinate identity.
