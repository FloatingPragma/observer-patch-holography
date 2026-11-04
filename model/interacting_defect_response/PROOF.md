# A quartic interaction on the common localized carrier

The source is the finite CAR construction in
[defect_current_bridge](../defect_current_bridge/PROOF.md). For every admitted
odd radial cube, its isometry `E=W tensor I45` embeds ninety annihilators into
the full field. The radial kernel obeys `hW=0`, `W*W=I2` and
`W_x*W_x=rho_x I2`. Every complete memory record has a quadratic Hamiltonian
which reduces this same CAR subalgebra. The native current probes act on fresh
auxiliary legs, so their transition probabilities are independent of the field.

## Interaction and exact inclusion

Let `P=EE*`, `N_W=dGamma(P)` and choose a real coefficient `lambda>=0`. Add

```math
V=\lambda(N_W-45)^2
```

to every record Hamiltonian. This is chosen simulator design. It uses the
source kernel, not an empirically fitted particle assignment. In the localized
annihilators `a_i`, its normal-ordered form is

```math
V=\lambda\left(2025-89\sum_{i=1}^{90}n_i
                  +2\sum_{i<j}n_i n_j\right),\qquad n_i=a_i^\dagger a_i.
```

The interaction is parity even. It commutes with every native internal current
because `P=WW* tensor I45`. The all-left charge
`Q=3(B-L)` assigns `(1,-1,-1,-3,3)` to `(Q,u^c,d^c,L,e^c)`.
The exact native matrices commute with this charge, and both components of
every charged Higgs tensor obey `QK+KQ^T=0`. Consequently every charged pairing
record and V conserve this charge. On the ninety-mode Fock space it lies in
`[-72,54]`. A neutral Majorana pairing has charge `-6` and does not obey this
conservation statement. No gauged B-L theory or anomaly cancellation is asserted.

The decomposition into localized and complementary CAR factors gives

```math
H_r^{\rm full}=H_{r,W}\otimes I+I\otimes H_{r,\perp},\qquad
H_{r,W}=H_{r,W}^{(2)}+\lambda(N-45)^2.
```

All terms are even, so the factors commute without a graded sign. Thus every
record preserves the localized algebra exactly. Tensoring a parity-even core
state with normalized trace on the complement gives the exact state inclusion
for every record. `N` is not conserved by the pairing term: it changes by two.
Neither commuting V with H nor replacing V by its expectation is used.

## A bound uniform in the radial cube

The ninety-mode pairing Hamiltonian with antisymmetric kernel Delta has norm
at most `90||Delta||/2`, by its skew singular-block decomposition. The declared
quark kernel has norm below one; a unit-gain lepton projector kernel has norm
one. The radial kinetic term vanishes on the kernel. If g is the chosen
feedback gain, then every record on the localized algebra satisfies

```math
\|H_{r,W}\|\le M_W:=45(1+g)+2025\lambda.
```

This bound does not depend on the cube radius. It applies to the actual g at
the chosen source value of P, without freezing the inherited numerical gain.

The complete memory transition is `T=(1-p)R+pS`; R resets to a fixed port and S
is stochastic. For zero-sum operator-valued memory blocks, its trace-norm
contraction is at most p. The record unitaries differ from identity by at most
`2 M_W h` as trace-norm channels. Iterating this contraction and the second-order
unitary remainder, with `h=t/n`, gives the same Doeblin averaging estimate as
the source memory theorem:

```math
\|\rho_{W,n}(t)-e^{-it\bar H_W}\rho_W(0)e^{it\bar H_W}\|_1
\le {4M_W t\over n(1-p)}
 +\left({8\over1-p}+4\right){M_W^2t^2\over n}.
```

The estimate includes arbitrary initial classical-quantum memory correlations
and a finite untouched reference. The first term is absent for stationary
product initial memory. The mean is the same native mass matrix plus V because
V occurs on every branch. All disabled records and all fine outcomes remain in
the transition law. The comparison is uniform in cube radius under the source
isometries, for this fixed ninety-mode algebra. It is not a uniform statement
about all spatial observables, volume densities or lattice-spacing refinement.

The `O(1/n)` estimate compares the unread core, summing every native outcome.
The complete history register has a separate bound. For fixed initial memory
label m and a purified core/reference state psi, let omega label the full
native history. Its probability `p(omega|m)` is independent of psi. Conditional
core evolution is unitary, so let `psi_omega=U_omega psi` and
`psi_bar=exp(-it Hbar) psi`. If delta is the unread trace-norm bound above,

```math
\sum_\omega p(\omega|m)
 [1-|\langle\psi_{\rm bar}|\psi_\omega\rangle|^2]\le\delta/2.
```

The pure-state trace-distance identity and Jensen's inequality then give

```math
\left\|\sum_\omega p(\omega|m)|\omega\rangle\langle\omega|
 \otimes\left(|\psi_\omega\rangle\langle\psi_\omega|
             -|\psi_{\rm bar}\rangle\langle\psi_{\rm bar}|\right)\right\|_1
\le\sqrt{2\delta}.
```

Partial trace proves the same result for mixed inputs and retained external
references. For an initially correlated classical-quantum memory, apply it to
each conditional initial state rho_m and retain m as a classical tag. Convexity
gives the same uniform bound; the ideal comparison keeps the actual
`p(m,omega)` and applies the mean unitary to each rho_m. It does not erase
initial memory/core correlations. The smaller stationary bound is valid for
stationary product initial memory. Tracing out the tag or applying any final
instrument preserves the estimate. This gives an `O(n^-1/2)` bound for the
complete history and final core state, capped by two. A symmetric pair of
opposite qubit phase records saturates `sqrt(2 delta)`, so the unread rate
cannot simply be assigned to the retained history. No uniform error for a rare
conditioned branch or an intervening field-measurement protocol is asserted.

The independent sixteen-state control uses all 512 fine records per old memory
input, compressed only by identical field actions. Its scalar projector probe
has a noncommuting quartic term. This checks the general averaging estimate;
the probe is not the full ninety-mode physical carrier.

## Spatial range and an onsite obstruction

This interaction is localized about the prescribed texture and is not strictly
onsite. For lattice distance L, the radial theorem gives

```math
\|(I-A_L)E\|^2\le\epsilon_L
 ={3/2\over62^2}\,32^{-2(L-1)},
```

where `A_L` projects onto `|x|_infinity<L`. Hence
`||P-A_L P A_L||<=2 sqrt(epsilon_L)`. Both operators have rank at most ninety,
so their difference has rank at most 180. With `N_L=dGamma(A_L P A_L)`,

```math
\delta_L:=\|N_W-N_L\|\le360\sqrt{\epsilon_L},\qquad
\|V-\lambda(N_L-45)^2\|\le\lambda\delta_L(90+\delta_L).
```

The second inequality uses a telescoping product bound. Thus the chosen quartic
can be approximated by finite-region operators with an explicit exponential
tail, uniformly in the containing cube. No translation-invariant many-defect
interaction or continuum limit is supplied by this estimate.

There is an exact obstruction to replacing it with an onsite Hubbard term.
On the two-particle sector, let
`D_g=sum_x g_x wedge^2(P_x)` be a local density interaction. The localized
two-particle isometry is `E_2=wedge^2(E)`. Since `E*P_x E=rho_x I90`,

```math
E_2^*D_gE_2=\left(\sum_xg_x\rho_x^2\right)I,
\quad
E_2^*D_g^2E_2=\left(\sum_xg_x^2\rho_x^2\right)I.
```

The exact leakage Gram is their difference after squaring the first scalar:

```math
E_2^*D_g(I-E_2E_2^*)D_gE_2
 =\left[\sum_xg_x^2\rho_x^2-
             \left(\sum_xg_x\rho_x^2\right)^2\right]I.
```

Write `c=sum_x rho_x^2`. The source has support on more than one site, so
`0<c<1`. Cauchy-Schwarz bounds the bracket below by
`(1-c)sum_x g_x^2 rho_x^2`. It is positive for every nonzero interaction on
the supported sites. For uniform g, the norm is `|g|sqrt(c-c^2)`.
The radius-one source gives approximately 0.02705 for unit g. This excludes
exact preservation of the source two-particle subspace by that specified
onsite density class. It does not exclude a dressed low-energy subspace,
other local interaction tensors, or different carrier designs.

## Exact connected dynamics and shifted spectral poles

In a charge-preserving singular basis, the mean neutral-Higgs Hamiltonian
contains 42 disjoint charged pairing blocks and six zero neutrino modes. Select
the two spin blocks of one charged-lepton singular value m. In each of the
other 40 blocks, fix exactly one occupied mode. Pair creation and annihilation
vanish on those states. Fix three of the six zero neutrino modes occupied.
These 43 spectators are invariant under the mean Hamiltonian and under V.
The selected four modes therefore carry the exact reducing restriction

```math
H_4=m(a_1^\dagger a_2^\dagger+a_2a_1
       +a_3^\dagger a_4^\dagger+a_4a_3)
           +\lambda(N_4-2)^2.
```

This sector is not the full core vacuum. It is generally not preserved by
individual record projectors before averaging. The full sixteen-state
characteristic polynomial, verified symbolically, is

```math
z^5(z-4\lambda)
 [z^2-4\lambda z-4m^2]
 [(z-\lambda)^2-m^2]^4.
```

For `m>0`, `lambda>=0`, its ground energy is
`E0=2lambda-2sqrt(lambda^2+m^2)`. The addition-plus-removal spectral measure
of a single chosen annihilator has two poles

```math
\omega_\pm=2\sqrt{\lambda^2+m^2}-\lambda\pm m,
\quad
w_\pm={1\mp m/\sqrt{\lambda^2+m^2}\over2}.
```

Both residues are positive when lambda is positive. At lambda zero, only the
pole m has weight. The interacting core therefore cannot retain the old
quadratic coefficients as its excitation spectrum by declaration. Its ground
cross-pair density cumulant is

```math
\langle n_1n_3\rangle-\langle n_1\rangle\langle n_3\rangle
 =-{\lambda\over4\sqrt{\lambda^2+m^2}}.
```

The cross-pair normal and anomalous two-point functions vanish. A nonzero
cumulant consequently proves a non-Gaussian state and connected interaction
within the charge-zero sector. Independent Jordan-Wigner matrices check all
sixteen energies and residues, including the three actual inherited source
lepton coefficients at the declared witness `lambda=1/128`.

## A faithful reference and its noncommuting response

On the same finite core Hamiltonian `H=bar H_W`, choose beta positive and
`rho=exp(-beta H)/Z`. It is faithful, with
`lambda_min(rho)>=exp(-2 beta ||H||)/2^90`. A chosen preparation service has
Kraus operators `K_ij=sqrt(p_i)|i><j|` in its spectral basis. Their Gram sum
is identity, and their complete channel sends every input state to rho.
This is an explicit global spectral reset, not a derivation from site C/E
repair or a claim of efficient local preparation. Tensoring with the Gibbs
state of the reducing complement gives the full mean-Hamiltonian Gibbs state.
That product is not asserted invariant under each finite record pulse.

The core modular Hamiltonian is `-log rho=beta H+log Z`. At fixed H, every
traceless state tangent obeys `delta S=beta Tr(H delta rho)`. For a genuine
Hamiltonian perturbation `H(s)=H+sX`, Duhamel differentiation instead gives

```math
{d\rho\over ds}\bigg|_0
 =-\beta\int_0^1\rho^t(X-\langle X\rangle)\rho^{1-t}\,dt.
```

This formula includes noncommuting X. The logarithmic-mean kernel in the
eigenbasis of rho is positive, giving the Bogoliubov-Kubo-Mori variance
`C_rho(X,X)>=0`. Consequently

```math
F'=\langle X\rangle,\qquad F''=-\beta C_\rho(X,X)\le0,
\qquad S'=\beta(E'-\langle X\rangle).
```

The sixteen-state witness uses `X=(N4-2)^2`, which does not commute with H.
Its derivative is checked independently by the matrix-exponential Frechet
derivative, and its free-energy curvature by centered differences. The
finite state/energy/modular identities thus refer to the actual nonquadratic
Hamiltonian. X is a coupling perturbation, not an identified metric stress.
Geometric modular flow, local conserved relativistic stress, an area entropy
normalization, and fixed-volume geometric entropy stationarity are not
consequences of these finite Gibbs identities. No Einstein equation is
derived from this construction.
