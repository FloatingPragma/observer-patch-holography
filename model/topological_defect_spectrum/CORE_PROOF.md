# A controlled interacting finite defect core

The source has its original two onsite fermion modes and centered contact
interaction. The quadratic mass profile is the supplied linear staggered
spin texture of `lattice.py`. Its preparation and selection are not derived
from that contact interaction. All claims in this proof concern the finite
27-site cube with coordinates -1,0,1 and parameters v=kappa=a=1.

## Exact one-body domain and the full Fock gap

After the explicit onsite phase of the lattice construction the one-body
operator is `h=sum_(j=1)^3 B_j tensor sigma_j`, with commuting real spatial
matrices B_j acting on the three coordinate factors. In one dimension,

```
B = [[1,-1/2,0],[-1/2,0,1/2],[0,1/2,-1]],
phi=(1,2,1)/sqrt(6),      B phi=0,
spec(B)={0,+sqrt(3/2),-sqrt(3/2)}.
```

Thus `h^2=sum_j B_j^2` has exactly two zero modes
`phi(x1)phi(x2)phi(x3) xi`, with arbitrary two-spin xi, and every other
one-body energy has magnitude at least `Delta=sqrt(3/2)>6/5`.
There are 26 negative and 26 positive levels. Fill all negative levels,
leave all positive levels empty, and allow both zero modes arbitrary
occupancy. After subtracting the fixed filled-sea energy, the complete
Fock Hamiltonian H0 has a four-dimensional zero eigenspace P and

```
H0 >= Delta (1-P).
```

This inequality concerns the full `2^54` carrier, by summing occupation
energies of the exact quadratic source. Constructing that full matrix is
unnecessary. The four states have particle numbers 26,27,27,28.

The original contact is

```
W=U sum_x (n_(x,up)-1/2)(n_(x,down)-1/2),  U>0,
||W||=w=27U/4.                                           (1)
```

This norm follows from the commuting onsite eigenvalues +/-U/4. The
contact is unchanged by the onsite scalar phase map.

## The exact sea contribution and compressed interaction

Let `rho_x=|phi(x1)phi(x2)phi(x3)|^2`. The inverse participation is

```
I=sum_x rho_x^2=[(1+16+1)/36]^3=1/8.                    (2)
```

Let R_- be the negative one-body spectral projection, with onsite
two-spin block R_-(x,x). Spatial inversion S anticommutes with h, and
the zero orbital is even. Consequently `S R_- S=R_+` and
`R_-+R_+=I-P0`. Weighted by the even rho_x, these identities give

```
sum_x rho_x R_-(x,x)=sum_x rho_x R_+(x,x)
                    =(1-I) I_2/2.                     (3)
```

This is a complete negative-sea identity, without dropping occupied modes.
Write the zero-mode annihilators as f1,f2. Wick expansion against the
filled sea gives for the local contact a constant, the quartic term
`U rho_x^2 n1 n2`, and the quadratic core matrix

```
U rho_x { [Tr R_-(x,x)-1/2]I_2-R_-(x,x) }.
```

Using (3), the exact four-state compression is therefore

```
PWP=C+g[n1 n2-(n1+n2)/2],   g=U I=U/8,                 (4)
C=U sum_x [det R_-(x,x)-Tr R_-(x,x)/2+1/4].
```

The empty and full core levels are C, while each singly occupied level
is `C-g/2`. Their first-order charged excitation gap is `U/16`.
The core interaction has a nonzero quartic coefficient. This is an
interacting finite-source statement, without assigning the levels to
observed particle species.

## Full-carrier spectral error and exact degeneracies

Assume `w<Delta/4` and use the conservative rational lower bound 6/5 for
Delta in every executable estimate. The four lowest eigenvalues of
`H0+W` lie in [-w,w]; its remaining spectrum is at least `Delta-w`.
For any low energy lambda, the complement block satisfies

```
Q(H0+W-lambda)Q >= Delta-2w,
F(lambda)=PWP-PWQ[Q(H0+W-lambda)Q]^(-1)QWP,
0<=PWP-F(lambda)<=epsilon P,
epsilon=w^2/(Delta-2w).                                  (5)
```

Schur-complement inertia, or the ordered min-max principle, puts the
ordered four low eigenvalues in `[a_j-epsilon,a_j]`, where a_j are the
ordered eigenvalues of PWP. The weaker symmetric enclosure +/-epsilon
also suffices. The same argument within each conserved particle-number
sector labels the two single-core and two charged levels.

The degeneracies are exact at finite U. Let `T_s=i sigma2 K` on spin;
it changes the sign of every sigma_j. The antiunitary `Theta=S T_s`
commutes with h, satisfies `Theta^2=-1` on one-body states and preserves
the contact. On Fock space its square is `(-1)^N`. The N=27 low sector
is therefore an exact Kramers doublet. This uses spatial inversion and
the supplied texture, without asserting ordinary spin time-reversal
invariance of that texture alone.

The onsite particle-hole transformation
`c_up -> c_down^dagger`, `c_down -> -c_up^dagger` also preserves the
transformed quadratic Hamiltonian: for `J_s=i sigma2`,
`J_s h^T J_s^dagger=-h`, and `Tr h=0`. It preserves each centered
onsite contact and sends N to 54-N. Hence the N=26 and N=28 low levels
are exactly equal. Conjugation by the phase map gives the corresponding
symmetries in the original source frame.

If `U/16>2epsilon`, the actual ground is the N=27 doublet, the lowest
charged states form the particle-hole-related pair, and their excitation
gap obeys

```
U/16-2epsilon <= gap_charged <= U/16+2epsilon.            (6)
```

At `U=1/65536`, w=27/262144 and `Delta_lower=6/5`, both hypotheses
hold by rational arithmetic and the lower endpoint is strictly positive.
Equation (5) retains all virtual states in the full Fock complement.

An independent three-site one-dimensional control keeps all 64 Fock
states. It has I=1/2 and projected gap U/4. Tests build its Jordan-Wigner
Hamiltonian, its filled-sea Slater states and the entire centered contact,
then compare the exact compressed matrix, ordered low levels, sector
labels and excitation gaps. A separate 54-dimensional one-body control
checks the complete cube's zero pair, gap and weighted sea matrix.

## Scope

This is a controlled interacting defect core in a finite supplied texture.
It does not derive the texture from a self-consistent source ground,
calibrate U or the clock, prove a thermodynamic defect gap, give an
interacting relativistic field limit, or assign a physical particle.
The continuum radial-hedgehog spectrum and the finite interacting cube
have distinct profiles and domains; their masses are not identified.
