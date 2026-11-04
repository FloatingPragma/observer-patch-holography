# Bulk mass shells and exterior channels of the bounded radial source

This construction uses the hopping and onsite mass channels of
`model/native_mass_channels` and the bounded radial operator of
`model/bounded_defect_index`. The one-body clock has its original sign and
energy origin. The supplied parameters are the lattice spacing `a>0`,
velocity `v>0`, and mass coefficients. Momentum is an angular wave number
in inverse model length. Evolution is `exp(-i H_model t_model)`.
The [calibration proof](CALIBRATION_PROOF.md) states the physical-unit map.

## Complete period-two Bloch operator

For lattice sites `n in Z^3`, let

```
h_(n,n+e_j) = i v sigma_j/(2a),
B_n = m0 (-1)^(n1+n2+n3) I2 + sum_j m_j (-1)^n_j sigma_j.
```

These are the existing four onsite channels, each anticommuting with every
kinetic generator. The three spin masses anticommute pairwise and commute
with the scalar mass. The profile is
constant in its four coefficients and periodic over two lattice sites.
Fold the eight microscopic momenta `a p+pi b`, `b in {0,1}^3`, into one
reduced zone. With the source convention `exp(i k.n)`, its sixteen-row
matrix is exactly

```
H_a(p) = (v/a) sum_j sin(a p_j) gamma_j + m0 beta0 + sum_j m_j beta_j,
gamma_j |b> = -(-1)^b_j |b> sigma_j,
beta0 |b> = |b xor111> I2,
beta_j |b> = |b xor e_j> sigma_j.                         (1)
```

No valley is discarded. The exact finite-source witness constructs all
128 Fourier columns of the periodic four-by-four-by-four source matrix,
including all four onsite masses, and checks (1) on every column. The
real-space recursion proves the formula at arbitrary momentum.

Write `r=|m|`. For `r>0`, the operator
`Q=beta0 (sum_j m_j beta_j)/r` is a Hermitian involution of trace zero.
It commutes with every `gamma_j` and with `beta0`. Since
`sum_j m_j beta_j=r beta0 Q`, its two rank-eight sectors have signed
mass `mu_s=m0+s r`, `s=+1,-1`. The source Clifford identities give

```
H_a(p)^2 |Q=s = [(v/a)^2 sum_j sin(a p_j)^2 + mu_s^2] I8,
E_a,s,±(p) = ±sqrt((v/a)^2 sum_j sin(a p_j)^2 + mu_s^2). (2)
```

The signed-energy projectors are
`P_s,±=(I+sQ)(I±H_a/E_a,s,+)/4` when `E_a,s,+>0`.
Their ranks are four. If `E_a,s,+=0`, the entire rank-eight Q sector is
the zero eigenspace; the expression dividing by energy is not used.
For `r=0` the spectrum has one positive and one negative rank-eight
level, or one rank-sixteen zero level. Coincident energies from different
sectors have their multiplicities added. These formulas retain the full
spectral resolution, including the tuned `m0=r` case with two massive
and two massless Dirac copies.

## Quantified relativistic continuum dispersion

At fixed physical-model momentum, (1) tends to

```
H(p) = v sum_j p_j gamma_j + m0 beta0 + sum_j m_j beta_j.
```

Within each Q sector, set `Gamma0=beta0` and
`Gammaj=beta0 gamma_j`. Their anticommutator is
`{Gamma_mu,Gamma_nu}=2 diag(1,-1,-1,-1)_mu,nu I8`.
The equation is the massive Dirac equation in the coordinate `x0=v t`,
with signed mass coefficient `mu_s/v`. Each sector is two four-row
Dirac representations. Its positive rest energy and inertial model mass
are `M_s=|mu_s|` and `M_s/v^2`. The continuous limiting dispersion is

```
E_s,±(p) = ±sqrt(v^2 |p|^2+M_s^2),
E_s,±(p)^2-v^2 |p|^2=M_s^2.                            (3)
```

For `|a p_j|<=1`, the elementary inequalities
`0<=x^2-sin(x)^2<=x^4/3` and `|sin x-x|<=|x|^3/6`
give the uniform signed shell defect

```
-(v^2 a^2/3) sum_j p_j^4
 <= E_a,s,±^2-v^2 |p|^2-M_s^2 <=0,                    (4)
|E_a,s,±-E_s,±| <= (v a^2/6) sum_j |p_j|^3.           (5)
```

The last bound follows from the reverse triangle inequality in the
four-vector `(v sin(a p)/a,mu_s)` and is uniform through zero mass and
zero momentum. If `M_s>0`, the stronger alternative
`v^2 a^2 sum_j p_j^4/(6 M_s)` follows by dividing (4) by
`E_a,s,+ + E_s,+ >=2 M_s`. For `|p_j|<=Lambda`, (4) is bounded by
`v^2 a^2 Lambda^4` and (5) by `v a^2 Lambda^3/2`.
The finite-spacing group velocity obeys `|grad_p E_a,s,+|<=v` wherever
it is differentiable: its squared numerator is
`(v^4/a^2) sum_j sin^2(a p_j) cos^2(a p_j)`, while its squared
denominator is at least `(v^2/a^2) sum_j sin^2(a p_j)`.
The equation and mass shell are Lorentz covariant in the continuum.
The lattice dispersion and a fixed momentum cutoff are not invariant
under all Lorentz boosts.

The energy calculator encloses all roots and sines by rational intervals.
The sine series uses adjacent alternating partial sums on `[-1,1]`;
integer square roots give outward dyadic root bounds. An irrational
length of the supplied spin-mass vector is retained by interval
arithmetic. Numerical eigensolvers serve only independent controls.

## The same bounded radial source

The supplied pinned radial operator has the same hopping and

```
B_n^rad = m sum_j [n_j/|n|] (-1)^n_j sigma_j  (n!=0),
B_0^rad = 0,             m>0.                          (6)
```

It is not periodic. Every nonzero lattice translation moves its unique
zero onsite block. If `T_z psi(n)=psi(n-z)`, then the origin block of
`T_z H_rad T_z* - H_rad` is `B_-z`, whose norm is `m`; hence

```
||T_z H_rad T_z* - H_rad|| >=m,             z!=0.        (7)
```

The translation-invariant hopping cancels. In particular an even-cell
translation does not preserve the pinned profile. A Fourier probe of a
localized charged state has that state's fixed excitation energy. It
does not obtain a moving-defect band from its probe label.

The existing radial theorem gives the essential spectrum

```
sigma_ess(H_rad) = [-sqrt(m^2+3v^2/a^2),-m]
                    union [m,sqrt(m^2+3v^2/a^2)].       (8)
```

Its compactness argument is explicit. Put `u(n)=n/|n|` for `n!=0`.
The staggering identity on a forward edge cancels the constant part of
the hopping-mass anticommutator and gives the block

```
{h0,B^rad}_(n,n+e_j)
 = i m v/(2a) sum_l (-1)^n_l [u_l(n)-u_l(n+e_j)] sigma_l sigma_j.
```

Its coefficients tend to zero at infinity. A bounded finite-range matrix
whose coefficients vanish at infinity is compact by finite-site
truncation. Since `(B^rad)^2=m^2 I` away from the origin,
`H_rad^2=h0^2+m^2 I+compact`. The Fourier range of
`h0^2=(v^2/a^2)sum_j sin^2(k_j) I2` is `[0,3v^2/a^2]`.
Spatial inversion anticommutes with `H_rad`, supplying both signs in
(8). This statement includes the original profile rather than replacing
its direction by a constant.

## Explicit compact exterior packets

The following construction supplies a momentum-resolved approximate
eigenvector family for (6). It describes bulk excitations arbitrarily
far from the pinned core; it assigns no kinetic energy to the core.

Choose positive integer `L`, positive even integer `R`, and anchor
`n0=(0,0,R)`. Write the support as `n=n0+2r+s`, where
`r_j=0,...,L-1` and `s_j in {0,1}`. Define the normalized cell envelope

```
f_L(r)=product_j sqrt(2/(L+1)) sin(pi(r_j+1)/(L+1)),
```

with zero extension. Let `u=(u_b)_b` be any unit sixteen-row eigenvector
of the exact frozen-direction matrix (1) at coefficients `(0,0,0,m)`,
with signed eigenvalue `E_a,±(p)`. Set

```
psi(n0+2r+s)
 = f_L(r) exp(i a p.(n0+2r+s))
   [sum_b (-1)^(b.s) u_b]/sqrt(8).                     (9)
```

The eight-by-eight Walsh transform is unitary, so `||psi||=1` exactly.
The finite coefficient vector contains every microscopic valley. Its
support has `(2L)^3` physical sites; all one-hop boundary output rows
are included when evaluating `H_rad psi`.

Let `H_z` be the periodic source with `B_n^z=m(-1)^n3 sigma3`.
The zero-extended one-dimensional sine vector obeys

```
||f_L(.+e_j)-f_L|| = 2 sin(pi/[2(L+1)]) <=pi/(L+1).
```

For example, its scalar product with its shift is
`cos(pi/(L+1))`, from the Dirichlet path adjacency eigenvalue;
its squared shift difference is therefore `2-2cos(pi/(L+1))`.
A one-site hop either preserves its cell coordinate or shifts that
coordinate once. In the commutator with the envelope, each shift has
norm at most the displayed difference: summing the squared cell-spin
weights on the shifted parity face is at most one. Its coefficient is
`v/(2a)`. Adding the six shifts bounds the kinetic envelope residual by
`3 pi v/[a(L+1)]`. The uncut periodic wave is an eigenvector, so

```
||(H_z-E_a,±)psi|| <=3 pi v/[a(L+1)].                  (10)
```

Every support point has `|n-n0|<=2sqrt(3)L`. For nonzero vectors,
`|x/|x|-y/|y||<=2|x-y|/|y|`. The three staggered Pauli matrices
anticommute, so the norm of their linear combination is exactly the
Euclidean coefficient norm. Consequently

```
||(B_rad-B_z)psi|| <=4sqrt(3)mL/R,
||(H_rad-E_a,±(p))psi||
 <=3 pi v/[a(L+1)]+4sqrt(3)mL/R.                      (11)
```

For fixed `a`, take `L->infinity`, `R/L->infinity`, with supports
escaping to infinity. The normalized packets converge weakly to zero
and their residuals vanish. They are Weyl sequences for every energy
(2) with pure spin mass `m`, on the actual radial operator.
This proves the stated asymptotic bulk channels. It does not construct
time-dependent scattering wave operators or a momentum operator
commuting with the entire radial Hamiltonian.

Combining (5) and (11) gives a direct continuum-shell residual bound:

```
||(H_rad-E_±(p))psi||
 <=3 pi v/[a(L+1)]+4sqrt(3)mL/R
   +min(v a^2 sum_j |p_j|^3/6, v^2 a^2 sum_j p_j^4/(6m)),
E_±(p)=±sqrt(m^2+v^2 |p|^2).                          (12)
```

For example `a=1/n`, `L=n^3`, `R=2n^8` makes (12) vanish at fixed
`p,v,m`. Both the physical packet width and its distance from the
defect grow. The finite bound calculator can certify much larger
supports than the explicitly executed numerical controls; those
large support counts are resource requirements, not executed arrays.

## Scope of the mass interpretation

Equations (2)--(5) are exact one-body and free quadratic field results
for a supplied periodic mass family. Equations (8)--(12) attach its
pure-spin bulk shell to exterior packets of the same pinned radial
one-body source. They retain negative branches and all sixteen folded
rows. A free filled-sea particle/hole reference gives positive free
excitation energies, but that state is not substituted for the
interacting radial ground.

The local contact interaction vanishes in the uncentered N=1 sector.
This fact does not identify the one-body exterior spectrum with the
charged poles above the interacting N=27 radial ground. In a centered
finite-volume N=1 sector the contact adds a constant `-U/2` relative
to the empty sector, which must be retained if that reference is used.
No such sector shift is made here. The finite strong-mass interacting
core theorem requires `m>sqrt(3)v/a`; fixed `m,v` violate this
condition at small spacing. Its charged-pole theorem is therefore not
continued along (12).

There is no dynamical texture coordinate, moving-defect rest mass,
physical species identification, calibrated physical time unit, or
interacting relativistic pole theorem in this calculation. The tuned
periodic massive-plus-massless case and the gapped radial mass case
are different source profiles. A localized protected zero mode and
a massless bulk are not asserted in one common regime.
