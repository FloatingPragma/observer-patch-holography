# Controlled relativistic charge dispersion on a specified running family

The thermodynamic Hubbard charge excitations have a controlled relativistic
scaling limit on the explicit coupling family already derived in
`mass_generation_audit`. This package proves a finite-parameter error bound
for that dispersion, rather than calling a charge gap a mass:

\[
 E\longrightarrow M\cosh\theta,\qquad
 p\longrightarrow (M/v)\sinh\theta,\qquad
 E^2-v^2p^2=M^2,
 \qquad m_{\rm rest}=M/v^2.
\]

For the specified trajectory through the original point
`a=1,c=v=1/8,t=1/16,U=1`, the limiting values are

\[
 M_0=\frac{e^{-\pi/8}}\pi=0.2149329913552644\ldots,
 \qquad m_0=64M_0=13.75571144673692\ldots .
\]

These are model energy and mass units for a **chosen, nonunique running
trajectory in one spatial dimension**. The original fixed-coupling source
does not select this trajectory or this mass. The package gives an explicit
family through the same starting point with any chosen positive limiting
mass scale. It does not calculate Standard Model masses, a three-dimensional
particle spectrum, or a native isolated-particle pole.

## Source law, reference and momentum

The owning `mass_generation_audit` receipt proves the exact onsite-frame
join from the original periodic Weyl/contact chain to the standard Hubbard
chain. The charge reference is its half-filled canonical ground sector.
Here the excitation generator is explicitly `K=H-(U/2)N`, with the irrelevant
constant `UL/4` omitted. Original-H charge addition/removal energies carry
the corresponding `+/-U/2` shift. This centering is a reference/reservoir
condition; it is not a silent change of the original `mu=0` ground.

The [primary exact dispersion](https://arxiv.org/html/cond-mat/9410043),
equations (2.1)–(2.3), is written for the attractive chain with hopping
`1/2`. The staggered particle-hole transformation

\[
 d_{x,\uparrow}\mapsto d_{x,\uparrow},\qquad
 d_{x,\downarrow}\mapsto(-1)^x d^\dagger_{x,\downarrow}
\]

preserves hopping and reverses the centered interaction. The attractive
spin doublet therefore becomes the repulsive charge doublet. Multiplying
the source normalization by `2t` gives attractive parameter `u=U/(4t)`.
`source.py` checks all these identities on every state of the 256-dimensional
four-site Fock space, in addition to the local proof for even bipartite
rings. The inherited original-frame join uses lengths `L=4m`.

Momentum also requires its reference convention. For ordinary Hubbard
translation `T`, the same exact witness gives
`P T P^dagger=-(-1)^N_down T`. Its constant phase cancels from momentum
differences; the staggered map shifts them by `pi Delta N_down` modulo
`2pi` on admissible whole states. The earlier onsite SU(2) frame conjugates
ordinary Hubbard translation to source translation combined with a uniform
internal rotation. The momentum below is measured relative to the chosen
low-energy charge branch, with these carrier offsets subtracted. It is not
identified silently with an unshifted raw source translation eigenvalue.

The even-ring SO(4) selection rules remain in force. An elementary
thermodynamic quasiparticle dispersion is not, by itself, a permitted
isolated finite-ring state or a single local-CAR pole. No ground preparation,
complete repair invariance, or native charge detector is executed here.

## Exact dispersion and finite-parameter control

Set `A=2pi t/U`, `theta=A sin k`, `|k|<pi/2`, and `v=2ta`. After the
explicit normalization and duality above, the exact dressed dispersion is

\[
 E=\frac{8t}\pi\sum_{\substack{m\ge1\\m\;\mathrm{odd}}}
       \frac{K_1(mA)}m\cosh(m\theta),\qquad
 vp=\frac{8t}\pi\sum_{\substack{m\ge1\\m\;\mathrm{odd}}}
       \frac{K_0(mA)}m\sinh(m\theta).
\]

The primary dispersion is an imported thermodynamic Bethe-ansatz result;
this receipt does not claim to reprove that theory. The finite bounds below
are derived independently and evaluated with exact rational arithmetic.

For `b(z)=sqrt(pi/(2z))exp(-z)`, the integral representation gives
`K_nu(z)/b(z)` as the expectation of `(1+S/(2z))^(nu-1/2)` under a
`Gamma(nu+1/2)` variable. Tangent inequalities yield

\[
 1-\frac1{8z}\le\frac{K_0(z)}{b(z)}\le1,
 \qquad 1\le\frac{K_1(z)}{b(z)}\le1+\frac3{8z}.
\]

With `A>=1`, `|theta|<=Theta<A` and `q=exp[-2(A-Theta)]`, the odd-mode
tail is bounded geometrically. In particular,
`cosh(m theta)<=cosh(theta)exp[(m-1)|theta|]` and
`|sinh(m theta)|<=m|sinh(theta)|exp[(m-1)|theta|]`.
On the specified family

\[
 a(A)=\frac{4v e^{-A}}{M\sqrt{2\pi A}},\qquad
 t=\frac{v}{2a},\qquad U=\frac{\pi v}{aA},
\]

the leading coefficient `(8t/pi)b(A)` is exactly `M`. Define

\[
 r_E=\frac{1+3/(8A)}{1-q}-1,\qquad
 r_p=\max\!\left\{\frac1{8A},\frac q{1-q}\right\}.
\]

Then uniformly on that rapidity interval,

\[
 |E/M-\cosh\theta|\le r_E\cosh\Theta,\qquad
 |vp/M-\sinh\theta|\le r_p\sinh\Theta,
\]
\[
 \left|\frac{E^2-v^2p^2}{M^2}-1\right|
 \le\cosh^2\Theta(2r_E+r_E^2)
    +\sinh^2\Theta(2r_p+r_p^2)\longrightarrow0.
\]

For `Theta=1`, the last upper bound is below `0.570,0.272,0.135,0.068,
0.034,0.017` at `A=4,8,16,32,64,128`, respectively. Thus this is uniform
dispersion and mass-shell control, not just a zero-momentum threshold limit.

`point()` additionally encloses particular finite-A dispersion values.
It integrates global alternating Taylor bounds for the same Gamma
expectation, with coefficient recurrence
`c_n/c_(n-1)=[4nu^2-(2n-1)^2]/(8nz)`. Adjacent finite truncations bound
the integral even when the infinite asymptotic series does not converge.
Eight odd terms and an explicit geometric remainder complete the enclosure.
All elementary functions use inherited rational interval primitives;
floating Bessel evaluations are independent tests only.

## What the numerical mass does and does not select

For `A0=pi/8,v=1/8,M0=exp(-pi/8)/pi`, exact substitution gives
`a(A0)=1,t(A0)=1/16,U(A0)=1`. The limiting rest mass of this **specified
trajectory** is therefore the exact model number `64 exp(-pi/8)/pi` above.
It differs from the starting lattice half-gap `0.3857758578690188...`.
Independent integral/Bessel diagnostics at that starting point show
`E^2-v^2p^2` changing from approximately `0.1488230` at `theta=0` to
`0.1763631` at `theta=0.3`; the starting dispersion is not being called
relativistic. These diagnostic decimals have an explicit floating admission,
separate from the exact rational scaling certificates.

The continuation is not unique. For any `M>0`, multiplying `a_M(A)` above by

\[
 f_M(A)=1+(M/M_0-1)e^{-(A-A_0)},\qquad A\ge A_0,
\]

while retaining `t=v/(2a)` and `U=pi v/(aA)` gives the very same starting
point. This factor stays positive and tends to one, so the limiting rest
energy is `M`. This exact counterfamily prevents inference of a unique
mass from the bare point alone. In contrast, the previous fixed-U
velocity-preserving continuum family has zero limiting charge mass.

Run `.venv/bin/python -W error -m pytest model/hubbard_relativistic_sector -q`
and `.venv/bin/python -W error -m model.hubbard_relativistic_sector --check`.
The receipt binds the owning source-law audit, new proofs, code, tests and
paper fragment; it reports exact rational bounds separately from diagnostics.
