# Pixel closure with the same charged-fermion coefficients

## Declared response and scale

Write `u=alpha_U`, `P` for pixel capacity and `E_*` for a declared energy unit.
The scale design is

```math
M_U/E_*=e^{-2\pi}P^{1/6},\qquad
v/E_*=P^{-1/2}e^{-2\pi/(\beta u)},\qquad
\alpha_i^{-1}=u^{-1}+\frac{b_i}{2\pi}\log(M_U/m_Z),
\quad m_Z/v=\sqrt{\pi(\alpha_2+3\alpha_1/5)}.
```

The pixel equation is the retained infinite heat-sum equation
`ell_SU2(4 pi^2 alpha2)+ell_SU3(4 pi^2 alpha3)=P/4`, with weights
`dim(R) exp(-t C2(R))` and observable `log dim(R)`.
For the minimal three-family, one-complex-doublet census, exact representation
indices give `b=(41/10,-19/6,-7)`. They are computed directly as
`Tr(C3)/8`, `Tr(C2)/3` and `Tr(Y^2)` from the actual fifteen-state
current identification operators. The complete joint charge/Casimir blocks
are checked against all fifteen native rows. The supersymmetric census gives
`b=(33/5,1,-3)` and requires its extra scalar and fermion fields. These
coefficients use the declared four-dimensional one-loop functional; finite
representation arithmetic alone is not a derivation of that functional.

The charged-lepton coefficients are exactly the conditional minimal complete
public response law realized by `model/mcpr_native_memory_bridge`:

```math
\begin{aligned}
\kappa&=\frac{u/50}{1+u/31+u^2/310},&
\chi&=\frac{u^2/512}{1-u/77},&
\zeta&=\frac{\alpha_P^2/21}{1-u/27-u^2/135},\\
\alpha_P&=(P-\phi)/\sqrt\pi,&
s_j&=[1+\sqrt2 e^{-\chi}\cos(2/9+\zeta+2\pi j/3)]^2,\\
t_j&=s_j\left[\frac{e^{3\kappa}}{2\,6^{14}\prod_k s_k}\right]^{1/3}.
\end{aligned}
```

The three ordered masses are `m_j=v t_j`. Quark coefficients are
`(6^-6,6^-3,1)/sqrt(2)` and `(6^-6,6^-4,6^-2)/sqrt(2)`. The same coefficients
therefore enter the field design, its memory realization and this response.
The architecture has target-informed provenance; its execution reads no measured
mass, measured coupling, GeV calibration, or stored pixel root.

For all three charged leptons and **all six** quarks, including top, define

```math
K_f=\frac{2N_cQ_f^2}{\pi}\int_0^1x(1-x)
\log\!\left[1+\frac{m_Z^2}{m_f^2}x(1-x)\right]dx,
\quad
A=\alpha_2^{-1}+\frac53\alpha_1^{-1}
 +\sum_{\ell}K_\ell+(1-3\alpha_3/\pi)\sum_qK_q.
```

The closed loop is `alpha -> P=phi+sqrt(pi) alpha -> (u,mZ,v,t_j) -> A -> 1/A`.
The quark multiplier is an inherited prescribed factor. This is a declared
charged-fermion response kernel. It supplies no charged-vector or charged-scalar
terms, supersymmetric threshold matching, confining hadronic spectral measure,
physical pole-mass convention, or native-current derivation of the Thomson limit.
Calling it the complete electromagnetic running of the Standard Model would be
incorrect.

## Exclusion of the minimal census with scale coefficient four

Assume `P>=1`, `M_U>=m_Z`, positive running couplings, the minimal census above,
no intermediate thresholds, and `beta=4`. Put `x=alpha2^-1`, `y=alpha3^-1`,
`L=log(M_U/m_Z)` and `r=m_Z/v`. Eliminating `u` and `L` gives

```math
0<y\le x,\quad
u^{-1}=\frac{42x-19y}{23},\quad
\frac{L}{2\pi}=\frac{6(x-y)}{23},\quad
\alpha_1^{-1}=\frac{333x-218y}{115},
```

The scale equation consequently reads

```math
\frac{\pi(18x+5y)}{46}
=2\pi-\frac23\log P+\log r,\qquad
r^2=\pi\left(\frac1x+\frac{69}{333x-218y}\right)
\le\frac{8\pi}{5x}.
```

Thus `F(x)=9 pi x/23-2 pi-(1/2)log[8 pi/(5x)]<=0`. But
`F'(x)=9 pi/23+1/(2x)>0` and the outward interval calculation proves
`F(6)>1`. Hence `x<6`, and `0<y<=x<6`; both nonabelian couplings exceed `1/6`.
For either heat sum its denominator is at least one. Its positive numerator
at any `t>=4 pi^2/6` is bounded above by the numerator at `4 pi^2/6`.
Direct outward summation plus geometric infinite-tail bounds gives a combined
upper bound below `1/80`. Pixel closure requires at least `1/4`, a contradiction.
This proves nonexistence on the entire stated domain, independently of the
outer readout kernel. It does not exclude another scale law, intermediate
thresholds, a different census, or OPH as a class.

## Two certified conditional fixed points

The first design retains the supersymmetric running ansatz and `beta=4` as a
proxy; its extra supersymmetric fields are not instantiated in the native
quadratic carrier. The second uses the minimal census and chooses `beta=7`,
the magnitude of its color coefficient. This scale choice is explicit;
no theorem identifies the Higgs scale with that color transmutation scale.

| Design | `P` | `A=alpha^-1` | `u` | `v/E_*` |
| --- | --- | --- | --- | --- |
| Supersymmetric running proxy, `beta=4` | 1.630971514159316 | 137.000994766342 | 0.041124260725897 | 2.01983668924449e-17 |
| Minimal census, color scale `beta=7` | 1.628533369611833 | 168.815082928453 | 0.019531387695090 | 8.61718526720080e-21 |

`closure.py` extends the source interval backend without modifying its historical
calculations. Each implicit root has verified endpoint signs and a sign-definite
partial derivative. The infinite heat tails are enclosed. Forward interval
automatic differentiation gives `|g'|<0.00734` and a strict self-map on each
reported alpha interval of radius `10^-14`. Banach's theorem proves a unique root
**in that interval**; no global uniqueness assertion is made. Every running
coupling and response denominator in the enclosure is positive.

`replay.py` independently solves three simultaneous equations in
`(1/u,log(M_U/m_Z),alpha)`, without importing the producer or its implicit-root
code. Direct quadrature of each fermion integral checks the algebraic kernel.
This is an independent numerical crosscheck, separate from the interval proof.

## Native memory on the full root enclosure

Normalize `lambda_j=t_j/sum(t)` and let `e2,e3` be their elementary symmetric
functions. The three projector weights are roots of

```math
z^3-z^2+\frac54e_2z-\frac{e_3}{d}=0,
\qquad d=\frac25-\frac{2\sqrt5}{25}.
```

For ordered roots `w<u_w<v_w`, assign weights `(u_w,v_w,w)` to the existing
port projectors `(P2,P3,P10)`. The complete C-duty inverse satisfies

```math
-v_w(v_w+w)p^2+[v_w(1-v_w+w)+(v_w+w)w]p-w=0,
\quad q=1-\frac{v_w}{1-p(v_w+w)},\quad g=\frac{\sum_jt_j}{p}.
```

`memory.py` uses rational interval arithmetic on the full certified mass boxes.
Three disjoint sign-changing cubic brackets and two disjoint quadratic brackets
exhaust the roots by degree. Both solutions obey `0<p<1/12`, `0<q<1/6`,
and the total mass coefficient is strictly below one. The generic exact
projector identities therefore give `spec(g Z_C)={t_j}` for every enclosed root,
with all disabled records retained. Full 512-record sums on twelve old-memory
inputs independently check each displayed center controller. The inherited
weak-pulse theorem applies with the enclosed `p,g`; it is a finite-volume theorem.

At the color-scale root the two center controllers are approximately
`(p,q,g)=(0.01915345877,0.05474578848,0.39849962002)` and
`(0.05998935985,0.01675156563,0.12723333039)`.
These controllers replace the rates attached to another pixel coordinate;
reusing those rates without recomputation would not realize the same mass law.
The resulting conditional finite operator is compatible with the generic
localized current/mass inclusion, but physical dispersion, vacuum, autonomous
defects, gauge response, and Einstein coarse-graining are not proved here.
