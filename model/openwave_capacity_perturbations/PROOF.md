# Perturbations and a canonical completion of the capacity history

Canonical reconstruction of holographic dark energy is an established
method; see Xin Zhang, [Reconstructing holographic quintessence](https://arxiv.org/abs/astro-ph/0604484).
The negative adiabatic-speed control also has precedent in Yun Soo Myung's
[Instability of holographic dark energy models](https://arxiv.org/abs/0706.3757).
The contribution here is the exact radiation-inclusive calculation for the
declared OPH comparison branch, its full-past bound and the explicit
Bernstein/Lean curvature certificate. General scalar reconstruction is not
claimed as a new principle or a native-source derivation.

Use the separately conserved ideal-radiation/matter/event-area background
and its supplied density amplitudes. Write x=sqrt(Omega_x),
r=rho_r/(rho_m+rho_r) and prime=d/dlog(a). The parent proves

```
x'=x(1-x²)(1+r+2x)/2,   r'=-r(1-r),   w_x=-(1+2x)/3.
```

## The adiabatic perfect-fluid control

Separate conservation gives
c_a²=dot(p)/dot(rho)=w-w'/[3(1+w)]. Exact substitution gives

```
c_a²=[-2-3x+3x²+2x³+r x(1+x)]/6.
```

This polynomial increases with r and is convex in x on the unit square.
For0<=z<=1 the supplied current fraction bounds x<5/6 and the conserved
radiation ratio gives r<1/1000. For all z>=1, monotonicity of the parent
branch and its certified z=1 value give x<11/20, while r<=1. Endpoint
evaluation of each convex polynomial therefore proves c_a²<-1/5 on the
entire past branch, including beyond the numerical integration endpoint.
No guessed perturbation speed enters this background identity.

An adiabatic barotropic perfect-fluid attachment identifies the local
rest-frame speed with c_a². Linearized local energy and momentum
conservation then give delta_rho_tt=c_a² Delta delta_rho. The complete
finite periodic comparison retains its constant mode and every positive
Laplacian mode; the latter have growing and decaying real exponents when
c_a²<0. A declared damped frozen-coefficient control illustrates that a
bounded Hubble friction cannot remove the high-wavenumber gradient sign.
It is not a computed CMB or structure-formation likelihood. This excludes
the stated adiabatic fluid attachment, not every capacity source or OPH.

The conserved Hubble-area control instead has
c_a²=4r/[3(3+r)] in[0,1/3], while retaining its separately failed background
comparison. Fixed capacity has dot(rho)=dot(p)=0; their ratio is undefined,
so assigning it a negative sound speed would be invalid.

## Constructive canonical scalar completion

There is a local canonical scalar action reproducing the declared background.
With reduced Planck mass M_Pl, define

```
phi'=M_Pl sqrt[2 Omega_x(1-sqrt(Omega_x))],   phi(1)=0,
V=rho_x(2+sqrt(Omega_x))/3,
dot(phi)²/2=rho_x(1-sqrt(Omega_x))/3.
```

The field increases strictly with a on0<Omega_x<1, so these expressions
define a single-valued parametric V(phi). They give exactly rho_x and p_x;
the original conservation equation implies the Klein-Gordon equation since
dot(phi)>0. The canonical rest-frame propagation speed is1. Its effective
nonadiabatic pressure coefficient1-c_a² is therefore positive, supplying
the stress response missing in the adiabatic fluid control.

Twice differentiating the potential with respect to phi, using the same
x' and r', gives the polynomial V''/H² in the receipt. Its complete
bivariate Bernstein coefficients of degrees(6,2) are all between2 and6.
The positive Bernstein basis sums to one, hence

```
2 <= V''/H² <= 6   for 0<=x<=1 and 0<=r<=1.
```

The executable reconstructs the polynomial from every coefficient by exact
rational arithmetic. `CapacityScalarCurvature.lean` checks both bounds
through explicit polynomial identities and positivity. The fixed-background
canonical quadratic kinetic, gradient and potential terms are consequently
positive. Gravitationally constrained cosmological perturbations and their
observational fits are not inferred from this local scalar statement.

In radiation domination, Omega_x proportional to a² gives
phi-phi_early proportional to a and V proportional to a^-2, hence an
inverse-square potential asymptotic. In the separate dust-dominated
small-Omega asymptotic, phi-phi_reference proportional to sqrt(a) gives an
inverse-quartic potential. These limits concern the stated respective
background regimes, not two simultaneous early-time limits.

This is a constructed scalar completion of a supplied background history.
It uses the parent's density amplitudes and horizon choice, with no extra
empirical fit. It does not derive a scalar particle from native observers,
prove its transfer from the finite capacity register, establish a quantum
ultraviolet completion, or validate a CMB/structure likelihood. Those source
and physical identifications remain explicit completion requirements.
