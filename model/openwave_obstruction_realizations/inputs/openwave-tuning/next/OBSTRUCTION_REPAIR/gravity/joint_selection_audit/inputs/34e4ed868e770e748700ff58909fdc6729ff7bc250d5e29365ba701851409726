# Higgs curvature coupling, stress and regional variation

The scalar in the conditional Standard Model action is the same complex
doublet H used here. Write X=H†H, use signature(-,+,+,+) and hbar=c=1,
and specify the curved action

```
S=integral sqrt(-g) [F(X) R + L_SM,min],
F(X)=F0-xi X,                 F0>0.
```

The metric, spin structure, curvature coupling xi and covariant extension
are specified inputs. L_SM,min contains the canonical gauge, Higgs,
fermion and Yukawa terms from the conditional field-action construction.
The package verifies the scalar/current/Yukawa tensor attachment and
classical variational identities. It constructs no native curved field
evolution, quantum Ward identity or physical entropy state.

## Metric variation and the common area coefficient

Holding matter fields fixed, integration by parts in the variation of the
Ricci scalar gives

```
delta(sqrt(-g) F R)
 = sqrt(-g) [F G_ab+(g_ab box-nabla_a nabla_b)F] delta g^ab
   + boundary divergence.
```

Consequently the metric equation of this action is

```
2F G_ab+2(g_ab box-nabla_a nabla_b)F=T_min,ab.             (1)
```

Equivalently 2F0 G=T_min+T_xi, with

```
T_xi,ab=2xi[X G_ab+(g_ab box-nabla_a nabla_b)X].          (2)
```

A constant vacuum term may be included in T_min; it is not set by this
calculation. An independent lapse/scale-factor Euler variation of the
homogeneous curved action checks both metric components of(1). The
curvature convention is R_ab=R^c_acb, with positive scalar curvature
for a round spatial sphere.

The curvature derivative is
E^abcd=(F/2)(g^ac g^bd-g^ad g^bc). A normalized Lorentzian binormal
epsilon has epsilon_ab epsilon^ab=-2. Its contraction gives the Wald
functional

```
S_W=-2pi integral_boundary E^abcd epsilon_ab epsilon_cd dA
   =4pi integral_boundary F dA.                          (3)
```

Thus eta(X)=4pi F=1/(4G0)-4pi xi X if F0=1/(16pi G0).
The coefficient is not constant for a varying Higgs norm. Positive F is
a required domain condition; its value and xi are not selected by the
finite gauge representation. Calling(3) a Wald functional does not
identify it with the entropy of a simulator state or prove an entropy
stationarity principle. The Noether-charge construction is due to
[Wald](https://arxiv.org/abs/gr-qc/9307038).

## Conservation with all source forces retained

The contracted Bianchi identity and the scalar commutator give

```
nabla^a T_xi,ab=-xi R nabla_b X.                         (4)
```

The improvement term is not conserved separately. To exhibit the other
terms, write the Yukawa interaction as V_Y=H†J+J†H, where J is the
Higgs-covariant even fermion bilinear constructed from the same Yukawa
tensors. Let D=partial+A with anti-Hermitian connection, F_ab=[D_a,D_b],
and define

```
E_H=D^2 H-V'(X)H-J-xi R H,
j_A^a=2 Re[(D^a H)† t_A H].
```

For the Higgs stress including V_Y and(2), direct differentiation yields

```
nabla^a T_H,ab
 =2 Re[E_H† D_b H]
  +2 Re[(D^a H)† F_ab H]-2 Re[H† D_b J].                (5)
```

The second term is the gauge force. Adding Yang-Mills stress replaces it
by F_ba^A E_A^a, where E_A^a=(D_c F^ca)_A/g_A^2-j_A^a,
with any other charged matter current retained. The last term is the
Yukawa exchange with the fermion sector. It cannot be discarded when J
varies. `ward.py` checks(4) and(5) directly in a curved rational metric,
with noncommuting weak connections and nonzero J, without imposing the
field equations. The even bilinear coefficients are formal commuting
coefficients for this scalar identity; the execution is not a fermionic
quantum state or a full Grassmann stress replay.

For the complete covariant classical action, varying all fields under a
compactly supported diffeomorphism gives its off-shell Noether identity.
On all matter equations, only the metric term remains:
integral sqrt(-g) T_ab nabla^a zeta^b=0. Integration by parts for arbitrary
zeta proves nabla^a T_ab=0. This derives total on-shell conservation from
the action's covariance. It does not assume that a scalar, fermion or
gauge sector is conserved by itself, and it is not a regulated quantum
Ward theorem. In the equivalent form(1),
nabla^a[2F G_ab+2(g_ab box-nabla_a nabla_b)F]=-R partial_b F;
the Higgs equation supplies the matching exchange.

## A finite first-variation identity on a ball

Let s be the variation parameter. The reference at s=0 is flat spacetime,
with constant Higgs norm and positive F*. On a spatial slice orthogonal
to a future unit timelike u, use radial normal coordinates x and a ball
B_ell. The first spatial metric variation is h_ij=partial_s g_ij|0,
with h(0)=partial h(0)=0 and h_ij(x)x^j=0. The latter is radial gauge;
coordinate radius is proper radial distance to first order. Write
f(x)=partial_s F(X_s(x))|0. For xi and nonzero reference w, any real
polynomial f is realized by the actual doublet variation
H_s=(0,w-s f/(2xi w)): partial_s(F0-xi H_s†H_s)|0=f.
At xi=0 only f=0 is allowed by this source coupling.

All following variations are first derivatives at s=0, not finite
differences at a fixed nonzero s. Products such as f h are second order
and do not enter. Exact first-variation formulas are

```
delta V=(1/2) integral_B tr h,
delta A=(1/2) integral_boundary tr h,
delta S_W=4pi[F* delta A+integral_boundary f].            (6)
```

The Noether generalized volume specialized to F R in four spacetime
dimensions is

```
W=integral_B (3F-F*)/(2F*) dV,
delta W=delta V+(3/(2F*)) integral_B f.                  (7)
```

This follows by contracting E^abcd with the slice normal and projector:
E^abcd u_a h_bc u_d=3F/2, while its reference coefficient is E0=F*/2.
The functional is the specialization of the generalized volume in
[Bueno, Min, Speranza and Visser](https://arxiv.org/abs/1612.04374).
Its physical use as an admissibility constraint is an additional premise.
It is distinct from ordinary fixed proper volume in the existing Einstein
interface.

Define the variation with the generalized-volume term removed by

```
B_ell=delta S_W-(8pi F*/ell) delta W.
```

For a scalar Taylor polynomial through degree three, spherical symmetry gives

```
integral_boundary f=4pi ell^2 f(0)+(2pi/3)ell^4 Delta_3 f(0),
integral_B f=(4pi/3)ell^3 f(0)+(2pi/15)ell^5 Delta_3 f(0).
```

Odd terms integrate to zero. In particular(7) cancels the ell^2 f(0)
term exactly. Among constraints delta V+alpha integral_B f, cancellation
for arbitrary constant f forces alpha=3/(2F*). The remaining scalar
coefficient is 16pi^2 ell^4 Delta_3 f/15.

The normal-coordinate curvature identity is
Delta_3 tr h(0)=-(2/3)delta R^(3)=-(4/3)delta G_uu.
Using the same moments in(6) gives

```
B_ell= -2pi C_ell [2F* delta G_uu-2 Delta_3 f(0)] + R_ell,
C_ell=4pi ell^4/15.                                    (8)
```

The bracket is exactly the timelike variation of the left side of(1):
(g_uu box-nabla_u nabla_u)f=-Delta_3 f at the flat reference.
The metric and scalar coefficients in(8) are therefore derived from the
same action as its Wald functional, rather than identifying a convenient
matrix with a desired stress tensor.

Assume the absolute fourth directional derivatives of tr h and f on
B_ell are bounded by M_h and M_f respectively, with Euclidean unit
directions and M_h,M_f>=0. Taylor's integral remainder is at most
M|x|^4/24. The sphere and ball integrals of this bound give

```
|R_ell| <= pi^2 ell^6 [(3/7)F* M_h+(20/21)M_f].          (9)
```

This is an explicit bound for the first variation. It bounds neither
finite-s quadratic corrections nor an ultraviolet quantum limit.
The executable polynomial moment calculator retains all degrees; for
degree at most three(8) is exact, and quartic/higher controls check the
remainder independently.

At ordinary fixed proper volume, the expression is instead

```
delta S_W|V
 =16pi^2 ell^2 f(0)
  -2pi C_ell[2F* delta G_uu-5 Delta_3 f(0)] + O(ell^6).  (10)
```

Fixing f(0)=0 alone does not remove the coefficient mismatch. A constant
Higgs norm throughout the allowed variations, xi=0, or quantitative
controls that make both unwanted contributions negligible are sufficient
specializations; a variable eta cannot simply be replaced by a constant.

## Conditional entropy-to-equation join

If the same regional state supplies a minimally coupled bulk first-law
variation delta S_bulk=2pi C_ell delta T_min,uu+e_bulk, and generalized
entropy is stationary at fixed W with error e_stat, then(8) yields

```
|delta[2F G_uu+2(g_uu box-nabla_u nabla_u)F-T_min,uu]|
 <= (|e_bulk|+|e_stat|+|R_ell|)/(2pi C_ell).              (11)
```

Probe completeness can use the existing ten-observer tensor theorem.
Equation(11) is a conditional linearized scalar-tensor equation. It
does not prove the bulk entropy identity, fixed-W stationarity, a native
spatial/refinement map or nonlinear field equations. Nonminimal boundary
and quantum entropy counterterms require a consistent split so that the
same coupling is not counted in both S_W and an improved bulk stress.
With f=0 and constant positive F*, it reduces to the constant-coefficient
linearized Einstein normalization eta=4pi F*. The finite responsive Gram
and regional parity entropy are not promoted to any of these hypotheses.
