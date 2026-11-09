# A nonzero finite response and its refinement boundary

## The actual source and the question being tested

Use exactly the [primary contract](../particle_common_contract/PROOF.md):
spacing `a=2^-n`, side `2^(2n+1)`, volume V, per-link word degree n+1,
all sixty CAR modes per cell, full native M15 links, and every Gauss
constraint. The selected polynomial density rho has ground excess at
most `epsilon=2^(-n-4)`, commutes with the actual H, and obeys

```
rho >= c P_G,       c=2^-N>0.
```

N is the exact factored floor exponent supplied by the
[physical boundary-state proof](../particle_physical_boundary_state/PROOF.md).
This binary density-floor exponent is distinct from a cosmological capacity
parameter written `N=ln D`; the two quantities are not identified.
The reference is its canonical conjugate with idle Hamiltonian zero.
No scalar energy shift is made in H or its controlled Cayley gates.

The [common field algebra](../particle_primary_field_algebra/PROOF.md)
gives a stationary positive-energy subsequential representation
of fixed bounded filtered probes. Its unit survives, but its nonconstant
probes may vanish in the selected representation. The result below proves
a nonzero finite primary response with explicit level-dependent times.
It does not assert a fixed-time or fixed-energy-band survival theorem.

## An exact plaquette commutator on a physical witness

For n>=1 the centered plaquette at target position zero is active in
the common algebra's authenticated menu. Its four links are distinct.
Let

```
A=P_ell ReTr_M15(U_p) P_ell/15,       ||A||<=1,
|b>=empty Fock tensor constant link functions.
```

The blank b is normalized and physical. It is used only to evaluate a
matrix element; it is not substituted for rho. The ten native and
conjugate irreducible character loops `L_(r,+),L_(r,-)` in the boundary
proof are normalized, mutually orthogonal, gauge invariant and
inside the product word cutoff. Their nonzero integer charges distinguish
the conjugate blocks. On this blank,

```
A|b>=(1/30) sum_(r,sign) |L_(r,sign)>.
```

The electric term has coefficient `e0/(2a)`, with e0=105. Each loop has
four links with the same trace-metric Casimir C_r. The actual generator
metric and source Casimir packet independently give

| Q6 | Dimension | C_r |
| --- | --- | --- |
| -4 | 3 | 4/5 |
| -3 | 2 | 9/20 |
| 1 | 6 | 21/20 |
| 2 | 3 | 7/10 |
| 6 | 1 | 3/10 |

Thus `sum_r C_r²=101/40`, and `H_E L_(r,sign)=2e0 C_r L_(r,sign)/a`.

Every term of the full matter Hamiltonian is number preserving and
annihilates the empty Fock factor, both before and after A. Its
commutator therefore annihilates b. For the magnetic term write its
parent multiplication operator as V and the parent plaquette
multiplication as A_parent. They commute before compression. Both
`Vb` and `A_parent b` lie in the product cutoff: every single plaquette
uses each of its four links only once. Hence

```
[P_ell V P_ell,P_ell A_parent P_ell]b
 =P_ell(V A_parent-A_parent V)b=0.
```

The common final compression is essential. Products of compressed
multiplication operators need not commute on other inputs. The proof
does not assume that the eleven-dimensional blank/loop span is an
invariant sector of the full Hamiltonian.

Since `H_E b=0`, these identities give an exact matrix element on the
actual full physical source:

```
[H,A]b=(7/a) sum_(r,sign) C_r L_(r,sign),
||[H,A]b||²=K/a²,              K=4949/20.                (1)
```

This is independent of box size and cutoff once the plaquette is
present. It includes the actual full-Fock and magnetic action through
their proved vanishing commutators on b, rather than through removal
of those terms from H.

## Original moments and a proved primary contrast

For the stationary original spectral measure

```
mu=sum_(i,j) p_i |A_ji|² delta_(E_j-E_i),
F(t)=Tr rho A alpha_t(A)=integral exp(it omega)dmu,
```

the second moment is `M2=Tr rho [H,A]*[H,A]`. The physical density
floor and (1) imply

```
M2 >= c ||[H,A]b||² = c K/a².                          (2)
```

The trace inequality uses `Tr_phys C*C>=||Cb||²` for C=[H,A]. It does
not use a population assigned to b in an energy eigenbasis.

Let `[L,U]` be the selected source's certified interval and W=U-L.
Every transfer frequency lies in [-W,W]. Choose the following distinct
experiment at each level:

```
T_n=t_n=1/W_n,
B_n=integral f_T(s) alpha_s(A)ds,
f_T(s)=(1-|s|/T)_+/T.
```

The inherited triangle theorem gives
`dnu=sinc(omega T/2)^4 dmu`. For `|omega|<=W`, the sine remainder gives
`sinc(omega T/2)>=23/24`. Set `r=(23/24)^4`. The cosine remainder on
`|omega t|<=1` gives `1-cos(omega t)>=(11/24)(omega t)²`. Consequently

```
M2_filtered >= r c K/a²,
D_n=F_B(0)-Re F_B(t_n)
    >= d_n=(11/24) r K c/(a² W²)>0.                    (3)
```

This is a source-derived lower contrast, not a supplied observed value
or a finite comparison-source transfer. The complete spatial response
matrix is not executed. Equations (1)-(3) certify its sign and a
potentially extremely small original lower weight.

## A positive finite band using polynomial passivity

On the physical space the selected density is a nonincreasing function
of energy: `p(E)` is proportional to
`[T_m(beta(U-E)/(2s))]^(2s)`. All polynomial coefficients are nonnegative,
`U-E>=0`, and beta,s>0. Therefore `p(E_i)>=p(E_j)` whenever E_i<E_j.
For a Hermitian A, `|A_ji|²=|A_ij|²`, so its original measure obeys

```
mu(S)>=mu(-S)       for every Borel S subset (0,infinity).
```

The triangle multiplier is even, so the same dominance holds for nu.
This passivity uses the actual polynomial density, without replacing
it by a Gibbs state or an exact ground vector.

For any `0<delta<=W`, total mass at most one and support [-W,W] give

```
M2_filtered <= delta² + W² nu({|omega|>=delta}),
nu([delta,W]) >= (M2_filtered-delta²)/(2W²).            (4)
```

Choose `delta_n=(8/a)2^(-ceil(N/2))`. Then
`delta_n²<=64c/a²`, and the source bounds give delta_n<W. Since rK>64,

```
nu([delta_n,W_n]) >= (rK-64)c/(2a²W²)>0.              (5)
```

Negative-weight subtraction from the general
[positive-window criterion](../particle_positive_window/PROOF.md)
is too weak for this tiny witness: its cooling allowance can exceed
the contrast. Equation (4) instead uses the additional proved passivity
and the exact finite spectral support. For fixed physical times and
fixed bands, `fixed_band_requirement` retains the existing criterion
and its actual supplied-record errors. Equation (5) does not discharge
that fixed-band premise.

## A complete finite apparatus that resolves the sign

Use the existing [filtered-probe apparatus](../particle_filtered_probe/PROOF.md)
for the authenticated plaquette contraction. Every selector, neutral
probe and Ramsey outcome and the original quantum reference is retained.
There is no division by success. The source interval has L<=0<=U, hence
`R_H=max(|L|,|U|)<=W`. With T=t=1/W, even quadrature count M and R Cayley
subdivisions, its complete original contrast comparison gives

```
|D_finite-D_n| <= 4/M + 5/(6R²),                      (6)
```

using exact declared algebraic preparation/apparatus services, with
zero additional compiler or preparation allowance. Approximate services
and empirical sampling must add their own errors. The square-root
Lüders instrument is a separate valid apparatus and is not substituted
for this Ramsey readout.

Write `d_n=C_n 2^-N`. Choose

```
M=2^m,  m=max(1,N+ceil(log2(16/C_n))),
R=2^r,  r=max(0,ceil((N+ceil(log2((10/3)/C_n)))/2)).
```

Each term of (6) is at most d_n/4. Thus the specified finite complete
apparatus has `D_finite>=d_n/2>0` on the actual selected state. The code
checks these exact inequalities using factored dyadics, so it never
allocates an integer with N bits or the enormous selector. This is a
finite computable construction under the inherited explicit services;
the full spatial gate tape and its sample counts have not been executed.

The finite filter also preserves a positive original band weight. In the
energy basis, the Hilbert-Schmidt vector `B sqrt(rho)` has coordinates
`B_ji sqrt(p_i)`. Project its `(j,i)` coordinate onto a transfer window S
exactly when `E_j-E_i` belongs to S. The squared projected norm is the
original measure nu_B(S). Both B and its finite contraction B_tilde
have weighted Hilbert-Schmidt norm at most one. Therefore, for every
transfer window,

```
|nu_Btilde(S)-nu_B(S)| <= 2||B_tilde-B|| <= 2e_B.
```

This is a frequency-difference projection, not a final-energy projector
on the source or an evolution of the idle reference. The precision
schedule above gives `2e_B<=d_n/4`. Thus the same band in (5) satisfies

```
nu_Btilde([delta_n,W_n])
 >= [(37/96)rK-32] c/(a²W²)>0.                        (6a)
```

This is the original spectral weight of the actually specified finite
filter operator. Its complete Ramsey readout retains the allowance in
(6). Neither this transfer interval nor the source spectral width is
identified with a rest-mass interval.

## Exact invariant sectors with vanishing selected weight

The total fermion number N_F commutes with the actual Hamiltonian,
product word projection and Gauss actions. Its empty and fully occupied
sectors are exact reducing sectors. On the empty Fock state every
`c_i* c_j` vanishes. On the filled state it vanishes for i!=j, and its
diagonal sum is the one-body trace. That trace is zero: beta_D has two
positive and two negative rows and all hopping blocks are spatially
off-diagonal. Thus both extreme sectors have exactly

```
H_extreme=H_E+H_B >= -P/(7a),                          (7)
```

where `P=3 side (side-1)²` is the actual number of interior plaquettes.
They retain the full native product group and every admitted link word.
They are not compact-simple factor sectors, and magnetic traces
couple the group factors through M15.

These Fock vectors are gauge singlets. The native representation is
`Lambda²(C5) direct-sum Lambda⁴(C5)` for a determinant-one fundamental
matrix; the exterior determinant identity gives `det R_M15=1`.
The empty vector is invariant and a filled native copy transforms by
this trivial determinant.

For a global ground upper bound, at every site fill both negative-beta
Dirac directions and all fifteen native internal states, and use constant
links. This is a physical trial vector by the same determinant identity.
Its onsite energy is `-90/a` per cell. Each hopping term changes the
fermion number of its two cells and therefore has zero expectation.
Its electric expectation is zero and every nontrivial plaquette character
has zero Haar expectation. It follows that

```
E_phys <= -90V/a,
Delta_extreme=(90V-P/7)/a>0.                           (8)
```

Let P_ext be the sum of the empty and fully filled sector projectors.
They are orthogonal and reducing. On the physical space,
`H-E_phys>=Delta_extreme P_ext`: use (7)-(8) on those sectors and
ground positivity on their complement. The selected mean-excess bound
therefore gives a bound on their COMBINED original weight,

```
Tr rho P_ext <= epsilon a/(90V-P/7)
             <= 7/(627*2^(8n+7)).                     (9)
```

For any contraction supported in these sectors, its total original Gram
spectral mass is at most (9). Conditioning on P_ext defines a different
state and removes the probability cost. An internal threshold of that
conditional pure-gauge experiment cannot be assigned surviving weight
in the primary preparation. A source operation that injects an extreme
sector instead encounters the extensive ground-relative separation in
(8), rather than an established finite rest mass. No extreme-sector
spectrum or controlled strong-coupling gap is computed here.

## The remaining primary noncollapse premise

The witness interval (5) widens towards both zero and infinite energy:
N is at least the full dimension exponent, delta_n->0 and W_n->infinity.
Its proved weight floor tends to zero. Both t_n and T_n tend to zero,
and their time-regularity cap 2/T_n diverges. Thus (3)-(6) do not enter
the prior fixed-filter, common strongly continuous algebra theorem as
a uniform nonzero witness. A vanishing lower bound does not prove that
the actual fixed-band response vanishes; it leaves that question open.

Shrinking a filter can retain weight whose energy grows as 1/a. For a
logical counterfamily take a two-level Hamiltonian with spacing 16/a,
the Hermitian exchange probe, and an excited population small enough
to obey the same cooling tolerance and a positive density floor. It is
passive and its commutator second moment grows as 1/a². At fixed T, the
filtered mass is bounded by `16/(T^4(16/a)^4)` and vanishes. At T=a/16
the multiplier is constant, but all nonzero energy escapes every
fixed compact band. These matrices are a counterexample to inference
from those inequalities, not a substitute for the primary Hamiltonian.

Normalizing rare records changes the measure and cannot restore their
original weight. Rescaling an operator by the reciprocal square root
of a vanishing probability also exits the uniform bounded-probe domain
and requires separate domain, moment and cooling estimates. The remaining
positive premise is a same-source, fixed-band original response with
uniform positive weight and the existing complete apparatus errors.
Rest-mass identification additionally needs momentum and a particle
spectral attachment. Source-action units alone do not select an absolute
SI scale; a declared comparison scale is allowed and retains its ancestry.
