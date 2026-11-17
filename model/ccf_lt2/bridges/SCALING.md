# A common uniform-chart principal speed with explicit cutoff bounds

Let `a,g,c>0` and `r>=0`. On a uniform cubic chart set

    g_E=g^2/a,      g_B=c^2/(g^2 a).

These are arguments of the existing source builder, not new field or matter
terms. Positivity preserves its finite Hermiticity and Gauss proofs. A
constant common lapse multiplies both source principal speeds by the same
factor; a dynamical/inhomogeneous geometry requires a separate analysis.

## Electromagnetic quadratic symbol

Write `d_j=exp(i a k_j)-1`. The oriented face derivative is `C(d)A=d cross A`.
Direct multiplication gives

    C(d)† C(d) = |d|^2 I - d d†.

It vanishes on the gradient direction `d`, and on `d†A=0` its eigenvalue is
`|d|^2=4 sum_j sin^2(a k_j/2)`. At `d=0` all modes have zero symbol.
The quadratic compact-field Hamiltonian is
`(g_E/2)|E|^2+(g_B/2) A†C†C A`, so canonical Hamilton equations yield

    omega_EM^2 = 4 g_E g_B sum_j sin^2(a k_j/2)
               = (4c^2/a^2) sum_j sin^2(a k_j/2).

This is an exact statement about the bulk quadratic classical symbol of the
actual compact plaquette potential. It is not an equality for the full
finite hard-wall quantum spectrum. The compact-to-quadratic approximation,
cutoff removal, quantum state and field readout remain separate requirements.

## Dirac and actual mass compatibility

The actual forward kinetic coefficient is `i c alpha_j/(2a)` and its adjoint
is backward. Fourier transformation therefore gives
`p_a,j=c sin(a k_j)/a`. The local and neighboring Wilson contributions give

    w_a(k)=(r/a) sum_j(1-cos(a k_j)),
    h_a(k)=M+w_a(k) beta+sum_j alpha_j p_a,j.

The parent receipt establishes the actual charged and neutral three-family
matrices `M` are Hermitian and anticommute with every actual `alpha_j`.
`beta` also anticommutes with them. Consequently

    h_a(k)^2 = (M+w_a(k) beta)^2 + |p_a(k)|^2 I.

This retains the actual matrix mass, both signs and Wilson term. It does not
replace the full source's interactions/field/geometry by this free symbol.
Even the squared mass matrix need not be a scalar across flavours.

Both principal slopes tend to `c`. Equivalently the required coefficient
condition is `a^2 g_E g_B=c^2`. With defaults `g_E=g_B=c=1`, the EM slope is
`a` and the Dirac slope is one; equality fails for every refined `a<1`.
The explicit family repairs this scale mismatch. Its value of `c` and
coupling `g` remain supplied parameters, not independent predictions.

## Uniform estimates

For every real `x`,

    |sin x-x| <= |x|^3/6,
    0 <= x^2-sin^2 x <= x^4/3,
    0 <= 1-cos x <= x^2/2.

The first and third follow by integrating the elementary cosine/sine bounds;
the second follows from `|x+sin x|<=2|x|` and the first. Thus, for Euclidean
`|k|<=K`, and using `sum k_j^4<=K^4`,

    0 <= c^2|k|^2-omega_EM^2 <= c^2 a^2 K^4/12,
    0 <= c^2|k|^2-|p_a|^2 <= c^2 a^2 K^4/3,
    |p_a-c k| <= c a^2 K^3/6,
    0 <= w_a <= r a K^2/2.

The vector bound uses `sum k_j^6<=K^6`; Clifford algebra gives the same norm
bound for `alpha.(p_a-c k)`. Hence the free Hamiltonian error obeys

    ||h_a-(M+c alpha.k)|| <= c a^2 K^3/6 + r a K^2/2.

For any explicit `M0>=||M||`,

    ||h_a^2-(M^2+c^2|k|^2 I)||
      <= M0 r a K^2 + a^2 K^4(r^2/4+c^2/3).

Indeed the difference is `w{M,beta}+w^2 I` minus the scalar kinetic deficit,
and `||{M,beta}||<=2M0`. This retains the Wilson mass cross term of order
`a`; it does not mislabel every mass error as order `a^2`. The bound holds
for each actual mass matrix using its own norm, without recalibrating it.
It is uniform on bounded momenta and does not cover the full Brillouin zone.

This is an analytic design only. No dedicated implementation or executed
scaling-family receipt is asserted.

## What remains

The supplied family closes a normalization prerequisite on a uniform chart.
It does not prove a joint physical field/matter/geometry refinement, identify
the interacting particle poles, establish a physical clock or rods, remove
Wilson species effects in interactions, or derive the coupling/velocity.
Open-box boundary effects are not diagonalized by the bulk Fourier symbol.
An irregular angular/radial mesh needs appropriate Hodge/volume weights or
a proved comparison; the uniform formula does not supply them. The finite
local-capacity locality theorem and this bulk speed match are separate
ingredients, not a completed Lorentz theorem.

Separate fixed-graph background and weak-coupling investigations leave the
specified default coefficients unchanged and do not establish a physical
common speed. Their independent experiments are outside this receipt; they
do not supply an executed spatial-refinement family.
