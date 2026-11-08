# Composite response on the original-default cube

The complete eight-site source at `c=1/8,U=1` has a nonzero composite
response residual in every original mode direction. In its unique full
centered ground, the exact residual satisfies

    (3/80) I_16 < R*R < (61/1000) I_16.                    (1)

The inequalities are certified over rational arithmetic from the actual
65536-state clock and its certified 4900-coordinate ground proposal. They
do not use a fitted energy, an aligned-spin reference, or a presumed
dispersion. The original16-mode response has three independent cyclic
levels, of total dimension48. These are response vectors on the
same source, not additional canonical fermions or identified particles.

## Source, ground and common clock

Use the open2x2x2 cube in `cube_eta_protection`, with interleaved hidden
CAR modes `(x,up),(x,down)`. Its exact integer action is `16H`, and

    K=H-N/2,    H=T+sum_x n_xup n_xdown,
    T=sum_<x,y>,s h_xy(a_xs* a_ys+a_ys* a_xs),
    h_xy=eta_xy/16,    eta_xy in {-1,1}.                   (2)

Every edge sign is read from the actual full source. The onsite frame
relates this signed Hubbard chart to the original Weyl/contact action.
The response statements transport by that unitary. A transported
experiment must transport its instrument too; the separate local-observer
theorem also treats the original Weyl action with its bare repair maps.

`centered_source_ground/gap.py` proves that the absolute ground ofK,
over all65536 states, is unique, lies in the complete balancedN8 block,
and has full excitation gap greater than1/128. Its integer rank-one
Gram certificate supplies the balanced excited-energy lower bound
`b=-1053/256`. The source proof also establishes hidden spin zero and
particle-hole invariance. The staggered eta generators commute withK.
Since eta raising or lowering changesN, uniqueness forces both to
annihilate this ground; it is eta invariant as well.

A fresh proposal `p` uses the same deterministic eigensolver initialization
and23-bit dyadic rounding as the gap owner. Its admission uses only the
actual integer source and exact values

    r=<p,Kp>/<p,p>,   v=||(K-r)p||^2/<p,p>,
    epsilon=v/(b-r)^2,
    r-v/(b-r) <= E0 <= r.                                 (3)

The inherited excited bound gives ground infidelity at mostepsilon.
Consequently the half trace distance between the normalized trial and
the exact full ground is at mostsqrt(epsilon). Integer arithmetic checks

    epsilon < (1/1920)^2.                                 (4)

No symmetry projection or normalization of selected measurement outcomes
is applied to the trial. The receipt retains the fixed rational bound(4),
rather than platform-dependent eigensolver residuals. For every bounded
operatorO, its trial and ground expectations differ by at most
`2delta||O||`, with `delta=1/1920`.

## Actual composite columns

Let `b_i=2(n_iopposite-1/2)a_i`. The complete source phase
`D=product_x exp(i pi n_xup n_xdown)` obeys `b_i=-D a_i D*`.
Each family is canonical; the mixed a/b anticommutators are nonzero
operators. The invariant-ground retarded embeddingsV0,V1, with the
conjugate reversed removal clock, are the orthogonal isometries of
`COMPOSITE_CAR.md`. Putalpha=1/2 and

    Lret V0=V0 h+alpha V1,
    Y=Lret V1-alpha V0,
    C=V1*Y,
    R=(I-V0V0*-V1V1*)Y,
    R*R=Y*Y-C^2.                                         (5)

The i-th column ofY is the retarded vector of `[T,b_i*]`. Therefore all
entries required in(5) are expectations of actual local even operators,

    C_ij=omega({b_i,[T,b_j*]}),
    (Y*Y)_ij=omega({[T,b_i*]*,[T,b_j*]}).                 (6)

The source code constructs the sparse integer matrices
`F_i=[16T,b_i*]` and checks on the full occupation space that
`[16K,b_i*]=F_i+8a_i*`. It does not substitute a projected spin clock.
EveryFi has maximum absolute row and column sum9, so

    ||[T,b_i*]|| <= 9/16.                                 (7)

The matrix inequality uses
`||Fi||^2<=||Fi||_1||Fi||_infinity`; sparse entries and both sums are
checked exactly. All16 original modes and both number-changing response
directions enter the calculation.

## Two independent readouts ofC

The spin/eta identity in `COMPOSITE_CAR.md` gives, on each actual edge,

    C_(x,s),(y,t)=delta_st h_xy omega(O_xy),
    O_xy=4(S_x.S_y+eta_x.eta_y),    ||O_xy||<=3,            (8)

and zero for onsite and nonedge entries. The two invariant local sectors
are disjoint, so the norm bound follows directly from the two-site spin
and eta singlet/triplet matrices. The source implementation constructs
Oxy in the full16-mode space, with the actual bipartite eta signs. It
compares its trial readout with the independent commutator readout(6),
including their exact ground-error allowances.

For every one of the12 edges, integer expectation arithmetic proves

    -157/100 < omega(O_xy) < -155/100.                     (9)

No equality between distinct edge expectations is needed. The signs and
normalizations ofC follow from the originalh. In particular this is a
source correlation coefficient, not a supplied Dirac lower block.

LetChat be assembled from the trial expectations in(8), with exact
structural zero entries. Since each cube vertex has degree3,

    ||C-Chat|| <= e_C=3(1/16)6delta=3/5120,
    ||C|| <= 3 c degree/2=9/16.                           (10)

LetGhat be the trial expectation matrix ofY*Y in(6). Each corresponding
operator has norm at most `2(9/16)^2`. Entrywise expectation bounds and
the symmetric row-sum norm give

    ||Y*Y-Ghat|| <= e_G=16*4delta*(9/16)^2=27/2560.        (11)

These estimates apply to arbitrary complex coefficient vectors in the
16-dimensional retarded mode space. No entrywise sign or trial-state
symmetry is assumed in the conversion to an operator bound.

## Positive residual and controlled response

Form `Rhat=Ghat-Chat^2` with rational matrix multiplication. Equations
(10)--(11) imply

    ||R*R-Rhat|| <= e_R=e_G+(9/8)e_C+e_C^2.               (12)

Exact Gershgorin row bounds verify

    min_i(Rhat_ii-sum_(j!=i)|Rhat_ij|)-e_R > 3/80,
    max_i(Rhat_ii+sum_(j!=i)|Rhat_ij|)+e_R < 61/1000.      (13)

This proves(1) without an eigensolver admission or a scalarity assumption
for the true residual. Deleting the residual or its removal contribution
fails the admitted complete matrix bounds. The exact dimer closure cannot
be imposed on the cube.

Since alpha is nonzero,V1 is generated byV0 andLretV0. The injective
16-columnR is orthogonal to both isometricV0,V1 ranges and is generated
by the next application ofLret. Thus

    dim span{V0,Lret V0,Lret^2 V0}=48.                    (14)

The exact fourth original moment differs from the two-block compression
`A=[[h,alpha I],[alpha I,C]]` by `alpha^2 R*R`. In particular,

    (3/320)I < M4-[I,0]A^4[I;0] < (61/4000)I.            (15)

The Feshbach identity and source-vector estimate in `COMPOSITE_CAR.md`
therefore have the concrete cube bounds, for eta=|Imz|>0,

    ||W*(z-Lret)^-1W-(z-A)^-1|| <= (61/1000)/eta^3,
    ||exp(-itLret)V0-W exp(-itA)[I;0]||^2
      <= (61/1000) min(|t|,t^2/4)^2,   W=[V0,V1].        (16)

These control retarded vectors, response correlations and smoothed
spectral scores. The signed retarded embedding does not by itself
implement a physical channel on arbitrary quantum inputs and references.
No native preparation of the full centered ground is executed here.
For each local original site/flavour probea_i, the addition and removal
weights are each1/2 and are not renormalized by(1)--(16). A general
normalized smearing has total signed weight1; its separate charge weights
need not be equal.

The cube has open boundaries. Spatially oscillating coefficients are
valid probe profiles; graph reflections or such profiles do not acquire
continuous momentum labels. Positive finite residuals establish neither
the existence nor the absence of a continuum pole. The source calculation
supplies no relativistic mass, extra CAR flavour or Standard Model label.
