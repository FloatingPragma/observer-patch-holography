# A local completion of the protected-flux Dirac route

This construction specifies a finite local program carrying the same
paired-flux link object through a four-component routing field. It uses
the twelve-word link domain of [the mobile source](PROOF.md), together with
declared routing memory and gates. It is not an invariant-subspace theorem
for the original occupation-controlled clock. Its one-object relativistic
window follows from the [Clifford routing proof](../spatial/ROUTE.md).

## The physical object and its code

Let `s` be the interior vertex-star word. The twelve admitted link words
form six two-element gauge orbits: one flat orbit `{0,s}` and five paired-
flux orbits `{g_f,g_f xor s}`, with `f=1,...,5`. With empty native matter
registers, the interior Gauss operator just exchanges the two words. Set

\[
 |[g]\rangle_+=(|g\rangle+|g\mathbin{\mathrm{xor}}s\rangle)/\sqrt2.
\]

These six vectors are orthonormal and Gauss positive. The flat orbit has
zero flipped faces. Every other orbit has one flipped face in the region
and its exterior compensator, hence two flipped faces on the full carrier.
Their regional charges are zero and one respectively. A paired-flux object
is globally neutral under the spherical parity identity; it is not an
isolated global monopole.

Introduce a neutral routing memory with an empty label and two occupied
labels. One logical lane has the eleven-dimensional code

\[
 \mathcal L=\mathbb C|0\rangle\oplus(\mathbb C^5\otimes\mathbb C^2),
 \quad |0\rangle\mapsto|[0]\rangle_+|\mathrm{empty}\rangle,
 \quad |f,\sigma\rangle\mapsto|[g_f]\rangle_+|\sigma\rangle.
\]

The routing memory is an architectural input. In particular, four spinor
components are not inferred from the number of faces or from a scalar
holonomy alone. The face syndromes, their gauge partners and the exterior
compensator are retained in the code, rather than replaced by labels with
no source encoding.

## Why two lanes are required in this completion

At each spatial cell take two independent lanes,

\[
 \mathcal H_x=\mathcal L_{x,+}\otimes\mathcal L_{x,-},\qquad
 \dim\mathcal H_x=121.
\]

The zero/one-object subspace has dimension `1+10+10=21`; its one-object
factor is `C5 tensor C4`. The remaining one hundred states have one object
in each lane. They are retained. Truncating the cell to dimension twenty-one
would discard states reached when two counterpropagating objects arrive at
the same cell. A one-object walk alone does not justify that truncation.

The two lane link registers are independent copies of the stated source
code. This is an explicit enlargement of local memory, rather than an
assertion that the original single link register already held two objects.
One can use a different many-object architecture, but it must also give a
unitary rule at collisions.

## Exact conditional translations by local swaps

Fix a periodic spatial direction `e_i`. First swap the two lanes at every
cell. Then, for every x, swap the lane factors `L_(x+e_i),+` and `L_x,-`.
Every factor appears in exactly one of the second-stage swaps. Therefore
those gates commute even when their bonds share a spatial cell. Both stages
are products of finite local unitaries; they preserve the exact object
number on the entire tensor-product Hilbert space.

Writing the initial lane contents as `(p_x,m_x)`, the two stages give

\[
 p'_x=p_{x-e_i},\qquad m'_x=m_{x+e_i}.
\]

Thus the plus lane moves right and the minus lane moves left. This identity
holds for empty and occupied lanes, including configurations with two
objects per cell. In the one-object wavefunction it reads

\[
 \psi'_x=P_+\psi_{x-e_i}+P_-\psi_{x+e_i},
 \qquad P_\pm=(I\pm Z)/2,
\]

where Z has two plus and two minus entries on the four routing components.

For each Clifford matrix `alpha_i`, choose a fixed unitary V_i with
`V_i Z V_i*=alpha_i`. Let any twenty-dimensional unitary W on the one-object
space act locally by

\[
 \widehat W=P_0\oplus W\oplus I_{100}
\]

on the zero-, one- and two-object sectors. This is a complete local unitary;
its action when two objects coincide is specified. Applying
`widehat(V_i)*`, the countertranslation, and `widehat(V_i)` realizes exactly
the directional shift with projectors `(I+/-alpha_i)/2` on one object.
No one-object amplitude is postselected or dropped.

The identity choice on the two-object sector is a declared collision rule.
It does not assert that the many-object evolution is an exterior lift of
the one-object operator. In particular, it does not establish fermionic
statistics. Tensor-factor swaps carry no CAR exchange sign by themselves.

`lane_collision.py` executes this completion on every one of the 121 local
basis states for an exact mass coin, including all one hundred coincident
states. It checks inverse, norm, and 729 row configurations. Two objects
arriving from opposite neighbors occupy both lanes of the middle cell with
probability one; the declared identity collision coin retains them, and
the inverse countertranslation recovers the incoming configuration. This
tests the multiple-object completion rather than truncating to one object.

## Lifting every gate to the Gauss-constrained source

Let E be the tensor product of the lane isometries above, restricted to the
finite support of a gate. For any complete logical gate W define on that
support's full finite physical Hilbert space

\[
 W_{\rm source}=EWE^*+(I-EE^*).
\]

Because `E*E=I`, this operator is unitary and its action on the code is
exactly intertwined by E. Every interior Gauss operator equals +1 on the
range of E. The orthogonal complement is invariant under that operator,
and the gate is identity there. Hence W_source commutes with each admitted
interior Gauss constraint. This proves a local gauge-compatible completion,
not merely preservation of a projected matrix element.

The gates transporting an occupied lane transport its entire link code,
including the carrier boundary data and the exterior compensator. They
change the old per-cell boundary conditions by an explicitly declared
interface operation. No claim of transporting flux through separately
fixed carrier boundaries follows. There is one interior constraint per
lane; boundary constraints or continuous gauge fields are not supplied by
the interior isometry alone.

## A face-count mass coin and its controlled free window

The exact full-carrier face-count operator equals two on every one-object
codeword. Choose the source law to use that count in its one-object mass
coin,

\[
 C_a=\exp(-ia\,\beta M),\qquad M=2I,
\]

extended as above on the other occupation sectors. The unit coefficient
between face count and coin angle, the Clifford memory and this gate are
construction choices. The resulting value `m=2` is an output of that
specified program, not a consequence of topology without an energy law.
In particular, the general two-object coin is the stated completion; it
is not silently replaced by a different additive mass operator.

With spacing a and tick duration a in units `c=1`, the complete one-object
tick has Fourier symbol

\[
 D_a(p)=e^{-2ia\beta}
        e^{-ia p_3\alpha_3}e^{-ia p_2\alpha_2}e^{-ia p_1\alpha_1}.
\]

The existing routing theorem gives, for `t=na`,

\[
 \|D_a(p)^n-e^{-it(\alpha\cdot p+2\beta)}\|
 \le ta\sum_{r<s}|h_rh_s|,\qquad h=(p_1,p_2,p_3,2).
\]

For `|p_i|<=P`, the bound is at most `ta(3P^2+6P)`. The limiting generator
squares to `(|p|^2+4)I`. It therefore has signed frequencies
`+/-sqrt(|p|^2+4)` and rest-frequency magnitude two on this window. The
[spin/current attachment](../physical_attachment/PROOF.md) gives its actual
spin-one-half rotation action and conserved one-object current. The finite
ordered tick has the stated continuum-approximation error. Exact Lorentz
covariance at arbitrary finite lattice momenta is not asserted; boosts also
transform time and the field's hypersurface data.

The position detector is the source projector onto the cell containing
the nonflat lane. On the one-object code its sum over internal labels is
`E*Q_x E=|x><x| tensor I20`. Translation of the same link code supplies
the conjugate finite Fourier readout. Thus the spatial variable in this
free dispersion is the position of the paired-flux object itself.

Finite eigenphases are quasienergies modulo `2pi`; a signed frequency branch
requires the small-spacing, bounded-momentum branch convention. A negative
frequency is not identified as an antiparticle by this construction. A
positive-energy many-body vacuum, physical statistics, chiral gauge matter,
an interacting relativistic limit and an operational conversion to GeV
require additional results. The finite tensor program is fully specified;
those physical conclusions are not inferred from its one-object limit.
