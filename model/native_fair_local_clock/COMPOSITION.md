# A prepared source with complete local readout

The preparation, repair and readout below are commands of one admitted
fair binary source programme. The two comparison channels are specified
on its quantum and record outputs. Neither comparison inserts a source
measurement or an independent archive into the executed programme.

## Reference preparation and oracle

Write $J_b$ for the owning eleven-E tree followed by $b$ double sweeps of
the thirty original-cell seams. Its output retains the eleven tree records
and the source, with an arbitrary untouched reference $R$; only the
$60b$ mixing-record values are traced by this observation map. The source
physically retains every original FIFO and exported record. The owning
preparation theorem gives a CPTP comparison $J_*$ with the same tree
archive and a sector-uniform conditional source reference, including its
empty/full coherent corner, such that

$$\tfrac12\|J_b-J_*\|_\diamond\leq\delta_b=2^{-b}. \tag{1}$$

The global source's exterior is part of the untouched reference during
this internal preparation. Subsequent source repairs may act on that
exterior. Equation (1) remains a bound on their joint input, and any common
complete continuation contracts it. The continuation must not consult the
mixing-record values marginalized in (1).

Let $D$ be the complete owning readout, including its initial-local-mask or
erasure flag, records, terminal disturbance and $R$. Its initial mask
refers to the occupations **at readout start**, after preparation. Let $O$
be the CPTP comparator that measures that same initial local mask first
and then runs the identical native readout instrument, retaining the true
label. This is the readout theorem's mathematical comparator. For the
actual postpreparation input its full-output error is

$$D_{\rm tr}(D J_b\rho,O J_b\rho)
 =\gamma=p_{\rm cov}+(1-p_{\rm cov})\eta_{\rm amb}. \tag{2}$$

This formula applies to $0<q<V$. At $q=0$ or $q=V$ the declared global
source has only one configuration, so the specification-aware decoder has
$\gamma=0$ without a coverage requirement.

Here $\eta_{\rm amb}$ is the actual probability of an unresolved local
empty/full input after any allowed global-number feasibility refinement.
Equation (2) retains the specified full history labels; further output
marginalization preserves the upper bound. It holds even when $R$ contains
the retained preparation tree archive and is correlated with the source.
The subsequent uniform seam selector is independent of those records.

Apply the same $O$ to both sides of (1) and use the triangle inequality:

$$D_{\rm tr}(D J_b\rho,O J_*\rho)\leq\delta_b+\gamma. \tag{3}$$

This is a complete source/reference and recorded-output comparison, not
a successful-branch normalization. Equation (2) also bounds the actual
readout against its *own* preparation oracle without the term $\delta_b$.
The two targets are distinct.

## Which input facts the composition needs

Every preparation seam lies inside the detector region $B$. Each complete
native branch preserves its occupation number $N_B$. Thus the preparation
preserves $\operatorname{Tr}(P_{\emptyset}^B+P_{\rm full}^B)\rho$ exactly,
for every input and reference. The ambiguity weight in (2) can be taken
from the actual initial input before preparation. No canonical state or
independence follows from this conservation identity.

There is a stronger statement for a declared canonical source marginal.
On the full four-dimensional endpoint factor the complete fair E channel
is unital: summing all 25 native raw slots gives $\Phi_E(I)=I$. Each
branch also conserves endpoint number. Tensor extension therefore gives
$\Phi_E(P_q)=P_q$ on every full global fixed-number sector. Any fixed E
word preserves $\tau_{V,q}=P_q/\binom Vq$. Its action on a correlated
extension can change correlations, but its unread source marginal is
exactly $\tau_{V,q}$. The conditional reference map in $J_*$ preserves the
same unread marginal. Neither statement factors the retained archive or
requires an independent external $R$.

Consequently a canonical source marginal, with an arbitrary correlated
reference extension, supplies the readout's hypergeometric ambiguity
weight. This use of a marginal is specific to (2). It does not provide
the *joint* product-reference estimate required by the separate spatial
Poisson comparison in `PROOF.md`.

## One finite canonical protocol

Take level six, $V=983040$, $q=V/3$ and $M=2580480$. Use $b=64$:
the preparation consists of $L=11+60b=3851$ fixed E attempts. Append one
uniform E command of $T=64M=165150720$ attempts. The literal programme has
3852 command nodes and $L+T=165154571$ attempts, with partition

$$\left(0,\frac{3851}{M},\frac{3851}{M}+64\right). \tag{4}$$

There is no zero-repeat or extra background command in the preparation
bin. The admitted finite rooted detector contains 348 ports in the
radius-four whole-cell ball, with all 922 incident seams authenticated
inside the radius-five window. Its tree has 347 edges. The native
resource calculation, including the retained-record bit bound, accepts
the combined programme.

The canonical ambiguity weight obeys
$\eta_{\rm amb}\leq(2/3)^{348}+(1/3)^{348}<10^{-60}$ and the uniform
coverage error is at most $347/2^{64}$. Equations (1)--(3) give

$$D_{\rm tr}(D J_{64}\rho,O J_*\rho)
 \leq\frac{348}{2^{64}}+\eta_{\rm amb}<10^{-16}. \tag{5}$$

This holds for any reference extension of the stated canonical source
marginal. The eleven preparation-tree records and every selected readout
record are retained. Only the preparation's mixing-record values are
marginalized for the comparison. The source has the full archive.
An arbitrary source marginal instead requires its actual ambiguity weight
in (2). A wholly empty/full local input can make that error large.

The 64-time-unit readout is the literal long protocol (4). Its bound (5)
uses readout coverage and the preparation contraction. It is not the
short-time Poisson/window estimate, nor does it claim a completed expanded
level-six host run.

## Fixed-length sparse-sector protocol

Public global number is part of the authenticated source specification.
For a twelve-port root, empty is feasible exactly when $q\leq V-12$,
and full is feasible exactly when $q\geq12$. If these possibilities do
not coexist, complete tree records determine the unique feasible constant
mask as well as every nonconstant mask. A fresh eleven-E root tree has
zero oracle error on that entire global sector. This uses no assumed
local number.

In particular take global $q=2$. All feasible local sectors $0,1,2$ are
retained, with arbitrary exterior and reference coherences. Append that
fresh eleven-E tree to the $b=64$ preparation. One admitted programme has

$$3851+11=3862\text{ fixed E attempts and command nodes},\qquad
D_{\rm tr}(D J_{64}\rho,O J_*\rho)\leq2^{-64}. \tag{6}$$

Every fresh readout record and both terminal quantum systems remain in
the output. Preparation mixing-record values have the same explicit
marginalization as (1). The appended readout does not consult them. The
combined two-tree protocol permits $b\leq67$ under the unchanged node
limit; the preparation-alone value $b=68$ cannot be appended to this tree
inside that limit.

The 3862-attempt instrument is independent of source capacity on its
addressed root factor wherever its original admission checks pass. Its
per-seam duration is $3862/M$, which decreases with capacity. The result
does not keep a positive per-seam duration or density fixed. In a canonical
global two-particle source, active local mass can vanish in the large
capacity limit; no uniform coherent signal or observed particle identity
is inferred. Fixed high-level cell words remain subject to the source's
unevaluated-carrier stop, and the complete finite-word theorem does not
override an expanded host's branch ceiling.
