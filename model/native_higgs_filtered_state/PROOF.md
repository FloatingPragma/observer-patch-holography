# Elementary and quantum-Higgs probes of the same filtered source state

## Source and exact energy selection

The source is the complete eight-cell native Hamiltonian
`h=H0+lambda T`, with lambda=2^-34 and ||T||<=360. The Yukawa, five-state
scalar and PW119 penalties have their source coefficients one. Every
native fermion occupation, scalar state and link state is retained.
The actual unique neutral ground is g at energy Eg. The full old/dual/port
reference has dimension 2^31+1 and its Hamiltonian is 4n on occupied
old/dual sectors, zero on the shared vacuum and even ports. These are
the source, ground and reference of the native abelian field join.

Use either radius2v or radius128v absolute-energy bin of
`model/native_charge_resolution`, where v=360 lambda. Subtract the
explicit reference cost four from its source-plus-reference bin to
define the **absolute native source** projector P=1_I(h). Its source
excitation enclosure [L,U] satisfies E-Eg in [L,U] for every E in I.
This enclosure is an alias for an energy bound; it does not define P.
The projector and its full complement are exact finite spectral services.
No phase-estimation implementation at these bin edges is assumed.

At the chosen source origin and spin component set

```
A_alpha=c_(alpha,ec)^dagger,
B_alpha=-sum_beta epsilon_(alpha,beta)
          (phi0^dagger c_(beta,10)+phi1^dagger c_(beta,11)).
```

The scalar matrices are the actual five-state operators. The native
Higgs spectral join proves the full-domain equation, reproduced by the
producer,

```
[h,A]=B+lambda[T,A],   ||A||=1, ||B||<=2, ||[T,A]||<=720.
```

Both probes are odd, color/weak singlets, with Y=1, 6Y=6 and raw lepton
number change -1. The actual trivial cubic ground character gives their
common origin G1 action. Their compensated opposite-charge channels
follow from the same-source antiunitary with conjugated reference.
The full scalar is never replaced by a number. In particular, the
nilpotency obstruction in `native_weak_current_interface` is compatible
with the filtered-state statement below.

## Ground filtered-state theorem

Write a=PAg, b=PBg, w=||a||^2>=w_lo, e=(L+U)/2 and d=(U-L)/2.
Using the actual ground equation and the fact that P commutes with h,

```
b-e a=(h-Eg-e)a-lambda P[T,A]g,
||b-e a|| <= d ||a||+720 lambda <= q ||a||,
q=d+720 lambda/sqrt(w_lo).
```

The exact rational arithmetic uses a downward square-root bound in the
denominator and rounds q upward. Both bins satisfy 0<q<e. Consequently
`||b||>=(e-q)||a||>0`. Project e a onto the line spanned by b. Its
orthogonal residual has norm at most q||a||, because b itself is a
point on that line. For normalized states this gives

```
1-|<a/||a||,b/||b||>|^2 <= (q/e)^2,
D(|a_hat><a_hat|,|b_hat><b_hat|) <= q/e.
```

Here D is half the trace norm and fidelity means squared Uhlmann
fidelity. This geometric bound is sharp for vectors of the stated
residual. It does not require dividing by e-q. The full operator norm
`||B-e A||` need not be small; only the spectral-bin ground response
has been bounded. Neither a canonical composite anticommutator nor a
constant scalar eigenbackground follows.

## Imperfect preparation, without a ground projector

The finite native preparation starts from the deterministic algebraic
trial. Let r0 be its certified excited-to-ground odds and epsilon_f the
mixed-sector amplitude of the finite agreement filter. Each binary
filter step obeys

```
r_next <= r^2+2 epsilon_f^2 r <= c r,  c=r0+2epsilon_f^2<1.
```

A depth d_p tree uses 2^d_p source trials and 2^d_p-1 actual filters.
All partners, spin references and filter outcomes are retained. Its
success probability and conditional main-source ground overlap obey

```
Q >= product_(j=0..d_p-1)(1+r0 c^j)^(-2^(d_p-j)),
1-<g,rho g> <= r0 c^d_p.
```

The first optional filter success appears in the j=0 factor. At depth
zero there is no filter, Q=1, and the unfiltered trial has its nonzero
state error. These formulas are the retained-tree result of
`interacting_block_propagation`; the underlying filters and raw clock
are the actual native preparation. They require no knowledge of Eg.

Keep a purification Psi of the complete conditional preparation.
Since g is pure, there is a normalized retained environment vector chi
such that, after a global phase choice,

```
||Psi-g tensor chi|| <= epsilon,
epsilon=sqrt(2 r0 c^d_p).
```

The same chi is used for both pump settings. It can include purification
registers for mixed initial data; tracing them cannot increase a bound.
No simultaneous exact P0 projection is applied. In the executed setting
d_p=3, all eight initial source trials and seven finite filters remain.

Let a_Psi=(P A tensor I)Psi and b_Psi=(P B tensor I)Psi. Set s=sqrt(w_lo).
Then

```
||a_Psi|| >= s-epsilon,
||b_Psi|| >= (e-q)s-2epsilon,
||b_Psi-e a_Psi|| <= q ||a_Psi||+(q+2+e)epsilon.
```

The last line follows by comparison to a tensor chi and b tensor chi,
using ||B-e A||<=2+e and ||a||<=||a_Psi||+epsilon. Thus, provided
epsilon<s and q_Psi<e,

```
q_Psi=q+(q+2+e)epsilon/(s-epsilon),
D(rho_A, rho_B) <= q_Psi/e,
F(rho_A,rho_B) >= 1-(q_Psi/e)^2.
```

The normalized successful states rho_A,rho_B can retain the entire
common preparation environment. Discarding any of it preserves the
inequalities. The same line-distance proof separately gives

```
D(rho_A,rho_A,ideal) <= epsilon/s,
D(rho_B,rho_B,ideal) <= 2epsilon/[(e-q)s].
```

These comparisons use the normalized successful purifications and then
contractivity, so no unnormalized norm error is silently reused after
postselection. Source preparation error is included before normalization.

## Complete pump, energy and subsequent-history records

Choose setting C=A or C=B/2 with classical probability one-half. Both
are contractions on the complete source. Attach their source gauge
charge to the even port with Y=-1. Retain a gauge-neutral fermionic
apparatus doublet of boundary charge -1, raw L=+1 and conjugate G1
action. A successful pump creates the selected apparatus fermion.
The combined operation preserves every Gauss constraint, total boundary
charge, raw lepton number and parity. The unused doublet mode and all
old/dual/port sectors remain in the carrier.

A Halmos dilation of C has the two maps C and sqrt(I-C^dagger C), with
the odd off-diagonal maps changing apparatus parity. On the occupied
pump branch measure P and I-P; on the unoccupied branch retain the
zero energy flag. The three nonzero Kraus maps on the declared input are

```
K_in=P C,  K_out=(I-P)C,  K_fail=sqrt(I-C^dagger C),
sum_j K_j^dagger K_j=I.
```

Their probabilities on Psi are their squared norms, with failure
`1-||C Psi||^2`. The full setting, apparatus-doublet and energy alphabet
has 2*4*2=16 records. Six can be nonzero on the specified initially empty
apparatus. All ten other records are retained as exact zeros. A regular
48-state binary-cubic reference lifts the complete apparatus and both
probes together. Every lifted Kraus map satisfies K_tilde J=J K, so
all probabilities and trace bounds are unchanged. These are joint,
reference-relative records. Source-only marginals are the corresponding
twirls, not the selected-frame state.

The source and full preparation environment give the success bounds

```
p_A,in >= (s-epsilon)^2,
p_B/2,in >= [(e-q)s-2epsilon]^2/4.
```

Multiplying each by Q gives its unconditional complete-preparation
success bound. Multiplying also by one-half gives the joint fair-setting
record. The preparation failure probability is 1-Q; each fine failure
label is retained, even where a displayed table coarsens it. The ideal
ground bounds are w_lo and ((e-q)s)^2/4. The producer reports both the
original unnormalized pump probabilities and their preparation-inclusive
versions. It does not rescale B/2 to equalize them.

Condition on successful preparation, pump and selected-bin readout in
each setting. Any common fixed subsequent detector policy, including
adaptive choices based on its own subsequent records and free waits
under the same h, defines a complete quantum channel into the joint
history register. Hence

```
TV(history_A,history_B) <= q_Psi/e.
```

This holds for arbitrary finite histories, with every detector complement
included, without a factor for history length. It does not compare
different policies that depend on the original A/B setting. Extra
conditioning on a rare subsequent event requires its own success bound.
The unconditioned A and B/2 laws need not be close: their pump success
probabilities differ. Retaining those flags is compatible with comparing
the normalized states they select.

For the verified same-source dressed weak density W=sum_a K_a j_a,
||W||=2, the fixed-time consequence is

```
|Tr[W(t)rho_A]-Tr[W(t)rho_B]| <= 4 q_Psi/e,
W(t)=exp(i h t) W exp(-i h t).
```

Its complete effects (I+W/2)/2 and (I-W/2)/2 are one possible history.
A gauge-singlet fermion can have nonzero dressed weak-current three-point
responses. The bound transfers a measured response between these two
filtered preparations; it does not calculate its value, turn variance
into scattering, or assert Standard Model V-A dynamics.

## Clock, control and capacity

The autonomous source clock is the original h, with the full reference
and retained source partners. For a finite desired unitary U over duration
tau, the explicitly declared programmable drive service chooses

```
L=i Log(exp(-i H_ext tau)U)/tau,
V=L-H_ext,
exp[-i(H_ext+V)tau]=exp(-i H_ext tau)U.
```

Use the Hermitian logarithm branch with phases in (-pi,pi]. The full
kinematic bound is `||H_ext||<=2^d_p(480+360lambda+120)`. All added
apparatus registers have zero Hamiltonian. Preparation, pump and
conditional energy readout each have duration one; a chosen bare wait
t>=0 gives elapsed time 3+t. The source evolves during every pulse,
including failed branches and partners. The common free factors commute
with the final pointer effects, and the post-filter source states evolve
by the same unitary. Thus the stated records and comparisons use this
physical schedule, with no source reset or silent pause.

Each pulse control has norm at most ||H_ext||+pi. Net bare-source energy
change is bounded by 2||H_ext||, separately from the gross switching bound
6(||H_ext||+pi). A realized drive error delta L_j changes the complete
purification by at most sum_j tau_j||delta L_j|| and the complete record
law by that amount. A successful branch with ideal norm at least r
then incurs normalized trace error at most eta/r when the total
purification error is eta<r. Conditional claims therefore require this
additional success-dependent allowance for approximate services.
Exact finite matrix functions are premises here; no efficient local
controller, passive C/E compilation or finite-clock implementation of
the bin projector is claimed.

With M the frozen finite agreement-filter pointer dimension, the retained
source/reference/pointer capacity in nats is

```
2^d_p [240 ln2+8 ln5+24 ln119+ln(2^31+1)+16 ln2]
 +(2^d_p-1) lnM +ln48+4ln2.
```

The last term counts the setting, pump doublet and energy flag. This
counts complete native kinematic carriers, all partners and records.
Programmable coefficient evaluation and its additional workspace are
separate apparatus premises.

## Typed free-spinor consequence of a refined shared energy

For either bin the source equation also gives
`Re<a,b>>=L w_lo-720lambda>0`. The sum of its spectral cross residues
is positive. Select the **actual absolute** native energy E_H,sharp in
the radius2v bin maximizing Re<A_up g,P_E B_up g>, with lower energy
breaking a tie. At least one maximizer is positive. Cubic Schur
orthogonality makes the cross matrix z_E I2, so this selected source
eigenspace contains an actual elementary G1 doublet with nonzero
composite projection. The conjugate doublet has the same native energy
under the proved source antiunitary and conjugated full reference.
Put m=E_H,sharp-Eg. Its strictly positive interval is the radius2v
excitation enclosure, not a measured mass input.

The free native-block Dirac construction is parametric in any such
selected positive doublet. Two copies of each charge encode vacuum,
both single occupations and the two-occupied singlet. Four copies give
all sixteen canonical fermion occupations. The producer replays the
exact exterior-action, parity and CAR identities. The constituent h,
ground, full reference and every source complement are unchanged. The
spectral matrix units **must be selected at this E_H,sharp**, and the
positive full-range kinetic matrix must be declared as
sqrt(m^2+L_a)-m. These choices define a separate free architecture;
they do not alter the frozen maximum-elementary-residue E* member.

On that declared reducing code the same proof gives H=dGamma(Omega_a)
and the positive/negative-frequency Dirac filters, complete detectors
and bounded-window free limit with this m. The inherited finite-volume
field error includes `T a^2(P+sqrt(3)h_p/2)^4/(24m)`, where h_p=2pi/L,
and is finite because the displayed lower bound on m is positive.
The literal source composite has only the proved vacuum-to-selected
doublet overlap. It does not equal the code field on all sixteen
occupations. No affine interaction, interacting continuum, observed
GeV mass or Standard Model weak chirality transfers to this free code.
