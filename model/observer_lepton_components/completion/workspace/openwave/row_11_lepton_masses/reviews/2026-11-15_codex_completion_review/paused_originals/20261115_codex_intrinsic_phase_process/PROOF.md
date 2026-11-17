# A positive native phase/reference process and its actual stopping law

Platform row11 remains untested. This packet constructs a new N+ memory
protocol on literal native currents. Its positive transfer and reconstructed
phase ratio are source results on that protocol. Identification of charged
lepton minima, a constant family interaction and masses is absent.

## 1. Native records and the continuing memory

Use the existing sixty ordered addresses of the twelve-port graph. A fresh
double at port2 is supplied for each draw. The observer records its address,
placement coin and native signed current j. Both placement coins give the
same current on this preparation. Exactly5 ordered addresses give j=1,
5 give j=-1, and50 give j=0. This is executed from the pinned native function;
the final loads are retained only as verifier data. No controller reads them.

For two independent draws X,Y put

    f=1[j_X=1],  g=1[j_Y!=0],  q=f XOR g.

All3600 address pairs remain. Thus

    Pr(f=1)=1/12, Pr(g=1)=1/6, Pr(fg=1)=1/72,
    p=Pr(q=1)=2/9.                                          (1)

This count itself was established previously. The new execution uses it in
a continuing memory process, without a coherent phase pulse or address
reflection. Four independent source draws supply the two phase counters in
each macrostep. Two more supply the pre/post oriented-memory refreshes.
Every source is prepared even if the associated update flag is zero. Spent
source cells and every fine record remain. Fresh-double production is a
declared preparation service; native transport of binary loads does not
recreate these double inputs.

A current refresh writes the two recorded bits `(r,o)`:

    r=abs(j),  o=1[j<0] if j!=0;
    o=the ordered-address orientation bit if j=0.

The four codes c=(0,0),(0,1),(1,0),(1,1) have probabilities

    pi=(5/12,5/12,1/12,1/12),
    J(c)=r(1-2o)=j.                                         (2)

The same polynomial defines the family object's oriented record in
`lepton_register_round/kick.py`. Here it is read from an actual current.
Both distinct J=0 codes are retained. This is a coordinate identity, not
an assertion that the entire original coherent family object has already
been implemented by this current-memory experiment.

The working state is x=(c,b0,b1,b2), with32 states. In a macrostep:

1. Independently refresh c with probability delta, otherwise hold it.
2. Flip b0 with probability epsilon/2. Flip b1 with probability epsilon p/2
   if J<=0, and b2 with probability epsilon p/2 if J=-1. The two p flags are
   fresh independent native XOR counts, and the three rare-update flags are
   independent native coin words.
3. Perform another independent current refresh of probability delta.

The rule is implemented with epsilon=2^-k, delta=epsilon^3. A native event
on a zero-load pair returns its supplied fair coin and current0. An all-zero
word of k+1 such coins has probability epsilon/2; a3k-bit all-zero word has
probability delta. Every other word has its original probability and remains
in the archive. No rare word is assigned probability one. Overwritten c,
the two refresh results and all unused count draws are archived as well.

The palindrome, the two count slots, the passive bits and their routing are
specified choices of this candidate. The construction has no claim that
all native observers must use this memory law.

## 2. One positive transfer and its reference

Use row-stochastic matrices: P_xy is the probability of y after x. Let T
refresh c and leave the three bits unchanged. On functions of c,

    T4=(1-delta)I+delta 1 pi^T.

Let X_a flip bit ba and P_a=(I-X_a)/2. Put A1=1[J<=0], A2=1[J=-1].
The counter-update kernel and whole macrostep are exactly

    D=(I-epsilon P0)(I-epsilon p A1 P1)(I-epsilon p A2 P2),
    P=T D T.                                                (3)

This is a derivation from the independent recorded update probabilities.
The right side is not installed as a separate event law. In the code it is
checked against a direct sum over both refresh results and the three flip
flags, on all1024 transitions.

The stationary law is rho(c,b)=pi(c)/8. Both T and D obey detailed balance
with rho. In the orthonormal representation

    S=diag(sqrt(rho)) P diag(1/sqrt(rho)),

T becomes `(1-delta)I+delta |sqrt(pi)><sqrt(pi)|`, tensored with the bit
identity. Its eigenvalues are1 and1-delta. D is a product of commuting
positive contractions, since each projector factor has its stated coefficient
between0 and1. Therefore

    m I <= S <= I,
    m=(1-delta)^2(1-epsilon)(1-epsilon p)^2 >0.              (4)

The inequalities follow by congruence: D>=dmin I implies T D T>=dmin T^2.
No approximate spectral diagonalization is used for positivity.

Every entry of P is positive. For an explicit path through one macrostep,
the first refresh can choose c=(1,1), where both phase counters are active;
the second can choose any target code. Each refresh contributes at least
delta/12. Every specified counter difference has probability at least
(epsilon/2)(epsilon p/2)^2. Hence

    P_xy >= delta^2 epsilon^3 p^2/1152.                      (5)

In particular, this includes transitions between different values of J.
Setting delta=0 before reconstruction would split that observation algebra;
the construction keeps delta>0 at every finite refinement.

The source reference observable is chi0(x)=(-1)^b0. Directly from (3),

    P chi0=(1-epsilon)chi0,
    E_rho[chi0(X0)chi0(Xn)]=(1-epsilon)^n.                  (6)

All currents and all unsuccessful rare-update words are included in this
correlation. The positive reference gap is

    g_epsilon=-log(1-epsilon).                              (7)

It is measured/reconstructed from the same process, independently of p.
Dividing the complete positive logarithm by g_epsilon fixes one common
reference-time scale. It does not select when a family analyzer must be read.

## 3. The full history algebra and the scope of reconstruction

The initial memory is public and every update is a function of the retained
native records. Its coordinate projectors Q_x are therefore observer fields.
Positivity (5) gives

    Q_x P Q_y=P_xy E_xy.                                    (8)

All1024 real matrix units are generated by these one-step words. This is
stronger than merely assigning an energy to three disconnected character
labels. The real matrix algebra has a canonical complexification. Selection
of that complexification as physical coherent observations remains an
interpretation, as in the earlier source reconstruction theorem.

For clarity the closed-return experiment has the separately disclosed uniform
preparation on32 memory states. It is a different permitted preparation of
the same process. Its L-step return probability is r_L=Tr(P^L)/32, with
nonreturn1-r_L retained. The original probability of x->y followed by L
further unrestricted steps ending at x is

    p_(L,x,y)=P_xy (P^L)_yx/32.

Consequently the scaled record statistic

    ell_L(E_xy)=p_(L,x,y)/(P_xy r_L)
               =(P^L)_yx/Tr(P^L)                           (9)

determines every matrix entry. In the weighted-coordinate algebra the adjoint
is A^sharp=diag(rho)^-1 A^T diag(rho). After the similarity in section2,

    ell_L(A)=Tr[S^L A_tilde]/Tr(S^L),
    A_tilde=diag(sqrt(rho)) A diag(1/sqrt(rho)).              (10)

Thus ell_L is positive on the corresponding positive effects, with density
S^L/Tr(S^L). Its value is a reconstructed statistic; division by the rare
path coefficient is not a successful preparation assigned unit probability.
Every original closed/nonclosed outcome is kept. The finite run checks (9)
at L=1 for every matrix unit.

The positive transfer also gives, for every finite collection of vectors v_n,

    sum_(m,n) <v_m,S^(m+n)v_n>
      = ||sum_n S^n v_n||^2 >=0.                            (11)

The same reversible Markov memory histories supply the usual reflected
cylinder forms, by folding their conditional future vectors; the positive
one-step form handles a reflection through a half-step. This establishes a
positive finite reconstruction for the specified memory-observable algebra.
It does not assert reflection positivity for every unreversed fine raw
source label or for the growing spent-source bank. Those records remain in
the complete original instrument; summing them to form a memory observable
is not deletion of their outcome probabilities. A physical common-clock
claim on a larger algebra would need that larger attachment separately.

The unique positive logarithm

    H_epsilon=-log S / g_epsilon                            (12)

then defines reconstructed unitary correlations by spectral analytic
continuation. No endpoint-logarithm 2pi ambiguity occurs: S is positive and
its real logarithm is fixed. This does not prove that a reduced classical
diagonal observation algebra has a nontrivial continuous automorphism, or
that coherent gates are being physically enacted between native attempts.
The reconstructed full matrix algebra and the directly committed operational
algebra are explicitly different objects.

The finite cost is visible. Equation(5) tends to zero as the refinement is
made; with delta=epsilon^3 its bound is order epsilon^9. Matrix-unit estimates
divide by these small coefficients and by r_L. Exact finite tables therefore
do not assert a bounded sample cost in the limit. A full empirical use needs
the corresponding accuracy/resource growth. The finite positivity theorem
and its operator limit do not use a zero-cost sampling assumption.

## 4. Every character sector and the actual conditional J operator

All8 bit characters chi_a(b)=(-1)^(a dot b) are retained. In character a,
the four-code transfer is

    P_a = T4 diag[d_a(c)] T4,
    d_a(c)=(1-epsilon)^a0 (1-epsilon p)^n_a(c),
    n_a(c)=a1 1[J(c)<=0]+a2 1[J(c)=-1].                    (13)

Its symmetric form uses the same stationary pi and the same reference gap.
No individual character receives a different time normalization.

Define h_epsilon=-log(1-epsilon p) and
r_epsilon=h_epsilon/g_epsilon. Without the small refresh perturbation,
the normalized logarithm is exactly

    H_D=P0+r_epsilon(A1 P1+A2 P2).                          (14)

The generator couples to the ACTUAL recorded J through A1,A2. It is not
obtained by attaching different names to three unrelated external rates.
For the fixed bit character a0=0,a1=a2=1, the active count is

    A1+A2=1-J.

Thus its limiting operator is p(1-J), acting on all four native record codes,
including both J=0 codes. With a common scalar reference removed its
reconstructed real-time action is exp(i p t J). The sign of complex
continuation fixes an orientation convention; conjugation reverses it.

The other characters have the following limiting phase energies, before the
reference bit adds a0 to every entry. Columns list J=+1, the two J=0 codes,
and J=-1 respectively:

| (a1,a2) | Four retained entries |
| --- | --- |
| (0,0) | 0,0,0,0 |
| (1,0) | 0,p,p,p |
| (0,1) | 0,0,0,p |
| (1,1) | 0,p,p,2p |

There are32 states, not three particle minima. At finite delta the J sectors
mix. In the `(a1,a2)=(0,0)` block the three nonconstant refresh modes have
gap `-2log(1-delta)/g_epsilon`, tending to0. They are not discarded.

The same J polynomial allows a CONDITIONAL family-object attachment:
tensor the fixed character corner with the untouched load/E/family factors
and identify its normalized coordinate vectors in the S representation with
the existing object's orthonormal `(r,o)` record basis. The weights pi are
not silently substituted for the original family's preparation amplitudes.
Then the limit acts as `exp(i p t J)` on that factor, up to a common phase.
This operator intertwiner is exact mathematics. To obtain the old family
analyzer as a measured response of this candidate, its preparation, U0,
analyzer and all complementary operations must be realized as observations
or interventions of this same history algebra. They are not supplied by the
coordinate identity. In particular, an old CP instrument certificate does
not execute this reconstructed coherent action.

Static operator insertions can be represented as linear combinations of (8),
with their supplied coefficients and the costs in section3. A family analyzer
whose additional factors are absent from the32-state source requires an
explicit extended source/algebra and preparation. Such a source join is an
outstanding construction, not silently included in this packet.

## 5. Uniform controlled refinement

In the symmetric representation ||T-I||=delta, ||D||<=1 and ||T||<=1, so

    ||T D T-D|| <=2delta.

For positive A,B>=mI, the resolvent integral for the logarithm gives
`||log A-log B||<=||A-B||/m`. Hence, with m from (4),

    ||H_epsilon-H_D|| <=2delta/(m g_epsilon)
                       <=2delta/(m epsilon).               (15)

The ratio has 0<r_epsilon<p. An exact integral identity is

    p g_epsilon-h_epsilon
      =p(1-p) integral_0^epsilon t/[(1-t)(1-pt)] dt.

Since g_epsilon>=epsilon, this proves

    0<p-r_epsilon
      <=p(1-p)epsilon/[2(1-epsilon)(1-p epsilon)].            (16)

Let H_infty=P0+p(A1P1+A2P2), on the WHOLE32-state space. The count operator
has norm2. Combining (15)-(16) yields

    ||H_epsilon-H_infty|| <= B_epsilon,
    B_epsilon=2delta/(m epsilon)
                 +p(1-p)epsilon/[(1-epsilon)(1-p epsilon)]. (17)

For delta=epsilon^3 this vanishes uniformly. Duhamel's formula gives

    ||exp(-it H_epsilon)-exp(-it H_infty)||
        <=|t| B_epsilon.                                   (18)

The corresponding unitary-channel trace/diamond bound is at most
`min(2,2|t|B_epsilon)`. These are bounds on the reconstructed finite
descriptions, not hardware pulse certificates. In the fixed phase character,
they compare with the actual conditional J operator above. The finite ratio
is not exactly p at nonzero epsilon and is not a result of rescaling the
family clock independently. At delta>0, h_epsilon by itself was introduced
as a gap of D; it cannot be called a pole of the whole P merely by that
definition. The following exact full-process attachment resolves this point.

### Post-run analytical pole attachment

This lemma was added during root's written review. It uses the frozen kernel
and adds no execution or check count. Set

    f(c)=1[c=(0,0)]−1[c=(0,1)],
    f1(c,b)=f(c)(−1)^b1.

Because the two J=0 codes have equal pi weight, pi(f)=0. Thus T f=(1−delta)f.
The counter update leaves f unchanged. On f's entire nonzero support, J=0
and the b1 counter is active with exactly the same probability. Consequently

    P f=lambda_r f,       lambda_r=(1−delta)^2,
    P f1=lambda_p f1,     lambda_p=(1−delta)^2(1−epsilon p).

Both fields have stationary squared norm5/6. Their ORIGINAL correlations,
without conditioning on J=0, are

    E_rho[f(X0)f(Xn)]=(5/6) lambda_r^n,
    E_rho[f1(X0)f1(Xn)]=(5/6) lambda_p^n.                    (18a)

All states outside the J=0 support give the retained value0; they are not
omitted trials. Together with the unit-weight reference correlation(6),
these three nontrivial full-process poles determine two gaps:

    h_epsilon=−log(lambda_p)+log(lambda_r),
    g_epsilon=−log(1−epsilon).

Their ratio is exactly r_epsilon. This is a same-process difference of
reconstructed energies, rather than a claim that h_epsilon is an unshifted
eigenvalue of the refresh-coupled transfer. In the count-HOLD control the
two poles coincide; under the count-one plant their relative gap equals
the unchanged reference gap. The fixed source correlation weights5/6 and1
make this calibration operationally explicit. The algebra proves it at
every allowed epsilon and delta, with no new numerical run.

The script encloses each logarithm with100 terms of
`-log(1-x)=sum_(n>=1)x^n/n`, with omitted tail bounded by
`x^101/[101(1-x)]`. All inequality gates use rational enclosures. Decimal
strings in the receipt are diagnostics of those exact bounds.

## 6. An actual reference stop has a distribution, not its mean

The b0 flip word contains k+1 native fair bits. Its first k bits being zero
is a reference MARK, of probability epsilon; the final fair bit decides the
flip. The mark is independent of the current/phase-memory evolution. Stop
on the first mark, retaining its whole word and the elapsed macrostep count N.
Then exactly

    Pr(N=n)=epsilon(1-epsilon)^(n-1), n>=1,
    Pr(N>L)=(1-epsilon)^L.                                 (19)

At the mark b0 has a fair flip. Before it b0 has not flipped. All other
source and phase-memory updates occur even on the terminal macrostep.

Let A be the16-state current/phase-memory transfer obtained by summing both
reference-bit outcomes of P. Its complete original stopped transfer is

    R_epsilon=epsilon A[I-(1-epsilon)A]^-1.                 (20)

The full final b0 map is the additional `(I+X0)/2`, and N remains a retained
record. Summing N in (20) describes a marginal, not an unnoticed timeout
removal. Finite-horizon duration/state laws through8 and the complete timeout
distribution are executed in the receipt.

Write A=exp(-g_epsilon H_phase,epsilon). By (17), its limit in (20) is

    R_epsilon -> (I+H_phase,infty)^-1.                       (21)

This follows uniformly on the finite dimension, using
g_epsilon/epsilon->1 and continuity of the resolvent; every denominator is
at least epsilon. On the fixed `(1,1)` phase character the limiting entries
are `1/[1+p(1-J)]`. This is an actual source-stopped Euclidean response. It
differs from evaluating exp(-H_phase) at the mean reference duration.

There is a related CONDITIONAL prediction if a future coherent family
implementation uses the reconstructed real-time phase at each retained
duration. Neglecting only the controlled refresh error, its phase is h_epsilon N.
The full, unconditioned harmonic is

    chi_(q,epsilon)
      =epsilon exp(i q h_epsilon)
        /[1-(1-epsilon)exp(i q h_epsilon)],
    chi_(q,epsilon) -> 1/(1-i q p).                         (22)

This is not exp(i q p). For q!=0 its limiting modulus is strictly below1.
The normalized duration g_epsilon N converges to an exponential law of mean1;
the mean does not specify its phase average. The expected error from (18)
is bounded by `(g_epsilon/epsilon)B_epsilon`, since E[N]=1/epsilon, so the
refresh approximation in this conditional average is controlled too.
The duration label is retained classically in this conditional family read;
no interference between distinct N branches is asserted. Their individual
common scalar phases therefore cancel in each complete branch effect.

For the inherited family expression

    M(theta)=aI+b(exp(i theta)R+exp(-i theta)R*),
    a=1/(3sqrt3), b=1/(3sqrt6), R^3=I,

the resulting conditional averaged effect retains all durations and is

    A_bar=(a^2+2b^2)I
       +2ab[chi_1 R+conj(chi_1)R*]
       +b^2[chi_2 R^2+conj(chi_2)(R*)^2],
    chi_q=1/(1-i q p).                                     (23)

Equation(23) is an analytical target on a still-unrealized family/source
join. No coherent family pulse or effect has been executed here. It shows
why replacing a native random dwell by one unit would change the prediction.
The unitary analytic continuation of the stopped transfer's own logarithm
is another operation again; it is not the mixture of real-time phases in
(22). Both distinct descriptions are explicitly separated.

A fixed observation at t=1 instead gives the conditional limit M(p), as in
the earlier proposed formula. Choosing that action window is a declared
protocol. The reference scale alone does not supply a native stop at t=1,
a continuing renewal rule that selects it, or a constant lepton mass. The
native first-mark stop supplied here gives (20), and under the extra coherent
family reading it gives (23), rather than M(p)^*M(p).

## 7. Controls and the precise advance

The count-source HOLD control replaces endpoint movement only in the four
count slots; all reference and oriented-current refresh slots remain native.
It sets both XOR counts to0. The phase bits never change, while (6) remains
unchanged. The original all-source/coarse matrix kernels of this control are
retained. The AND protocol instead has p=1/72. The forced-count1 plant has
p=1 and is detected by the same phase-update read. In the frozen J=-1
middle-update setting, native b1 flip probability is epsilon p/2; HOLD is0,
AND epsilon/144, and the plant epsilon/2. These are original probabilities.

The unit character identity and the bit dimension alone carry no physical
clause. The source dependence is in the current flags, their measured p,
the actual transition law and the failing HOLD comparison. The one-reference
ratio, positive reconstruction and conditional J operator constitute a new
source-controlled route to an oriented phase description. It fixes a
frequency ratio within this specified protocol and removes the endpoint
branch ambiguity for its positive transfer. It does not prove a universal
beta=1, identify a physical charged lepton, or select the family analyzer's
action window. No row gate or OpenWave award follows.

The frozen finite run passed37 checks, executing960 literal native columns
and480 separately counted changed-HOLD calls. It reconstructs four full
32-state kernels, all8 character blocks at each refinement,4096 history
matrix-unit words, complete first-mark/timeout distributions through8, and
exact logarithm enclosures. At epsilon=1/64 the ratio enclosure is near
0.2208652231697421 and the whole-generator bound(17) is0.003252541649632.
The first run took0.787540s. Complete producer replay reproduced the receipt
bytes. These are finite implementation checks of the displayed derivation;
they do not execute a coherent family analyzer or certify a physical row.
