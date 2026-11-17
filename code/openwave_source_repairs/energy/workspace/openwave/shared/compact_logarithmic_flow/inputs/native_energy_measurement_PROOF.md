# A complete finite energy measurement using native source services

The fair twelve-port source fixes a positive transfer. Its energy can be
measured by a complete controlled instrument on the same sixty-six-dimensional
carrier, with errors bounded independently of spectral degeneracies. This
construction uses the declared proper-code coherent services. Reading the
repair source specification does not cause passive repair to issue those
services.

## The same fair source and the measured operator

Use binary total occupation two, the thirty original seams, independent
uniform selection, fresh fair endpoint coins, repair probability one and the
fixed C/E mixture. In the original source configuration basis,

    P=I-G/60,       G=sum_e(I-S_e),       (2/3)I<=P<=I,
    H=-60 log(P).                                                        (1)

The counted convention is t=N/60. Each entry of P is obtained from the actual
native endpoint rule; the exact integer matrix 60P is symmetric, stochastic
and has diagonal at least 50. The Gershgorin bound gives the lower bound in
(1). Its basis is checked equal, entry by entry, to the existing coherent
sixty-six-configuration chart. Changing fair to weighted repair probabilities
does not change this chart. The weighted source's eigenvalues are not reused.

The fixed rational proposal supplies columns v_j and transfer values lambda_j
on the 2^-44 grid. A numerical eigensolver is only their proposer. Exact
rational arithmetic checks

    ||V*V-I||<=g=10^-9,    ||PV-V Lambda||<=r=10^-9,
    2/3<=lambda_j<=1.                                                    (2)

Absolute row and column sums bound the operator norms. With the mathematical
polar basis W=V(V*V)^(-1/2), the established finite polar estimate gives

    ||P-W Lambda W*||<=r+2g,
    ||H-H_tilde||<=epsilon_H=60(r+2g)/(2/3)=27/10^8,
    H_tilde=-60 W log(Lambda) W*.                                       (3)

The logarithm bound follows by integrating the positive resolvents. No polar
gate or numerical eigensolver output is executed as an ideal operation.
For completeness, g<1 makes V invertible. Its singular values lie between
sqrt(1-g) and sqrt(1+g), so ||V-W||<=g from
|sqrt(x)-1|<=|x-1|. Since ||P|| and ||Lambda|| are at most one,
||PW-W Lambda||<=||PV-V Lambda||+2||V-W||<=r+2g. Multiplication by W*
proves the first bound in (3). This is the same polar comparison used by
the [controlled transfer proof](../native_transfer_control/PROOF.md).
The actual compiler uses u_j=v_j/||v_j||, for which ||u_j-w_j||<=2g and
||P_(u_j)-P_(w_j)||<=4g. These estimates require no eigenvalue gap.

## One complete destructive filter

Use the existing nine-bit work register and the source embedding C. Let p be
one existing padding address and set

    z=(C u+|p>)/sqrt(2),       U=I-2|z><z|.                              (4)

The amplitude tree followed by phases prepares z from literal work address
zero. The existing rank-one phase word at angle pi implements U. Exporting
the source chart gives the source-success map and padding row

    C*UC=I-P_u,             <p|UC=-<u|.                                (5)

All other ideal padding rows vanish but remain in the complete alphabet.
The padding flag is an observed destructive detection of energy label j;
it retains the untouched reference and terminates the experiment. It does
not leave a source eigenstate available for further use. The source-success
branch retains (I-P_u)rho(I-P_u) and proceeds to the next stage.

This is the survival-zero case of the
[native weak filter](../native_coexistent_effects/PROOF.md). Its correcting
phase word acts as the identity on the successful complement and is omitted.
Equation (5) proves that omission, including the state update. No projector,
postselection, extra ancilla or free probability draw is inserted.

## Sequential energy readout with every outcome retained

For j=1,...,66 put A_0=I and A_j=(I-P_(u_j))A_(j-1). The jth energy effect
and final quantum residual effect are

    E_j=A_(j-1)* P_(u_j) A_(j-1),       E_res=A_66* A_66.

Their sum is I by telescoping. The complete instrument has reference-valued
terminal blocks <u_j|A_(j-1) rho A_(j-1)*|u_j> and a residual source block
A_66 rho A_66*. All effects are positive; no survivor is renormalized.
Degenerate energies retain their separate measured basis labels. Coarsening
those labels is permitted only after the complete instrument is defined.

Replacing u_j by w_j gives the complete destructive spectral measurement of
H_tilde; its residual is zero because the w_j form a full orthonormal basis.
For one stage a Stinespring isometry consists of the detection row <u| and
the successful block I-P_u on orthogonal flags. The norm difference from the
w stage is at most ||u-w||+||P_u-P_w||<=6g. The half-diamond distance between
the complete channels is therefore at most 6g. Composition, with identity
continuation on earlier terminal flags, gives

    d_diamond(instrument_u,instrument_W)<=396g.                         (6)

This is a complete instrument bound for arbitrary entangled source/reference
inputs. It compares destructive detections; it asserts no nondestructive
Luders state update.

## Exact words, numerical gates and finite resources

Each stage is the original admitted `phase_word(u_j,1/2,pi_grid)` with 58,392
pulses. There are 66 separately admitted source words, totaling 3,853,872 pulses.
Every original pulse and its negative or zero-angle request is retained.
Each word obeys the 65,536-node and 2^24-byte description ceilings. The
controller advances on source success and records a terminal stage otherwise.
It needs a seven-bit stage counter and the original sixteen-bit position
counter per word. The literal pulse definitions can be generated and executed
one word at a time from the bounded proposal; all words are not stored in one
monolithic repair program.

This is a bounded composition of declared source services. It is not an
admitted single repair AST or single coherent word. No claim is made that
the passive repair grammar compiles this controller automatically. The
generic coherent evaluator retains its six-bit/256-pulse RESOURCE_STOP.
The direct nine-bit compiler evaluates the success prefixes. All first
processor failures are defined by the original pulse, complement label and
successful prefix; each export keeps all 446 padding rows behind its one
observed flag. The final source residual is retained. There is no expansion
of every complete history.

For N=58392, pulse gate h=2*10^-11 and angle grid Q=2^41, the existing exact
stored-array certificate gives, per word,

    x=Nh,      e=x/(1-x)+23/Q,       f=N h^2/(1-x)^2,
    delta=[e(2+e)+f]/2,             gain=(1-x)^(-2)+f.                  (7)

The 23 grid units include both vector preparation errors and the pi rounding.
Every emitted primitive must pass its directed interval check. The packet
performs these checks for all 66 words. A rejected primitive rejects the
receipt. Combining (6) and (7) gives

    complete half-diamond error <1/12000,
    complete stored-map gain <1001/1000.                               (8)

The ideal service instrument is trace preserving. Stored numerical maps have
the bounded gain in (8); they are not declared exactly trace preserving.
Floating success-prefix propagation is an independent diagnostic. It is not
the exact-product error proof or a hardware error model.

## Preparation, clock and finite energy accuracy

An optional preparation uses eleven original fair E tree commands. Every
live fine branch is 2^(-11/2)|y><x| and the records decode x and y. Conditional
on y, existing X routes move its source address to one chosen source
configuration. They use at most 18 pulses and preserve every tree record.
The ideal prepared source factors from the archive and arbitrary untouched
reference; the archive can remain correlated with that reference. Zero raw
words remain zero maps. The code checks the fair raw weights directly; it
does not substitute the weighted detector's branch probabilities.

The full controlled experiment has eleven repair attempts, at most 3,853,890
coherent pulses, one repair FIFO slot and ten exported tree records. The
complete source-service flags and first-failure labels are retained outside
that repair FIFO; no finite FIFO claim is made for storing the full service
archive. The counter pair (repair attempts, service pulse ordinal) orders the
experiment. It does not convert pulse durations into repair attempts or
laboratory seconds. The measured label H is fixed by (1), independently of
the duration needed to implement its detector. The optional route error has
its separate small complete-channel allowance.

Each energy label E_j=-60 log(lambda_j) has an exact rational interval from
the 96-term positive logarithm series. For a normalized eigenstate of the
actual H with energy E, (3) gives

    sum_j |<w_j,psi>|^2 (E_j-E)^2<=epsilon_H^2.

For any positive margin d, the probability of |E_j-E|>=d is at most
(epsilon_H/d)^2. Add 396g for the ideal compiled ray instrument. For the
stored bounded-gain map, all failure and residual flags can be counted as
bad labels and their unnormalized score mass is bounded by

    (epsilon_H/d)^2 + 2/12000.                                       (9)

At d=1/1000 the nearby spectral term is 7.29*10^-8. No spectral isolation or
relativistic dispersion is inferred from this finite detector accuracy.

The construction measures a source-derived counted transfer energy through
existing coherent services. A growing all-number coherent loading map,
autonomous field dynamics, physical time normalization, stress tensor and
the identification of rest masses require further source-to-field results.

## The source-prepared transfer reference and modular clock

The existing [native phase experiment](../native_phase_experiment/README.md)
prepares the count-derived ideal reference rho=P/Tr(P)=P/56. Its independent
rational implementation supplies the same sixty-six configuration coordinates
and every entry of the same fair P. The attachment checks those equalities
against (1), and checks the owning preparation receipt with its independent
verifier. No second preparation law or imported target spectrum is introduced.

Since P>=2I/3, this reference has rho>=I/84. Functional calculus gives

    K=-log(rho)=H/60+log(56)I,
    rho^(i s) A rho^(-i s)=U_H(s/60) A U_H(s/60)*,
    U_H(t)=exp(-i t H).                                                (10)

Thus the modular action of this prepared reference and the controlled transfer
action are exactly the same finite operator family under the displayed
parameter conversion. The reference is stationary for U_H. It is not
stationary for unread native C/E repair; the owning phase experiment measures
that distinction. Equation (10) does not identify a service pulse duration
with a native attempt or a laboratory second.

The preparation uses the declared twelve-bit occupation load and export, and
the measurement uses the declared nine-bit source chart. Their sixty-six
configuration basis equality authenticates the state dictionary; the proper-code
interfaces remain operational premises. The finite selector retains its timeout
(23/128)^16 and all failures. Its numerical full-density admission is a
2*10^-10 Frobenius bound on the represented arrays, without certified outward
rounding. Equation (10) is an exact statement about the ideal source-derived
reference. This attachment does not promote the preparation's numerical
admission to an arbitrary-reference channel bound or hide timeout by normalization.

This reference differs from the diagonal tree-weighted preparation used below.
The native occupation exchange observable on ports zero and one is a direct
sum of two-state swaps and zeros, hence has norm one. The independent exact
source calculation gives Tr(rho A)=1/168, whereas every diagonal occupation
reference has Tr(sigma_tree A)=0. Consequently

    d_trace(rho,sigma_tree)>=1/336.

Regional modular compatibility proved for the tree-weighted reference cannot
be combined with (10) as properties of one reference without an additional
map. Neither construction by itself supplies all geometric stress premises.

## A complete repair-prepared energy experiment

The [retained-record preparation](../native_noncentral_reference/PROOF.md)
supplies a second, entirely native preparation. On the fair total-two carrier,
eleven tree attempts and sixty ordered E attempts give an exact acceptance
flag with probability a=7771/65536, independent of the input. Conditional on
each retained tree record and that flag, after marginalizing the sixty mixing
record values, the accepted source is

    sigma=diag(W(x)/409),     W(x)=product_(i occupied in x)(1+depth_i).

The depths come from the actual rooted native tree. No chosen transfer ground
is substituted for this density. The input reference can be arbitrary and
entangled. The source factors from that reference on the stated accepted
observation. The actual archive retains all seventy-one records.

Append the same energy instrument only on acceptance. Its coherent requests
depend on the acceptance flag, not on the marginalized mixing values.
Rejection terminates with its actual source/reference block. The complete
reported outcome law is

    Pr(REJECT)=1-a,
    Pr(j)=a Tr[E_j sigma],       Pr(RESIDUAL)=a Tr[E_res sigma].        (11)

Every term is nonnegative and the sum is one. No rejection, energy label or
residual is removed. The code evaluates (11) exactly. Writing the integer
proposal ray as v_j, N_j=v_j^T v_j, and A_(j-1)=B/D with integer B,D, its
energy numerator is

    Tr[E_j sigma]=sum_i W_i (v_j^T B)_i^2/(409 N_j D^2).

The exact residual update is B<-N_j B-v_j(v_j^T B), D<-N_j D. This avoids
normalizing square roots, numerical matrix products and a sampled polar basis.
The final residual is evaluated from B^T B with the same exact weights.

The fair repair program is admitted with seventy-one nodes/attempts, one FIFO
slot and seventy exported records. Its expanded-host step ceiling remains
sixty-four; the complete kernel and path-rank postprocessing are the owner's
execution domain. Coherent measurement uses the separate finite service
composition above, with at most 3,853,872 pulses on the accepted path. The
native attempt counter is seventy-one at that path's measurement start;
the pulse ordinal is a separate counter. No physical duration map is assumed.

Because acceptance is input-independent and the rejected instrument is
unchanged, (8) gives full unnormalized half-diamond error below a/12000 for
this composition, on the stated observation that marginalizes mixing values.
It does not bound a comparison which uses those values for extra feedback.

For a sampled coarse histogram, group every processor failure into one label
while retaining the sixty-six energy labels, residual and preparation
rejection. This gives k=69 reported labels. For L independent normalized
ideal repetitions, Cauchy--Schwarz and the multinomial variance identity give

    E[d_TV(p_hat,p)^2]<=k/(4L),
    Pr(d_TV(p_hat,p)>=s)<=k/(4L s^2).                                 (12)

Thus L>=k/(4p s^2) suffices for error s and failure probability p. This is a
bound on the coarse histogram, not on its entire fine archive. It assumes
independent fresh trials. The exact ideal instrument has normalized trial
semantics; a normalized physical device realizing the rounded arrays is a
separate apparatus claim. No displayed number of repetitions is asserted
to have been run.
