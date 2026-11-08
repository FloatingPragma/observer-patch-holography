# Native boundary Bell resources and their charge boundary

Fix two adjacent cells of the literal fair binary federation. Let a,b be
their actual boundary-seam endpoints, with one distinct spare port u,v
inside each respective cell. Each logical qubit uses one occupied rail:

    |0>_A=|0_a 1_u>, |1>_A=|1_a 0_u>,
    |0>_B=|0_b 1_v>, |1>_B=|1_b 0_v>.

Other occupied ports are spectators. The two active rail pairs contain
two native occupations in total. The declared local occupation-wire service
implements the existing 22-pulse beam, or its inverse, to prepare |+> on
each known one-occupation pair. No cross-cell coherent word is assumed.

## One exact native attempt produces the Bell flag

The native E partition retains incoming contrast d=n_a-n_b, together with
coin c and current j. For each coin, raw address zero has d=0 and j=0.
Its actual Kraus operator is

    K_(c,raw0)=P_even/sqrt(2),
    P_even=|00><00|+|11><11|.                                (1)

Both even endpoint states are unchanged and have the same raw record and
phase. In particular their coherence |00><11| is preserved. Both odd
endpoint inputs are excluded. This is an exact one-attempt parity filter;
there is no limiting repetition or exponentially small odd leakage.
C's no-current record has a different projector and admits one odd
orientation for each coin, so it cannot replace E in (1).

Apply (1) to the product of the two local |+> states. Each retained coin
branch has probability 1/4 and normalized logical output

    |Phi+>=(|00>_logical+|11>_logical)/sqrt(2).

The total Bell-success mass is exactly 1/2. Four other nonzero E branches
each have mass 1/8 and a fully known basis output. Every hold, C, reserved
and zero E address remains in the complete 25-slot native instrument.
The code executes that original instrument on an actual federation seam
and also on the full encoded/reference Choi operator. It checks (1)
including arbitrary encoded input/reference coherences. This establishes
the exact ideal coupled operation and its complete native output; no
successful flag is silently renormalized into an unconditional resource.

For arbitrary source input, the existing eleven-E tree may be run in each
cell first. On each nonconstant branch it supplies the terminal basis
configuration. If the boundary endpoint has occupation zero, some spare
port has occupation one, and conversely. Thus a suitable rail pair exists
in every nonconstant cell. The controller chooses its beam orientation
from that committed configuration. Constant-cell branches retain their
state and an inconclusive label. If a is the initial mass of branches in
which both cells are nonconstant, unconditional Bell mass is a/2. The
construction does not assume a positive capacity-independent lower bound
on a. The tests use a declared fixed pair with one occupation in each;
the all-sector tree attachment is the owning preparation theorem.

## Bounded retries use native restoration, without reset

The four rejected branches of a Bell attempt are known basis states. With
total mass 1/4 the boundary current is zero and each rail pair retains one
occupation. With total mass 1/4 it moved one occupation, giving the two
rail pairs totals zero and two. In the latter case repeat the same native
E seam until its recorded current reverses the move, or until r attempts.
The first successful reversal occurs with probability 1/2 per trial.
Every unsuccessful trial has zero current and leaves that basis state
unchanged. Thus restoration times out with conditional probability 2^-r.
The timeout retains the actual terminal state and its complete history.

Once restored, the pair is a known basis state. Two local beams prepare
the next product |++> resource. This uses no state reset, undisclosed input
inspection or fresh source occupation. For at most R Bell rounds put

    c=2^(-r-2),       f=1/2-c.

The complete conditional law given two eligible cells is

    p_Bell=(1/2) sum_(j=0..R-1) f^j,
    p_restoration_timeout=c sum_(j=0..R-2) f^j,
    p_round_exhaustion=(1/2) f^(R-1).                         (2)

These masses sum exactly to one. There is no restoration after a final
unsuccessful Bell round. All finite coin/raw histories behind these coarse
labels remain part of the native instrument.

Worst-case resource counts are 22 initial tree attempts when needed,
R+(R-1)r boundary native attempts, 44R local proper-code pulses and twelve
reused local occupation wires per cell. FIFO capacity one with
`export_oldest` retains all native records. Round, restoration, native
attempt and source-pulse counters are separate. The local beam words are
individually admitted source services; the record-controlled composition
is a declared bounded external controller, not an unimplemented adaptive
node inserted into the native grammar. Expanded execution limits continue
to apply. No pulse-to-native-attempt duration conversion is inferred.

The local word owner supplies complete primitive, first-failure and gain
bounds. Its ideal beam preserves local number, while rounded stored work
maps can have a small component outside the fixed global-number sector.
The federation interpreter explicitly rejects such input masks. Its carrier
schema has one exact global total, with no direct-sum alternative. The
existing all-occupation 4096-state law is a twelve-port component and does
not amend that growing carrier admission. Therefore the rounded numerical
beam output is not fed to the fixed-q interpreter here. An exact source
return interface or an authorized complete-carrier extension is a distinct
premise. No free number projector, illegal-state clipping or source-class
change is used to conceal it. The ideal coherent-service coupling, supplied
and bound by the occupation-wire construction, is sufficient for (1)-(2).

## Encoded entanglement does not transport native charge coherently

Let N_c denote total occupation in cell c. Suppose an initial state is
block diagonal in the complete vector of cell totals. Every allowed local
number-preserving control commutes with every N_c. For a complete native
branch on a seam from c to d with recorded current j,

    N_c K_h=K_h(N_c-j),     N_d K_h=K_h(N_d+j),
    N_e K_h=K_h N_e for other cells.                         (3)

Hence K_h rho K_h* remains cell-number block diagonal whenever rho is.
This holds on each complete record branch, so adaptive classical selection,
stopping, coarse-graining and retaining or forgetting records preserve it.
The Bell state above has a fixed total in each cell and obeys this
invariant. Number-preserving local Bell measurements and logical
teleportation, if supplied, also preserve it.

The actual bare transfer reference P/TrP violates this invariant. Each
crossing exchange connects two configurations x,y with different cell
totals, and

    <x|rho|y>=1/(2M TrP) > 0.

The norm-one observable |x><y|+|y><x| has expectation 1/(M TrP) in rho and
zero in every reachable cell-number block-diagonal state. Their half trace
distance is therefore at least 1/(2M TrP). This excludes exact preparation
by the specified service grammar. The lower bound can shrink with capacity;
it is not a uniform obstruction to all local limiting observations.

The uniform global-q Dicke ground has an even stronger separation. For
q=2 and at least two twelve-port cells, the largest cell-number block has
dimension 144, with one occupation in each of two specified cells. If
D=binom(V,2), the normalized Dicke vector assigns that block mass 144/D.
For any normalized cell-number block-diagonal sigma,

    <Omega|sigma|Omega> <=144/D,
    (1/2)||sigma-|Omega><Omega|||_1 >=1-144/D.                (4)

To prove the first inequality, decompose Omega into its mutually
orthogonal cell-number blocks. Each rank-one compressed projector has norm
equal to its block mass, bounded by 144/D; sum against the positive blocks
of sigma. The second inequality measures the ground-state projector.
At V=240 the fidelity ceiling is 6/1195 and the distance lower bound is
1189/1195. The result concerns this local-service/committed-current grammar,
not every admitted coherent coupling or every alternative state.

Thus the exact Bell construction supplies logical entanglement across an
actual source seam. It neither derives a bare coherent charge-transfer gate
nor prepares the growing transfer reference or Dicke vacuum. Local Bell
analyzers, protected storage, delivery schedules and their complete costs
require explicit constructions before a universal logical network is
claimed. No three-dimensional physical geometry follows from connectivity.
