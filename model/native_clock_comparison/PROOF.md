# Comparing the working clock with record transfer

## Source and proposed identification

The source is the biased native H66 member in
[`native_record_response`](../native_record_response/PROOF.md): two binary
occupations on twelve ports, weights (25,25,8,...,8), thirty uniformly selected
seams and the complete C/E instrument. Let P be its column-stochastic population
transition, pi its faithful stationary distribution, and Pi=diag(pi). The
working reference is Pi. Its modular clock is K=-log(Pi).

Correctly reflected scalar record histories give a different Hilbert
representation, with

    S = Pi^(-1/2) P Pi^(1/2),    H = -log(S),    Omega = sqrt(pi).

The positive transfer theorem proves (296/495)I <= S <= I. Its ground direction
is Omega; it does not prepare the working mixed density Pi. The following
comparison tests the specific proposal that the record representation is
identified with the working carrier in its existing configuration basis, with
the same state Pi and accessible measurements, and that exp(-itH) supplies its
clock. This is a hypothesis under test, not an executed source instruction.
Neither a physical time calibration nor a unitary continuation protocol is
inferred from positive transfer.

## An actual measurement separates the two proposals

Use the complete Y instrument on the native probe (1,7), as constructed in
[`native_clock_observables`](../native_clock_observables/PROOF.md). For every
spectator k, let i={1,k} and j={7,k}. There are ten pairs. The Hermitian
observable is

    Y = sum_k (-i |i><j| + i |j><i|).

It has eigenvalues -1,0,1 on the full carrier; the zero eigenspace has dimension
46. Its native implementation uses the declared coherent analyzer followed by
the original E repair. The score is its incoming contrast, including all coin
and raw records. Every abort is assigned score zero. No odd-sector conditioning
or removal of the zero outcomes occurs.

The modular proposal leaves Pi fixed for every real parameter s:

    exp(-isK) Pi exp(isK) = Pi,    Tr(Y Pi) = 0.                  (1)

For the proposed transfer continuation define

    f(t) = Tr[Y exp(-itH) Pi exp(itH)].                          (2)

The matrices H and Pi are real; Y is purely imaginary. Complex conjugation
therefore gives f(-t)=-f(t), so f(0)=f''(0)=0. Differentiation gives

    f'(0) = -2 sum_(i,j) (pi_i-pi_j) H_ij.                     (3)

The sum runs over the ten oriented probe pairs, not all ordered matrix
indices. This sign uses exp(-itH). Reversing the candidate clock orientation reverses
the mean and preserves the nonzero separation.

All probe pairs have pi_i/pi_j=25/8. Define the rational polynomial

    F_m(P) = sum_(n=1)^m (I-P)^n/n,
    r = 199/495,    R_m = r^(m+1)/[(m+1)(1-r)].

The existing transfer theorem bounds ||H-F_m(S)|| by R_m. Similarity gives
F_m(S)_ij=sqrt(pi_j/pi_i) F_m(P)_ij. Consequently the polynomial value of
(3) is sqrt(2) times the exact rational

    q_m = -(4/5) sum_(i,j) (pi_i-pi_j) F_m(P)_ij.

Its derivative error is at most 2R_m, by
|Tr(Pi [Y,H-F_m(S)])| <= 2||Y|| ||H-F_m(S)||. For m=10 the executable rational
calculation gives

    q_10 = 6842377389660770935808120212145413
           / 1351444261278934581551654400000000000.

Bounding sqrt(2) between 1414213562373095/10^15 and
1414213562373096/10^15, with both squared inequalities checked exactly, proves

    f'(0) > 7/1000.                                            (4)

The scalar series with the same remainder proves
||H|| <= -log(296/495) < 103/200. The third derivative of (2) is bounded by
(2||H||)^3, since Pi is a density and ||Y||=1. Taylor's theorem, including
the exact vanishing second derivative, thus gives

    f(1/10) > (7/1000)(1/10) - (4/3)(103/200)^3(1/10)^3
             = 3107273/6000000000 > 1/2000.                    (5)

This is a finite-time bound, not an inference from a nonzero infinitesimal
commutator. An arbitrary rescaling or offset of the native modular parameter
does not change (1) on this reference.

## Complete implementation allowance

Let M be the ideal complete Y instrument and Mtilde its specified numerical
implementation. The owning receipt retains the full 512-address arrays, all 25 native
raw slots, address-export aborts and the bound on every first proper-code
failure. Its exact uniform coefficients imply

    ||Mtilde-M||_diamond/2 < delta = 1/40000000,
    ||Mtilde||_diamond < G = 20000001/20000000.                 (6)

The comparison checks the receipt's source bindings and independently
recomputes these inequalities from its fixed returned-array acceptance gate,
stored-primitive multiplication allowance and first-failure bound. No
platform-dependent residual is serialized or polar normalization introduced.
The owning compilation proof and acceptance checks are the evidence for the
implemented analyzer; no repeated full compilation is asserted by this
comparison's receipt.

The native preparation theorem, after 1600 unread attempts, bounds the reduced
source/reference half trace error by epsilon=542159/10^12. It compares the
source to Pi with the corresponding untouched reference state. This statement
traces the preparation archive. Conditioning the later protocol on that
archive is outside its domain.

For a candidate clock service C whose complete map is within half diamond
error eta of the relevant ideal unitary, composition bounds the output error
by

    E = G(epsilon+eta)+delta.                                  (7)

If clock failures are present, retain their labels, terminate those branches
and give them score zero; measurement acts only on the success branch. Its
extension by identity on terminal failure blocks has norm at most G, so the
same bound holds. The compared ideal maps have zero failure blocks. The
clock-service allowance is a conditional acceptance budget; no implementation
of H, or arbitrary-parameter compilation of K, is supplied here.

The score operator on the complete output has norm at most one. Equation (5)
therefore gives an ideal complete-output half trace separation greater than
1/4000. If both candidate services have eta<=1/100000, the triangle inequality
gives

    D(output_H, output_K) > 1/4000 - 2E > 1/5000.               (8)

An arbitrary untouched reference is allowed. All actual numerical blocks
remain unnormalized. Their complete classical-quantum half trace separation
is the quantity in (8); it is not called total variation of normalized
probability distributions without an exact channel realization. Dropping
the terminal working carrier keeps the ideal score witness, so the same
lower-bound argument applies to the complete records/reference output.

Thus a proposed service accurate to the stated budget cannot realize both
clock identifications with this reference and readout. This acceptance test
does not constitute an experimental detection of a transfer unitary that the
source never executes.

## Reference-preserving encodings

The derivative (4) implies [H,Pi] != 0. Suppose a unitary encoding W preserves
the specified reference, W Pi W†=Pi, and identifies H with a rescaled modular
clock up to an energy offset: W H W†=aK+bI. The right-hand side commutes with
Pi, which would force [H,Pi]=0. Such an encoding is impossible. This conclusion
needs no assumption that W separately preserves Y.

The same stationary-reference test applies to a proposed dilation. If a joint
reference tau is invariant under a joint clock and its working marginal is
Pi, that marginal stays Pi. It cannot also equal exp(-itH)Pi exp(itH) for all
t under the tested identification. A different reference, a nontrivial
readout/encoding map, or a clock that changes the preparation is a different
proposal and needs its own complete-record comparison. These alternatives
are not excluded by (8).

## What the bare repair records can identify

Every actual complete C/E branch is a scalar partial injection in the
configuration basis. Composing branches and arbitrary diagonal phase
operators preserves that form. For a complete raw history h its map can be
written

    K_h = sqrt(w_h) sum_(i in D_h) exp(i phi_i) |f_h(i)><i|,

where f_h is injective. For any initial source/reference operator rho, tracing
the terminal working carrier gives

    Tr_source[(K_h tensor I) rho (K_h† tensor I)]
       = w_h sum_(i in D_h) rho_(ii).                          (9)

All phases cancel. This holds branch by branch, including zero, hold and
reserved outcomes. It extends to finite adaptive policies using the retained
records and to their FIFO/export/stop records by induction on each prefix.
The working carrier must be traced for (9); a subsequent coherent Y measurement
need not be phase blind.

Consequently bare repair-record/reference experiments cannot distinguish
interleaved identity phases from arbitrary diagonal candidate clocks. This
does not forbid reading the declared simulator or clock specifications:
service identifiers, source code and clock-specification records can
distinguish those proposals. Equation (9) concerns the repair outcomes and
untouched reference under a common declared record policy. The coherent
analyzer in (5) supplies phase-sensitive measurement data; its primitives are
explicit source inputs, not consequences of the blindness statement.
