# A common operational clock for the interacting finite source

The carrier is the complete 384-state all-number modular CAR source in
`next/C10/joint_source.py`. Its supplied Hamiltonian is the same interacting
matrix with mass 0.23, electric coefficient 0.37 and face coefficient 0.61.
The preparation is its actual global ground reference followed by all four
branches of its fair-0/1 charged reservoir instrument. The reservoir labels
and their subnormalized input and final density matrices are retained.
No fixed-number or field-free reference is substituted.

The construction below defines a finite quantum process with a common
operational parameter. Hamiltonian duration, the event rate and the service
schedule are declared apparatus inputs. The result supplies no equality
between a source duration, a native attempt counter and physical proper time.

## Complete marked instrument

Let `K_r` be the actual C10 native Kraus operator, including the uniform
seam probability 1/3. The MIX law has thirty marks, each containing seam,
law, coin and record value. The underlying integer repair gives its current
`j_r`. For C marks the record is the current. For E marks it is the contrast
`d`, with `j_r=[d-(2 coin-1)|d|]/2`. The complete operators obey

    sum_r K_r* K_r = I.

Each mark costs one native commit pulse and one native attempt, including
`j_r=0`. A zero-current Kraus operator can change coherence through its
measurement record. It is not an identity controller idle. The exact finite
carrier has no rejected-boundary outcome; every admitted basis event has
a target. No unimplemented reset or detector success branch is assumed.

Choose a native rate `lambda>0`, an auxiliary uniformization rate
`Lambda>=lambda` and a duration `T>=0`. A Poisson process of rate Lambda
supplies ordered times `0<t1<...<tm<T`. Between consecutive times apply the
actual coherent service `U(delta)=exp(-i H delta)`. At each time apply the
complete instrument

    M_idle = sqrt(1-lambda/Lambda) I,
    M_r = sqrt(lambda/Lambda) K_r.

For a complete word `b=(b1,...,bm)`, put

    A_b = U(T-tm) M_bm U(tm-t_(m-1)) ... M_b1 U(t1).

Its density-valued instrument is

    exp(-Lambda T) Lambda^m A_b rho A_b* dt1 ... dtm.

This formula retains every mark, including controller idles and zero-current
attempts. For each fixed time sequence, summing all words preserves trace
because every M instrument is complete. Integrating the simplex gives the
Poisson factor `exp(-Lambda T)(Lambda T)^m/m!`; summing over m gives one.
Consequently the construction is a complete positive instrument on the
same preparation, without a success-conditioned source replacement.

The executable `trajectory()` implements this formula for a specified word.
Its conditional state normalizations are ordinary measurement conditioning;
the product of every selected probability is retained in the full path
density. Its complete alternatives at each visited node sum to one. The
complete output is calculated by the generator below, rather than inferred
from that single illustrative word.

## Exact generator and common resource relation

Define `J(rho)=sum_r K_r rho K_r*`. The infinitesimal time-ordered expansion
gives the exact unmarked generator

    L(rho) = -i[H,rho] + lambda (J(rho)-rho).

It is independent of Lambda. This equality requires no commutation between
the coherent and native maps. Replacing the random waits by their mean,
or using a finite product of the two channels, generally changes the map.
The independent reference checks this distinction on the actual preparation.

Let N count native attempts, M count nonzero native currents, and I count
controller idles. Write `J=J_zero+J_move`. The density generating function is

    G(z,x,y;T)=E[z^N x^M y^I rho_record]
             = exp(T L(z,x,y)) rho,

    L(z,x,y)=-i ad_H + lambda z (J_zero+x J_move)
             - lambda Id + (Lambda-lambda)(y-1) Id.

For `x=1`, tracing this equation closes exactly:

    Tr G(z,1,y;T)
       = exp(T[lambda(z-1)+(Lambda-lambda)(y-1)]) Tr rho.

Thus N and I are independent Poisson counts with means `lambda T` and
`(Lambda-lambda)T`. This count independence does not make the source state
independent of N or its native marks. M is state-dependent, and the code
calculates its marked density instead of substituting N for M.

For each retained history, the source-service durations sum exactly to T.
Under the declared segmentation the resource identities are

    native commit pulses = N,
    coherent interval services = N+I+1.

The means are `lambda T` and `Lambda T+1`. Combining adjacent coherent
intervals separated only by an identity idle reduces the service-call count
to `N+1` while preserving the source output. Thus this count is a specified
implementation cost. A native pulse implementation of each U(delta), its
cost and a physical pulse duration are not supplied. In particular, native
commit pulses here do not mean the coherent M6 source pulses of the separate
192-CAR rest-clock construction.

The source-energy expectation obeys

    d Tr(H rho_t)/dt = lambda Tr(H[J(rho_t)-rho_t]).

The Hamiltonian intervals contribute zero to this derivative. The native
instrument can exchange energy with its measurement apparatus. No energy
conservation or physical work cost for that apparatus is inferred.

## Complete infinite-count complement

The count-zero and count-one density maps are exactly

    rho_0(T)=exp(-lambda T) U(T) rho U(T)*,
    rho_1(T)=lambda exp(-lambda T)
             integral_0^T U(T-s) J(U(s)rho U(s)*) U(T-s)* ds.

The retained density `rho_ge2=exp(TL)rho-rho_0-rho_1` is the sum of all
instrument paths with at least two native attempts and is positive. Its
trace is `1-exp(-lambda T)(1+lambda T)` for normalized input. This class is
not a halted simulation with a finite attempt cap. It is an exact coarse
outcome of the unlimited process; all its source-state weight is retained.
If a consumer uses only the zero/one-attempt outputs, this class is its
explicit inconclusive complement and cannot be discarded or renormalized.
The receipt stores all twelve reservoir-label/count-class joint states.

The calculation uses number-sector compression because the prepared input,
H and every K preserve number blocks. The implementation rejects input
with cross-number coherence for this solver. Its `evolve()`, `native()` and
`trajectory()` APIs act on the full 384-dimensional space, including those
coherences. The instrument theorem itself holds on the full space.

## Falsifying clock interventions

For `T=0.37`, `lambda=0.8`, `Lambda=1.6`, the exact count probabilities are
approximately 0.7437874280, 0.2201610787 and 0.0360514933 for zero, one and
at least two native attempts. The zero-moving-attempt probability is
0.9304557347. Counting only nonzero current loses a substantial part of
the native pulse budget.

On the positive-probability zero-attempt branch, the conditional state is
`U(T)rho U(T)*`. It changes a terminal actual native-mark probability by
0.0019959540 relative to the same preparation at duration zero. Therefore
a deterministic identity `T=c N` cannot represent all retained histories,
and replacing the zero-attempt branch by an unchanged source is falsified
by the actual C10 native readout.

Reducing only the native rate by a factor a gives

    L_(rate only) = -i ad_H + a lambda (J-Id).

Changing a common time parameter by a gives `a L`, which also multiplies
the coherent commutator by a. They agree universally only under the
appropriate vanishing commutator condition. The executed half-rate control
differs from common half-time evolution by 0.0963135872 in Frobenius norm.
Scaling both H and lambda by a gives the exact identity
`L_(aH,a lambda)=a L_(H,lambda)`, which the execution verifies.

Adding recorded controller idles through any Lambda greater than lambda
leaves the source output unchanged. This is the precise uniformization
identity. It does not contradict the finite-attempt modular rescaling
obstruction in `routes/clock_transfer/native_clock_transfer.py`.

## Source transfer boundary and reproducibility

The passive two-time rest-clock theorem uses a distinct longitudinal
192-CAR source and explicitly supplies its duration. The routed-flux join
uses a different encoded carrier, routing tick and supplied weak
Hamiltonian service. Neither theorem identifies its clock with the C10
attempt counter. This adapter consumes their stated scope without importing
their rest frequencies, continuum bounds or source-pulse compiler.

Run the producer and its separate numerical reference with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 common_clock.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 reference_review.py

`COMMON_CLOCK.json` pins the source and retained density archive.
`REFERENCE_REVIEW.json` reconstructs all 2,304 native event columns from
integer repair and direct CAR parity, then compares a direct matrix ODE
and an ordered-time integral with the sparse counting generator. It rejects
no-op deletion, complement deletion or renormalization, and a mean-wait
substitution. It is a separate implementation by the same author and does
not count as an external scientific audit. No acceptance configuration or
official OpenWave row is changed.
