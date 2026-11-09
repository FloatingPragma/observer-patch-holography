# Complete odd-probe energy filters

## 1. Actual source and the available odd response

Use the [primary construction](../particle_common_contract/PROOF.md), with
its unchanged finite Wilson/Gauss Hamiltonian \(H_n\), selected stationary
polynomial density \(\rho_n\), idle conjugate reference and source clock.
Write \(a_n=2^{-n}\) and \(K_n=H_n-E_{0,n}\ge0\). Ground subtraction
appears only in bounds. The
[selected odd creation proof](../particle_primary_selected_odd_response/PROOF.md)
provides the local contractions

\[
 A_{x,\alpha b}=c^\dagger_{x,L_1,\alpha}
 c^\dagger_{x,L_2,\alpha}c^\dagger_{x,e,b},\quad\alpha,b=1,2,
 \qquad [N_q,A]=\Delta_q A,\quad\Delta=(0,2,0,0,1).
\]

Each is gauge invariant and fermion odd. Its loaded odd apparatus memory
supplies the opposite native increment. The complete creation instrument
has success A and failure \(D=I-A^\dagger A\), with a physical even
unitary extension. Both outcomes, the apparatus and the original reference
are retained. Average the four spin choices and all cells through a
retained classical address. No source input is replaced by a trial.

Its original energy-transfer measure \(\mu_n\) has mass at most one and

\[
 \int(\omega-9/a_n)^2d\mu_n\le9/a_n^2,\qquad
 \mu_n([1/(2a_n),35/(2a_n)])\ge1/50.                \tag{1}
\]

At member zero the second-moment cap is \(9/2\) and
\(\mu_0([3,15])>23/32\). These are actual interacting source bounds;
they do not isolate a spectral line. The pointwise full-CAR residual
bound is \(\|[H_n,A]-9A/a_n\|\le3/a_n\) on every member, including
the boundary. The sharp member-zero value is also used by the producer.

## 2. Complete filtering has an exact selector sum rule

Choose times \(t_j\), normalized selector amplitudes \(v_j\), a unitary
selector matrix R whose first column is v, and scalar phases \(z_j\) of
unit modulus. Conjugate the complete loaded-memory unitary by the same
source clock at each selected time, multiply by \(z_j\), and uncompute
the selector. If p labels memory success or failure, its Kraus operators
are

\[
 K_{r,p}=\sum_j\overline{R_{jr}}v_jz_j\alpha_{t_j}(M_p),
 \qquad M_0=A,\quad M_1=D.                              \tag{2}
\]

This is one finite complete instrument, since every controlled block is
unitary. All selector values r and both memory values p are retained.
On a transition from energy \(E_i\) to \(E_f\), the creation amplitude
is multiplied by
\(\sum_j\overline{R_{jr}}v_jz_je^{it_j(E_f-E_i)}\).
Unitarity of R gives the exact Parseval identity

\[
 \sum_r\left|\sum_j\overline{R_{jr}}v_jz_j
                   e^{it_j(E_f-E_i)}\right|^2=1,
 \qquad \sum_r\mu_{K_{r,0}}=\mu_A.                      \tag{3}
\]

The same statement holds separately for D. It holds with exact finite
Cayley source clocks: each clock is a unitary function of the original H,
so the transition-dependent factors retain modulus one even when they are
not exponentials of a linear energy difference. Scalar Cayley analyzer
phases also retain modulus one. Thus (3) is exact for the specified finite
algebraic implementation. Selected source stationarity is sufficient;
joint stationarity of its retained purification is not assumed.

Both energy clocks may be appended. Their effects are functions of H, so
the same identity holds pointwise in the two clock labels on this
stationary marginal. Summing the selector preserves the original complete
creation/energy law. The coherent selector and source/reference output
retain additional information; arbitrary phase-sensitive quantum outputs
are not asserted identical to the unfiltered instrument.

A successful low-pass blank can be a useful observable. Equation (3)
does not forbid one: it proves that discarded blank weight is routed to
other retained selector outcomes, with unchanged source transition
energies. No probability is divided by the blank success probability.
A different complete dilation of the same blank effect need not have
the selector-summed law (3).

## 3. A populated finite filter with explicit implementation errors

Use the positive triangle \(f_T(t)=(1-|t|/T)_+/T\) and set

\[
 T_n=a_n/6,\qquad\nu_n=9/a_n,\qquad
 B_n=\int f_{T_n}(t)e^{-i\nu_nt}\alpha_t(A)\,dt.           \tag{4}
\]

The frequency \(\nu_n\) selects the source-derived onsite transfer
scale. It is an analyzer phase, with no subtraction from H or the recorded
energy. The original measure of B is

\[
 d\mu_B(\omega)=\operatorname{sinc}^4
             ((\omega-\nu_n)T_n/2)\,d\mu_A(\omega).
                                                               \tag{5}
\]

On the band in (1), the sinc argument is at most \(17/24\). Since
\(\sin x/x\ge1-x^2/6\) for \(|x|\le1\), the original blank weight
there is at least

\[
             \frac1{50}\left(\frac{3167}{3456}\right)^4.
                                                               \tag{6}
\]

Use M=1024 equal cells on \([-T_n,T_n]\), midpoint nodes
\(t_j=-T_n+(2j+1)T_n/M\), and exact positive weights
\(w_j=(2/M)(1-|t_j|/T_n)\). Their sum is one. The derivative of the
demodulated integrand in (4) has norm at most \(3/a_n\) by the actual
centered residual. Hence the blank quadrature error is at most

\[
        \delta_q\le(3/a_n)T_n/M=1/2048.                  \tag{7}
\]

No large-source Hamiltonian norm enters this quadrature estimate. It
does enter implementation of the source clocks. For Hermitian H of norm
at most W, use the same finite algebraic Cayley unitary
\(C_r(t;H)=[(I-itH/(2r))(I+itH/(2r))^{-1}]^r\).
The inherited bound is
\(\|C_r(t;H)-e^{-itH}\|\le|t|^3W^3/(12r^2)\).
The producer chooses integer subdivisions so each source echo error is
at most \(1/65536\), and the scalar analyzer error at most \(1/32768\).
There are two echoes. Thus the exact-node unitary error is at most
\(1/16384\), independent of the number of selector labels, and the blank
operator error to (4) is at most \(9/16384\). Its original event error
is at most twice that number.

The complete two-energy-clock plan uses resolution \(1/(8a_n)\), with
each one-clock tail and implementation budget \(1/16384\). Its combined
event error is at most \(1/4096\). Subtracting these errors from (6)
proves original probability at least \(1/80\) for the joint event
consisting of memory success, blank selector, and reported transfer in

\[
                         [1/(4a_n),71/(4a_n)].          \tag{8}
\]

At member zero, the narrower exact band \([3,15]\) has sinc factor at
least \((23/24)^4\). Its reported band \([11/4,61/4]\) therefore has
original blank probability greater than \(3/5\), including the same
complete error budgets. The receipt retains source-specific sharper
bounds. There are \(8VM\) complete filter records before the two clocks.
All failures, scalar phases and finite subdivisions are specified.
The large-source histogram and the compiled native tape are not executed.

This is a populated finite-capacity spectral readout. Its frequency band
scales as \(1/a_n\); demodulation does not turn that band into finite
physical mass or prove a low-energy particle branch.

## 4. What an actual zero-centered filter can certify

The same moment, without an assumed gap, gives for \(a_nR<9\)

\[
 \mu_n([-R,R])\le\frac9{(9-a_nR)^2}.                    \tag{9}
\]

A zero-centered triangle of fixed half-width T has multiplier
\(\operatorname{sinc}^4(\omega T/2)\le\min(1,16/(\omega T)^4)\).
Split its success measure at
\(R_n=\max(4,2^{\lfloor n/2\rfloor})\). The original blank probability obeys

\[
 p_{n,T}\le\min\left(1,
       \frac9{(9-a_nR_n)^2}+\frac{16}{T^4R_n^4}\right),
 \qquad \limsup_n p_{n,T}\le1/9.                        \tag{10}
\]

At T=1 the actual member-zero moment improves the bound to
\(p_{0,1}\le97/400\). For every \(n\ge1\), \(a_nR_n\le2\) and
\(R_n\ge4\), so \(p_{n,1}\le9/49+1/16<1/4\).
Thus this zero-centered fixed-width blank has original probability below
one quarter on every member.

The contribution from the certified escaping band alone is suppressed
by at most \(\min(1,(4a_n/T)^4)\). Neither estimate supplies a positive
lower bound on the blank's fixed-energy weight. A lower-frequency
component is permitted. Equations (3) and (10) distinguish the complete
filter's retained energy records from the possibly useful filtered field.

## 5. Number-phase conjugation and growing dressing depth

For any static real number profile, let
\(Z=\exp(i\sum_y\theta_y(N_y-30I))\). On the actual source,

\[
                 ZA_{x,\alpha b}Z^\dagger
                    =e^{3i\theta_x}A_{x,\alpha b}.       \tag{11}
\]

Number phases commute with the cutoff and Gauss projections. Arbitrarily
many such phases at one source time combine into another Z. Every local
success effect and the complementary failure are unchanged. Thus the
complete classical-address creation instrument, including its source and
reference outputs after recording the outcomes, is exactly unchanged.
An unmeasured coherent pointer instead receives a known diagonal phase.
A coherent spatial analyzer can detect relative address phases; its
instrument is not included in this invariance statement. Time-separated
phase operations with interacting source evolution are also a different
operation.

There is a quantitative bound for the latter possibility. Let sigma be
the original unnormalized odd-creation success output, including retained
registers, with \(\operatorname{Tr}\sigma\le1\). Equation (1) implies
\(\operatorname{Tr}\sigma 1_{[L_n,\infty)}(K_n)\ge p\), where
\(L_n=1/(2a_n)\) and \(p=1/50\), because initial excitation is
nonnegative. Any word U containing k actual critical-density Weyls or
their adjoints, with arbitrary source waits, obeys the inverse
[energy-form bound](../particle_primary_current_tail/PROOF.md)

\[
 \|K_n^{1/2}U^\dagger v\|\le\|K_n^{1/2}v\|+D\|v\|,
 \qquad D=45k/2.
\]

For \(P=1_{[0,R]}(K_n)\), \(Q=1_{[L_n,\infty)}(K_n)\), this yields
\(\|QU^\dagger P\|\le(\sqrt R+D)/\sqrt{L_n}\).
Purify sigma as a vector of norm at most one, split its transformed vector
with P and I-P, and apply the triangle inequality to its original Q part.
The original tail probability q consequently satisfies

\[
 q:=\operatorname{Tr}(U\sigma U^\dagger1_{(R,\infty)}(K_n))
 \ge[\sqrt p-(\sqrt R+D)/\sqrt{L_n}]_+^2.                 \tag{12}
\]

All projectors act trivially on the retained registers. The same estimate
holds for address-controlled words with the same maximum k. To obtain
\(q\le\eta<p\) therefore requires

\[
 k\ge\frac2{45}
 [\sqrt{L_n}(\sqrt p-\sqrt\eta)-\sqrt R]_+ .             \tag{13}
\]

In particular, at R=1 and \(\eta=1/200\),
\(k\ge[(a_n^{-1/2}-20)/450]_+\). The certified necessary integers
are at least three at level20 and 2331 at level40. More generally,
\(k_n=o(a_n^{-1/2})\) retains limiting tail at least \(1/50\) above
every fixed R, so it cannot make these full outputs energy tight.
This is a necessary count of the specified actual Weyl operations,
not a lower bound for arbitrary native gates, a sufficiency theorem for
large k, or an exclusion of another dressed low-energy observable.

The current source constructions supply no positive original fixed-energy
odd residue at the interacting endpoint. Obtaining one requires control
of that channel's source spectral projector and matrix elements, or an
equivalent source correlation bound. A small-coupling series, an auxiliary
prepared vector or the populated moving band (8) does not supply that
estimate. Physical dispersion and a GeV rest-mass identification are
separate from the finite filter proved here.
