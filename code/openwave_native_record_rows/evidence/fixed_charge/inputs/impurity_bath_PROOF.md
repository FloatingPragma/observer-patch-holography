# Interacting native impurity and neutral-rung bath

## Full source and its actual reference

Add a binary impurity port to the previous rung source. In addition to N
internal-rung macro addresses and N−1 whole-rung transports, add two native
addresses connecting the impurity to the two endpoints of core rung 0.
Both use the same passive counter rule z ^= coin XOR n_first(before).
Every primitive call, cap failure and partial swap remains. There are
W=2N+1 equally weighted macro addresses, with no supplied coupling strength.
The auxiliary zero pair is an invariant preparation, not a hidden reset.
All 2^(2N+1) load configurations and both counter values remain.

Every internal/impurity transfer is (I+Swap tensor X)/2. Each capped transport
is the complete positive mixture in the copied rung-bath proof. Hence the
full transfer is positive, selfadjoint and stochastic at every finite cap.
In the counter-odd sector, in per-address Poisson units H=W(I−P_odd), its
ideal limit is

  H = sum_r Pi_triplet,r + Pi_triplet,(s,a0) + Pi_triplet,(s,b0)
          + (1/2)sum_r(I−R_(r,r+1)).

Here Pi_triplet,ij=(I+Swap_ij)/2=Si·Sj+3/4. Thus the two impurity terms are
S_s·T_0+3/2, where T_0=S_a0+S_b0. The constant is 3/2, not 1/2. It is fixed
by the complete native counts and cannot be removed without naming the
reference.

The three bonds of the core triangle sum to
(1/2)S_core²+9/8. Their minimum is 3/2 in total spin 1/2; their spin-3/2
value is 3. All remaining rung and transport terms are positive. The state
of an arbitrary impurity spin times singlets on every rung attains 3/2:
T_0 annihilates a rung singlet and every whole-rung swap fixes this reference.
Therefore 3/2 is the exact full ideal-source ground energy, with the stated
impurity doublet among its ground states. We reconstruct excitation energy
K=H−3I/2 with respect to this same-source reference. The finite capped
source is kept separate; its partial swaps need not fix that reference.

## Full one-triplet sectors, charges and products

The spin labels below are the exact internal exchange-SU(2) representation.
They are not, by themselves, an attachment to emergent spacetime spin.

The ideal source preserves the number of triplet rungs. Whole-rung swaps
transport their positions; S_s·T_0 preserves the core rung's singlet/triplet
Casimir. With one triplet and the impurity spin, the entire sector separates
into total spin J=3/2 and J=1/2. On a path of N rungs the restrictions are

  K_(3/2) = I + L_path/2 + (1/2)|0><0|,
  K_(1/2) = I + L_path/2 − |0><0|.

This follows from S_s·T = [J(J+1)−3/4−2]/2, equal to 1/2 or −1. It is
proved as a reducing-subspace identity in the full source matrix, before
examining any particular product channel. Other total-load sectors and
multi-triplet states remain in the full source. Finite cap leakage is
reported rather than discarded.

Use the explicit initial core state with J=3/2, M=1/2:

  h = sqrt(2/3)|up,t0> + sqrt(1/3)|down,t+>
    = (|110>+|101>+|011>)/sqrt(3),

where the first bit is the impurity and the remaining two are rung 0.
All other rungs start in singlets. This is the higher core eigenmode of
its actual isolated native triangle. A lower core mode is |up,s>, with
the same core load 2, and core energy 3/2. Define core relative charge as
n_s+n_a0+n_b0−1 and bath-rung charge as n_a+n_b−1. The complete initial
core charge is +1. When the triplet is outside the core, the complete
product outcomes are

  probability weight 2/3: |up,s> at core, outgoing t0 (bath charge 0);
  probability weight 1/3: |down,s> at core, outgoing t+ (bath charge +1).

The second core then has charge 0, and total relative charge remains +1.
Both channels are kept. The outgoing neutral mode is individually neutral,
not an oppositely charged pair grouped together. A source-defined same-charge
lower core plus neutral ejectum therefore occurs with asymptotic weight
2/3 on this branch. The charged alternative remains 1/3; the packet does
not claim the full physical muon channel or discard unwanted charge outcomes.

The other prepared core mode in the same total charge sector is
b=(|110>+|101>−2|011>)/sqrt(6), with J=1/2. Its product weights are reversed:
1/3 neutral same-charge core and 2/3 charged alternative. Its bound part
below is essential, not an omitted preparation failure.

These vectors are explicit reconstructed source functionals. The packet
executes the full native event census and reconstructs their probabilities;
it does not claim that coherent vectors were prepared in classical native
working loads. Their charge and energy assignments follow from the same
complete source and declared primitive charge operators.

## Radiating branch and retained bound branch

On the half-line, K_(3/2) has diagonal 2 everywhere and offdiagonal −1/2.
The unitary sine transform has generalized eigenfunctions
sqrt(2/pi) sin[k(r+1)] and energy 2−cos k. The state initially at the core
has amplitude

  u_r(t)=exp(−2it) i^r [J_r(t)+J_(r+2)(t)].

This follows by writing 2 sin k sin((r+1)k)=cos(rk)−cos((r+2)k) and
integrating exp(it cos k). The full sine transform is unitary, so total
probability is 1. In particular

  P_core(t)=[J0(t)+J2(t)]²=[2 J1(t)/t]²,

with its continuous value 1 at t=0. Fixed finite regions empty as t grows.
An elementary proof splits the oscillatory integral into intervals of
width delta around 0 and pi and integrates by parts in the complement;
with delta=t^(-1/3), every fixed-r amplitude tends to zero. No exponential
lifetime is assumed. Neutral escaped probability tends to 2/3, and charged
escaped probability tends to 1/3. The reference/lower core remains a
stationary exact ground without outgoing triplets.

For K_(1/2), endpoint diagonal is 1/2 and bulk diagonal 2. There is one
normalizable bound mode

  E_b=1/3, phi_r=sqrt(8/9)(1/3)^r.

Direct substitution proves the eigenvalue and normalization. There are no
additional bound modes: outside the band [1,3] a square-summable bulk solution
is C r^n with |r|<1; the endpoint equation forces r=1/3. At threshold r=±1
is not square summable.

For completeness the endpoint resolvent of K_(1/2), in the convention
G(z)=<0|(z−K)^−1|0>, is

  G(z)=2/[z+1+sqrt((z−2)^2−1)],

where the square root is asymptotic to z−2 at infinity. Its pole has residue
8/9; its continuous density is

  rho(E)=sqrt(1−(E−2)^2)/[pi(3E−1)], 1<E<3.

The denominator is strictly positive on the band. Finite-site matrix elements
have the same integrable density times fixed spectral polynomials. Splitting
an integral into compact interior and small endpoint tails, then integrating
by parts on the interior, proves decay of its continuous component. Thus
endpoint survival tends to (8/9)^2=64/81. Probability in r=0,...,R tends to
(8/9)[1−9^(−R−1)]. Total bound probability is 8/9; the remaining 1/9 escapes.
The source therefore exhibits both a radiating higher mode and a bound mode.

## Finite source, limits, control and costs

Every capped transport differs in norm by at most 2*2^(-B) from its ideal
version, so ||H_N,B−H_N||<=2(N−1)2^(-B). This comparison applies on the full
source; no unproved positivity of a projected cap process is used. At fixed
t, unitary state error is at most 2|t|(N−1)2^(-B) and a norm-one effect error
at most twice that. Reference subtraction is the fixed ideal reference;
finite-cap ground/reference shifts are bounded by the same source error.

For the radiating branch, subtracting 2I gives norm at most 1 in finite and
half-line volume; powers from the core agree through order N−1. Thus the
full-vector finite-volume error is at most 2|t|^N/N!, and any probability
error at most 4|t|^N/N!. For the bound branch the same proof uses norm at
most 2, giving 2(2|t|)^N/N! in norm and 4(2|t|)^N/N! for probabilities.
These bounds are intentionally conservative. Fix a finite observation window,
choose N for volume error, then B for full-source cap error, and only then
take the local-relaxation limit. Fixed finite chains have recurrences.

The changed source disables only transport by executing the same actual
auxiliary native decision and holding. It keeps both impurity edges and
all macro marks. Each prepared core spin eigenmode then stays at position 0:
survival is exactly 1 and no external neutral or charged product escapes.
The planted outside-core mode is detected by the same source projector.

Per-address reconstructed time t is mean Wt macro marks in the total-attempt
Poisson convention. Ideal expected primitive calls are [N+2+3(N−1)]t;
finite transport costs 3−3delta+delta² with delta=2^-B. All random durations,
abort/partial records and chosen preparation data remain. This is an actual
native interaction/bath attachment, but not a physical μ/e/ν identification,
a source-derived generation register, an SI clock or a fitted decay lifetime.
