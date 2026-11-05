# A complete adaptive observation law for the compact interacting source

## Fixed microscopic law and common quantum carrier

Use the same compact opposite-degree texture, four-dimensional spatial band
and 45 internal annihilator slots as
[`compact_interacting_response`](../compact_interacting_response/PROOF.md).
There are 180 localized CAR annihilators. Every native memory record has
the same compact kinetic operator A and the same interaction

```
V=lambda (N_W-90)^2,             lambda=1/128.
```

Its quadratic pairing kernel is the actual record-dependent native family
kernel K_r. The localized and complementary CAR factors reduce every
record Hamiltonian. The memory is the existing classical 12-state process,
with 128 native raw labels and four independent binary-filter labels at
each step. All 512 labels, including zero-probability outcomes, are part of
the complete instrument. Probabilities and memory destinations do not
depend on the field state. The supplied coherent pulses and filter services
are retained as specified apparatus.

The forward pixel-memory design has feedback gain g<1 and reset
contraction p<1/12. For M=64 and half-separation one, the spatial-band
bound is r=7/240. Every record's full localized Fock generator therefore
obeys the inherited analytic bound

```
||H_r,W|| <= M_W =90r+90(1+g)+8100 lambda <= 7869/32.
```

Its mean is the same interacting compact generator Hbar_W, including A,
the mean K_r and V. This does not replace V by its expectation or assume
commutation with the pairing term. The pulse parameter is the original
Hamiltonian time; no physical duration or particle mass calibration enters.

## The complete one-interval channel

Let an interval have duration t>=0 and n>=1 steps, each of duration t/n.
The input is any classical-quantum memory/field state, with any finite
untouched quantum reference. A fixed incoming memory value m and complete
native record string omega have probability q(omega|m), ending memory
m(omega) and conditional field unitary U_omega. The native channel retains
m, omega and m(omega) as classical labels. The ideal channel uses exactly
the same q and classical labels, and replaces U_omega by exp(-it Hbar_W).
It does not replace the memory process by independent stationary samples.

The existing memory averaging proof gives, for every initial memory value
and arbitrary field/reference input, the unhalved unread trace bound

```
delta(t,n)=4 M_W t/[n(1-p)]
          +(8/(1-p)+4) M_W^2 t^2/n.
```

The first term is absent for an independently supplied stationary-product
initial preparation. A stationary memory marginal alone is insufficient.
The bound comes from the reset decomposition T=(1-p)R+pS, its contraction
on zero-sum operator-valued memory blocks, and the bounded-unitary
second-order remainder. The inherited proof is uniform in the size of
an untouched reference and in the compact source's exterior volume.

For clarity, the upgrade from the unread bound to complete records can be
proved without a measurement assumption. Purify the conditional input and
put psi_omega=U_omega psi, psi_bar=exp(-it Hbar_W)psi. The trace bound and
the projection onto psi_bar give

```
sum_omega q(omega|m) [1-|<psi_bar,psi_omega>|^2] <= delta/2.
```

The pure-state trace-distance identity and Jensen's inequality imply

```
sum_omega q(omega|m)
  || |psi_omega><psi_omega|-|psi_bar><psi_bar| ||_1
<= sqrt(2 delta).
```

The left side is precisely the trace norm with the entire classical
archive retained. Convexity extends this to correlated classical memory
and mixed field/reference inputs, retaining the initial memory tag.
Discarding that extra tag is contractive. Thus the complete channels differ
by at most

```
epsilon(t,n)=min(2,sqrt(2 delta(t,n)))
```

uniformly on the admitted classical-memory inputs and their arbitrary
quantum-reference extensions. This is a completely bounded comparison on
that classical-quantum domain; it makes no assertion about coherent
superpositions of classical memory values.

## Adaptive interventions and actual continuing memory

Choose m successive evolution intervals, with duration upper bounds t_j
and step counts n_j. Between intervals apply an arbitrary complete quantum
instrument to the localized field and any auxiliary quantum observer
register. It may depend on every retained native and detector record.
Its classical outcome and all failure outcomes are retained. Conditional
quantum feedback is included in the instrument. The native memory label
continues according to the actual transition law; it is not reset by this
comparison. The same declared instruments and feedback policy are applied
to the ideal mean-dynamics experiment.

The instruments must be defined on all their admitted inputs and form
trace-preserving complete channels after retaining their outcomes. No
postselection is inserted. Durations may be chosen from earlier records
provided they obey the fixed t_j and n_j bounds. The memory transition and
record Hamiltonian family between interventions are the stated source law.
Changing that law or intervening after each individual native operation
requires its own schedule and resulting error budget.

Replace the actual interval channels by their ideal counterparts one at
a time, keeping the rest of the experiment unchanged. At each replacement,
its incoming memory, field, auxiliary observer and earlier classical
archive form an arbitrary allowed classical-quantum state. The
one-interval bound therefore applies even when a preceding detector has
correlated them. A record-dependent complete instrument is a channel on
this joint classical-quantum algebra, so trace distance contracts through
every subsequent instrument and interval. The triangle inequality yields

```
||rho_actual,all_records - rho_mean,all_records||_1
<= min(2, sum_(j=1)^m epsilon(t_j,n_j)).             (1)
```

The archive in (1) includes each native record, detector result, feedback
choice and the final quantum field/reference state. In particular every
joint detector probability distribution has total-variation error at most
half this bound. Every bounded final observable has expectation error at
most its operator norm times this bound. No small-probability conditioned
branch receives an error guarantee independent of its success probability.

This is a multi-time comparison of complete controlled observable laws.
The field algebra, initial-state inclusion, interacting generator, memory
and clock are the same at each interval. The exterior-volume-uniform
constant is inherited from the same localized algebra. Arbitrarily many
external sites are not silently promoted to additional controlled field
observables, and no lattice-spacing or Lorentz limit is taken.

## Explicit precision and intervention-count schedules

Let a requested complete-state accuracy be 0<epsilon<=2. Allocate
epsilon_j=epsilon/m and define

```
B_j=4 M_W t_j/(1-p)+(8/(1-p)+4) M_W^2 t_j^2,
n_j=max(1,ceil(2 B_j/epsilon_j^2)).                (2)
```

Then each squared interval error is at most epsilon_j^2, proving (1)
with error at most epsilon. These are rational integer calculations for
rational input bounds. Only an explicitly prepared first stationary
product may omit its transient term. No later interval acquires that
premise from an earlier unconditioned memory marginal.

For the actual compact bound M_W<=7869/32, p<=1/12, three intervals of
duration 1/3 and epsilon=1/100, (2) requires 46,370,150,349 native steps
in total. This is a finite resource bound, not a claimed executed tape or
an optimized algorithm. The complete archive is defined by its finite
instrument factors; an exponentially large branch array is unnecessary
for specifying its exact law.

If a total of N steps is divided equally among m equal intervals in a
fixed total time T, with N divisible by m, squaring the sum in (1) gives
the upper bound, before capping at four,

```
8 M_W T m^2/[N(1-p)]
 +2(8/(1-p)+4) M_W^2 T^2 m/N.                    (3)
```

Hence m=o(sqrt(N)) is sufficient for this displayed estimate to vanish
at fixed T,M_W,p. Fixed finite intervention schedules are included. The
estimate gives no uniformly accurate limit for an observer applying
unrestricted feedback at every microscopic step.

## Literal native record control and noncommuting detectors

The independent execution uses the original native 128-record C/E law
and independent filter probabilities 3/8 on both binary inputs. The exact
effective event factors are p=1/32 and q=1/16. All 12 old memory inputs
and all 512 fine labels are retained. Exact summation reproduces the
source memory recursion and verifies its reset decomposition and stationary
distribution. This is a specified member of the same memory grammar;
these rational factors are not claimed to be the forward pixel-memory
design parameters.

Its complete four-mode probe is

```
H_r=e1(n1-n2)+e2(n3-n4)
    +(1/4+(1/2) enabled_r (P_new)00)
        (a1^* a2^*+a2 a1+a3^* a4^*+a4 a3)
    +(1/2)(N4-2)^2,
e1=1/64, e2=3/128.
```

The supplied kinetic, pairing and quartic terms act on all sixteen Fock
rows. The quartic and pairing do not commute. A direct triangle bound is
453/128<4 for every generator, independent of memory. The initial state
has two classically correlated memory blocks and an entangled two-row
untouched reference, with total field parity even.

Two intervals of duration 1/32 contain n=1,2,4 native steps each. Between
them an occupation instrument is followed by an even pair-rotation pulse chosen
from its result and the last native enabling flag. The final detector is
a pair-rotated occupation instrument, which does not commute with the
first. The mean-dynamics comparison uses exactly the same continuing
memory transitions and feedback policy. Every final detector pair has
positive probability; no successful outcome is normalized separately.

Fine records are grouped only when they have exactly the same new memory
and field operation. The factor table retains every constituent label and
its exact probability, including zeros. Since the feedback uses only data
constant on each such class, independent classical splitting recovers the
full fine-label law at every history. For identical quantum states the
trace norms add with the label probabilities, so this factorization
preserves the complete-history trace distance exactly. The n=4 execution
retains 52,488 class-history rows on a 32-row field/reference carrier.
It does not replace the archive by the final memory or an unread average.

The executed full-history errors decrease under these three refinements
and obey (1). Stored decimal bins summarize floating diagnostics; the
analytic error bound and rational resource schedule are the proof.
The four-mode probe is not an invariant slice of every 180-mode family
pulse and is not substituted for the compact source in that theorem.
Native compilation of the additional detectors and feedback is not
claimed by this execution.

Two exact countercontrols protect the scope. The next-step law conditional
on memory port two differs from the stationary law by total variation
56025/67055104; replacing the continuing memory by its stationary marginal
changes the record law. Separately the equiprobable pure states
`(4/5)|0> +/- (3/5)|1>`, compared with `|0>`, have unread trace distance
18/25 and retained-record distance 6/5=sqrt(2*(18/25)). Thus the unread
O(1/n) error cannot be assigned unchanged to the complete archive.

The construction closes the stated intervening-observation comparison for
this finite interacting source. It supplies no physical species, common
laboratory energy scale, gauge dressing, interacting chiral spatial
continuum, geometric modular flow or Einstein realization.
