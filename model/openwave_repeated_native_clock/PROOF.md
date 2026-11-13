# Two complete native phase cycles with a finite reference and reset battery

## Native carrier and exact collision

The inherited Gauss-dressed charged pair restricts to
`N=diag(0,0,1,1)` and `A=|3><0|`, with `[N,A]=A`. These operators are
restrictions of the actual native Fock Hamiltonian and its bounded partial
isometry. The native orthogonal complement is retained. The input lies in the
reducing carrier; its complement has probability zero. Inactive occupations1
and2 are also present in every displayed internal occupation instrument.

For a reference lowering operator b, set `B=A b+A* b*`. Then B is Hermitian,
`B^3=B` and `[N+n_reference,B]=0`. For real c,s with `c^2+s^2=1`,
`U=I+(c-1)B^2-i s B` is unitary and conserves the resonant source-plus-reference
energy. Its action is exact on every spectator resource coordinate.

Start the source in its vacuum, a preparation reference P and readout references
R1,R2 in `|+>`, and a reset battery B0 in its vacuum. All four resources have two
levels. The ideal sequence is preparation P, free source flight, readout R1,
full source-battery exchange, preparation using the SAME P, second free flight,
and readout R2. Each terminal instrument contains all64 native/resource
occupations and the orthogonal native-complement effect.

## Reset, depletion and both free frequencies

After preparation, the source coherence is `i s/2`, P has population `c^2/2`
and coherence `c/2`. Operations on source, R1 and battery do not change P's
marginal. A full exchange with an initially empty battery sends
`alpha|0,0>+beta|3,0>` to `alpha|0,0>-i beta|0,1>`, including arbitrary
correlations with P and R1. The source is reset on this active input; the
battery retains its state, energy and correlations. No erasure occurs.

Repeating preparation with the depleted P gives source population
`s^2 c^2/2` and coherence `i s c/2`. After both preparations P has population
`c^4/2` and coherence `c^2/2`. Under the isolated native Hamiltonian mN,

    <A>_1(t1)=i s exp(i m t1)/2,
    <A>_2(t2)=i s c exp(i m t2)/2.

Thus both ideal free frequencies are exactly m, independently of the coupling
and preparation duration. The second contrast is reduced by c. For a readout
exchange with coordinates cr,sr, the complete excited-reference probabilities
are

    p1=cr^2/2+sr^2 s^2/2-cr sr s cos(m t1)/2,
    p2=cr^2/2+sr^2 s^2 c^2/2-cr sr s c cos(m t2)/2.

The certificate executes the full4-by4 quarter-phase grid and changed
preparation coordinates, with every64-outcome probability exact. Twirling all
three initial references removes the two source coherences and the tested
phase dependence of the complete terminal distribution. These statements do
not identify a phase-sensitive quadrature with a zero-resource detector.

## Gap schedule and work

A reference held at zero gap retains its phase relative to the declared
controller. This is a resource. Initially only source, P and battery have gap
m, so the mean energy is m/2. Turn P's gap off after preparation1: work
`-m c^2/2`. Turn R1 on for readout1: work `m/2`. Turn P on for preparation2:
work `m c^2/2`, then off: `-m c^4/2`. Turn R2 on for readout2: work `m/2`.
R1 and battery retain their gaps. Total gap work is
`m(1-c^4/2)` and final energy is `m(3-c^4)/2`. Every exchange conserves its
resonant energy. The exchange switching work sums to zero because B is
conserved during its collision and commutes with the resonant free energy.
The final occupation measurement commutes with the remaining free energies.
These are identities for the stated collision schedule, whose timing is
implemented approximately by the separate finite controller below.

## Static controller and complete finite outcomes

Use an open512-position clock, hopping J=16 and scale L=4096. The complete
executed dimension is `512*4*2^4=32768`. The clock packet is the literal integer
vector in controller.py, normalized by the exact square root of its squared
norm, multiplied by `i^x`, supported on positions4 through36. It is a declared
finite resource, not a target-dependent mass fit.

The time-independent Hamiltonian divided by L contains the positive clock
Laplacian, `(m/L)N`, and position-dependent reference gaps. P is on for x<64
and256<=x<272; R1 for x>=160; R2 for x>=400; battery everywhere. Exchange
regions are preparation1[48,64), readout1[160,176), reset[208,224),
preparation2[256,272), readout2[400,416). Both preparation strengths are the
same g, either2 or3. Readout strengths are2, reset22/7. The latter is a
declared rational pulse strength, not an asserted exact full reset.

The full Hamiltonian obeys `H/L>=-max(|g|,2,22/7)`. All regions are disjoint;
the other summands are positive. A scalar shift gives a positive Hamiltonian
without changing records. The bounded exchanges preserve the native reducing
carrier and its complement. This construction does not truncate the native
Fock spectrum by declaring other sectors nonexistent.

All six native-mass/strength cases and a fully twirled three-reference control
are propagated. Every final clock/native/resource occupation probability is
retained, including reflected packets, failed resets and unread clocks.
The reset diagnostic conditions on clock positions228 through251 and reports
both the condition probability and all excluded outcomes. The first reference
is a quantum record until the terminal instrument; its overwrite bound is
reported. All clock kinetic, native, reference, battery and exchange energies,
and the work of the final complete occupation pinching, are recorded. The
battery is neither emptied nor its entropy discarded.

## Finite propagation enclosure and arithmetic assumptions

Subtract the scalar2J from H/L for propagation. Its exact real symmetric
matrix has at most four nonzeros in a row and column: two clock hops, one
diagonal and one exchange. For `0<m<=4096`, g in{2,3}, the absolute row and
column sums are below40. The centered operator norm is therefore at most40.
The largest diagonal contains at most four unit source gaps. Converting its
exact rational value and22/7 to binary64 gives an operator error at most
`beta=2^-46`, by the row/column norm inequality. Integer hops are exact.

Use216 steps of dt=1/16 and the degree28 Taylor polynomial. With z=40dt,
its exact one-step tail is at most `exp(z) z^29/29!`. The exponential upper
bound is obtained by the inherited rational series at z/4 and raised to the
fourth power. The sparse matrix-vector product has at most four real-times-
complex products and three complex additions per output; the subsequent
complex scaling and accumulation require fewer than64 real operations per
coordinate. Under correctly rounded IEEE binary64, their normwise error is
bounded by `gamma=2^-47`, using the absolute row/column sum bound. The
polynomial recurrence and summation have at most `k=2*28+2` such stages.
Their accumulated error is bounded by

    k*gamma/(1-k*gamma) * exp(z) * (1+2^-20).

The last factor bounds `exp(gamma*z)`. Add the Taylor tail, `dt*beta`, and
an absolute underflow allowance2^-900, then round outward to a2^-60 grid to
obtain delta. No overflow is possible in these bounded vectors. Normalization
of the declared integer packet and reference amplitudes has normwise error
at most `e0=2^-43`. For step n the state error is bounded by
`E_n=(e0+n*delta)/(1-n*delta)` because the error recurrence is dominated by
`(1+delta)E+delta`. The same Frobenius estimate applies to the eight-column
ensemble of the twirled input. A complete POVM's total variation in the
unhalved l1 convention is at most `(2+E_n)E_n`.

These statements assume the stated arithmetic behavior; they are not a
formal verification of NumPy or SciPy. Sparse storage does not approximate
or discard entries. Norms and expectation sums contain at most262144 complex
coordinates. The explicit2^-28 reduction allowance exceeds the worst linear
binary64 summation error for fewer than2^22 real operations. Every projected
term norm is increased by a relative2^-28 and an absolute2^-48. Products of
these upper bounds are rounded upward to a2^-40 grid with an additional grid
unit, which covers their bounded binary64 multiplication/addition errors.
Exact rational receipt fields use Fraction arithmetic. Closed-energy bounds
include the same reduction allowance and the full uncentered norm.

For times inside a Taylor step, the sum of norms of the projected endpoint
Taylor terms bounds the projected amplitude at every intermediate time.
Add the incoming state error and twice delta for local numerical and tail
errors. Integrating these positive bounds yields genuine interval-in-time
contact and boundary-current bounds, not an endpoint sampling assumption.

## Free-flight boundary current, actual frequency and residual bias

For a clock projector P onto an exchange-free corridor, all reference-gap
terms commute with A and P. Therefore

    d<P A>/dt = i m<P A> + i L<[K,P] A>.

The boundary term is at most `2 L J` times the sum of the two adjacent-cell
amplitude products. The two windows use the entire corridors[64,160) and
[272,400). Let D be its integrated scaled-time bound and a a lower bound on
the initial coherence amplitude. If a>D, demodulated coherence moves by at
most D and stays nonzero. Its phase differs from the native phase by at most
`D/(a-D)`. Dividing by the native time interval gives the stated conservative
frequency bound. Every coherent case also verifies
`m*Delta_t + D/(a-D) < pi/2` with a rational lower bound on pi.
Consequently the entire continuous phase change lies in the principal half
plane: the endpoint atan has no unrecorded winding. Conditioning on clock membership divides the coherence by
a positive real probability and does not change its phase; all failures are
retained in the complete outcome table.

A separate, tighter enclosure describes the actual endpoint phase increment.
The numerical coherences have complex-disk errors from the state and reduction
bounds. Their product is enclosed, its positive real part is checked, and
atan is bounded with an exact alternating rational series, range reduction,
and Machin's formula for pi. The frequency interval is rounded outward to a
2^-32 grid. Floating exp and atan2 are used only for diagnostics and redundant
consistency checks, not to construct this interval.

The finite controller does not exactly implement the ideal free-flight
collision schedule. Its second-window frequency has a resolved finite bias
that depends on g; the full intervals and boundary bounds retain that result.
The exact ideal frequency theorem must not be reported as exact frequency
independence of this complete detector. Nonzero second visibility and small
reset failure establish a finite repeated resource protocol, not a sustained
self-starting clock, a reference-replenishment theorem, physical seconds, a
native twelve-port process implementation, or an OpenWave row validation.
