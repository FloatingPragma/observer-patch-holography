# An accessible charged threshold of the native finite source

## Source and question

The Hamiltonian is the unchanged208-state native quantum-flux source in
`particle_observer_outputs.field_route.exact_rotor(6)`. It retains all16
occupations of four CAR modes and all13 flux labels. Its clock is
`U(t)=exp(-itH)`, with the same energy convention as
`particle_ground_spectral_records`. This calculation asks whether its
lowest charged spectral thresholds have nonzero original CAR detector
weight. A sector minimum with zero matrix element would be invisible to
that detector.

This source is a finite directional control. Its energies are not GeV
masses, and its charge labels do not identify observed particles. The
spatial Wilson family, its Gauss-constrained ground preparation, the
208-state thermal preparations, and this pure ground approximant are
distinct experiments. No state or clock is transferred between them.

## Isolating the source spectral windows

The inherited rational LDL calculation encloses the minimum of each
occupation sector. The calculator uses48-bit intervals, and independently
counts the complete spectrum with Sturm's theorem on each real rational
reducing block. These blocks partition all208 source states.

Square-free factorization restores algebraic multiplicity: Sympy's
`count_roots` counts distinct roots, so each square-free factor's count
is multiplied by its multiplicity. All interval endpoints are checked
to be off the spectrum. For a Hermitian matrix, this count equals the
spectral projector rank.

Exactly one source eigenstate has energy below-2. It lies in N=2, and
the inherited interval contains exactly this eigenvalue E0. In each of
N=1 and N=3, exactly two eigenstates lie below-3/5, and both lie inside
the inherited lowest-energy interval `[l_c,u_c]`. The source factor
polynomials are

```
p_ground(x)=x^5-54x^4+958x^3-5328x²-5816x+23184,
p_charge(x)=x^3-35x²+286x+341.
```

Exact polynomial division attaches these factors to their source blocks.
Each displayed interval contains one root of its factor. Hence the
lowest charged energy E_c is a single algebraic energy with multiplicity2
in each adjacent-charge sector. Its complement has energy at least-3/5.
The transfer `E_c-E0` is approximately1.107713183 source units.

## A prepared source ground approximant

On the exact ground-containing block, start from its first basis vector
and apply `(H+9I/4)^(-1)` eight times. Clear rational denominators and
common integer factors after each multiplication; this changes no ray.
The resulting integer vector v defines the normalized pure state
`rho_v=|v><v|/<v,v>`, restored to the native phase chart.

This is a deterministic finite preparation. A unitary taking a declared
ready vector to `v/sqrt(<v,v>)` has an algebraic Householder construction;
the vector is specified before preparation. Inverse iteration here
calculates that vector and is not a successful-measurement filter. The
preparation supplies no discarded trials or empirical mass inputs.

Let `[l_0,u_0]` enclose E0, let `c_0=-2`, and calculate the exact
Rayleigh energy `E_v=<v,Hv>/<v,v>`. The unique ground projector P0 obeys

```
1-Tr(rho_v P0) <= (E_v-E0)/(c_0-E0)
               <= (E_v-l_0)/(c_0-u_0).
```

Taking an outward rational square-root enclosure gives the half-trace
distance delta between rho_v and the exact pure ground rho_0. The
certificate has `delta<1.065e-7`. Both states lie in N=2. The trial is
not asserted to be stationary. Complete-channel contractivity transfers
all its detector probabilities to the stationary ground within delta.
There is no assumption about convergence of inverse iteration: the exact
Rayleigh calculation admits the resulting trial independently.

## A complete compensated CAR instrument

The inherited eight branches are

```
K_(+,i)=c_i†/2,       K_(-,i)=c_i/2,       i=0,...,3.
```

Their effects sum to the identity on the entire208-dimensional source.
The addition and removal effects are `I-N/4` and `N/4`; each is exactly
`I/2` on the N=2 input sector. References carry opposite charge and odd
parity, so the complete source/reference map is charge conserving and
even. All four modes, both signs and all quantum references are retained.
Mode labels specify this apparatus's native component basis; they are
not asserted to measure physical spin.

The reference's idle energy is declared to be zero. The source energy
transferred by K belongs to the intervention's controller account.
Changing reference energy changes the total apparatus energy account;
it does not change the source threshold being calculated here.

## An exact rational energy-window apparatus

After sign s, the apparatus targets occupation `2+s`. On that complete
sector define, with `w=1/16`,

```
D=(H-l_c I)/w,
A=(I+iD)^(-1),        B=iD(I+iD)^(-1),
T=A†A=(I+D²)^(-1),    B†B=I-T.
```

All matrix entries are rational complex numbers. A, B and their adjoints
commute, so the pointer collision

```
U_pointer = [ A  -B† ]
            [ B   A† ]
```

is exactly unitary. Outside the target occupation sector use A=0 and
B=I. This extends the apparatus to every source state: other charges
produce a retained failure rather than an inadmissible input.

Repeat four collisions, stopping subsequent collisions after the first
failure and retaining its time label. A five-state pointer has one
continuing state and four failure states. Each collision rotates the
continuing state with its corresponding failure state and acts as the
identity on the other stored records. Its complete Kraus list is

```
A^4, B, BA, BA², BA³.
```

Their effects telescope to the identity. The original CAR branch and
these five pointer outcomes give40 records. The calculator evaluates
every probability exactly from the rational effects
`T^4, I-T, T(I-T), T²(I-T), T³(I-T)`. Their raw sum is one, and the
sum for each charge sign is1/2. First-failure histories and their
quantum references are retained. No division by a success probability
occurs.

The apparatus is a finite, source-derived controller for an energy
observable of the same H. It is not claimed to be a local spacetime
detector or an efficient native gate sequence. The separate Fourier
construction below implements energy readout through the declared
source time evolution, with its error and observation time explicit.

## Bounding the actual lowest-window residue

The success kernel is

```
f(E)=[w²/((E-l_c)²+w²)]^4.
```

Inside the lowest source window it is at least
`a=[w²/((u_c-l_c)²+w²)]^4`. On the target-sector complement it is at
most `b=[w²/((-3/5-l_c)²+w²)]^4`. These statements use the exact
spectral counts and the monotonicity of f above l_c.

Let p_v be the original probability of window success for one sign,
summed over all four CAR modes. Let W be the corresponding exact-ground
weight of the true lowest spectral projector. Both ground and trial
states lie in N=2, where each sign's entire effect is `I/2`. Therefore
the pass probability error is at most `delta/2`. The original charge
branch has weight1/2, so

```
p_v-delta/2-b/2 <= W <= (p_v+delta/2)/a,
0 <= W <= 1/2.
```

No eigenvector residue is taken from a floating eigensolver. Exact
inverse iteration, integer square-root bounds and rational detector
weights certify a strictly positive W close to0.498558 for each sign.
The JSON receipt encloses the exact rational bounds outward.

The inherited antiunitary charge conjugation preserves H, reverses G
and exchanges the complete addition/removal frame. Uniqueness implies
that it preserves the ground ray. It also exchanges the two lowest
charged spectral projectors. Consequently the two W values are exactly
equal. Intersecting their independent rigorous intervals tightens their
common enclosure. The equality is a finite source symmetry statement,
without a physical antiparticle assignment.

## A finite source-clock energy-readout plan

Prepare a uniform M-state clock pointer, apply controlled powers of the
same `exp(-i tau H)`, Fourier-transform the pointer and retain every
clock label. The inherited positive Fourier kernel gives a complete
energy POVM. Here the full source has `|E|<=22`.

The declared plan uses `M=131072`, `tau=1/16`, energy resolution1/64
and a reported window of radius1/32 about `(l_c+u_c)/2`. The true
lowest charged window lies inside the reported interval eroded by the
resolution. Its enlarged interval remains below-3/5, so its complement
carries no additional charged source line. The alias bound is zero for
the full source support. `spectral_readout.kernels.real_energy_leakage`
supplies the finite filter error.

Each clock step uses32768 Cayley substeps of the same H. If its operator
error is e, the largest controlled power has error at most `(M-1)e`.
Since one sign has original weight1/2, its event error is bounded by

```
(Fourier_leakage + complete_clock_error)/2 + delta/2.
```

The certificate combines this with the ground residue enclosure and
gives a positive original probability interval for each charge's
reported energy window. The exact rational40-record apparatus is
evaluated; the full131072-bin Fourier histogram is a bounded protocol,
not a claimed executed histogram. The maximum source observation time
is `(M-1)tau`.

Subtracting the certified source ground-energy interval translates the
final-energy readout into excitation energy, with the ground interval's
width retained. It neither inserts a chemical potential nor selects a
different clock. A relativistic dispersion relation, spatial localization,
physical spin, scattering and measured mass remain separate source
requirements.
