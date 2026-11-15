# Native autonomous N16 experiment on the compact harmonic fiber

This packet implements a short-time experiment for the proposed autonomous
source `native-local12-autonomous-CAR-N16-twisted-fiber-v6-proposed`, frozen
in autonomous_car_vacuum df34a04c. It is an implementation theorem for that
declared source, domain and preparation. It does not identify this compact
fiber with a normalizable state of the original noncompact harmonic field,
prepare the interacting ground, establish a particle mass, or award a row.
The clock and all couplings remain those of the source: h=1/8, m=1, a=1/2,
b=3. The selected cell duration is u=1/64, hence physical duration t=1/512.

The initial source state is the exact product Omega that fills orbitals
4x+1 and4x+3 at every cell, the actual190-mode coexact oscillator vacuum,
and harmonic Fourier momentum p_phi=0. Omega is the free comparison
preparation; it is not the interacting ground. The N16 source proves its
nonzero interacting residual. Any untouched quantum reference is retained.
All operator and complete-instrument estimates below are stable under
tensoring that reference. The state-specific field estimate is not claimed
uniformly over arbitrary initial excited field states.

## The actual additive CAR generator

Use the 32 ordered occupation modes, alpha_a=sigma_a tensor X,
beta=I tensor Z and the actual positive integer routes from frozen field25.
For each positive cell/axis route x to y and each spin s, exactly one t has
alpha_a[t,s] nonzero. In cell units its Hermitian hopping block is

    B_l = [i alpha_a[t,s] c_yt^dagger c_xs W_J + h.c.]/2.

There are24 positive routes times4 spins, hence96 blocks. W_J is a unitary
on the single shared field. It multiplies by exp(i J_co.q_co) and shifts
p_phi by e_a. The exact quotient identity J^T F=e_a is verified for every
block. There is no separate field assigned to each fermion. The block norm
is1/2, including all spectators and Pauli-blocked occupations. The source
kinetic norm is therefore at most48 in cell units.

Let i be the destination mode and D_i=(1-n_i)+n_i W_J. Since W_J commutes
with all CAR operators and the origin differs from i,

    D_i (c_i^dagger c_j) D_i^dagger = c_i^dagger c_j W_J.

Consequently exp(-i du B_l) is D_i times the bare CAR exponential times
D_i^dagger. The chronological native word is D_i^dagger, bare hopping, D_i.
Codex's frozen ca3199f5 compiler supplies the exact Jordan-Wigner, number,
and commuting Pauli expansion functions. Their byte-pinned definitions are
loaded without its different species/Floquet source. For eta=i alpha[t,s],
eta real uses the X hopping quadrature with coefficient eta; eta imaginary
uses the Y quadrature -i c_i^dagger c_j+i c_j^dagger c_i with coefficient
-Im(eta). Both signs and both quadratures are checked.

This uses additive dGamma of the one-body Hermitian generator. It does not
use (Gamma(T)^dagger-Gamma(T))/(2i). The distinction has an executed control:
for T=exp(i pi/6)X on two orbitals, the correct additive generator on double
occupation is its one-body trace0, whereas the latter Gamma expression is
sin(pi/3). No exterior-power route unitary is substituted for this sum.

## Exact finite compact flux and coexact code

The finite code has32 matter bits,190 oscillator number registers0..R,
three even cyclic flux registers p=-M,...,M-1, and a clock plus two pointer
work bits. M is a power of2. Every binary oscillator padding column has zero
field and displacement generator. It is retained, never projected away.
The selected R95 uses7 bits per oscillator and labels96..127 are padding.
The selected M16 uses5 bits per harmonic coordinate. Total memory is1380
binary registers. The physical finite field subspace has dimension
96^190 times32^3; the complete code defines all128^190 times32^3 columns.

The actual kappa0 sector obeys p_phi mod2=r(n), where r_a(n) is the parity
of the occupied cell coordinates. A hop changes r_a and p_a together.
The cyclic wrap changes p by an odd number because2M is even, so every
complete hopping block preserves the twist, including at the boundary.
Mass, Coulomb and free field preserve it. Intermediate D_i controls can
enter another parity sector; the full code retains those sectors and the
conjugation returns to the correct one. There is no parity postselection.

Free harmonic energy is exactly |p_phi|^2/48. The cyclic shift is implemented
as a reversible increment controlled on n_i. The binary high-to-low ripple
uses a target flip controlled on all lower bits1 (increment) or0 (decrement).
Each multiply controlled flip is exp[i pi P(X-I)/2]; both the X and scalar
terms are compiled. This preserves conditional phases on every column.

For any history of k<=K actual directed kicks starting at p_phi=0, each
coordinate has absolute value at most k. If M>K, no wrap occurs in that
history, and the finite and infinite Fourier shifts and all intervening
free harmonic phases agree exactly. Thus the harmonic comparison error
through orderK is zero. Higher histories remain covered by the Dyson
remainder. This is not a Gaussian/Hermite truncation; no old harmonic shear
or tail estimate is imported. The actual long finite product may wrap;
its difference from the finite Hamiltonian exponential is separately bounded.

Coexact W is exp(i P_R J_co.q_co P_R). The190 coordinates commute, and the
matrix is the identity on every corresponding padding column. Free coexact
energy is the full diagonal omega_j(n+1/2) for n<=R, zero on padding.
The exact physical curvature/frame and vacuum discrepancy are paid below;
the approximate spectral frame is never called the exact physical frame.

## Longitudinal charge and mass

Keep the source Gauss reference G=2 at every cell root. With P_c the frozen
rational longitudinal vector for delta_0-delta_c and i=(c,s),

    P(n)=sum_i n_i P_c - 2 sum_c P_c,
    C(n)=||P(n)||^2/4.

The native polynomial includes the constant ||2 sum P_c||^2/4, the linear
coefficient (P_c.P_c/4-sum_d P_c.P_d)n_i, and pair coefficients
(P_c.P_d/2)n_i n_j for i<j. It therefore implements the G=2 source,
not the root-background convention in the inherited species example.
Every coefficient is rational. The diagonal polynomial gives C(Omega)=0
exactly and is checked on98 additional masks including N0 and N32.
Its Pauli coefficient l1 norm13771/924 bounds its operator norm on the
entire code, including occupations outside N16. This also avoids using a
fixed-N16 Coulomb bound as a bound on all Fock columns. The source sector
remains N16. The diagonal mass is h sum_i beta_i n_i. All free terms commute.

## Physical source to finite Hamiltonian

Expand the exact source and finite Hamiltonian in interaction-picture
Dyson histories of the bounded kinetic sum, using their respective free
field, mass and Coulomb operators. Each directed CAR term has coefficient
modulus1/2. A nonzero history from Omega remains a single occupation label;
at each prefix at most16 occupied modes times6 directed neighbors can move.
Pauli-blocked histories vanish in both sources. The absolute sum of weights
at order k is therefore bounded by (48u)^k/k!. Fermion signs, spin phases,
and exact Coulomb phases remain in each history; taking absolute values
only enlarges the comparison estimate. Both unitary-final Duhamel
remainders are bounded together by2(48u)^(K+1)/(K+1)!.

Use the independently reviewed orderwise coefficient certificate11f5772b:
all48 actual coexact displacement coefficient squared norms are bounded by
the exact upward Gamma^2=7115711258619106807/9223372036854775808. All full
integer currents have norm at most sqrt5. Thus at order k use
lambda=k^2 Gamma^2 and S=k sqrt5, for arbitrary ordered free intervals
of total duration at mostu. Geometry and coexact boundary estimates from
the field25 source transfer apply to each such shared-field history;
particle number changes the history weight, not those single-field bounds.

For lambda<=R+1 put

    p=exp(-lambda) lambda^(R+1)/(R+1)!,
    z=p(R+2)/(R+2-lambda),
    e_co=sqrt(2[z+190 lambda(R+1)p]).

The local Poisson endpoint bound requires the displayed monotonicity
condition. When it fails, use the complete-history bound2 instead.
The executable encloses exp(-lambda) by the reciprocal of a finite positive
Taylor sum, and all subsequent square roots by outward rational intervals.
It uses no floating value as a certificate endpoint.

Let eta be the frozen actual-curvature error, w0 the lower bound
sqrt[(3/2)(min Lambda-eta)], v=1/(4w0), and
delta=frame_residual/min Lambda+polar_Gram_error. The inherited actual-frame
bound for each history, with this packet's actual total duration, is

    e_geom = (3u eta v/2) sqrt[(4lambda+190)^2+16lambda+380]
             + delta S sqrt[v(4lambda+1)].

The complete history difference is at most min(2,e_co+e_geom); the harmonic
term is zero for k<=K. Both nonzero history outputs are normalized after
the field embedding, so2 is a bound on their complete vector difference.
It is not a normalization of a projected successful probability. If the
Poisson condition fails the complete difference is simply bounded by2.
Pay the independently certified actual/comparison vacuum difference once
outside the full evolution. The total physical comparison is

    e_phys <= e_vac + 2(48u)^(K+1)/(K+1)!
                + sum_(k=0)^K (48u)^k/k! min(2,e_co,k+e_geom,k).

For u1/64,K12,M16,R95 the exact receipt gives e_phys<1.921837e-6.
The embedding retains the actual field and all quantum-reference outputs.
It proves the specified initial-state experiment, not a uniform diamond
bound over arbitrary field preparations or a ground-preparation theorem.

## Finite Hamiltonian to a native word

On the complete finite code let H0 be the free diagonal and let B_1,...,B_96
be the hopping blocks. A conservative full-code bound is

    X0=190*11*(R+1/2)+3M^2/48+13771/924+32h.

The frozen frequencies are below11. Since ||B_l||<=1/2, commutator norms
give ||[H0,B_l]||<=X0 and ||[B_l,B_j]||<=1/2. The first-order Lie product,
chronologically free then the96 blocks, repeatedN times has error at most

    (u^2/N)(48X0+1140).

This is a finite-code bound after the physical comparison. The N1
energy-dependent splitting theorem is not silently extended to N16.
Its proposed source has a different current norm and initial domain.

Each controlled W is made from all190 adjacent-level pair products with
s subdivisions. Their operator error is at most Gamma^2 R(R+1)/(4s).
There are192N controlled W or inverse uses. A midpoint norm check on all
24 positive routes at the selected precision ensures that the same Gamma
bound covers the implemented midpoints. The inverse reverses the adjacent
pair order and all signs, using exactly the same rounded coefficients.
The bare CAR exponential's two Pauli terms commute, so it needs no further
splitting. Neither compact flux nor its free energy needs a product error.

The conservative coefficient/root rounding bound is

    [192N*190*(2R*1024+2R)+u*190*(R+1/2)] 2^(-bits).

All4560 selected actual coefficient intervals are checked to have radius
at most1024*2^(-bits). Compact increments, harmonic energies, mass and
Coulomb use exact rational or pi-multiple angles. The primitive-response
error uses the frozen full six-level12-control norm bound for every pulse,
not an ideal gate substituted for a native event.

For k=bit_length(R), a controlled adjacent pair has at most2^k Pauli terms,
each using at most8(k+1)+1 pulses. One controlled W is therefore bounded by

    190*s*R*2^k*[8(k+1)+1]
       + sum_(j=0)^(b_h-1) 2^(j+2)[8(j+2)+1].

The free word is priced for every padded diagonal term, all harmonic
columns, every Coulomb Pauli term and all32 mass terms. Each bare32-mode
CAR hopping costs at most514 pulses. The total includes a conservative2000
pulses for preparation, the common clock, and all33 terminal pointers.
Reset event sequences are distinct adaptive native services and are counted
separately. flow.py checks every actual bare block against514 and short
pair/flux words against their formula. The full program is a finite lazy
iterator; its enormous selected execution has not been performed.

## Complete source, clock and apparatus records

The native preparation first executes conclusive blank resets of all1378
source/clock bits. It flips the16 beta-negative occupation bits and the top
bit of each flux register, encoding p_phi=0 as labelM. It prepares clock|+>.
The actual reset retains every fine native letter, success location and
terminal timeout, together with the remaining reference. Conditional on
a conclusive blank outcome the prepared code state is exact. The reference
and classical reset transcript remain part of the output, not discarded.
The ideal comparison uses that same reset instrument/transcript and the
specified physical preparation; timeout branches are included in the bound.
No arbitrary input-reference correlations are assumed erased from records.

The program runs the finite words, advances the independent common clock
by exp[-i Z h u/16], then reads all32 matter occupations and the clock.
For X or Y clock measurement the declared analyzer is included. Each
pointer resets two work registers, applies the actual frozen native pointer
word, and executes the literal C or E event law with its fair coin. All
fine labels, coin values, output code, and remaining source/reference state
are retained. Sequential measurement backaction is part of the instrument.

Every controlled native response retains its first full-source export
complement and stops that branch. For a frozen pulse tape e_1,...,e_P the
failure map at j is its literal complement row times the preceding proper
code success maps, tensor identity on untouched registers. This formula
specifies every quantum-reference failure map, not just its trace. Terminal
pointer complements and all reset timeouts are likewise retained. There is
no successful-record renormalization or return-ancilla postselection.

There are1378+66=1444 reset register uses. Two uniform blocks per reset
give timeout union upper442225/1152921504606846976<3.836e-13. The uniform
conditional reset theorem, not independence of failed resets, justifies
this union bound. Each proper unitary error is a full response/isometry
error including export complements. Telescoping these isometries and
contractivity under the complete subsequent apparatus bounds the total
variation distance of every complete classical record law by the sum of
the physical, finite splitting, pair, algebraic, native and reset errors.
The same construction controls the corresponding retained-reference
instrument; scalar Born probabilities are not the only validated outputs.

The selected exact complete-record upper in BUDGET.json is below
3.921837e-6. It uses2,339,691,255 segments,790174877262418599815 pair
subdivisions,256-bit certified controls and a native pulse upper
53306912682581980540567117995463547586845. The huge program is formal,
finite and explicitly bounded. No claim that these resources were used
in a full physical193-coordinate execution is made.

## Executed scope and independence

controls.py executes four local native words on all128 input columns:
three CAR bits (origin, spectator parity, destination), one actual
coexact coordinate with cap2 and retained padding3, and a four-value
compact flux coordinate. It tests both hopping quadratures, a sign reversal,
and subdivisions1 and4. The Y/sign-1 word uses the first actual source
block; the X and opposite-sign cases are algebraic controls using that
same retained field coefficient, not additional full physical routes.
The entire Jordan-Wigner spectator parity is
represented by one bit for this local control; this is not a native
enumeration of all2^32 matter columns. The algebra proves that reduction
for each bilinear, and checks the phases of all96 actual source blocks.

clock_control.py also executes the common t1/512 clock's X and Y words
on all2 columns and reads each with both raw C/E laws. Their unconditional
contrasts agree with cos(t/8) and sin(t/8); complete failures remain in each
law. This is the same clock phase and analyzer used by the full iterator.

The raw native responses, success maps, complete sum of first-failure
Grams, every per-input first-failure probability and every fine C/E pointer
record are saved. The frozen tape plus explicit prefix formula above
reconstructs every failure Kraus map. LOCAL_MAPS.npz stores all128-column
success maps and aggregate failure Grams without storing a large redundant
archive of prefix matrices. The pointer itself is executed on all8 columns.
All-column native/logical residuals are below3e-14; the exact-W residual
decreases from0.00124924 to0.000312138 as pair subdivisions go1 to4.
These are finite diagnostic residuals, not the full physical budget.

The sign-erased Jordan-Wigner and omitted-W controls differ from the target
by nonzero operator distances. The Gamma(T) substitution fails separately.
Padding, flux wraps, CAR number conservation, exact G=2 Coulomb and the
physical clock remain explicit. Tagged full-flow checks exercise the actual
outer iterator with short replacement segment markers; they are labeled
as control-flow checks, not unexecuted native bodies counted as executions.

The N16 source, Codex CAR compiler, field25 native compiler, harmonic
quotient, orderwise field comparison and their opposing reviews are pinned.
The author of this packet also authored the quotient and reviewed several
inherited components. Their separate opposing-review manifests identify
the inherited boundary. This N16 composition has no producer
self-approval. The independent review 4f647b3c approves the component implementation on this domain. No physical-row criterion follows from that component alone.

The independent replay in `autonomous_car_native_review.py` consumes the frozen sibling-workspace source inputs named in its receipt. That dependency closure is not a self-contained public export.
