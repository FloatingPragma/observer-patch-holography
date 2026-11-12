# Prepared-domain averaging of the complete native current

## Source and preparation

The source, all270 CAR coordinates,31464 native quartic words, filled negative
free bands, vacant zero bands and opposite-charge reference are exactly those
of `openwave_current_averaging`. Only the declared contact coefficient G varies.
Let P span its normalized parent and receiver. Both have free energy1. Exact
full-current free-clock averaging preserves P, with matrix

    P Vbar P = [[173/2,1],[1,85]].

Consequently the averaged receiver probability is
(16/25)sin²(5Gt/4). This package compares that derived orbit with the complete
unfiltered evolution on P tensor an arbitrary reference. No invariant claim
for P under the unfiltered source and no full-input unitary approximation are
made. Every evolved complementary state is retained.

## Exact columns and reference contractions

Each original CAR coordinate is transformed by the source's exact
Bogoliubov coefficients. Coefficients are polynomials in x_f,y_f and sqrt2,
with x_f²+y_f²=1. Polynomial reduction performs all cancellations exactly;
source intervals then enclose the remaining algebraic coefficients.
The complete parent and receiver current columns contain4236 and4087 nonzero
entries, respectively. Removing their two resonant outputs leaves actual
nonzero energy differences of modulus at least0.00481278138. No minimum over
unreachable current transitions is substituted for this column calculation.

The producer also normally orders every transformed word, using
 a_i a_j*=delta_ij−a_j* a_i, retaining every contraction. Thus

    V = 84 I + sum_(a+b=2 or4) W_ab,
    W_ab = sum_(|I|=a,|J|=b) C_ab(I,J) a_I* a_J.

The scalar84 and all quadratic terms are part of the actual source. The
coefficient matrices may be rectangular. Their operator norms are bounded
by the smaller of their Frobenius bound and the square root of the product
of maximum absolute row and column sums. Rational coefficient enclosures and
outward square roots certify every bound. Free-clock conjugation multiplies
rows and columns by unit phases and preserves these bounds.

## Number-sector norm lemma

For a normalized n-particle vector, the ordered-annihilation map satisfies

    sum_(|J|=b) ||a_J psi||² = binomial(n,b).

Its adjoint creation map into k=n−b+a particles has norm
sqrt(binomial(k,a)). Hence

    ||W_ab restricted to N=n||
       <= sqrt(binomial(n,b)binomial(n−b+a,a)) ||C_ab||.

These identities concern all coordinates, with every spectator included.
They follow directly from CAR occupation counting. They do not truncate the
dynamics. The maximum particle number needed to bound the current times the
oscillatory primitive is11; the actual evolved state is unrestricted.

## Prepared-domain primitive and Duhamel bound

Set c=343/4, Vc=V−cI and B=P(Vbar−cI)P. Then ||B||=5/4. Subtracting this
scalar changes only a common phase. Let

    A(t) = integral_0^t (V_I(s)−Vbar) P ds.

For each nonresonant column entry v_rp of frequency omega_r, its coefficient
is v_rp(exp(i omega_r t)−1)/(i omega_r). Orthogonality of Fock words and the
Hilbert--Schmidt bound give, in each output number sector n,

    ||Pi_n A(t)|| <= A_n
       = 2 sqrt(sum_(p in P,r in n outside P) |v_rp|²/|omega_r|²).

Only n=3,5,7 occur. Combining them gives A=sqrt(sum_n A_n²). Applying the
number-sector lemma to each A_n and grouping equal output numbers yields

    ||Vc_I(t) A(t)|| <= K.

The source certificate gives K below2672350.169 and A below3269. Duhamel's
identity between the full interaction-picture unitary and the averaged
unitary is integrated by parts on A(t). The initial boundary vanishes;
unitarity on the left and invariance of P on the right give

    ||(U_I(t)−exp(−iG Vbar t))P||
       <= G A + G² |t| (K+5A/4).

This estimate does not use a norm bound for V on the whole Fock space inside
the integral. The left unitary remains the complete unfiltered propagator.
Tensoring an untouched reference changes none of these operator bounds.
Using2pi/5<44/35 gives a first-peak error coefficient below3367930.443.

## Complete records and controls

At G=2^-30 the first-peak receiver probability differs from16/25 by less
than0.00628. This coefficient is2^42 times2^-72. At G=2^-44 the probability
error is below4e-7. These are alternative declared source coefficients, not
an inferred physical weak strength. G=1/64 has no nontrivial guarantee from
this bound and is retained as a failed approximation control.

The complete instrument records parent, receiver and the entire orthogonal
source complement. If delta is the state error above, each displayed event
error is at most2delta and the complementary probability is at mostdelta².
There is no conditional renormalization. The exact two-state records are
executed; the large complete-source comparison is analytical.
Free energy commutes with these projectors. The contact contribution bounds
the absolute mean readout work by2G times41952; the interaction energy and
external opposite charge are not dropped.

The result is a controlled finite coherent conversion. It does not identify
physical muon decay, a neutrino spectrum, an irreversible limit or a continuum
rate. The apparatus and source-time convention retain their specified roles.
