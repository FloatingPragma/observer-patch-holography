# Exact first rest-mass correction and a compact-shell obstruction

This audit uses D3's actual electromagnetic source and its free filled
negative-band reference. That preparation is not the interacting ground
state. The result is about the strict order-e² correction above that
reference, conditional on placing the actual shared lepton-register
coefficients into its bare mass slots. No native family attachment or
interacting pole is asserted.

## 1. The scalar reduction is exact

Write s² = Σa sin²ka, E = sqrt(μ²+s²), and

    Bj(k) = Σa |βa(k,j)|² cos²(ka/2) ≥ 0.

At external momentum zero the two route vertices coincide. Projecting onto
the β=+1 block and averaging its two spin components gives

    Splus,j = −Bj(E−μ) / [2E(E+ωj−μ)],
    Spair,j =  Bj(E+μ) / [2E(E+ωj+μ)].

Their sum is exactly

    Splus,j + Spair,j = μ Bj ωj / [E(s²+2Eωj+ωj²)].

The projected instantaneous term is μ v(k)/(2E); the seagull vanishes at
rest. Consequently the fractional correction is

    F(μ) = ∫dk/(2π)^3 {v(k)/(2E)
             + Σj Bj ωj/[E(s²+2Eωj+ωj²)]}.

The same identity holds for every positive quadrature sum. The script
verifies the rational algebra symbolically and checks the actual Clifford
matrix contractions on 27 cached mass/momentum combinations; the maximum
absolute mismatch was 1.82e−12. This is a cross-check, not a numerical proof
of the exact identity. The matrix derivation supplies that proof.

Every displayed term is nonnegative and decreases as μ increases. For the
second term the denominator is D(E)=E(s²+2Eω+ω²), with
D'(E)=s²+4Eω+ω²>0. Thus F(μe)−F(μμ)>0. No field-mode contribution can cancel
a lower bound obtained from the Coulomb term alone. The rest spin average
is the quantity explicitly audited here; no unproved full pole or
nonperturbative degeneracy assumption enters the bound.

## 2. A source-specific bound that needs no momentum grid

The LOCAL12 Laplacian has 30 internal links and two links along each of the
three spatial directions. For u=(1,...,1), the internal links contribute
zero and

    u†L(k)u = 4 Σa(1−cos ka) ≤ 2 |k|².

Cauchy-Schwarz in the positive L(k) inner product gives
(L(k)⁻¹)00 ≥ |e0†u|²/(u†L(k)u). Since the actual electric coefficient is
a=1/2, the actual root Coulomb kernel satisfies

    v(k) ≥ 1/(4|k|²).

Let μe=h me and μμ=h mμ. The complete certified register-effect intervals,
not just floating midpoints, obey me<1/12500 and mμ>3/250. On the compact
spherical shell h/250≤r≤h/125, which excludes every singular corner and
lies inside the Brillouin zone for all three declared h values, s²≤r².
The difference 1/sqrt(μe²+x)−1/sqrt(μμ²+x) decreases with x. Hence throughout
the shell it is greater than

    [1000000/8001 − 625/9] / h.

The rational square inequalities justifying these two square-root bounds
are checked exactly in run.py. Integrating the Coulomb bound over this
shell and using π<22/7 gives the exact, quadrature-independent result

    F(μe) − F(μμ) > 2765/1966976
                       = 0.0014057111017114596.

This exceeds 1405 ppm, uniformly at h=1/8,1/16,1/32. The strict order-e²
relative correction to mμ/me is F(μμ)−F(μe), so it has the needed sign but
already exceeds the approximately 9.83 ppm asymptotic discrepancy by more
than two orders of magnitude. The direct bare-register transplant into
this leading filled-reference correction therefore does not produce the
desired small repair. It overshoots before other positive contributions
are included. The finite register's discrepancy is also much smaller than
this bound. No fitted counterterm is introduced.

This is **not** an all-orders exclusion. Higher-order remainders are not
bounded, the interacting ground state differs from the reference, and a
physical common-clock/native family attachment remains unproved.

## 3. The existing tiny-mass grid cannot be used as a prediction

D3's cached grid refines the zero-momentum cell but includes the seven
other sine-doubler corners as individual coarse samples. For the tiny
transplanted electron mass their point weights give artificial 1/μ terms;
the corresponding continuum neighborhood is integrable. At h=1/8 the
diagnostic electron F is 244.262, of which 244.008 comes from those seven
points. Dropping them gives 0.254677, but that is also not a certified
integral: the missing central-cell error and unresolved neighborhoods have
not been bounded. The exact compact-shell bound above is independent of
all these quadrature values.

The script retains all four cache prefixes, all three h values, the D3
mass-one reference, all three actual register coefficients, coupling-off
controls, and omitted-pair controls. Omitting the pair term destroys the
positive scalar combination. No apparent grid convergence or omitted
sample is promoted to a physical lepton mass.

The remaining meaningful task is an interacting-vacuum response with a
verified family attachment and controlled remainder or nonperturbative
calculation. Retuning the original phase, round count, or an empirical
counterterm would not supply that missing derivation.
