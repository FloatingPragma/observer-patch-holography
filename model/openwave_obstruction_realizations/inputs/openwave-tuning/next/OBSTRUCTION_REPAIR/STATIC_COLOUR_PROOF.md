# Actual static probes and screening on the signed vacuum

The previous bare string cost was not an interacting potential. This
calculation adds an explicit external fundamental/antifundamental pair to
the unchanged four-site signed-hopping source at fixed fermion number6.
It retains their nine colour states, the whole fermion Cartan carrier,
and the total colour Gauss constraint. It is a new declared static-probe
experiment; there is no native apparatus or physical unit transfer claim.

For a prefix, define the actual U(3) currents J_ab as the fermion currents,
plus E_ab for the external fundamental when present, minus E_ba for the
antifundamental when present. If n is their signed U(1) trace, the electric
Casimir is

    C² = (3 sum_ab J_ab J_ba − n²)/6.

All prefix representations are retained. The complete total Cartan-zero
carrier has dimension1224 at N=6. Every source component commutes exactly
with the total Casimir. A positive100 C_total² penalty isolates the
physical sector: a nonsinglet has C²≥4/3, while the kinetic and mass terms
are bounded below by−1.5−6m. A local singlet trial is below8+6m throughout
the three declared couplings. Thus no nonsinglet can be the lowest state
of the penalized matrix. This proves sector isolation, not the accuracy
of a particular numerical ground-state search.

The returned low eigenvalues are numerical Ritz diagnostics, with ground
residual below2.0e−12 and total Casimir below6.6e−17. No certified ground
interval is claimed. The complete zero-fermion control exactly reproduces
the bare4κd/3 cost to numerical precision. At coincident probes the
original vacuum plus an external singlet supplies the reference. The full
coincident carrier also includes the adjoint pair, which is retained.

Measured finite potentials E_with_probes−E_original_vacuum, in source units:

| κ | d=1 | d=2 | d=3 | Bare value at d=3 |
| --- | ---: | ---: | ---: | ---: |
| 1/2 | 0.06038729 | 0.10827975 | 0.05948282 | 2 |
| 1 | 0.02906153 | 0.05081240 | 0.02658196 | 4 |
| 2 | 0.01427365 | 0.02421110 | 0.01256034 | 8 |

The responding fermions strongly screen the probes. None of these short,
boundary-sensitive sequences supplies a linear-potential regime. The
large bare coefficient is therefore not a physical string-tension
prediction on this source. Screening itself is not an exclusion of colour
confinement; a physical theory with dynamical quarks can break strings.
The required intermediate-distance regime and physical scale remain open.

An independent exact variational witness confirms a mechanism for the
screening. At the fundamental's cell put two fermions in its colour
antitriplet; at the antifundamental's cell put one fermion in its colour
fundamental; put the remaining three fermions in a local colour singlet
on an unused cell. `colour_screening_bound.py` constructs all nine integer
amplitudes, checks every prefix generator exactly, and verifies zero
electric field and zero hopping expectation. The trial energy is−2m,0,
−2m for separations1,2,3 respectively. With the independently certified
original vacuum lower bound E0_low, this gives

    V(d) ≤ E_trial(d) − E0_low.

The bound does not grow with the bare string length cost on this fixed
box. It is exact and independent of the Ritz calculations. It does not
identify a physical nucleon, an infinite-volume potential, or a mass/unit
map. A larger geometric and interacting source must carry those tests.
