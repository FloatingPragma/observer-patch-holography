# Colour records after one unit of free flight

Prepare the published finite timing-filter output from the sharp face (0,1,5) in owner zero. The filter admits 0.199529 of launches. Open the full federation, including every outward bridge. On the admitted source, read the three fitted triplet receivers of the original cell after a wait of one exchange unit. Every other position and failed triplet admission counts as zero for the joint singlet probability. The receivers do not select a face or a connected triple before testing colour. Source timing, the wait, the receiver frames and the joint singlet analyzer are supplied. [Source and measurement settings](../frames.json) · [Finite flight calculation](flight_check.py) · [Results for all three boards](flight_check.json)

Write S for the effect of triplet admission followed by its antisymmetric singlet and T for triplet admission alone. The expectation of W = S − T/3 is positive only when the admitted triplets pass the biseparability bound. The following probabilities are per filter admission. The intervals account for every triple outside a four-owner numerical simulation. They are deterministic propagation estimates under the stated arithmetic tolerance.

| Owners | Native joint singlet, lower bound | Exchange-only joint singlet, upper bound | Native W, lower bound | Exchange-only W, lower bound |
| --- | ---: | ---: | ---: | ---: |
| 20 | 0.27482 | 0.16098 | 0.04196 | 0.01547 |
| 80 | 0.27526 | 0.15509 | 0.04240 | 0.02138 |
| 320 | 0.27526 | 0.15509 | 0.04240 | 0.02138 |

The contrast control opens the same ports on the same prepared source and drops only the contrast term during flight. The prepared colour correlation survives in both laws at this wait. The higher singlet yield per filter admission on the native law distinguishes their flight. A separate receipt panel gives both laws at 1/16 of an exchange unit. The result at these appointments uses a local colour read. The all-time result concerns spatial connection. [Source and fragmentation proof](PROOF.md)

## Exterior error

Keep the home owner and all its neighbours. Let P project onto configurations with all three loads inside this set. The compressed generator A = PHP includes the diagonal cost of every native hop out of the set. Let B = (I−P)HP. Each configuration reached in one outward hop has exactly one load outside. Its port can return along only one external seam, so this configuration has a unique inside parent. Thus B*B is diagonal. Its element on an inside triple is one quarter the number of available outward hops, at most three.

For a normalized source ψ inside, the distance between full and compressed evolutions at time t is bounded by the integral of the residual norm f(s) = ‖B exp(−iAs)ψ‖ from zero to t, by variation of constants and unitarity. The spectrum of A lies in [−6κ,18]. The three-load exchange graph has degree at most 18, and the occupied-edge potential is −2κ times zero through three. Subtracting the spectral midpoint leaves f unchanged and gives it Lipschitz constant at most half the spectral width times ‖B‖. On a time grid of spacing h, integration of the nearest-endpoint Lipschitz bounds gives the trapezoid estimate plus t h times this constant divided by four.

The receivers use the nearest isometries to their published rounded frames. The program uses h=1/256 and degree-ten Taylor steps for the centred Hermitian matrix. If its centred norm is at most r, one substep has norm error at most exp(rh)(rh)¹¹/11!. For m substeps and per-step bound d, the accumulated error is at most exp(md)−1. The enclosure includes this error at the endpoint and in the residual integral. The independently reconstructed finite-filter source is compared, up to a common phase, with the source that has the dyadic spectral residual certificate. The enclosure adds their vector discrepancy, the certified calculation’s finite-source error budget (about 3.02×10⁻⁶ including its scalar tolerance), and 10⁻⁸ for the floating-point propagation. The latter is a numerical tolerance without an interval-arithmetic guarantee. [Source error calculation](measure.py) The effects S, T and W each have spectral diameter at most one, so the resulting state-distance bound also bounds each reported probability or witness. The independent short-time argument uses the initial variance and the fitted source overlap.

From the RER root, replay with `OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.candidate_laws.interval_echo.openwave.row_13_baryon.open_release.flight_check`.
