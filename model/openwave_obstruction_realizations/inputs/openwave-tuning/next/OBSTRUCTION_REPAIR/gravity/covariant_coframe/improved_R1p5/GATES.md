# Fourth-order representation diagnostic, before execution

The source changes explicitly to the standard fourth-order differences
proved in ../NEXT_CANDIDATES.md and implemented in ../improved_source.py.
R/h1.5, radius4.5 and L60 match the completed second-order run exactly.
No coefficient is selected from entropy or gravity results.

The existing gates are unchanged: preserve the raw zero mode and all
derivative stencils; require decreasing uniform-sigma residual, at most10%
of the nonuniform q1 signal before interpreting its sign. Do not subtract
the residual from the preparation kernel. The separation-aggregation trace
must first agree with the full finite trace and independent dense control.

Passing this gate would validate only improved representation of the
geometric observer response. The new flat source needs its own eta, modular
and finite-state tests. The fixed R/ell ratio, mean-energy condition, lapse,
native source selection and physical unit map are not discharged.
