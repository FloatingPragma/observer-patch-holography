# Independent full-source platform review

The review reconstructs the exact native Nambu symbol independently of the portable exporter, checks its full Clifford algebra, and verifies the Runge–Kutta contraction polynomial. It then builds the complete periodic sparse action on both cube sizes and evolves all three charged preparations with a separate matrix-exponential routine. Every source-coordinate probability is compared with the saved Fourier result. Preparation failure, source complements, integration bounds, numerical drift and the changed-sign control retain their stated domains.

Run `python3 -m model.openwave_platform_review --write` or `--check` with one BLAS thread. The review pins its own code and all target sources and records; it does not modify the reviewed package.
