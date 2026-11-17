# Intermediate matched representation refinement

This run inherits the unchanged gates in ../refinement_R2/GATES.md. It
replaces the resource-stopped run, not a scientific negative result.

R/h=1.5, nominal radius=4.5, L=60 and envelope width=9 preserve the same
nominal physical R, ball radius/R, qR and envelope/R as R1/r3/L40. The digital
region is the389-cell ball of squared radius20.25, giving1556 modes. Each
regional complex matrix is about39MB before eigensolver copies/workspace;
Fourier spinor arrays are about55MB each. Multiple arrays are accounted for
in the aggregate working set; this is not a claim that the entire process
uses only one matrix.

The raw uniform-sigma residual must decrease and be at most10% of the
nonuniform q1 signal before the latter's sign is interpreted. No subtraction
is applied to the preparation kernel. All derivative stencils and q0 remain.
R/ell is unchanged and is not tested by this representation refinement.
