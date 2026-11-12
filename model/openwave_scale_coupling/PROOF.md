# Scale dependence measured from the declared colour ensemble

## Source and input ancestry

The five-level seam holonomy cost restricted to colour is exactly

    1-Re tr_5 diag(U,I_2)/5 = (3/5)(1-Re tr_3 U/3).

The declared four-dimensional periodic Gibbs construction with cost multiplier
10 therefore has the SU(3) Wilson weight at beta=6. The bare coupling in this
normalization is g0^2=6/beta=1. The lattice action and thermal state are declared
architecture. Neither beta nor a physical lattice spacing is derived from a
narrower repair law, and no experimental alpha_s value is supplied.

The existing native colour matrices give eight fundamental generators with
tr(T_a T_b)=delta_ab/2 and sum_a T_a^2=(4/3)I. The producer recomputes this
normalization and the seam-cost identity. The continuum SU(3) link ensemble
is the declared colour construction, distinct from the finite PW119 particle
Hamiltonian. No Wilson matter operator or quenched-hadron identification is used.

## Loop scheme and its finite geometry normalization

Let W(R,T) be the ensemble mean of Re tr_3 U(C)/3 for a rectangular loop,
averaged over all sites and both orders of all six coordinate planes. Set
W(R,0)=W(0,T)=1. The retained observables are

    chi_R = -log[W(R,R) W(R-1,R-1)/W(R,R-1)^2],
    g_C^2(R;N) = chi_R/k_R(N).

This defines a finite Wilson-loop coupling scheme at lattice length R-1/2.
The numerical ratio requires positive ensemble loop means; a nonpositive
estimate is reported as undefined rather than postselected away. Symmetry
justifies averaging the two orders of an asymmetric rectangle. The logarithm
is taken after ensemble averaging, not on each configuration.

The coefficient k_R removes the free finite-box geometry. On the N^4 periodic
grid let hat(k)^2=4 sum_mu sin^2(k_mu/2), with the constant mode removed.
For a rectangular conserved link current J_C, the exact finite Gaussian
covariance gives the quadratic small-coupling coefficient

    f(R,T) = C_F/(2N^4) sum_(k!=0) sum_mu |J_C,mu(k)|^2/hat(k)^2,
    k_R = f(R,R)+f(R-1,R-1)-2 f(R,R-1),  C_F=4/3.

For a loop in the first two axes, writing x=exp(ik_0), y=exp(ik_1),

    J_0=(1+x+...+x^(R-1))(1-y^T),
    J_1=(1+y+...+y^(T-1))(x^R-1).

Their discrete divergence vanishes. In particular k_1=(1-N^-4)/3 exactly.
The quadratic coefficient is verified against an independent position-space
inverse Laplacian and explicit oriented link currents. A Gaussian control
W=exp(-g^2 f) gives g_C^2=g^2 at every measured R. Without k_R the geometry
alone produces a false scale trend. This normalization is a defined free
covariance comparison, not a claim about the nonperturbative zero-mode
expansion or a derived continuum renormalization prescription.

The interacting ensemble means are measured using the existing SU(3)
heat-bath/overrelaxation code. No one-loop beta function, interpolation of
experimental values, or target running curve enters this calculation.

## Executed ensembles and statistical scope

Four independent productions retain 400 loop records each: cold/hot starts
on sides six and eight, all at beta=6. Each measurement follows two compound
sweeps, each compound sweep consisting of one heat bath and two
overrelaxations. Warm-ups are respectively 200,300,300,400 compound sweeps;
seeds are 3101,3102,3103,3104. All five rectangle measurements and the
independent plaquette measurement are retained at every sampling time,
together with the final link configuration and production metadata.

Errors use common-block deletions across all five loops so correlations
survive the nonlinear logarithm and scale differences. Block sizes
5,10,20,40,80, first/second halves, independent starts and both volumes are
reported. These are Monte Carlo diagnostics, not rigorous mixing guarantees
or certified confidence intervals. No sample is removed on the basis of its
loop value. Gauge transformations of the final configurations preserve every
loop; the separately implemented plaquette agrees at every record.

At N=8, combining the two chains with blocks of twenty records gives
g_C^2 approximately 1.56186, 2.73542, 4.13672 at lengths 0.5,1.5,2.5,
with blocked standard errors approximately 0.00050,0.00444,0.03837.
The corresponding N=6 values are 1.55976,2.75157,5.43421. Thus the coupling
weakens toward smaller measured lengths on both volumes under one unchanged
bare law. The large N=6 versus N=8 difference at the longest loop explicitly
exposes finite-volume sensitivity even after free-geometry normalization.

This is an own-ensemble finite scale-dependent coupling. It does not determine
a continuum beta function, a physical GeV scale, the observed QCD running
curve, or an onset at an identified particle core. A continuum step-scaling
study would require matched renormalized couplings and a family of cutoffs;
varying the bare beta alone would not establish that result.

## Reproduction

`--write` and `--check` regenerate the analysis from every raw loop record.
The separate production command generates the Monte Carlo records themselves;
it does not merely replay stored summaries. The producer binds the code,
proof, raw files and dependency versions. Bitwise raw-array reproducibility
uses the declared Python/NumPy/Numba stack; elapsed runtime and archive bytes
are not deterministic outputs of a stochastic-law replay. The physical
statistics and numerical comparison tolerances are separate from file custody.
