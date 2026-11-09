# Original boundary coherence, record work and regional information

## 1. One actual source and one induced regional cut

Use the [primary contract](../particle_common_contract/PROOF.md) and
[induced cut](../particle_primary_geometric_cut/PROOF.md), with
R={0,e1}. All matter in R and its one interior link are inside. Crossing
links, other matter, the complete apparatus and the idle canonical
reference remain present. Source entropies and mutual information below
refer to specified marginals of that retained experiment.

The selected source state is the actual positive polynomial rho=f(H) on
the physical Gauss carrier. Its authenticated floor is rho>=c P_G, with
c=2^(-N). The actual polynomial modular derivative satisfies

    beta_low <= K'(E) <= beta,  K(E)=-log f(E),
    beta_low=beta(1-delta)>0.                              (1)

These are the exact source parameters and rational delta from
[the modular-work calculation](../particle_primary_gravity_premises/PROOF.md).
No Gibbs density, trial vector, physical temperature or proper clock
replaces this state. Since f(E)>=c on its spectral enclosure,

    -f'(E)=f(E)K'(E)>=c beta_low.                        (2)

The finite regional restriction, support and Schur decomposition are
supplied by the actual cut, using the same finite restriction theorem.
The distinction between an algebra with center and an extended regional
factor follows the entropy conventions discussed by
[Casini, Huerta and Rosabal](https://arxiv.org/abs/1312.1183).
All estimates here are finite identities for this specified source.

## 2. Exact mutual information with matching boundary sectors

The physical space decomposes as

    H_phys = direct_sum_lambda
       M_R,lambda tensor M_out,lambda tensor span(Omega_lambda),

where Omega_lambda is the normalized invariant pair of dimension d_lambda.
Let Pi_lambda be the actual central projectors. Set

    D(rho)=sum_lambda Pi_lambda rho Pi_lambda,
    p_lambda=Tr(Pi_lambda rho),
    Pi_lambda rho Pi_lambda=p_lambda tau_lambda tensor
                              |Omega_lambda><Omega_lambda|,
    S_edge=sum_lambda p_lambda log d_lambda.

The multiplicity density tau_lambda can have arbitrary correlations.
Coherences between different global lambda sectors are retained in rho.
The source partial traces annihilate those coherences and give

    S_R=H(p)+sum_lambda p_lambda[S(tau_R,lambda)+log d_lambda],
    S_out=H(p)+sum_lambda p_lambda[S(tau_out,lambda)+log d_lambda],
    S(D(rho))=H(p)+sum_lambda p_lambda S(tau_lambda).       (3)

Define the mutual information of the gauge-invariant regional algebras,
with the inherited unweighted block traces, by

    I_alg=H(p)+sum_lambda p_lambda I(tau_lambda).

The ambient source mutual information and center coherence obey

    C_center=D(rho||D(rho))=S(D(rho))-S(rho)>=0,
    I_ambient=I_alg+2 S_edge+C_center.                    (4)

The relative-entropy identity holds because log D(rho) is block diagonal.
Equation (4) includes all center weights, multiplicity correlations and
global center coherence. In particular I_ambient>=H(p)+2 S_edge.
The code checks the native five representation dimensions3,2,6,3,1 and
independent controls with nontrivial multiplicities and cross-sector
coherence. These controls verify the matrix identity; they are not the
selected spatial density.

Representation factors are not local gauge-invariant matrix algebras.
Neither 2 S_edge nor C_center is asserted to be distillable information
under arbitrary allowed local gauge-invariant operations. Joint
operations crossing the cut have a different observable domain.
The reference-inclusive quantity I(R:out,reference) is also different:
for the canonical pure input it equals2 S_R, not (4).

## 3. The actual magnetic action forces global center coherence

Let b be empty matter with constant links. It is a normalized physical
source vector, used solely for a matrix element. On the elementary
plaquette with one link in R, let L_(r,+),L_(r,-) be the ten native and
conjugate character loops. They are orthonormal, physical, and inside
every selected word cutoff. The
[native loop proof](../particle_primary_survival/PROOF.md) gives

    A_p b=(1/30) sum_(r,sign) L_(r,sign),
    H_magnetic,p=-(kappa/a) A_p,  kappa=1/7.              (5)

Matter terms annihilate b; its electric energy is zero. Different
plaquettes have zero overlap with these ten character states by Haar
integration on a link occurring in only one loop. Consequently the full
compressed physical action has the exact matrix elements

    <L_(r,sign),H b>=-gamma,
    gamma=kappa/(30a)=1/(210a).                          (6)

This is an identity of the full source, not an assertion that the eleven
vectors span an invariant subspace. Both endpoints lie in the actual
cutoff and Gauss space, so compression preserves this matrix element.

The blank has trivial boundary representation. Each loop has a distinct
nonzero hypercharge pair (+/-q_r,-/+q_r). Thus the full central projector
Pi_0 onto the trivial boundary representation fixes b and annihilates
these ten loops. Orthogonality and self-adjointness give

    ||[H,Pi_0]||_HS² >= 20 gamma².                       (7)

In an energy eigenbasis, (2) implies

    ||[rho,Pi_0]||_HS²
       >= (c beta_low)² ||[H,Pi_0]||_HS².               (8)

Let X=rho-D(rho). Its off-diagonal sector blocks contain those of
[rho,Pi_0], so ||X||_HS²>=||[rho,Pi_0]||_HS². For density matrices on
the finite physical support,

    D(rho||D(rho)) >= (1/2)||rho-D(rho)||_HS².            (9)

One direct proof applies the matrix-log Hessian along D(rho)+tX.
Its eigenvalues are at most one, so every divided difference of log is
at least one. The second derivative of the relative entropy is at
least ||X||_HS². Its value and first derivative at zero vanish;
integrating twice proves (9). The source floor makes the path positive.
This argument requires no approximation to the large density matrix.
Combining (7)-(9) gives the strictly positive original-state bound

    C_center >= 10 (c beta_low gamma)² > 0.              (10)

Thus the complete classical center record actually disturbs this global
source. Its regional marginal remains unchanged, as (3) shows. A
coherent record isometry preserves the full pure source-reference-record
entropy; classical dephasing and marginal entropy have distinct values.

## 4. Source work of the complete center record

The nonselective source output is D(rho). Its mean source work is
W=Tr H(D(rho)-rho). In a common energy basis of H and rho, put
S_ij=sum_lambda |(Pi_lambda)_ij|². This matrix is symmetric and doubly
stochastic. Therefore

    W=(1/2)sum_(i,j)(p_i-p_j)(E_j-E_i) S_ij
      >= (c beta_low/2)sum_lambda ||[H,Pi_lambda]||_HS².  (11)

The ten transitions in (6), their adjoints and their two sector
projectors contribute at least40 gamma² to the last sum. Hence

    W >= 20 c beta_low gamma² > 0.                      (12)

The full source-energy change is supplied by the declared measurement
apparatus; it is not a free nondisturbing readout. A complete two-energy
experiment has this mean because the input source density commutes H.
Its initial record can change reference correlations; those outputs
remain retained and are not identified with an undisturbed purification.

Pinching is a finite convex mixture of sector-phase unitaries. The
actual-state passivity and slope bounds in (1) apply to every such
unitary and give the exact entropy/work comparison

    beta_low W <= C_center+D(D(rho)||rho) <= beta W.       (13)

Here the middle expression is Tr[(D(rho)-rho)K]. The entropy increase
of the source is C_center; it is not the entropy of all retained records
and reference. The independent matrix controls check (11)-(13) for
strict positive polynomial states with noncommuting projectors.

## 5. Explicit information floors and the geometric boundary

The blank and ten character sectors each have source probability at
least c. Coarsening all other sectors into one of these labels gives
eleven probabilities at least c. Shannon entropy is concave, so its
minimum on this truncated simplex is attained at probabilities
(1-10c,c,...,c). Consequently

    H(p) >= -(1-10c)log(1-10c)-10c log c
          >= 10 N c log2,  c=2^(-N).                   (14)

The [actual representation bound](../particle_primary_geometric_cut/PROOF.md)
also gives S_edge>=4c log108. Thus

    I_alg >= 10N c log2,
    I_ambient >= 10N c log2+8c log108
                   +10(c beta_low/(210a))².            (15)

All exponents, slopes and scales come from the same member. These
finite lower bounds deteriorate with capacity. No regional probability
distribution has been replaced by equal weights. The ambient and
regional-algebra information differ by actual source terms, so an
Einstein entropy functional must specify which one it uses and account
for its retained records. Equations (4), (10), (12) and (15) supply no
proper area, universal Newton coefficient, modular diamond kernel or
matched small-ball error. They do not prove a continuum gravity limit.
