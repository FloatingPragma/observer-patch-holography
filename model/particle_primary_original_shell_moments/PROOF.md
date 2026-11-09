# Original density spectral moments and finite spatial energy contrast

The unchanged primary member-zero state has an exact four-shell first
spectral moment in its native density sector. Its neighboring shell is
strictly negative, its two non-neighboring shells vanish, and its original
Walsh record first moments are in the exact proportions 0:1:2:3. A complete finite
two-clock apparatus resolves the corresponding original mean contrast.
The spatially constant field is elastic, whereas every nonconstant Walsh
class has inelastic spectral weight. These are even density observables;
the odd G1/H particle-creation kernels have separate matrix elements.

## 1. Same selected source and complete density analyzer

Fix member zero of the [primary contract](../particle_common_contract/PROOF.md),
with a=1, V=8, twelve interior hopping bonds, all60 CAR modes per cell,
actual quantum links and Gauss projection, selected polynomial density rho,
idle reference and original source Hamiltonian H. Put

    n_x=N_x-30I,     N_x=sum_(all60 modes at x) c* c.                    (1)

Each n_x is Hermitian, gauge invariant, even, preserves every native
occupation and satisfies ||n_x||<=30. The eight n_x commute. The complete
local instrument has success n_x/30 and failure sqrt(I-(n_x/30)^2).
A coherent uniform address and complete Walsh rotation of success give

    K_q = (1/(30V)) sum_x (-1)^(q dot x) n_x,   q in {0,1}³,
    L_x = sqrt(I-(n_x/30)^2)/sqrt(V).                                  (2)

All eight successes and all eight original-address failures are retained.
Walsh unitarity proves
sum_q K_q* K_q+sum_x L_x* L_x=I on the full physical source. This is one
sixteen-record instrument, not eight normalized experiments. The source
input is rho; no prepared anchor is substituted. Apparatus operations are
even and preserve source native numbers. The source, address, clock and
reference outputs remain in the complete experiment.

For these scalar fields define the original spectral kernel

    M_xy(B)=sum_(i,k) 1_(E_k-E_i in B)
               Tr[rho Pi_i n_x Pi_k n_y Pi_i].                         (3)

Stationarity and the actual proper cube action are supplied by the source.
The [four-shell reduction](../particle_primary_coherent_odd_response/KERNEL_PROOF.md)
uses the same24 spatial rotations: M_xy=m_d for d=Hamming(x,y), with
real signed off-diagonal measures. Here there is no internal-component
sum and no odd creation increment. The full-menu normalization in (2)
is different from that of the trilinear apparatus.

## 2. Exact original off-diagonal first moments

The [native density Ward identity](../particle_primary_energy_tightness/PROOF.md#6-selected-hopping-energy-and-a-complete-density-sum-rule)
uses the literal Wilson matrices and exact CAR bilinear relations.
Electric, magnetic and onsite terms commute with n_x. Writing T_xy for
the complete Hermitian hopping bond, it gives on the full source

    [n_x,[H,n_y]] = T_xy                     if x,y are neighbors;
    [n_x,[H,n_y]] = 0                        if distinct non-neighbors;
    [n_x,[H,n_x]] = -sum_(y neighbor x) T_xy.                         (4)

The product link projector and Gauss projector commute with n_x, so the
same identity holds after the actual source compression. No independent
link-coordinate average replaces rho. The receipt reuses the exact
source-Wilson and quantum-link control, including all native blocks.

Let M_xy^(1)=integral omega dM_xy(omega). Source stationarity gives

    M_xy^(1)=Tr[rho n_x[H,n_y]],
    Tr[rho[n_x,[H,n_y]]]=M_xy^(1)+M_yx^(1).                           (5)

The spectral matrix is Hermitian. The actual cube group exchanges every
ordered vertex pair, so its component-scalar entries are real and
M_xy^(1)=M_yx^(1). All twelve interior bonds form one orbit; let

    t=Tr(rho T_xy)=Tr(rho H_hop)/12,     tau=-t.                       (6)

Equations (4)-(6) compute the complete shell first moments:

    (m0^(1),m1^(1),m2^(1),m3^(1))=(3tau/2,-tau/2,0,0),
    M^(1)=(tau/2)(3I-Adj_cube).                                      (7)

The same-source selected-energy inequality establishes tau>0. In detail,

    H=-720I+3D+H_electric+H_hop+H_magnetic,
    D>=0, H_electric>=0, H_magnetic>=-6/7,
    Tr(rho H)<=U0+1/16.

Here U0 is the authenticated physical ground trial upper bound; its vector
is only used for an operator-energy inequality. Consequently

    tau >= gamma := (-U0-1/16-720-6/7)/12 >0.2373.                    (8)

The owning density sum-rule producer supplies the exact outward rational
gamma. Every native Wilson bond has parent full-Fock norm at most
2 sum d=30: its one-body hopping singular values are one in2d columns
per native block. Compression does not enlarge the norm. Thus tau<=30.
The receipt gives the corresponding signed intervals for all four entries
in (7), including a strictly negative original neighboring first moment.
These are expectations in rho, not moments in a trial state.

## 3. Exact Walsh moments and a genuine energy-shape distinction

The normalized Walsh vectors f_q(x)=(-1)^(q dot x)/sqrt(V) diagonalize
3I-Adj_cube with eigenvalue2s for s=|q|. Therefore

    integral omega dmu_q(omega)=s tau/(900V)=s tau/7200,
    multiplicities(s=0,1,2,3)=(1,3,3,1).                             (9)

The factor900V follows from K_q=(sum_x f_q(x)n_x)/(30sqrt(V)). It
differs from the divisors6 and48 of the odd trilinear instrument. Summing
all eight first moments gives tau/600. An independent exact eight-state
cube control evaluates all64 spectral entries directly in an energy
eigenbasis and checks (7)-(9). Its density is not substituted for rho.

The constant row is K_0=N_total/(30V)-I. Since [H,N_total]=0, its
ideal transfer measure is supported exactly at zero. It is not an absent
row: the authenticated physical floor rho>=cP_G gives

    Tr(rho K_0²)>=c>0,                                               (10)

using the physical empty-Fock, constant-link vector on which K_0=-1.
The code binds the exact factored c from the same polynomial preparation.
For each s>0, (8)-(9) give a positive original signed first moment and
hence positive spectral weight at strictly positive transfer. Thus its
normalized ideal energy shape differs from the constant row's point
measure at zero. The bounds do not establish different normalized shapes
among s=1,2,3; differing intensities can contribute to their moment ratios.
No spin, conserved translation momentum, relativistic shell or particle
mass is assigned to these density waves.

## 4. Complete finite clocks retain the original mean contrast

Use two copies of the same [source clock](../particle_primary_joint_spectroscopy/PROOF.md)
before and after (2), retaining every energy-bin, success and failure
record. Let R bound ||H|| and use its alias-free tick tau_clock=1/(4R).
The signed bin energies have absolute value at most pi/tau_clock<14R;
their reported difference is bounded by28R. True transfers lie in[-2R,2R].
If each exact clock has tail at most b outside resolution eta, the ideal
reported-versus-true first-moment error is at most

    2eta + 30R(2b).

If the complete implemented two-clock law has half-trace distance at
most e from that ideal finite-clock law, an event's signed energy
variable has range at most56R, giving additional error at most56R e.
The conservative combined bound used in the receipt is

    E_mean = 2eta+56R(2b+e).                                        (11)

This applies to each unnormalized success-row moment; all other records
are assigned zero in that event statistic. There is no division by its
success probability. It is uniform in the untouched reference.

Put mu_star=gamma/7200. Choose eta=mu_star/16 and each tail and one-clock
implementation allowance mu_star/(1024R). The actual Cayley/Fejer planner
supplies the full finite clock size, substeps and attained errors. Since
2b+e<=4mu_star/(1024R),

    E_mean <= 11mu_star/32 < mu_star/2.                              (12)

Every nonconstant reported row therefore has positive original signed
mean, at least s mu_star-E_mean. Adjacent nonzero shell means differ by
at least mu_star-2E_mean>0. These are expected statistics of the complete
implemented source instrument; a full histogram or a finite-shot estimate
is not claimed.

For q=0 the reported signed mean is exactly zero, rather than bounded
only by (11). K_0 commutes with both identical finite clock kernels, so
the joint law of initial and final clock bins is symmetric under their
exchange. Their difference has zero mean, including the finite Cayley
implementation. This argument does not assume that an initial energy
record leaves the source-reference purification unchanged.

The complete experiment therefore resolves an original spatial energy
contrast on a specified finite interacting source. Its even density
sector and normalization differ from the odd G1/H creation observables.
The calculation supplies no trilinear shell value, isolated pole,
refinement-uniform residue, continuum field theorem or GeV mass.

## 5. Strict hopping sign and original shells in every native block

Let r be one of the five native representations, with dimensions
(3,2,6,3,1) and Q6=(-4,-3,1,2,6), respectively. Write T_r for its
complete twelve-bond hopping operator. The source cube is bipartite, so
the physical, even Hermitian unitary

    U_r=(-1)^(sum_(black x) N_(r,x))

preserves every source term other than T_r and reverses each of its bonds:

    U_r H U_r=H-2T_r,      ||[H,U_r]||_HS²=4||T_r||_HS².             (13)

The product link cutoff commutes with U_r; the native number at a cell
commutes with its Gauss action. Thus (13) holds on the full physical
source. It is not an assertion about an invariant compressed trial space.
The Hilbert-Schmidt norms use the entire finite physical Hilbert space.

For the actual selected density rho=f(H), the
[polynomial modular derivative](../particle_primary_gravity_premises/PROOF.md)
obeys -f'(E)>=c beta_low on the source spectral interval. Here beta_low>0
is its authenticated logarithmic slope, and the same factored c used in
(10) bounds f from below throughout that interval. Its scalar bound uses
the polynomial normalization and source dimension; no Gibbs replacement
is made. For any physical Hermitian unitary U, the exact spectral sum is

    W(U)=Tr[rho(UHU-H)]
        =(1/2) sum_(i,j) (f(E_i)-f(E_j))(E_j-E_i)|U_ij|²
        >=(c beta_low/2)||[H,U]||_HS².                              (14)

The filled-beta, constant-link anchor vector Omega is physical.
It is used only as an operator-norm witness. The
[actual native Haar/CAR first-action calculation](../particle_primary_sector_threshold/PROOF.md)
gives, separately in every native block,

    <Omega,T_r Omega>=0,       ||T_r Omega||²=12d_r.

All three actual Wilson matrices are evaluated; every one of the twelve
bond packets has squared norm d_r. Different bonds have orthogonal
nonconstant link functions under the source Haar inner product, and the
degree-one first actions are admitted by the fixed link cutoff. With
chi=T_r Omega/sqrt(12d_r), the two orthogonal physical basis vectors
Omega and chi contribute both conjugate matrix entries to the full norm:

    ||T_r||_HS²>=24d_r,
    tau_r:=-Tr(rho T_r)/12 >=2d_r c beta_low >0.                    (15)

The receipt executes the native first-action calculation with zero
occupation increment and verifies the twelve actual bipartite bonds.
An independent full-Fock polynomial-state control verifies (13)-(14).
The extremely small lower bound in (15) is retained in factored form.
It is an exact source consequence, not a useful finite-shot precision.

For n_(r,x)=N_(r,x)-2d_r I, ||n_(r,x)||<=2d_r, the identical Ward proof
gives M_r^(1)=tau_r(3I-Adj_cube)/2. A complete sixteen-record analyzer
with success operators sum_x (-1)^(q dot x)n_(r,x)/(2d_r V) has original
first moments

    integral omega dmu_(r,q)=s tau_r/(4d_r² V),
    s=|q|,       s=1 lower bound=c beta_low/(2d_r V).                (16)

Each native density channel therefore has a strictly nonconstant
original shell measure. These even, occupation-preserving probes do not
create the fermions associated with the native internal weights. The
total-density clock bound in Section4 does not resolve the much smaller
individual bounds (16), and no such precision is claimed.

## 6. Complete density-event weight and conditional signed mean

Let p_x count the30 positive-beta particles and h_x the30 negative-beta
holes at x. Then n_x=p_x-h_x and D_x=p_x+h_x are commuting count
operators, with |n_x|<=30 and |n_x|<=D_x. Consequently, on the full Fock
space and its physical restriction,

    n_x²<=30D_x,
    S:=sum_q K_q* K_q=sum_x n_x²/(900V)<=D/(30V).                   (17)

The source-selected [occupation bound](../particle_primary_occupation_bound/PROOF.md)
supplies <D><4.992. All eight density successes therefore have combined
original probability p_success<=D_upper/(30V)<0.0208 at member zero.
The initial finite clock is a complete function-of-H instrument and
rho commutes with H, so its unrecorded source marginal remains rho;
the final complete clock is trace preserving. The same total success
probability applies with both actual clocks, including their failures.
No assertion of unchanged source-reference purification is needed.

The original sum of signed means is tau/600. Since the constant row has
exact reported mean zero, the total reported sum differs from its ideal
value by at most7E_mean. A separate conditional statistic can therefore
be bounded without replacing any original record weight:

    E[omega | density success] >= (gamma/600)/(D_upper/(30V)),
    E[reported omega | density success]
       >= (gamma/600-7E_mean)/(D_upper/(30V)).                       (18)

The receipt gives bounds greater than0.019 and0.015, respectively, in
the source energy convention. Conditioning is explicit in (18). These
signed transfer means are neither a sharp spectral line nor a rest mass.
The inequalities concern the same even density apparatus as (2).
