# Original energy-current records from a finite local control

The selected physical source has a strictly nonzero original response of
its conserved local energy current to an explicit finite control pulse.
The pulse, complete readout and probability bound use the same selected
polynomial state, idle reference and Hamiltonian clock. A source-relative
evolution derivative supplies the integrated energy generator. Identifying
this control with a proper lapse or deriving a Lorentz stress tensor requires
additional geometric data.

## 1. The fixed source and its energy generator

Use the [primary contract](../particle_common_contract/PROOF.md) and the
[cube-covariant energy allocation](../particle_primary_covariant_stress/PROOF.md).
On physical Gauss support, H=sum_x h'_x. At the corner x=0 put

    A=h'_0,
    B=i[H,A]=-sum_y J'_0y,    J'_xy=i[h'_x,h'_y].             (1)

Thus B is the actual net incoming energy current, not an independently
chosen stress observable. The selected rho=p(H) obeys Tr(rho B)=0.
The code reconstructs ||H||<=R and ||A||<=K from the source and its
actual corner terms. Set C=2RK and G=R+K, so ||B||<=C and ||H+A||<=G.
All terms and instruments preserve physical Gauss support. For a scalar
control f_x(t), use H_f=sum_x f_x h'_x, e_x=f_x h'_x and
J_xy[f]=f_x f_y J'_xy. The exact controlled balance is

    d<e_x>/dt + sum_y <J_xy[f]> = dot(f_x)<h'_x>.            (1a)

The applied-work term must be retained. Integrating it into the endpoint
energy increment supplies the hypotheses of the finite
[energy Ward identities](../../Lean/QFT/CoarseEnergyWard.lean), including
arbitrary spatial coarsening. Differentiating these finite identities at
f_x=1+lambda g(t) delta_x0 preserves the explicit control-work term.
No homogeneous Ward identity for the unperturbed current is asserted
during a pulse, and no limiting stress tensor is assumed.

## 2. A physical witness for the full current

Let Omega fill all negative-beta rows at every cell, leave positive rows
empty and take constant links. Each site is an even gauge singlet. At
x=0 move one charge-six fermion from its negative-beta row to the matched
positive-beta row. The normalized Slater vector phi preserves every
native occupation and every Gauss constraint. It is an analysis witness,
not a replacement preparation.

Project [H,A]phi to states with exactly one nonconstant native link word,
with the matter profile reached by one hop. Only the commutators of onsite
and electric terms with hopping contribute to this projection. Two hops
preserve the matter parity of one bipartite sublattice, whereas one hop
reverses it. Magnetic multiplication leaves at least three other
nonconstant link words when multiplied by a one-link hop. Pure magnetic
or electric/magnetic terms retain the original matter profile. These
selection statements apply after product-cutoff compression: a projector
cannot turn a nonconstant word on an untouched link into a constant.

For a bond from the corner to a neighbor, its share in A is1/2. If the
hop changes the two onsite energies by dM_x,dM_y and adds electric energy
E, the projected commutator coefficient multiplying the hopping amplitude
is exactly

    (dM_x+dM_y+E)/2 - (dM_x+E/2) = (dM_y-dM_x)/2.          (2)

All other native blocks give zero in (2). In the charge-six block, moving
the positive particle away or filling the negative hole from the neighbor
gives dM_x=-3/a and dM_y=3/a. The diagonal Wilson coefficients have
magnitude1/(2a). Each of the three corner bonds therefore supplies two
orthogonal physical one-word components with magnitude3/(2a^2). Haar norm
of the one-dimensional character is one. Consequently

    ||[H,A]phi||^2 >= ||P_one[H,A]phi||^2 =27/(2a^4).       (3)

The code enumerates all120 modes in a two-cell chart, all five native
representations and all three literal Wilson matrices. Filled spectator
cells have even parity, so this chart embeds with the same signs at every
selected side and cutoff degree. A separate68-dimensional charge-six
Gauss control evaluates the full commutator as matrices and confirms9/2
per axis at a=1. That control is not substituted for the full source.
The one-word output in (3) is admitted whenever the source degree is at
least one. No invariant subspace or complete histogram is asserted.

## 3. An original-state positive current susceptibility

Write rho>=c P_G, c=2^(-N), using the authenticated physical density floor.
The exact polynomial logarithmic derivative obeys
-beta<=p'(E)/p(E)<=-beta_low<0 on the source spectral interval. In an
energy basis, ordering E_i<E_j gives

    p_i-p_j >= c beta_low (E_j-E_i).

Thus the actual original double commutator has the lower bound

    q=Tr rho [A,[H,A]]
      =2 sum_(i<j)(p_i-p_j)(E_j-E_i)|A_ij|^2
      >=c beta_low ||[H,A]||_HS^2
      >=c beta_low 27/(2a^4)=:q0>0.                        (4)

Degenerate energies contribute zero. The Hilbert-Schmidt norm is only
used as an inequality on the finite physical carrier; no large matrix
is allocated. Equation(4) transfers a positive operator norm estimate to
the original state. It transfers no signed expectation of the witness.
The original susceptibility is capacity dependent; no uniform lower
bound is inferred from the tiny source faithfulness floor.

## 4. An exact finite pulse and a complete current detector

For a source duration tau, apply the declared local scalar-energy control
H+A, then return the generator to H. This is the corner setting f_0=2,
f_x=1 elsewhere in the existing H_f=sum_x f_x h'_x coupling. The source
program, initial density, reference and duration convention are unchanged;
the control is an intervention, not a new state-selection law or a proved
proper metric. A complete reference-inclusive implementation uses the
exact finite Cayley unitary

    Q_tau=(I+i tau(H+A)/2)^(-1)(I-i tau(H+A)/2).             (5)

Measure both effects (I+B/C)/2 and (I-B/C)/2, with their positive-root
Kraus operators and numeric labels +C,-C. The intervention and both
outcomes, including their retained quantum outputs, are part of the
experiment. No success conditioning is used. In the no-pulse experiment
the probabilities are exactly1/2 by stationarity.

For the ideal e^(-it(H+A)) pulse, the mean b(t) obeys b(0)=0,
b'(0)=i Tr rho[A,B]=-q, and |b''(t)|<=4G^2 C. Hence

    b(t)<=-t q0+2t^2 G^2 C.                               (6)

Functional calculus gives ||Q_t-e^(-it(H+A))||<=t^3G^3/12.
Its current-mean error is at most C t^3G^3/6. Choose the explicitly
source-derived finite duration

    tau=q0/(16G^2 C).                                    (7)

The receipts verify q0<1 and R,K,C,G>=1 in the declared source units.
The quadratic term in (6) is tau q0/8. The Cayley error divided by
tau q0 is q0/(1536 G C)<=1/1536. Therefore

    Tr(Q_tau rho Q_tau* B) <=-tau q0/2,
    TV(binary pulse law,binary no-pulse law)
                              >=tau q0/(4C)=:d0>0.       (8)

All tiny numbers are stored as rational multipliers times powers of two;
none is rounded to zero. The exact unitary and square-root services are
finite apparatus specifications, not a claimed execution of the large
source matrix. The current record has no energy-transfer or rest-mass
interpretation.

The [source modular-work bound](../particle_primary_gravity_premises/PROOF.md)
applies to this exact Q_tau. Data processing to the complete binary
record and Pinsker's inequality give

    beta W(Q_tau)>=D(Q_tau rho Q_tau* ||rho)>=2d0^2,
    W(Q_tau)>=2d0^2/beta>0.                              (9)

This is cyclic work relative to the unchanged H. Source marginal relative
entropy is used, not the generally singular relative entropy between two
retained pure reference states.

## 5. Source-relative evolution and its geometric boundary

For a bounded real pulse g(t) on[0,T], let U_lambda solve the exact finite
Schrodinger equation with H+lambda g(t)A and U_lambda(0)=I. Define
W_lambda(T)=e^(iHT)U_lambda(T). Duhamel's formula on the actual carrier gives

    d_lambda W_lambda(T)|_0=-i int_0^T g(t) alpha_t(A) dt,
    d_lambda[W_lambda(T)* O W_lambda(T)]|_0
                       =i[int_0^T g(t) alpha_t(A)dt,O].    (10)

This is a finite source-relative response, with the original H on both
sides. Functional calculus keeps the reference idle. The integrated
energy in (10), its full boundary current, the applied-work term in (1a)
and their coarse sums instantiate the finite Ward relation. The existing local commutator bounds control
responses of remote observables. They supply no Lorentz boost or proper
clock identification.

The locally covariant field-theory notion of metric response in
[Brunetti, Fredenhagen and Verch](https://arxiv.org/abs/math-ph/0112041)
requires spacetime embeddings, a time-slice structure and an actual metric
perturbation. Equation(10) supplies the finite dynamical derivative on
this selected source. It does not identify the control with that geometric
relative Cauchy evolution, supply the proper-area law, or prove Einstein
dynamics. The finite [path time-slice interface](../../Lean/QFT/PathTimeSliceInterface.lean)
is also a separate source domain; no unchanged formal inheritance claim
is made for its full geometric interpretation.
