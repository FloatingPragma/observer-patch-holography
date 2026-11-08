# Native two-body contact response

The fair twelve-port native source has a primitive two-body response below
the threshold of its centered independent-pair comparison. Complete native
repair records observe this mode with strictly positive weight. Connecting
the cell to the native federation produces a capacity-independent leakage
in the same exchange-action convention. These are finite correlation and
attachment results. A localized particle, scattering states, a relativistic
mass shell and a physical clock require separate constructions.

## Source, reference and response

Use the actual thirty icosahedral seams, binary total occupation two,
independent uniform seam selection, a fresh fair placement coin and repair
probability one. The declared C/E mixture has equal weights. The occupation
transfer is the same under C and E. The working reference is the uniform
66-state mixture `tau=I66/66`. This is the fair local source admitted by
`native_program.interpreter.program(biased=False)`, rather than a change to
the biased source's parameters without declaration.

For a graph with V ports and M seams let `G1` be its graph Laplacian and

    G2 = sum_e (I-Swap_e),       P2 = I-G2/(2M).                (1)

The swaps act on unordered distinct pairs. Equation (1) follows by summing
the literal fair placement branches, including unsuccessful movement.
`native_effective_fields` derives the hard-core expression and contact
coefficient; `native_correlation_sectors` derives its primitive decomposition;
`native_pair_response` derives the complete-record comparison to independent
pairs killed upon collision. Those results retain their source hypotheses.

Let D send pair coefficients to their sums at each port. The primitive
two-body space is `W2=ker D`, with orthogonal projector

    Q=I-D*(DD*)^-1 D,       DD*=(V-2)I+11*,

and dimension `binomial(V,2)-V`. It reduces G2. The fixed actual probe uses
seams `(0,1)` and `(2,4)`:

    phi=(n0-n1)(n2-n4),       h01,24 in W2,       ||h||²=4.

Two initial E reads score their *outgoing* contrasts `d-2j`; two final E
reads score incoming contrasts d. With n unread uniform attempts between
them, their complete scored expectation is

    C(n)=(2/33) <h,P2^n h>/4.                                 (2)

Each E command retains all 25 raw slots; each pair retains 625 slots,
including 589 zero slots. A uniform unread attempt has 721 slots. The
original `interpreter._outcomes` is replayed to check (2), with all source
weights and backaction. Signed summaries are linear functionals of this
complete instrument, rather than replacements for its retained branches.
The observer stores one FIFO record and exports earlier records. Eight
attempts require at most seven exported records. The expanded interpreter's
resource-stop outcome is explicitly tested; exact factorization of scored
records does not assert that an exponentially expanded run completes.

The exact local one-body gap is `(5-sqrt(5))/60>1/30` per attempt. The
existing exclusion-gap and polynomial-preparation theorem therefore gives
480 unconditioned warmup attempts for normalized signed-score error at most
`33/65536`, uniformly over arbitrary initial source/reference states in
this fixed occupation sector. Warmup records are retained and ignored in
this comparison. This is a signed reference-valued score bound, not a
full-state or preparation-archive independence statement. Warmup, lag four
and all four probe reads fit an admitted 488-attempt program.

## The independent-pair comparison and its contact equation

Embed unordered-pair coefficients isometrically into the symmetric tensor
square by

    J e_{ij}=(e_i tensor e_j+e_j tensor e_i)/sqrt(2), i<j.

The full independent generator is `L=G1 tensor I+I tensor G1`.
Compression retains its diagonal loss into coincidence states:

    J*LJ=G2+2A,                                               (3)

where A is the diagonal projector onto adjacent pairs. Deleting 2A from
the right side gives the actual exchange generator. Thus its coefficient
is fixed by the native source; no contact coupling is fitted or selected
after inspecting the spectrum. The compressed independent process retains
collision as an absorbing comparison outcome, not a postselected native
outcome.

For h in W2, every row sum of Jh vanishes. Each tensor factor is therefore
orthogonal to the constant vector. If gamma1 is the one-body spectral gap,

    H=G2|W2,    H0=H+2QAQ|W2,    H0 >= 2 gamma1 I.             (4)

This comparison is between operators on a specified coefficient space.
It does not identify the occupation variables with physical Bose or Fermi
particles. In particular `2 gamma1` is a finite variational threshold,
not an infinite-volume scattering threshold.

Let E_A insert the adjacent-pair coordinates and U=Q E_A. For z outside
the spectra involved, put `R0(z)=(H0-z)^-1`. Direct multiplication gives

    (H-z)^-1 = R0 + 2R0 U (I-2U*R0U)^-1 U*R0.                 (5)

For real z below `2 gamma1`, H0-z is positive. An eigenvector of H at z
has nonzero contact component: otherwise (4) would contradict its energy.
Equation `(H0-z)f=2UU*f` then identifies it with the kernel of
`I-2U*R0U`, preserving multiplicity. Conversely every vector in that kernel
produces an H eigenvector by R0U. This proves the finite contact criterion,
including degenerate eigenvalues. The inverse API rejects both reference
and correlated poles. It makes no claim about boundary values of an
infinite-volume resolvent or an on-shell scattering matrix.

## Exact visible pole

Exact integer matrix arithmetic gives

    det(x-G1)=x (x-6)^5 (x²-10x+20)^3,
    det(x-H)=(x-12)^5 (x-10)^11 (x-8)^5
      * (x²-20x+88) (x²-20x+90)^4 (x²-20x+98)^4
      * (x³-30x²+280x-776)^5.                                (6)

The 66-state characteristic polynomial is the product of these two
polynomials. The code also checks Q, its commutation with G2, and its
rank 54 independently. The lowest root E of

    c(x)=x³-30x²+280x-776

lies in `(5,51/10)`, while its other roots lie in `(11,12)` and `(13,14)`.
Every other primitive factor has its lowest root above six. Since
`sqrt(5)<9/4`,

    2 gamma1=2(5-sqrt(5))>11/2,    2 gamma1-E>2/5.             (7)

The mode has multiplicity five. The 12 Krylov columns
`h,G2 h,...,G2^11 h` are independent. Their annihilator d(x) is the product
of the seven distinct factors in the second line of (6). Thus the chosen
native probe observes every one of its 12 distinct primitive poles.

For `m_j=<h,G2^j h>/4`, the exact scalar response is

    R(z)=<h,(z-G2)^-1 h>/4=N(z)/d(z),
    N_j=sum_(k=0..j) d_k m_(j-k), j=0,...,11,                 (8)

where d_k are the descending coefficients of the monic degree-12 d.
The inverse in (8) is on the cyclic primitive space. It can be finite at
an eigenvalue outside that space, including zero and the one-body gap.
The code checks the polynomial identity, all moments and `gcd(N,d)=1`.
Self-adjointness makes every residue a positive spectral weight. At E,
exact rational interval evaluation of `N(E)/d'(E)` gives

    E approximately 5.01827276928,
    1/100 < w_E approximately 0.0194379829330 < 1/50.         (9)

The decimals illustrate certified rational intervals. The raw observable
weight is `(2/33)w_E`, and the actual attempt-transfer eigenvalue is
`1-E/60`. Its logarithmic decay rate and E are distinct clock conventions.
Neither is a measured particle mass. No success conditioning, independent
oscillator assumption or fitted spectral list enters (8).

## What connection to the federation does

Embed a pair coefficient supported inside one original cell into the
connected fair federation. Each cell has three distinct boundary ports,
each joined to a distinct outside neighbor. The all-level combinatorial
attachment uses the routing postconditions in `native_spatial_scaling`;
literal constructor and host evaluation retain their own level limits.
The implementation checks every cell at evaluated levels zero, one and two.

Let P_cell select pairs wholly inside the cell and let f be a normalized
vector in the 15-dimensional sum of the three cubic eigenspaces in (6).
For boundary-port set B the outside images of distinct boundary moves are
orthogonal pair coordinates, so

    ||(I-P_cell)G_federation J_cell f||²
       = sum_(i<j) |f_ij|² |{i,j} intersect B| > 1/100.       (10)

This inequality is exact. An integer basis N of `ker c(G2)` satisfies
`QN=N`. For each of the 220 possible three-port sets, all 15 leading
principal minors of

    N* [100 diag(|pair intersect B|)-I] N

are positive by fraction-free integer elimination. Sylvester's criterion
proves (10) for the entire cubic subspace, not only the displayed probe or
sampled floating eigenvectors. The connected graph's lower primitive map
also annihilates J_cell f, so its occupation and correlation-degree domain
have not changed.

For a normalized local eigenvector at E, the residual of its unchanged
embedding is at least 1/10 for every successful federation attachment.
This disproves that particular isolated-cell eigenvector identification
as a growing-source eigenvector. It does not exclude a dressed or spatially
extended eigenvector, or a resonance with a separately proved limit.

On q=2 and maximum port degree six, `0<=G_federation<=24I`. The scalar
Taylor remainder for its heat contraction is at most `288 t²`. Therefore

    ||exp(-tG_federation)J_cell f-J_cell exp(-tE)f||
       >= max(0,t/10-288t²),       t>=0.                     (11)

Projecting outside the cell removes the entire reference exponential.
Only the source exponential's remainder is needed in this lower bound;
no remainder for the reference exponential is discarded.

The existing all-sector native endpoint bound gives a further subtraction
`3/(5n)` for `P_federation^n`, at `t=n/(2M)`, n>=1. This is an operator
comparison in the transfer Hilbert space, rather than an executed coherent
evolution or a complete-record channel-distance bound. At fixed positive
t small enough, its positive lower bound persists as M grows. Increasing
capacity does not repair the unchanged cell embedding in this clock.

## Conditions for a persistent observable sector

A sequence of isolated-cell poles alone supplies no growing-source particle.
A candidate on the same connected source must specify normalized dressed
vectors or an isometry, a common counted-time convention, and residuals
`||(G_l-E_l)f_l||<=eta_l`. The spectral theorem gives

    mu_f_l({|x-E_l|>=delta_l}) <= eta_l²/delta_l².              (12)

Thus concentration requires `eta_l/delta_l ->0`, with an actual separated
spectral window, nonvanishing accessible residue and a proved localization
or propagation map. Equation (12) alone does not prove any of those inputs,
a simple eigenvalue, or an infinite-volume mass shell. The code's budget
helper consumes independently proved inputs, rather than certifying them.

The unnormalized readout cannot be omitted from this criterion. In a q=2
V-port reference the same four-E contrast has variance
`v_V=4/binomial(V,2)=8/[V(V-1)]`. A normalized weight of order one can have
a raw probability-score weight vanishing as V^-2. For independent complete
trials with scores in [-1,1], Hoeffding gives normalized error epsilon with
failure probability alpha once

    samples >= 2 log(2/alpha)/(v_V² epsilon²).                (13)

This sufficient budget grows as V^4 at fixed epsilon and alpha. Independent
repreparation is an explicit experimental contract; the bound is not applied
to correlated successive records as if they were independent shots. Native
warmup error, implementation error and sampling error must be retained with
their own normalization. The finite record and source-action calculation
proves the correlated mode (7)--(9), while (10)--(13) specify why it does not
by itself establish a capacity-stable physical particle.
