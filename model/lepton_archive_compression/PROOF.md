# Exact outcome-conditioned lepton archive compression

## Source and complete observation domain

The source is `model/lepton_cycle/instrument.py`, with its actual sixty
native address records on each of two independently prepared cells. The
family is C3 and can be entangled with an arbitrary untouched reference R.
The phase is theta=2/9. At n=1,2,4 rounds, each kick has angle 1/(7n), so
the total declared kick parameter is 1/7. Each round uses a fresh independent
source pair and the original four committed classical returns. All 4^n
words, with their original weights, are retained. Compression is applied
to the exported source after the original instrument, not to a selected
success channel.

The source archive occupies C3600 per round, embedded in C4096. The
compression domain consists of all states produced by this prepared
instrument for arbitrary family/reference input. It is not the whole
space of arbitrary unprepared C4096 inputs. The constructed transformations
extend to unitaries on their whole working spaces; their unused coordinates
are fixed blank only on the admitted span.

## A literal fixed local archive basis

The native record partition has fifty blank, five positive and five
negative addresses. Write their normalized characteristic vectors as
b,p,m, and define d=(p+m)/sqrt(2), c=(p-m)/sqrt(2). Let

    s=sqrt(5/6), t=sqrt(1/6), u=s b+t d, v=t b-s d.

Then u,v,c are orthonormal, u is the actual uniform native preparation,
and the right-source analyzer occupies span{u,v}. In the original outcome
order, the following columns are orthonormal archive bases:

    yes_yes: u tensor u;
    x_no:    (v tensor u, i c tensor u);
    y_no:    u tensor v;
    both_no: (v tensor v, i c tensor v).

Thus the original sixty-by-sixty arrays lie in dimensions 1,2,1,2. This
follows directly from the three native analyzer masks and the two uniform
reflections. The executable also compares every native source coordinate
and every family column with this expansion, rather than relying on a
numerically fitted low-rank approximation.

Set q_k=theta+2*pi*k/3, C_k=cos(q_k), S_k=sin(q_k),
lambda_k=sqrt(5/108)+2sqrt(5/216)C_k and z=exp(-i/(7n)).
The compact environment columns are

    a_yy(k) = 1+(z-1)lambda_k²,
    a_xn(k) = (z-1)lambda_k A_k,
    a_yn(k) = (z-1)lambda_k (5-sqrt(2)C_k)/(6sqrt(3)),
    a_nn(k) = (z-1)lambda_k B_k,

where the real two-vectors are

    A_k = ((5sqrt(6)C_k-sqrt(3))/18, sqrt(5)S_k/3),
    B_k = (-sqrt(15)(1+sqrt(2)C_k)/18, -S_k/3).

All three lambda_k have certified strictly positive rational intervals.
The y_no scalar is bounded away from zero because 5>sqrt(2). The yes_yes
modulus is at least 1-2 max(lambda_k²)>4/9. The declared nonzero kick
angles make z-1 nonzero. These diagonal family scalars therefore do not
change any environment rank.

## Exact ranks, including small history directions

For each of A and B, all three two-column determinants are nonzero. The
code prints their exact expressions in c=cos(theta), s=sin(theta), reduces
them modulo c²+s²=1, and encloses each in rational intervals excluding
zero. Certified sine Taylor remainders and integer-square-root enclosures
supply these intervals. This proves pairwise projective independence of
the three family columns in both rank-two outcomes.

For a history h, let m be its number of x_no or both_no outcomes. Its
environment columns are tensor products of the local columns. If m=0
their rank is one, and if m=1 it is two. For m>=2, select two rank-two
factors. Suppose sum_k alpha_k x_k tensor y_k=0. A covector annihilating
x_0 but neither x_1 nor x_2 leaves a combination of the independent
y_1,y_2. Hence alpha_1=alpha_2=0, and then alpha_0=0. Additional nonzero
tensor factors cannot reduce independence: apply dual functionals to an
independent factor to isolate each coefficient. The exact rank is
therefore three for every history with m>=2.

Consequently the rank-one, rank-two and rank-three history counts are
(2,2,0), (4,8,4) and (16,64,176) for one, two and four rounds. This is an
algebraic result, independent of the relative probability of a history or
the electron-family amplitude. A displayed floating eigenvalue-threshold
countercontrol misclassifies small but nonzero history directions; it is
never used to choose the compression rank.

For fixed history h, the reduced family channel is a Schur multiplier by
the conjugate Gram matrix G_h of its environment columns. Its Choi matrix
on span{|k,k>} is G_h transpose, so its Kraus rank is rank(G_h). Any
isometric environment realization preserving this whole outcome channel
needs at least that many quantum dimensions. The construction attains
this bound. The statement fixes the original classical history and does
not permit extra measured sublabels to hide a discarded quantum direction.

## Streaming isometries and retained references

Column normalization is used only to construct a well-conditioned
coordinate basis. Its removed complex family scalars are retained in the
subnormalized environment columns and original outcome probabilities. No
physical state is normalized or conditioned on success by this numerical
basis choice.

Rank-one outcomes multiply these scalars and add a known one-dimensional
archive factor. Rank-two outcomes append a normalized real two-vector.
The first such vector needs one archive qubit. A second is compressed by
a complete real orthogonal transformation on C4, with its three active
directions retained in a two-qubit register. Further rank-two outcomes use
a complete real orthogonal transformation on C8. Its active outputs occupy
coordinates 0,2,4, leaving the incoming working qubit exactly blank while
retaining the two archive qubits.

Each orthogonal transformation is constructed from the certified rank and
the actual compact source columns, completed on every unused coordinate,
and checked on every input column. There are four distinct two-factor,
eight three-factor and sixteen four-factor transformations. Rank-one
outcomes between these factors do not change the geometry. Their positions
and amplitudes remain in the original classical history and family scalar.
The full compression controller thus needs at most 28 small transformations
for the bounded four-round domain, in addition to the local source service.

For each history the receipt retains compact columns C_h and a decoder
E_h with E_h†E_h=I and E_h C_h equal to the complete tensor product of
local rank-two columns. The raw C3600 archive is reconstructed by applying
the stored local source-basis isometries at the original history positions,
including the known rank-one factors. This gives the entire vector, not
only its Gram matrix or diagonal probability law. The original ordered
history specifies the inverse placement; its order is not erased.

The identity on all three family input columns gives the identity on
every |k><l| operator. Tensoring with I_R proves equality for any untouched
quantum reference. In particular, every original subnormalized
family/reference output and the complete archive are related by the
declared isometry. The executable separately checks every history Gram
against the actual native source arrays, the complete trace-preserving
sum, and relative reconstruction errors on each small family column.

Two invalid replacements serve as controls. Collapsing rank-two or
rank-three histories to one environment vector preserves the original
classical weights but changes reference-sensitive off-diagonal entries.
Mapping nonorthogonal family-conditioned archive states to orthogonal
labels also changes their Gram matrix and cannot be an isometry. Their
complete reference-output trace distances are reported without dividing
out an outcome probability.

## Actual compression services and resource accounting

Three real Householder reflections map the original sixty-address frame
(u,v,c) to the first three coordinates of C64, acting unitarily on every
other native and padding coordinate. The left register additionally
phases its third coordinate by -i, to map i c to its real compact code.
The right register uses the same real frame. Conditional on the retained
return bits, a four-state permutation maps the two nonreturn coordinates
to one qubit, and an X gate maps the known right nonreturn coordinate to
blank. Every other address bit is blank on the admitted source span.

The Householder reflection for a real normalized q is implemented by a
native preparation word V and a blank phase:

    V exp(i*pi |0><0|) V† = I-2|q><q|.

Actual native preparation columns are replayed for all three reflections;
one complete Householder word is also replayed literally on all 64 input
columns. The all-input conjugation identity supplies the remaining large
block evaluation. The complete phase and routing words are replayed on
all their columns. All 28 real C4/C8 streaming matrices are compiled by
the existing native orthogonal compiler and replayed on every column.
Each primitive uses the same six-level, twelve-port source and supplied
coherent proper-code loading/export service as the bounded service compiler.

Exact target identities concern the ideal real-angle synthesis recipe.
The emitted binary64 angles approximate that recipe. Native replay errors
include angle representation, synthesis and matrix arithmetic; they are
diagnostics, not certified tolerance or physical gate accuracy. The raw
executed columns are never renormalized. Pulse and coherent-service counts
include the compression itself, with conditional worst-case totals for
every complete history. Original source preparation, measurement and
archive transport are additional services.

At four rounds the original padded archive uses 48 qubits. The compressed
quantum archive uses at most two qubits, alongside all eight original
classical record bits. Known blank working coordinates may be recycled by
the constructed unitary on this admitted span; an unknown source state
is not erased. Classical program storage, routing and analog control have
their stated supplied status. No equality between information dimension,
physical energy or thermodynamic work is claimed, and no autonomous reset
of the native source is derived.
