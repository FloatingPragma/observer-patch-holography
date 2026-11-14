# Finite current attachment to the common native walk

This package extends the single proposed `native-local12-current-v0-proposed`
source by coherent modular current memory. It does not install a second
Hamiltonian, thermal scheduler, Coulomb energy, or physical time calibration.
The geometry, reflection, native event, record recoder, and transport are the
executed DYN01 program in `../native_walk_bridge`. The extension and its joint
preparation are declared architecture. No official OpenWave status is changed.

## Current and Gauss constraint

Orient the actual 288 seams of the 96-port size-two graph from their lower to
their higher numbered endpoint. Let B have entries +1 at the lower endpoint
and -1 at the higher. For a one-load position x, write n=delta_x. A source
transport x to y on seam e has the literal native current

    j = 1[x=a_e] - 1[y=a_e].

The unchanged local coin reflection R acts on the event register. Extend the
actual reversible transport by

    S_F |x,event,F> = |y,event,F+j epsilon_e>.

Here every edge field takes values in Z/6Z. Source or event padding has j=0;
unused three-bit field labels 6 and 7 are fixed by the controlled addition.
Reversing the same matching sends y back to x and changes j to -j. Thus S_F
is an involution on every basis column, including the declared complements.
On the admitted modular field domain,

    delta_y + B(F+j epsilon_e) = delta_x + BF  (mod 6).

This proves unitary current transport and Gauss preservation for every F,
every Gauss sector, and any untouched quantum reference. It does not assume
that the input is a classical position or a single field configuration.
`CURRENT_WALK.json` independently checks all 1,152 active source/event atoms,
all 2,304 C/E terminal native labels, and the local modular identity at
moduli 3, 5, and 6. Modulus 6 is the shared implementation in this package.

The native event register remains coherent during repeated W_F=S_F R.
Measuring every j would define another instrument. The inverse record recoder
still supplies the original C/E labels at the declared terminal intervention;
the reused coherent coin is not a fresh uniform attempted native event.

## A normalized affine Gauss preparation

Fix G=delta_0. Choose the deterministic spanning tree in `geometry()` and
let P_x be the signed unit flow along its root-to-x path. Then
BP_x=G-delta_x. There are 288-96+1=193 chords. For chord e=(a,b), define

    C_e = epsilon_e + P_a - P_b.

BC=0, and the rows of C at the chords form the identity. Every integral
divergence-free flow is uniquely C z: subtract the chord values, leaving a
divergence-free flow on a tree; leaf elimination makes that remainder zero.
This argument also works over Z/kZ for composite k, with no division.
Consequently the finite modular fibre has exactly 6^193 elements and the map

    J_G |x,event> = 6^(-193/2) sum_z |x,event,P_x+Cz>

is an isometry. This is a coherent, position-correlated preparation. It is
not a product of the position state and an independent definite edge tape.
For every source atom, P_x+j epsilon_e-P_y=C h. Here h is zero on tree edges
and is j in the traversed chord coordinate otherwise. Translation permutes
the finite sum, so

    S_F J_G = J_G S,   R J_G = J_G R,   W_F J_G = J_G W.

The source-word initial cell preparation is mapped through this same J_G.
All logical observable laws transfer through the isometry. Physical source
positions have the same probabilities. Tracing out the physical edge field
does **not** recover the bare source density matrix: fibres at distinct x
have different divergence and are orthogonal. Source off-diagonal operators
must be dressed by the corresponding field translation. This distinction is
essential even in the zero-character preparation.

There is no normalizable uniform state on the infinite integer fibre
Z^193. The exact isometry above is finite and modular. For an integer cube
packet [-K,K]^193, a shift in one coordinate has overlap 1-1/(2K+1), hence
state-vector difference sqrt(2/(2K+1)). Operator norms telescope to the
T-step comparison bound min(2,T sqrt(2/(2K+1))). The same proof is valid with
an arbitrary reference. A finite implementation of that integer packet
requires a margin of at least T in every coordinate to avoid wrap during T
steps. It is a separately specified approximate preparation, not the finite
exact Z6 state.

A modular current register aliases winding 0 with winding 6. It is therefore
not an unbounded integer history. In the admitted known-G one-load sector,
regional charge G_R-n_R lies in {-1,0,1}; labels 5,0,1 decode it injectively.
That bounded charge statement does not reconstruct every integer current
representative or permit unlimited histories to be unwrapped without
additional records.

## Field-sensitive records in the same evolution

Replace the uniform cycle coefficients by a finite character
exp(2 pi i alpha.z/6). Under a current displacement h, the source transport
acquires exp(-2 pi i alpha.h/6). This follows by changing the summation index;
it does not add a new phase gate to W_F. Character zero gives the original
logical walk. A nonzero character is a different field preparation under the
same augmented dynamics.

The executable retains one explicit chord coordinate and puts the other 192
in their exact uniform states. The retained-array dimension is 128*16*8,
including source, event, and field padding. This is an invariant compression
of the finite Gauss sector. It is not an execution of the full 6^288 field
array and not a finite-field approximation to an integer uniform state.

All six characters are propagated with the actual compiled R and S for eight
development steps. A definite retained cycle coordinate yields the same
position law as the equal mixture of the six character laws. Both this law
and the character-one law differ from the zero-character law. The receipt
retains complete padded position probabilities and original norms without
success normalization. A reference Gram check uses arbitrary amplitudes on
all source, event, and retained-field padding labels.

These are field-dependent interference records and reciprocal conjugate
current shifts. The character is conserved by this evolution. No autonomous
field energy, propagating transverse mode, Maxwell limit, or force law is
proved here. Such an extension must amend this same common program and its
state, rather than import the earlier thermal feedback comparison.

## Native circuit and complete-outcome scope

`compiled_retained_increment()` compiles the two controlled six-cycles using
the existing phase-exact M6 `flip` and `ccx` source words and clean reversible
control ladders. All 16,384 source/retained-field clean columns, including
labels 6 and 7, are executed exactly as a gate word. The inverse word is
checked. Every gate has a complete unitary action on dirty scratch, but the
advertised current action requires its declared clean scratch sector.

`CURRENT_PROGRAM.json` supplies this word, the native primitive pulse words,
and all 288 full-field block addresses. For the physical edge representation,
apply each edge's two oriented controlled additions before the unchanged S.
They commute because source controls are unchanged and field targets are
distinct. The full register allocation is the existing 23 work bits, 864
edge-field bits, and 11 extra clean ladder bits. The block generator emits
the finite native circuit; a 2^898 amplitude array is not executed.

Preparing the affine code is also a finite reversible circuit: source-word
prepare uniform6 on the 193 chord registers, accumulate the tree rows of C
from those values, then add P_x controlled by source position. Each modular
addition is a finite permutation on admitted labels with a retained unitary
complement. `field_preparation()` executes the three-bit uniform6 source word.
The full 898-cell initialization and the complete affine-preparation pulse
word are not executed in this receipt. That resource and complete apparatus
instrument remain an explicit integration obligation; an abstract isometry
is not misreported as an already replayed full laboratory preparation.

All M6 proper-code load/export complements and native reset timeout branches
from the common construction remain part of the complete instrument. Extra
register reset failures must similarly retain their tags and final carrier
state. This package executes a bounded prepared-code-sector current bridge;
it does not silently condition away failure, claim an aggregate success
probability for 898 resets, or identify service pulse counts with the
coherent walk clock. No reserved acceptance configurations are run.
