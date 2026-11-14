# Many-body constructions from the local walk

The fixed one-particle operator is W=SR on the complete padded
128-by-16 position/event space. R reflects the uniform twelve active event
slots. S holds the even slots and transports the odd slot of each of six
perfect matchings. The two constructions below are typed amendments of this
operator. Their architecture and native-law scopes differ.

## Binary ports with a shared event

Place a binary occupation n_x at each of the 96 physical ports and retain
one four-bit event register. For event 2d+1, exchange the occupations on all
48 seams of matching d. For even events and event padding, hold every port.
Call this involution S_coll and set W_coll=S_coll R. It preserves the number
of occupied ports. Its one-load restriction is exactly W on active source
positions, with the same complete event padding. Invalid position addresses
of the one-particle encoding have no interpretation as binary-port loads;
they remain separate code complements of the original representation.

For a seam (a,b), let sigma be the shared move bit. Set the original native
placement coin to c=sigma xor 1[n_a=1,n_b=0]. The literal binary balancing
law produces the exchanged pair when sigma=1 and the original pair when
sigma=0. The actual coherent native-event permutation followed by its
record recoder returns that pair, sigma and a clean record. This holds on
all four binary occupations and both native laws. A terminal inverse
recoder restores every original seam, coin and C/E label.

The original coins on different seams are correlated and depend on the
occupations through reversible control. The scheduler performs an entire
matching layer. It is therefore a declared coherent scheduler, without an
independent-attempt interpretation. The macro program uses phase-exact
native controlled swaps. Its complete local dirty-scratch action is unitary;
the stated occupation action requires clean scratch. Native preparation,
reset failures and proper-code export outcomes accompany a full execution.
No reset probability is silently discarded by the finite gate check.

If two occupied ports are nonadjacent, both move exactly when sigma=1.
With a uniform event their two movement indicators have law
P(0,0)=P(1,1)=1/2. Each marginal moves with probability 1/2, but their
covariance is 1/4, and the joint law has total variation distance 1/2 from
the product of its marginals. This persists at arbitrary separation in the
graph family. The shared event is a long-range prepared resource. This
calculation refutes cluster factorization for that preparation; it does not
by itself assert superluminal signaling.

## Unrestricted integer balancing needs its record

For fixed total two on one seam, each of (2,0), (1,1) and (0,2) is repaired
to (1,1) for either placement coin. The six input/coin columns therefore
give only two distinct output/coin columns if the current is erased. Keeping
the current distinguishes all six. No unitary can perform that erasure on
the declared pair-plus-coin state, since it would merge orthogonal inputs.
The binary construction is an invariant-domain construction, rather than
an extension of the unrestricted integer balancing law with erased memory.

## Local fermionic mode amendment

Associate one fermionic occupation mode with each padded one-particle
address (x,event). The physical sector uses twelve modes at every physical
port. Choose the exterior algebra and its occupation basis explicitly.
For a one-particle unitary U, Gamma(U) is its exterior-power action, with
Gamma(U) fixing the empty vacuum. Define

    W_Fock = Gamma(S) Gamma(R) = Gamma(W).

The one-particle restriction on all 2048 addresses is exactly the original
W. Every particle-number sector is unitary, and its k-particle amplitudes
are determinants of k-by-k minors of W. In particular the two-particle
amplitude from occupied modes c<d to a<b is

    W_ac W_bd - W_ad W_bc.

The code evaluates the same source callback on two orthonormal orbitals and
retains the complete 128-by-128 position-pair law, including coincident ports
and padding. It makes no projection onto the 191-dimensional cyclic sector.
The complete one-particle flat and quasienergy branches are inherited.

R is implemented locally within each sixteen-mode site. Eleven real Givens
rotations map the first mode to the uniform twelve-mode vector. Conjugating
occupation phases pi on the other fifteen modes gives Gamma(R), including
its vacuum and many-particle phases. Each Givens rotation has the native
Pauli-word form

    exp[-i theta X_a Z_(a+1)...Z_(b-1) Y_b/2]
    exp[+i theta Y_a Z_(a+1)...Z_(b-1) X_b/2].

The two generators commute. Native primitive words and parity accumulation
give a finite compiler. An occupation phase on b followed by the Givens
rotation at theta=pi/2 gives the exterior lift of a mode swap. The executable
checks its full vacuum, single, double and spectator action for several
separations using actual M6 pulse responses. The finite local identities
extend tensorially to arbitrary external references.

The exterior lift of a remote mode swap includes occupied-spectator parity.
It also assigns a minus sign to a doubly occupied pair. Those phases are
supplied by the chosen statistics. Ordinary swaps and the bare native
binary recoder do not supply them. The even CAR operators are onsite or
adjacent in the routed graph; their Jordan-Wigner implementation on encoded
native qubits can contain long strings. Distributed native timing and
physical locality require the separate clock/implementation contract.

This amendment changes the carrier state from one integer load per port to
multiple binary event modes per port. The aggregate physical port load can
range from zero to twelve. Its general aggregate dynamics differs from native
integer balancing. For example, occupy the uniform local mode and an
orthogonal superposition of two holding modes at one port. After one W_Fock
step the probability that both particles remain at that port is 1/2.
A literal addressed native event on (2,0) leaves one unit at each endpoint
with probability one. The aggregate discrepancy is an executable control,
not a transferred native-law claim.

The CAR choice does not establish spin, an observed particle species or the
spin-statistics relation. A symmetric Fock lift is another architecture;
it has different coincidence amplitudes. Neither statistics follows from
the one-particle W alone.

The same circuit has a negative-current one-hole sector. Fill all 1152
physical modes and use the basis c_j|filled>, including its occupation-order
phase. The active one-particle determinant is one: S has 288 transpositions,
and R contributes determinant (-1)^(11*96). Therefore Gamma(W) fixes the
filled state and acts on its holes by the complex conjugate matrix W. The
ideal W is real, so this is the same real one-hole matrix. Its current is
the negative of the corresponding particle current. The executable checks
all 1152 hole transport columns and their complement phases.

For a reference hole at port zero and a mobile hole at x, choose
G=12*ones-delta_0, n=12*ones-delta_x and particular field -P_x. These satisfy
the same Gauss law. This supplied filled state is not identified as an
energy ground state or a physical antiparticle vacuum. Empty-reference
particles and filled-reference holes are separate preparations. A simultaneous
particle-hole force experiment requires one common partially filled reference.

## One shared field and decoder domain

For collective binary transport, the current on each selected seam is
j_e=n_a-n_b. For mode transport it is the occupation difference of the
corresponding two odd event modes. Let n_x denote the sum of the physical
event-mode occupations at port x. Both transports obey

    n' + B(F+j) = n + BF.

The local mode reflection preserves each n_x, so it preserves the Gauss
constraint without changing the field. In a fixed total-Q sector choose
G=Q delta_0 and the inherited tree flows P_x. Put P_n=sum_x n_x P_x.
Then BP_n=G-n. The same 193 integral cycle columns C give the finite fibre

    J_Q |n,internal> = m^(-193/2)
        sum_(z in (Z/mZ)^193) |n,internal,P_n+Cz>.

Translation of z proves S_F J_Q=J_Q S and R J_Q=J_Q R. Exterior signs are
the same on both sides and do not affect the translation proof. The proof
works for every integer modulus m>=2, including composite moduli. It does
not assume a normalizable uniform state on an infinite integer fibre.
The full many-body affine preparation and reset instrument is a separate
finite-circuit integration obligation.

With known Q and reference G_R, a regional charge lies in
[-G_R,Q-G_R]. This interval has width Q. One regional residue modulo six
therefore determines the charge when Q<=5. The reference-independent
symmetric bound [-Q,Q] only establishes that conclusion for Q<=2.
At Q=6, charges zero and six alias.

Pointwise Gauss readout has a different domain. Binary occupations zero and
one are distinguishable modulo six for every total Q. In the mode amendment,
port occupations zero, six and twelve have the same residue. A modulus
m>=13 distinguishes all physical local occupations zero through twelve
when the reference is known. A reference-independent signed interval
[-12,12] instead requires m>=25. Enlarging the modulus changes the field
register and its finite preparation; the tree/cycle proof transfers, but a
native apparatus implementation must be regenerated for that register.
The modulus-six, thirteen and sixteen current-addition blocks are executed
as complete native reversible circuits. For modulus sixteen, descending
binary increment/decrement gates supply the addition directly and retain
all field labels. This is a circuit simplification for a declared register;
it does not identify theories at different moduli.

For a general declared reference with sum_x G_x=Q, the particular field is
sum_x(n_x-G_x)P_x. Its divergence is G-n. A filled reference state or a
different charge sector therefore requires its own attached field preparation.
The empty Fock vacuum and an occupied negative-quasienergy reference cannot
silently use the same Gauss background.

No autonomous field energy, force, physical clock or continuum claim follows
from these many-body extensions. Candidate selection and field-register
amendments require a common source version used by every downstream result.
