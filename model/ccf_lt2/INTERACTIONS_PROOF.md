# Exact interaction structure and scope

`interactions.py` returns complete local native CAR polynomials. The source
assembler adds them to the complete retained matter source. Structural
identities and every control below are executed by the full verification
command and recorded in [RESULTS.json](RESULTS.json).

## Algebra and preserved quantities

For a normal CAR word `W`, a diagonal number charge
`D=sum_i d_i c_i^dagger c_i` satisfies

    [D,W] = (sum_creators d_i - sum_annihilators d_i) W.

Consequently checking this coefficient on every nonzero combined word
proves the complete operator commutator. Electric weights are the actual
native `3Q`; baryon weights are +1 on Q coordinates and -1 on their
conjugates, measuring `3B`; Dirac lepton weights are +1 on the left
doublet and -1 on its conjugates. Every supplied delta preserves these
three quantities. Every word is supported at one vertex, so electric
neutrality is local. Electric rotor operators commute with these deltas,
and all local CCF Gauss operators are preserved without adding a Wilson
route. This would not justify omitting routes for a nonlocal variant.

Each current is Hermitian because both `bar_sigma_mu` and `T_a` are
Hermitian and act on separate factors. The rational family frame is
exactly orthogonal, so the weak doublet has the same current algebra after
rotation. Real signed sums of current squares are Hermitian, even though
the Minkowski contraction is not a positive sum. Both orders and all
constant/lower-degree terms arising from CAR normal ordering remain.

The leptonic current contains only native lepton coordinates; the quark
current and family shifts contain only colourful coordinates. Thus no
quark-lepton cross term is generated. The parent semileptonic contact is
neither replaced nor doubled by these additions.

## Complete quark Nambu frame

Pair a Q annihilator b with the opposite-spin native conjugate a. With
`r=epsilon_spin*eta_colour` equal to +1 or -1, set

    psi = (b+r a^dagger)/sqrt(2),
    chi = (a-r b^dagger)/sqrt(2).

These obey the full CAR and invert to

    b = (psi-r chi^dagger)/sqrt(2),
    a = (chi+r psi^dagger)/sqrt(2).

The map therefore changes coordinates on the complete colourful Fock
algebra. It does not discard the right/complementary fields or occupation
sectors. The left current acts on psi but the Hamiltonian and reference
retain the complete carrier. The colour signs `(1,-1,1)` express the
epsilon intertwiner in the actual ordered exterior basis. The replay
checks that intertwiner against the exact colour generators read from
the bound native artifact.

Under the parent hadron dictionary, b and a are physical-candidate quark
and antiquark annihilators. A square of left currents therefore contains
terms of the form

    (b_u^dagger bar_sigma_mu b_s)
    (b_d^dagger bar_sigma_mu epsilon a_tilde_u^dagger)
    + adjoint,

with their full spin/colour sums and family-rotation coefficients. Such a
word preserves Q and B but changes total colourful occupation by two and
changes the strange-family occupation. This directly removes the two
selection rules of the imported colourful source. Its coefficient and
nonzero action must be read from the complete polynomial, rather than
installed as a separate selected channel.

The same exact reasoning uses the lepton source's charged/neutral flavour
Nambu frame: native left annihilators have particle-annihilation and
antiparticle-creation parts. The unfiltered pure-lepton contact can thus
map one muon-flavour excitation to an electron-flavour excitation and a
neutral particle-antiparticle pair. The exact flavour CAR frame does not
replace the original source's rounded neutral tensor or assert that these
outgoing neutral flavours are stationary mass eigenstates.

## Extensive lower bound

For each fixed weak/spin current the one-body matrix has eigenvalues
`+1/2` and `-1/2` in equal numbers. There are 12 left leptonic modes at a
vertex and 36 left quark modes. Second quantization therefore gives
`||J_L||<=3` and `||J_Q||<=9`, including their full Nambu complements.
The family rotation does not change those norms. Bounding all twelve
current squares by their norms yields, for V matter vertices,

    H_LL >= -216 |g_L| V I,
    H_QQ >= -1944 |g_Q| V I.

The declared family shifts are nonnegative. Thus these fixed-capacity
local additions preserve an extensive lower bound whenever the assembled
parent/spatial Hamiltonian has one. A stronger bound or a proof that the
native empty vacuum reduces the source is not asserted. The ground state and low-energy
particle spectrum can change.

## Why the checks are not a decay result

The finite occupation probes test matrix elements of the entire operator.
They are not exact Lambda, proton, pion or dressed muon eigenstates. A
nonzero matrix element only removes a selection-rule prohibition. It does
not prove an open energy-conserving shell, outgoing separation, an
exponential survival regime or a physical lifetime.

The full source requires its own ground, threshold, dispersion, complete
outcome and physical-unit calculations before a decay claim.

The added family shifts and current strengths are rational construction
choices. They do not derive measured quark masses, CKM data, the Fermi
constant or any OpenWave numerical target. All original calibrations
remain disclosed. Coherent quantum memory is required for the specified
amplitudes; the diagonal conserved-charge bookkeeping alone does not
establish those amplitudes.

## Executed structural result

The `interactions` result in [RESULTS.json](RESULTS.json) contains
204 leptonic-contact monomials, 96 family-splitting monomials, and 47,041
nonleptonic monomials. Every word passes local charge, baryon and lepton
conservation; Hermiticity and the absence of additional quark-lepton cross
terms pass exactly. The latter count comprises one combined scalar and
47,040 quartic words. There are 29,400 colourful-number-changing words and
37,656 colourful-family-changing words.

The exact checks cover all 5,184 ordered pairs of the full 72-mode quark
Nambu map and all 576 pairs of the 24-mode lepton flavour map. Full native
and canonical expressions of the leptonic contact agree. Its action on
the specified one-muon occupation retains 91 nonzero occupation outputs;
two requested electron/neutral-pair outputs each have amplitude -1/4096.
The nonleptonic contact's action on the declared uds occupation retains
1,451 nonzero outputs; five light five-occupation outputs have amplitude
of magnitude 3/25600. These are Hamiltonian matrix elements in source
units, not time-evolved probabilities or exclusive hadronic amplitudes.

The identity-family-rotation control removes the specified strange
conversion; zero contact and number-conserving truncations remove their
respective vertices. A missing adjoint and an altered charge assignment
are rejected. Setting every addition to zero returns exactly zero delta.

Those are component controls. The full enlarged source may produce
additional virtual channels: the nonleptonic Nambu contact need not leave
the native colourful vacuum reducing. In particular, deleting the explicit
leptonic contact is not a proof that the complete model has no
quark-mediated leptonic conversion. This is another reason to recompute
the full vacuum and physical channels.

