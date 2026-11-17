# Four-site vacuum and interactions from hopping

The next even chain uses the same signed hopping t=1/8, raw 432-register
mass and kappa=1/2,1,2. No new value is selected using a hadron comparison.
There are twelve CAR orbitals and 4096 kinematic states. Zero exterior
SU(3) flux allows exactly the complete singlet sectors N=0,3,6,9,12.
Their independently computed exact dimensions are 1,20,50,20,1.

For N=3k the Cartan-zero occupation space has binomial(4,k)^3 states.
The positive total colour Casimir has an exact rational nullspace R. Its
columns are retained with their full Gram R^T R; the calculation does not
mistake a nonorthogonal basis for a unitary embedding. Every mass, hopping
and prefix-electric component maps this complete kernel into itself.
The prefix electric operator is constructed from the actual CAR currents:

    6 C_prefix² = 3 sum_ab E_ab E_ba - N_prefix²,
    E_ab = sum_(x in prefix) c_xa† c_xb.

This includes every allowed intermediate colour representation of the
physical tree field. It is not a fundamental-representation truncation.

The half-filled N=6 sector has the unique lowest source energy for all
three members. Every other physical baryon-number sector is included and
its ground is higher. The complete matrices yield rational-inertia ground
intervals; a rounded rational mass is used only computationally and its
uniform 12|delta m| error is added outward. The N=3 and N=9 excitation
intervals coincide. Conditional electron-anchored gaps are:

| kappa | Two sites, GeV | Four sites, GeV |
|---|---:|---:|
| 1/2 | 0.780724218 | 0.541545301–0.541545366 |
| 1 | 0.353111464 | 0.271997971–0.271998035 |
| 2 | 0.163155173 | 0.140454860–0.140454924 |

The volume shift is material. None of these finite spectra establishes an
infinite-volume particle mass or the physical proton value. The same
local epsilon creator and remover have nonzero numerical ground-line
residues, approximately 0.055–0.522 across these members. Every spectral
weight and the unsuccessful probe effect are retained. These displayed
probe residues are diagnostics; no exact positive residue interval is
claimed by the ground-energy certificate alone.

## Interaction computed from the source

At zero staggered mass and zero hopping, let P contain every occupation
pattern with either zero or three quarks at each site. All its prefix
colour fields vanish. One literal quark hop across an occupied-empty edge
has electric cost D=4 kappa/3. Direct full CAR multiplication gives

    P K² P = 3 diag(number of occupied-empty edges),
    P K³ P = -6 (nearest-neighbour baryon transfer),

where K is the hopping matrix without t. Global colour constraints are
verified on all these images. Thus expansion of the exact Schur operator
about this degenerate manifold gives

    H_eff^(2) = -9t²/(4 kappa) sum_edges 1_occupied-empty,
    H_eff^(3) = -27t³/(8 kappa²) sum_edges(B_x† B_y+h.c.).

`induced_interaction.py` checks the identities on the two-site and all
three nontrivial four-site baryon sectors. These coefficients are derived
from the source hopping and colour cost; they are not independent inputs
or fits. The degree of freedom is baryon occupancy, not established
physical spin. The terms are perturbative coefficients, and no omitted
remainder is declared zero. The full finite spectra above use the exact
Hamiltonian rather than this truncated expansion.

## Conditional tension and units

The bare tree electric cost of a fundamental string is 4 kappa/3 per cell.
The free staggered signal slope is 2t=1/4. If this slope is identified with
physical c, the same energy unit E0=8.474166... GeV gives

    a = hbar c/(E0 2t),
    sigma_bare = E0² (2t)(4 kappa/3)/(hbar c).

Using the existing conversion hbar c=0.1973269804 GeV fm gives a≈0.09314 fm
and sigma_bare≈121.3 kappa GeV/fm. The common spatial Dirac velocity map has
not been proved, and this bare cost is not an extracted interacting
potential: string breaking and fermion screening remain relevant. It is
a conditional diagnostic showing that a different mass clock does not
automatically repair row22. It is neither a physical confinement result
nor a no-go for the theory.
