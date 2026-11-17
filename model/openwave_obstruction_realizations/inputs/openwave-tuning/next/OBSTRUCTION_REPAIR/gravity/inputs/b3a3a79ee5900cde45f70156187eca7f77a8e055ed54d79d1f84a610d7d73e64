# X1f stage 2: the same kernel on the four-component Dirac source (premise 3 where premises 1 and 2 hold)

Worker X1f, 2026-11-14. **Premise 3 is not produced on the Dirac source.** The entanglement area of a ball does not
respond at second order to a local rescaling of the terms of the law: there is no Poisson equation and no 1/r
solution. No OpenWave status is claimed; every official cell is untested.

Files: `stage2_lib.py` (source, response engine, planar reduction), `stage2_run.py` (parts in its header; each
writes `stage2_part_<name>.json`, merged into `stage2_results.json`), `stage2_verify.py` (17 of 17 checks, each
paired with an input on which the same test fails; `stage2_verify_log.txt`), logs `stage2_log_*.txt`. Stage-1
files are unchanged. Units: cell side 1, one step of the walk 1, ħ = 1; energies are levels of G = sin H;
c = √(μ cot μ) cells per step; entropies in nats; η = 1.49 per cell face (X1c).

## 0. Result in brief

1. Matter term: positive, first law holds, the same within ± 6 percent for every low-band pair; no number term.
2. Lapse (rest and hopping scaled together): the response of the ball entropy has no isotropic second-order
   term at the scale of the matter term (remainder −2 × 10⁻⁴ to −1.1 × 10⁻³ of Q/δE, inside an indefinite form
   that is negative along x and z and positive along y); it begins at fourth order and follows the
   vacuum-energy susceptibility, χ_NN(q)/q² → 2 × 10⁻⁵ (scalar source: 0.349).
3. Relative factor (hopping against rest, the spatial factor at fixed lapse): responds only through the mass,
   with a term at q = 0 (0.0054, 0.0195, 0.0557 for μ = 0.05, 0.1, 0.2; 27 cells) and second-order terms of the same
   small size.
4. Area identity of the matched limit: every ratio is between −1.2 × 10⁻³ and +1.5 × 10⁻³ where +1 is needed.
5. The equilibrium condition therefore fixes the relative factor locally (a stretch where the excitation is)
   and leaves the lapse free up to fourth order. The two-factor ratio is not found.

## 1. Existing results (additions to the table of `SOLUTION.md`)

| Needed | Existing result (path, statement, domain) | Action | Remaining difference |
| --- | --- | --- | --- |
| the Dirac step | the other team's `openwave-tuning/common_dynamics/spin/dirac_proposal.py::bloch(p, h, m)`, sha256 78b5b44b…da42a ("This file declares new dynamics. It does not identify them with W=S R") | reuse with an adapter: the five-line formula is copied and compared (5 × 10⁻¹⁶); the file is not imported | none |
| the state, its entropy per face, its modular generator | X1c (`blockers/X1c_tower_premises_on_walk/SOLUTION.md`, sections 5 and 6): principal quasi-energies in (−π, 0) filled; η = 1.49 ± 0.01 per face; half space: K = s (x + ½) sin H with residuals 6 × 10⁻⁸ to 1.3 × 10⁻⁴, s c/2π = 1.0002 at μ = 0.05; balls: slope 0.96 to 0.99, residual 9 to 15 percent; gapless at (±π/2, ±π/2, ±π/2) for every mass | reuse directly: same state (negative branch of sin H); slab entropy 1.49649 per face reproduced | handles and responses are new |
| KMS and collars | X1e (`blockers/X1e_cap_theorem_hypotheses/VERDICT.md`): 2π-KMS defect 8.7 × 10⁻⁶ on the Dirac half space; collars fail (eight gapless points) | cited | none |
| the two-factor target | audit item 7 (`openwave-tuning/next/CHAIN4/gravity_packet_reviews/FOR_FABLE.md`): "Rest scales by N and spatial kinetic terms by N²"; X1g: one-factor register N^(1/2) H N^(1/2), bending one half | reuse as the frame of the dictionary (section 2.3) | the handles of this source are defined here |
| the composition | `model/matched_einstein_limit/PROOF.md`, equation (1) | reuse directly | residuals for one observer, first variations (section 4.6) |
| the free Weyl source of the repository | `model/weyl_geometric_modular`, `model/lattice_modular_limit` (U = 0 Weyl source, eight valleys): geometric modular group for double cones; "not a derivation of … a dynamical metric" | cited | a different lattice source; no response to geometry is computed there |

Missing (searched "lapse", "coframe", "scale factor" together with "entropy" in the `.md` files of `model/` and
of the push): a response of a regional entropy to a lapse or coframe handle on a Dirac or Weyl lattice source.

## 2. Source, state, handles, dictionary (item 1)

2.1 Source and state. One step U(k) = S(k) e^(−iμβ), S = e^(−ik_z α_z) e^(−ik_y α_y) e^(−ik_x α_x), α_j = σ_j ⊗ X,
β = 1 ⊗ Z. The filled set is the negative branch of the strictly local Hermitian generator

    G = sin H = i(U − U†)/2 = cos μ · G_hop + sin μ · G_rest,   G_hop = i(S − S†)/2,   G_rest = (Sβ + βS†)/2.

Facts used: G(k + π e_a) = −G(k) for every axis (the sign (−1)^(x_a) anticommutes with G: filling one half,
no number term); the levels at the eight valleys k_a ∈ {0, π} are ± sin μ; G vanishes at the eight reflection
points (±π/2, ±π/2, ±π/2) for every μ, with velocity matrix σ_x ⊗ 1 along x there: two Weyl cones of one
chirality at each point, which no translation-invariant on-site term gaps. Every step moves by (±1, ±1, ±1):
the lattice is four decoupled body-centred sublattices, and cells that do not differ by an odd step in all
three directions are uncorrelated. Balls of 1 and 7 cells are maximally mixed (4 ln 2 per cell) and respond to
nothing; the balls used have 27, 93 and 257 cells (r² ≤ 3, 8, 16; 108 to 1028 modes).

Treatment of the gapless points. With them: blocks of side 16 to 48 with a half twist in every direction (side
divisible by 4) contain none of the 16 points; their cones are present (smallest level 0.30 on 16³, 0.21 on
24³ at μ = 0.2). Without them, three ways: (i) the decomposition of the response by the momentum of the virtual
transition, leaving out the neighbourhoods of the reflection points or of the valleys (section 4.7); (ii) a
variant source with a staggered term i m_s (−1)^(x+y+z) sin k_x sin k_y sin k_z, which gaps the reflection points
by m_s and vanishes cubically at the valleys, on a dense 8³ block; (iii) the mass scan (the valleys are gapped
by sin μ, the reflection points never).

2.2 Handles (declared). Two cell fields ν_h, ν_r:

    G(ν_h, ν_r) = ½{1 + ν_h, cos μ · G_hop} + ½{1 + ν_r, sin μ · G_rest}

(each hop of the generator between two cells is scaled by the mean of their factors). Uniform ν_h = ν_r is an
exact null direction (G → (1 + ν)G; measured 10⁻¹⁶ to 10⁻¹⁸). The staggered profiles (−1)^(x_a), (−1)^(x_a + x_b) and
(−1)^(x+y+z) are exact null directions of the one-factor handle as well (seven zeros on 16³). Any first-order
change δG acts on the state through its filled–empty block only, and that block is produced by a kick of the
step, U → U e^(−iδA) with δA_ij = 2 δG_ij/(e^(−iE_i) + e^(iE_j)), whose denominator does not vanish between a filled and
an empty mode. The mass-coin angle itself, μ → μ(1 + ν), is the combination ν_r = (μ/tan μ) ν, ν_h = −μ tan μ · ν.

2.3 Oriented dictionary. At each valley G ≈ ± [(1 + ν_h) α·q + (1 + ν_r) μβ]; a static metric
ds² = −N²dt² + e^(2σ)dx² gives E² = N²(μ₀² + c²e^(−2σ)q²). Hence rest term × N, hopping term × N e^(−σ):

    δN = ν_r,    δσ = ν_r − ν_h.

| Family | ν_h, ν_r | δN, δσ | G_uu (first order) |
| --- | --- | --- | --- |
| lapse at fixed coframe (one factor; X1g's register) | n, n | n, 0 | 0 |
| coframe at fixed lapse | −σ, 0 | 0, σ | −2∇²σ |
| two-factor target | 2Φ, Φ | Φ, −Φ | 2∇²Φ |
| rest only | 0, n | n, n | −2∇²n |

Orientation: ν_h > 0 is faster hopping, ν_r > 0 a faster rest phase. A positive mass needs Φ < 0: slower clock
and slower hopping.

2.4 Fixed proper volume. The proper volume of a set of cells is Σ e^(3σ) = Σ (1 + 3(ν_r − ν_h)). The lapse family
leaves it unchanged; the other families change it. The balls here are fixed sets of cells; a comparison at
fixed proper volume adds −2S₀⟨δσ⟩ over the ball by the area law (the term of stage 1; not executed here).

## 3. Controls

| Control | Result | Failing value |
| --- | --- | --- |
| step against the other team's formula | 4 × 10⁻¹⁶ | 1.94 with the factors in the order x·y·z |
| engine against dense kernel and finite differences (6³ dense, 27 cells) | hop 0.16746971 / 0.16746972, rest 0.00494197 / 0.00494197; engine against kernel 4 × 10⁻¹⁶ | the rest kernel against the hop difference: 0.163 apart |
| first law on this source | point: 9.664172 against 9.664168; number-neutral Gaussian pair: 1.423223 against 1.423223 (8³: 9.616261 / 9.616258; 1.268345 / 1.268344) | amplitude 10⁻⁵: 9.663839 |
| the scalar source of stage 1 through the same engine | ζ(q₁) = 1.722635, neutral Q/δE = 18.33566, number term −4.3387: the stage-1 values in all digits | |
| a non-neutral excitation | Dirac: number term below 1.4 × 10⁻¹⁰ for every mode (one particle and its staggered partner hole change the entropy alike) | scalar source: −4.3387 (seven cells) |
| no excitation | zero source, zero field | with the Gaussian pair the field components reach 2.3 |
| uniform one-factor lapse | Z_N(0) = 10⁻¹⁶ to 10⁻¹⁸ | hop factor alone: +0.0027 to +0.138 |

## 4. Results (μ = 0.05, Compton length 20 cells, unless stated)

4.1 Matter term (48³). Q = Σ_x δS_m(R_x) for a number-neutral excitation of energy δE.

| Ball (cells) | S (per face) | Q/δE, low band near the valleys | near the reflection points | all low-band pairs, min to max | Gaussian pairs 1, 1.5, 2 | point | continuum boost value (2π/c)(4πℓ⁴/15) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 27 | 56.209 (1.041) | 31.77 | 30.89 | 29.88 to 32.47 | 30.78, 30.80, 30.83 | 30.76 | 63.2 |
| 93 | 135.151 (0.901) | 263.3 | 248.3 | 232.8 to 263.3 | 241.2, 243.9, 245.6 | 238.2 | 328.6 |
| 257 | 273.735 (0.931) | 1192.6 | 1113.4 | 1052.4 to 1192.6 | 1086.3, 1096.3, 1103.7 | 1075.0 | 1274.3 |

The ratio to the boost value is 0.50, 0.80, 0.94: the matter side approaches the boost form with the size of
the ball (scalar source: 0.18 to 0.85). No eigenvalue of the three balls is unresolved (smallest 3 × 10⁻⁵).

4.2 Lapse response. δS(R_x) = Σ_q Z_N(q) n̂_q e^(iq·x), Z_N = Z_h + Z_r.

| Ball (cells) | −Z_N/q² at the smallest momentum of 48³ (q = 0.131), axes x, y, z | isotropic mean at q = 0.131, 0.196, 0.262, 0.393 | extrapolated to q = 0 (two ways) | as a fraction of Q/δE |
| --- | --- | --- | --- | --- |
| 27 | −0.0375, +0.0408, −0.0375 | −0.0114, +0.0021, +0.0289, +0.0763 | −0.022 or −0.007 | −7.0 × 10⁻⁴ or −2.2 × 10⁻⁴ |
| 93 | −0.338, +0.096, −0.338 | −0.193, −0.084, +0.125, +0.480 | −0.280 or −0.145 | −1.1 × 10⁻³ or −5.5 × 10⁻⁴ |
| 257 | −0.912, +0.237, −0.912 | −0.529, −0.046, +0.853, +2.262 | −0.915 or −0.217 | −7.7 × 10⁻⁴ or −1.8 × 10⁻⁴ |

(The two extrapolations use the first momenta of 48³ and 32³, and their second momenta; the first momentum of
a block carries an offset of about −7/L² for 27 cells.) The second-order form is indefinite: negative along x
and z, positive along y (the step is ordered z·y·x); the difference between the axes, 0.078, 0.43, 1.15, is
1.3 to 11 times the isotropic part. From q = 0.39 on the response has one sign on all axes and grows as
q^3.0 to q^4.3 along an axis (between q = 0.39 and 0.59 on 32³). The static susceptibility of the vacuum energy to the same handle is
χ_NN(q)/q² = 0.00094, 0.00172, 0.00270 at |n|² = 1, 2, 3 of 32³, extrapolating to 2 × 10⁻⁵ (0.0012 at μ = 0.5), and
−Z_N/((Q/δE) χ_NN) = 0.49 to 0.94, 0.40 to 0.72, 0.50 to 0.63 for the three balls at |n|² ≥ 2 of 32³:
the entropy response is the first law applied to the vacuum-polarisation energy, within a factor two.
Full zone (27 cells, 16³): the one-factor response is negative at 4026 momenta, positive at 62 (from
|n|² = 57 on) and exactly zero at the seven staggered momenta.

4.3 Relative factor. Z_h(0) = −Z_r(0) = c; Z_r(q) = −c + ζ′q² at small q.

| μ | c, 27 cells | ζ′, 27 cells | c, 93 cells | ζ′, 93 cells | χ_hh(0) |
| --- | --- | --- | --- | --- | --- |
| 0.05 | 0.0054 | 0.0115 | 0.0413 | 0.130 | 0.0029 |
| 0.1 | 0.0195 | 0.0373 | 0.1456 | 0.415 | |
| 0.2 | 0.0557 | 0.0849 | 0.3979 | 0.932 | |
| 0.5 | | | | | 0.207 |

(32³.) Raising the hopping against the rest term lengthens the Compton length in cells and raises the entropy
of a ball (c > 0). Both coefficients fall with the mass. On 48³: c = 0.0055, 0.0420, 0.1377 for 27, 93, 257 cells.

4.4 Which factor the condition fixes (item 3). In the variables (n, σ) the response is Z_N n − Z_h σ. The lapse
n enters with no term at q = 0, with the indefinite second-order remainder of 4.2 and with the fourth-order
term. The relative factor σ enters with −c at q = 0. One ball family therefore fixes σ locally: σ = (Q/δE)/Z_h
per unit matter amplitude, +8213 (27 cells) and +7978 (93 cells) at |n|² = 3 of 48³: a stretch (longer proper
step, slower hopping against rest) where the excitation is, of size ∝ 1/c. The lapse at fixed σ is the
undetermined combination. Two ball families at one momentum determine both formally: σ/n = +6.6 and +2.1
(27 and 93 cells, |n|² = 3 and 6), +0.22 and +0.39 (93 and 257 cells), with condition numbers 10 to 44: no stable
value, and not the two-factor ratio −1. In the matched limit the same structure appears: the rest-frame
equation G_uu = 8πG T_uu contains the spatial factor only.

4.5 Planar reduction (item 4). Slab of 16 layers in a ring of 64, 16² transverse momenta, half twists: entropy
1.49649 per face (X1c: 1.49 ± 0.01). Lapse kernel per layer and face, counted from a cut: +0.04254, −0.03821,
−0.00023, −0.00310, −0.00023, −0.00054, −0.00010, −0.00014 on either side; total 10⁻¹⁷; mirror symmetric within
3 × 10⁻⁵; deep inside −0.00025 per layer. Its second moment about the cut is −0.30: per face
δS = −0.15 ∂_n²n for a profile along the normal. On a closed surface this normal response is cancelled by the
response of the other faces to the tangential derivative to about one percent (18 faces × 0.15 against
0.0375 along x for 27 cells; the isotropic remainder of 4.2). For half spaces at every position the condition is of fourth order, τ n′′′′ = s ρ with τ = 0.15 and
s = 2π/c (X1c). A large-ball or half-space limit of a Poisson strength does not exist on this source. Hop
factor alone: total +0.00093 per face (the mass term).

4.6 Area identity (item 5) and the matched limit. Needed: entanglement contrast = η × (−C_ℓ G_uu), that is a
second-order coefficient 2ηC_ℓ per unit ∇² of the spatial factor (29.9, 155.8, 604.1 for the three balls) and
none for the lapse. Found at the smallest momentum of 48³, as fractions of that scale:

| Ball (cells) | lapse (0 needed) | coframe at fixed lapse (+1 needed) | rest only (+1 needed) | two-factor (+1 needed) |
| --- | --- | --- | --- | --- |
| 27 | −3.8 × 10⁻⁴ | −7.9 × 10⁻⁵ | −4.6 × 10⁻⁴ | +3.0 × 10⁻⁴ |
| 93 | −1.2 × 10⁻³ | +2.5 × 10⁻⁴ | −9.9 × 10⁻⁴ | +1.5 × 10⁻³ |
| 257 | −8.8 × 10⁻⁴ | −3.5 × 10⁻⁴ | −1.2 × 10⁻³ | +5.3 × 10⁻⁴ |

In the lapse channel the entanglement area behaves as a covariant area would (no response, to 10⁻³). In the
spatial channel it does not respond either, where the identity needs the full coefficient. Residuals of
`model/matched_einstein_limit/PROOF.md` on this source (one observer, first variations): e_bulk/(2πC_ℓT_uu) = 0.50,
0.20, 0.06 for 27, 93, 257 cells, falling; e_area/(C_ℓ|G_uu|) = 1 within 1.5 × 10⁻³ for every family with a
spatial factor; e_stationarity cannot be brought to zero by a field entering at second order, only by the
local term c of the relative factor, which differs between ball families by up to a factor 1.9.

4.7 With and without the gapless points. Decomposition (27 cells, 48³): isotropic −Z_N/q² at the smallest
momentum −0.0114 with all transitions, −0.0159 without those within 0.4 or 0.8 of a reflection point, −0.0147
and −0.0145 without those near a valley; at the second momentum +0.0289, +0.0116 to +0.0069, +0.0138 to +0.0093.
The fourth-order part comes in equal shares from the valleys and from the reflection points; the second-order
remainder comes from the rest of the zone. Variant with the reflection points gapped (m_s = 0.3, dense 8³): the
uniform lapse stays a null direction (3 × 10⁻¹⁶, even and odd centres alike); ζ(q₁) = 0.135 → 0.107 (27 cells),
0.183 → 0.146 (33 cells); c = 0.00271 → 0.00283. An 8³ block does not resolve the cones (smallest level 0.53
with and without); the variant was not built in Bloch form (eight cells per unit cell).

4.8 Strength against 1/(4η). A_chain = 1/(4ηc) = 0.168. No coupling of that kind exists here: the stiffness
that a Poisson equation needs is absent to 10⁻³. If the isotropic remainder of 4.2 alone is read as a
stiffness, the lapse would be n = −|A| δE/r with |A| = 75 to 440, the attractive sign and 450 to 2600 times
the chain's value; the anisotropic part of the same form is larger and has the other sign along y, the two
extrapolations differ by a factor 2 to 4, and the scaling with the ball is not established. This is not a
result.

## 5. Why the Dirac vacuum does not respond

A linear lapse profile n = g·x is the lattice boost Σ (g·x) G_x. X1c measured that the entanglement Hamiltonian
of a half space of this state is the lattice half-boost; for a pure state the entanglement Hamiltonians of the
two sides act equally on it, so the full boost annihilates the vacuum and a linear lapse changes nothing at
first order. In a conformal theory the isotropic quadratic profile Σ |x|² G_x is the generator of special
conformal transformations and annihilates the vacuum as well; the massless valleys and cones of this source
obey that to 10⁻³. What remains at second order is the traceless anisotropic part and terms ∝ μ. On the scalar
source of stage 1 the vacuum is invariant under neither (c₂ = 0.349), and its repulsive Poisson term was that
failure. Premise 2 and a second-order response of the entanglement area to a rescaling of the law on a fixed
complex exclude each other. In the matched limit the same statement reads: the lapse does not enter G_uu, and
the area responds only if the proper length of the cutoff does; on a fixed complex the modes at the cutoff are
massless and see the lapse and the spatial factor only through their ratio.

## 6. The route, link by link

| Link | Statement | Status | Evidence |
| --- | --- | --- | --- |
| D1 | The state is X1c's and has the area law | reused; reproduced | slab 1.49649 per face (W13) |
| D2 | The response engine is exact on this source | executed with failing controls | W5, W6, W7 |
| D3 | The matter term is positive and of the boost form | executed | 4.1; first law (W6, W8) |
| D4 | The ball entropy responds to the lapse at second order | **broken** | 4.2; χ_NN/q² → 2 × 10⁻⁵ (W9, W10, W11) |
| D5 | It responds to the spatial factor at second order with the coefficient 2ηC_ℓ | **broken** | 4.6: at most 1.5 × 10⁻³ of it, mixed signs |
| D6 | The equilibrium solution is a 1/r field | **broken** | fourth order in the lapse; local in the relative factor (4.4, 4.5) |
| D7 | The solution has the two-factor ratio | **broken** | σ/n = +0.2 to +6.6, no stable value |
| D8 | The result depends on the reflection points | not found | 4.7 (decomposition; 8³ variant) |

## 7. What was put in

The other team's Dirac step (a declared replacement dynamics) and its mass μ = 0.05, 0.1, 0.2, 0.5; Fermi
statistics and the principal-phase filling (X1c); the generator sin H as the operator whose terms are scaled;
the two handles with the mean-of-two-cells rule; the reading of the handles as lapse and coframe; twisted
blocks of 6³ to 48³ cells; balls of 27, 93, 257 cells around a cell; number-neutral excitations (a particle and
the hole of its staggered partner; Gaussian pairs; the point mixture); the eigenvalue floor 10⁻¹³; for the
variant the staggered term and m_s = 0.3; for comparison X1c's η = 1.49 and the continuum boost coefficient.
Nothing is fitted; no measured constant enters.

## 8. Not computed

- Anisotropic hopping factors (one factor per axis: shear of the coframe). The planar kernel shows that a face
  responds to the normal derivative with −0.15 per face; a direction-dependent handle is not constrained by the
  two symmetries above and may respond at second order. It is not the Newtonian (conformal) sector.
- Balls beyond 257 cells; blocks beyond 48³; the limit q → 0 of the isotropic remainder (two extrapolations
  that differ by a factor 2 to 4).
- The variant without reflection points beyond 8³; a symmetrised step (which would remove the anisotropy).
- Boosted observers; time-dependent and nonlinear responses; the Wilson primary.
- Balls resized to fixed proper volume.
- A wave packet on this source (stage 1 ran it on the scalar source; X1g ran the one-factor register).

## 9. Reproduction

    cd openwave-push/blockers/X1f_entanglement_area_response
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
    PY=../../../reverse-engineering-reality/.venv/bin/python
    for p in symbol scalar "dense 0.2 0 6 0.25" "dense 0.05 0 8 0.5" "dense 0.05 0.3 8 0.5" \
             "ball 0.05 16 3" "ball 0.05 24 3" "ball 0.05 32 3" "ball 0.05 40 3" "ball 0.05 48 3" \
             "ball 0.05 32 8" "ball 0.05 48 8" "ball 0.05 32 16" "ball 0.05 48 16" \
             "ball 0.1 32 3" "ball 0.1 32 8" "ball 0.2 32 3" "ball 0.2 32 8" \
             "chi 0.05 32" "chi 0.5 32" "refl 0.05 48 3" "slab 0.05 16 64 16" "full 0.05 16 3" summary merge; do
        $PY stage2_run.py $p; done; $PY stage2_verify.py

CPU per part: 45 s (slab), 39 s and 38 s (dense 8³), 37 s (full zone), 14 s or less (every ball), 10 s (χ);
one script at a time; at most 2048 modes.
