# Protected objects and source excitation families

The eleven protected-object entries retain their distinct construction domains.
They include conserved charges, protected memory, two kinds of winding and an
exclusion control. They are not eleven established topological particle species.
The sixteen excitation identifiers form a separate list. Their recorded
`protected_defect_id` fields are null: the defect-to-excitation attachment is unproved.

## Protected-object families

### Exchange-cover memory

**Exchange-orientation bit** (`ordered_exchange_charge`). A coherent orientation bit carries two protected charge sectors under the specified ordered-cover repair.

Domain: H66 with a coherent orientation bit; complete ordered-cover repair.

Protected object: Q=I66 tensor X, two sharp eigenvalue sectors.

Topology: Z2 character of the ordered exchange cover.

[Source proof](../../model/defect_protection/README.md).

### Spin, symmetry and memory sectors

**Spin-rotation charge** (`selected_spin_charge`). A local spin operation creates two opposite conserved internal-charge channels. This charge is not identified with electric charge.

Domain: actual eight-site selected spin cube, J=1, B=1/4; inherited two-flavor repair.

Protected object: integer U1 spin charge Q; global X parity is exp(i pi Q); local Z creates two sharp Q=+/-1 channels.

Topology: global internal U1 symmetry, not a winding number or identified electric charge.

[Source proof](../../model/defect_localization_readout/PROOF.md).

**One-particle-per-site spin memory** (`single_occupancy_spin_memory`). Repair preserves the spin information in a fixed-occupation code. Identity repair supplies no winding invariant.

Domain: eight-site one-particle-per-site spin code with complete C/E repair.

Protected object: full M256 operator algebra under repair alone.

Topology: no topological label supplied by identity repair.

[Source proof](../../model/defect_localization_readout/PROOF.md).

**Canonical spin fields** (`joint_positive_field_canonical_sector`). A specified spin family admits two canonical field pairs with one positive-energy clock. Its physical particle dictionary is unassigned.

Domain: protected spin interaction H_B=J sum(S_x.S_y-1/4)-B sum epsilon_x S_x^1; supplied p=n^-2 and B=J n^-5; positive-B ground and retained H_B clock.

Protected object: two regular canonical pairs (Q3,P2) and (Q2,P3), with common positive-energy time and full internal U1 on their bounded histories.

Topology: internal source spin rotation, no spatial winding or identified electromagnetic charge.

[Source proof](../../model/spin_joint_field_limit/PROOF.md).

**Shared spin subsystem** (`common_hidden_spin_factor`). A qutrit spin factor survives both the stated field dynamics and repair, with 105 multiplicity states.

Domain: same six-site Weyl/contact field and complete charged-block repair.

Protected object: at N=4, S=1, a qutrit with 105 multiplicity states, total dimension315.

Topology: internal spin factor, no topological boundary invariant.

[Source proof](../../model/common_spin_protection/README.md).

**Eta-pairing sector** (`cube_eta_range`). A rank-1820 symmetry sector is protected by the specified dynamics and repair. Its label is internal eta symmetry.

Domain: original eight-site Weyl/contact field with the stated four-dimer repair.

Protected object: rank1820 image of eta+ from N4 to N6; eta>=2 sector.

Topology: eta symmetry representation, no spatial winding invariant.

[Source proof](../../model/cube_eta_protection/README.md).

### Occupation and current responses

**Mobile vacancy with a spin record** (`one_hole_spin_readout`). One missing occupation can be located, while recorded transport carries the remaining spin information.

Domain: full two-flavour CAR repair on H_h: one empty site and one fermion at every other site; dimension V 2^(V-1).

Protected object: hole-position projectors and internal spin information carried by the recorded signed transport; conditional block dimension 2^(V-1).

Topology: one-unit occupation deficit, without a topological invariant or physical particle assignment.

[Source proof](../../model/source_family_compatibility/REPAIR_TRAJECTORY.md).

**Conserved-number response family** (`centered_full_CAR_charge_response`). Particle addition, removal and spin probes share an interacting fermion algebra with conserved source occupation.

Domain: fixed-coupling original Weyl/contact action with bare full-domain repair; invariant centered full-ground subsequence; the signed Hubbard chart is a separate direct local instance.

Protected object: odd addition/removal responses and even spin observables on the same full CAR algebra; total occupation is conserved under source propagation and repair.

Topology: source U1 charge grading, without an assigned winding invariant or Standard Model charge.

[Source proof](../../model/source_family_compatibility/CHARGED_SPECTRAL_BANDS.md).

### Spatial winding textures

**Opposite winding cores** (`compact_opposite_texture_cores`). A prescribed two-core texture has shell degrees +1 and -1 and four localized band modes. Spontaneous production from repair is not established.

Domain: prescribed two-core texture in the declared cubic two-spin field.

Protected object: opposite shell degrees +1 and -1, constant exterior, four localized band modes.

Topology: genuine shell winding of the prescribed order-parameter field.

[Source proof](../../model/compact_defect_pair/PROOF.md).

### Momentum-space winding

**Oppositely wound Weyl nodes** (`kinetic_Weyl_nodes_and_free_fields`). Eight kinetic nodes have four positive and four negative windings. A stated scaling family has a conditional free-field map; this does not identify eight observed massless species.

Domain: original cubic Weyl symbol -c sum sigma_j sin(k_j); centered full-ground family c=v/a and 0<=U=o(c^(2/7)).

Protected object: eight isolated kinetic nodes, four of winding +1 and four of winding -1; conditional free vacuum map to eight Weyl fields.

Topology: momentum-space kinetic winding with explicit boundary-gap stability; period-two onsite mass channels couple opposite winding sectors.

[Source proof](../../model/source_family_compatibility/KINETIC_NODE_MAP.md).

### Exclusion control

**Bare source has no nontrivial reducing charge** (`bare_native_source`). The full branch commutant is scalar in the specified bare source-only domain. This entry is an exclusion, not a protected defect.

Domain: complete all-seam C/E on H66 or H78; source-only charges.

Protected object: only scalar branch-commuting charges.

Topology: no protected source-only reducing charge in this domain.

[Source proof](../../model/completion/defect_algebra.py).

## Source excitation families

Names below describe source components or conditional analogies. Repeated
electron, muon and tau comparison labels do not create three generations.

### Internal matter weights

**Up-type doublet component** (`up_left`). Exact internal weights; physical Lorentz spin and flavor are unassigned.

**Down-type doublet component** (`down_left`). Exact internal weights; physical Lorentz spin and flavor are unassigned.

**Conjugate up-type singlet** (`up_right_conjugate`). Exact internal weights; physical Lorentz spin and flavor are unassigned.

**Conjugate down-type singlet** (`down_right_conjugate`). Exact internal weights; physical Lorentz spin and flavor are unassigned.

**Neutral lepton-doublet component** (`neutrino_left`). Exact internal weights; physical Lorentz spin and flavor are unassigned.

**Charged lepton-doublet component** (`charged_lepton_left`). Exact internal weights; physical Lorentz spin and flavor are unassigned.

**Conjugate charged-lepton singlet** (`charged_lepton_right_conjugate`). Exact internal weights; physical Lorentz spin and flavor are unassigned.

### Charged endpoint and pair records

**Charged regional endpoint** (`particle_primary_regional_fermion:regional_endpoint_irrep`). charged-lepton-singlet-type component; conjugate endpoint has opposite charge

**Compensated pair probe** (`particle_primary_regional_fermion:internal_pair_probe_irreps`). pair apparatus with retained charged endpoints; no single-photon assignment

### Neutral weak composites

**Weak-pair trion (direct)** (`particle_primary_composite_spin:weak_pair_charge_six:0,2,0,0,1`). neutral lepton-composite or sterile-like channel from LL e^c; not an identified active neutrino

**Weak-pair trion (conjugate)** (`particle_primary_composite_spin:weak_pair_charge_six:0,-2,0,0,-1`). conjugate neutral-lepton-composite or sterile-like channel; physical antiparticle identity unproved

### Neutral color composites

**Color-singlet trion (direct)** (`particle_primary_composite_spin:color_triple:1,0,0,2,0`). neutral-antibaryon (antineutron-type) composite channel from u^c(d^c)^2, with the opposite increment the neutron-type conjugate

**Color-singlet trion (conjugate)** (`particle_primary_composite_spin:color_triple:-1,0,0,-2,0`). neutral-baryon (neutron-type) conjugate channel; physical baryon identity unproved

### Neutral even composites

**Color determinant with singlet hole (direct)** (`particle_primary_composite_spin:color_determinant_charge_six_hole:0,0,0,3,-1`). neutral even composite channel; no selected named-particle hypothesis

**Color determinant with singlet hole (conjugate)** (`particle_primary_composite_spin:color_determinant_charge_six_hole:0,0,0,-3,1`). neutral even composite channel; no selected named-particle hypothesis

### Gauge-sector response

**Isolated odd-parity cube-scalar cluster** (`particle_primary_reduced_sector:isolated_A1u_cluster`). finite neutral collective gauge line; no named photon, gluon, Z or Higgs assignment

## Early-universe comparison and source time

The hot Big-Bang comparison depends on epoch and temperature. Around nucleosynthesis,
photons, electrons and positrons, neutrinos, protons and neutrons enter the thermal
history. At higher temperatures the Standard Model describes a quark-gluon plasma
with leptons and antiparticles; the populated degrees of freedom depend on temperature
and equilibration. This is not an observed inventory of the first instant.
See [PDG cosmology](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-bbang-cosmology.pdf)
and [CERN plasma experiments](https://home.cern/science/physics/heavy-ions-and-quark-gluon-plasma/).

The primary source state is stationary under its declared Hamiltonian. For a fixed
source effect E, [rho,H]=0 gives Tr(rho exp(iHt) E exp(-iHt)) = Tr(rho E).
This statement concerns the source marginal. Retained-reference and conditional
multi-time records can change. A first detector click is not a proof of a first
particle-production event. The refinement index is not cosmological time.

A formation history requires an initial preparation, a physical time/temperature
dictionary and complete production, annihilation and decay records, including
measurement backaction. Catalogue order, mass order and first-click order supply
none of these choices. See the [source-state proof](../particle_common_contract/PROOF.md)
and [candidate confirmation steps](CONFIRMATION_STEPS.md).
