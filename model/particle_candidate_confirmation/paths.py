"""Per-object particle attachments within the common field dependency map."""
from copy import deepcopy


PATHS = {
    'bare_native_source': {
        'result': 'The complete source-only branch commutant is scalar in the stated H66/H78 domain.',
        'classification': 'excluded_in_stated_domain',
        'first_required': 'Specify a different declared memory or field extension before proposing a nontrivial protected charge.',
        'remaining': ['Prove protection on the enlarged full process, including every outside and failure outcome.',
                      'Apply the source-operator, spectral and physical tests to that explicitly different construction.'],
        'shortcut': 'Do not search the excluded scalar commutant for a hidden particle label.',
    },
    'ordered_exchange_charge': {
        'result': 'The coherent ordered-cover bit gives two sharp repair-protected sectors.',
        'classification': 'protected_memory',
        'first_required': 'Attach the bit to a spatially accessible source excitation and a specified energy generator.',
        'remaining': ['Prove that its local creation and transport retain the charge under the complete repair and field process.',
                      'Calculate the original joint energy/momentum law and test a particle threshold or branch.',
                      'Identify physical exchange, currents, interactions and the mass observable.'],
        'shortcut': 'The exchange-cover character supplies a source label; its conservation does not supply a dispersion relation.',
    },
    'selected_spin_charge': {
        'result': 'The selected eight-site spin source has local creation, two sharp internal charge sectors, a positive finite energy response and a controlled detector lifetime.',
        'classification': 'finite_created_response',
        'first_required': 'Resolve the created charge channel in operational momentum on its selected state and clock.',
        'remaining': ['Separate the moving response from other spin channels and control volume, preparation and cutoff errors.',
                      'Prove the physical Lorentz/current dictionary and interaction records.',
                      'Extract a rest threshold in that dictionary; its finite selection-energy gap alone is insufficient.'],
        'shortcut': 'Reuse the existing local preparation and lifetime proof; no second proof of a nonzero finite response is needed.',
    },
    'single_occupancy_spin_memory': {
        'result': 'Complete repair acts identically on the one-fermion-per-site spin code.',
        'classification': 'repair_protected_algebra',
        'first_required': 'Choose a physical excitation within the code and bind its Hamiltonian, state and clock, including any leakage under evolution.',
        'remaining': ['Use the complete same-process record law to control leakage and controller energy.',
                      'Calculate a populated moving sector, its physical spin/current attachment and mass definition.'],
        'shortcut': 'An invariant information algebra is not a catalogue of its energy-carrying excitations.',
    },
    'joint_positive_field_canonical_sector': {
        'result': 'The specified positive-field spin scaling gives two regular canonical pairs on one ground-state family with a positive-energy clock.',
        'classification': 'constructed_field_family',
        'first_required': 'Determine the joint energy/momentum spectral measure and stiffness of these actual canonical fields.',
        'remaining': ['Control the selected-field scaling and identify whether the response is massless, massive or nonrelativistic.',
                      'Prove physical Lorentz transformations and interactions of the populated sectors.',
                      'Identify any observable species and its mass without replacing this state or clock.'],
        'shortcut': 'Canonical commutation and positive energy are supplied; a numerical stiffness or relativistic mass shell is not.',
    },
    'common_hidden_spin_factor': {
        'result': 'The same chain field and repair preserve a spin factor with H=h tensor I. Complete finite spectra and candidate probes are available in their declared reference.',
        'classification': 'protected_internal_factor',
        'first_required': 'Select an operator that changes the energy-carrying multiplicity or number sector, rather than only rotating the protected spin.',
        'remaining': ['Control energy/work during repair: conservation of the spin factor does not imply conservation of the joined Hamiltonian.',
                      'Establish spatial propagation and a controlled higher-dimensional extension, respecting the proved frame-alignment condition.',
                      'Identify Lorentz spin and physical currents before extracting a particle mass.'],
        'shortcut': 'Every operator I tensor A on the isolated protected spin factor commutes with H and has only zero ideal energy transfer.',
    },
    'cube_eta_range': {
        'result': 'The eta>=2 range is exactly protected; H_6 V=V(H_4+U) and the prepared fixed-number local probe has eight interacting spectral lines.',
        'classification': 'protected_model_spectroscopy',
        'first_required': 'Attach the existing local spectroscopy to a moving branch, including complete detector and evolution implementation errors.',
        'remaining': ['Keep one reference and clock: the eta ladder shifts H energy by U while changing occupation by two.',
                      'Control growing-volume dispersion and physical currents/spin, including repair and compensating occupations.',
                      'Extract the branch mass; U or U/2 is not thereby an elementary rest mass.'],
        'shortcut': 'Reuse the exact eta intertwiner and local resolvent; the model spectral lines need no guessed particle names.',
    },
    'compact_opposite_texture_cores': {
        'result': 'A specified opposite-winding texture has finite energy and four isolated fermion modes with an exterior-volume bound.',
        'classification': 'specified_texture_and_bound_modes',
        'first_required': 'Specify and solve the dynamical creation and motion of the texture itself in the declared source.',
        'remaining': ['Prove stability under the admitted dynamics and retain the opposite core and boundary energy.',
                      'Distinguish the energy and inertia of the texture from the occupations and eigenenergies of its bound fermion modes.',
                      'Calculate moving dispersion, quantum statistics and physical interactions of the selected excitation.'],
        'shortcut': 'The existing Feshbach isolation bounds control fermion modes in the supplied background; they do not derive an autonomous texture or its mass.',
    },
    'one_hole_spin_readout': {
        'result': 'Recorded repair transports a vacancy and its spin information in the one-hole sector.',
        'classification': 'finite_transport_records',
        'first_required': 'Bind an energy generator and original joint energy/momentum response to that same recorded transport.',
        'remaining': ['Retain all holds, collisions and controller work while proving a stable propagating channel.',
                      'Identify its physical charge, Lorentz spin and interaction records, then its rest threshold.'],
        'shortcut': 'A position record and a conserved occupation deficit do not fix an energy-frequency relation.',
    },
    'centered_full_CAR_charge_response': {
        'result': 'Odd addition/removal and even spin responses share the interacting full fermion algebra and conserved source occupation.',
        'classification': 'interacting_response_family',
        'first_required': 'Prove positive original weight in a fixed energy/momentum branch or charged threshold of this centered family.',
        'remaining': ['Control dressed charged sectors and the actual capacity dependence instead of substituting another ground state.',
                      'Attach physical Lorentz/current transformations and interaction vertices.',
                      'Calculate a controlled threshold, resonance or matched mass in a common clock.'],
        'shortcut': 'Reuse the full algebra and current Ward identities; algebraic existence alone supplies no positive particle residue.',
    },
    'kinetic_Weyl_nodes_and_free_fields': {
        'result': 'Eight winding nodes map to eight massless Weyl fields on the stated weak-interaction scaling family, with a controlled kinetic error.',
        'classification': 'conditional_free_field_map',
        'first_required': 'Prove the physical identification and interaction content of the actual desired family, retaining all opposite-chirality nodes.',
        'remaining': ['For a massive branch, determine the mass-generating coefficients or order from that declared dynamics and state.',
                      'Establish populated interacting spectral branches, physical currents and the fate of mirror modes.',
                      'Control refinement and map the resulting thresholds to observed species.'],
        'shortcut': 'The free limit and supplied mass-channel formulas are available; the default interacting coupling lies outside the proved free-limit hypotheses.',
    },
}

REUSE = [
    {'path': 'Lean/Dynamics/ProtectedCharge.lean',
     'result': 'chargeExpectation_preserved_iff_dual_fixed',
     'premise': 'The actual finite channel and its expectation dual obey the stated pairing.',
     'conclusion': 'All-state conservation is equivalent to a dual fixed observable.',
     'boundary': 'No physical charge, energy generator or particle identity follows.'},
    {'path': 'Lean/QFT/GroundRecordTail.lean',
     'result': 'downwardWeight_le_excess_div',
     'premise': 'The actual nonnegative excitation energies, complete contractive records and ground-excess bound.',
     'conclusion': 'A controlled probability tail for records that lower energy by a fixed threshold.',
     'boundary': 'The positive-energy field construction and a populated particle branch require their own additional premises.'},
    {'path': 'Lean/QFT/PeriodicMomentum.lean',
     'result': 'polynomial_conserved',
     'premise': 'The actual Hamiltonian commutes with the declared translation matrix.',
     'conclusion': 'Every polynomial of that translation matrix is conserved.',
     'boundary': 'A finite open box does not acquire periodic translation symmetry from this theorem.'},
    {'path': 'Lean/Geometry/MassShellKinematics.lean',
     'result': 'fourMomentum_time_sq_eq_mass_sq_add_spatial',
     'premise': 'A declared mass parameter and frame in the Lorentz module.',
     'conclusion': 'The kinematic mass-shell identity for the scaled frame vector.',
     'boundary': 'The theorem does not calculate mass or identify the vector with the source energy/momentum readout.'},
]


def build(dictionary):
    rows = dictionary['protected_catalogue_rows']
    if set(PATHS) != {r['id'] for r in rows}:
        raise ValueError('per-defect path coverage differs')
    return {
        'scope': 'Specializations of the existing six confirmation tests and ten-node dependency map.',
        'route_rule': 'Finish the stated construction directly, or prove a transfer to the primary source before combining their records. Architectural uniqueness is not a prerequisite.',
        'rows': [{'source_id': r['id'], 'source_domain': r['source_domain'],
                  'evidence': r['evidence'], **deepcopy(PATHS[r['id']]),
                  'physical_particle_id': None, 'physical_rest_mass': None} for r in rows],
        'finite_theorem_reuse': deepcopy(REUSE),
        'source_energy_shortcuts': {
            'protected_factor': {'identity': 'H=h tensor I; A_spin=I tensor A',
                'consequence': '[H,A_spin]=0; ideal transfer support is {0}',
                'source': 'model/common_spin_protection/README.md',
                'physical_mass_zero_inferred': False},
            'eta_intertwiner': {'identity': 'H_6 V=V(H_4+U); V^dagger V=I_on_N4',
                'consequence': 'V maps H_4 energy E to H_6 energy E+U',
                'occupation_change': 2, 'source': 'model/cube_eta_protection/README.md',
                'elementary_mass_U_or_U_over_two_inferred': False},
        },
    }


def render(value):
    data = value['protected_particle_paths']
    lines = ['# From each protected object to a particle mass', '', data['scope'], '',
             data['route_rule'], '',
             'Protection, finite response, a populated field branch and a named physical species ',
             'are separate assertions. The [six confirmation tests](CONFIRMATION_STEPS.md) ',
             'supply the common physical criteria. The paths below retain the source of every ',
             'protected-object entry; none is assigned a physical rest mass.', '']
    for row in data['rows']:
        lines += ['## ' + row['source_id'], '', row['result'], '',
                  '**Construction:** ' + row['source_domain'] + '.', '',
                  '**First required attachment:** ' + row['first_required'], '',
                  '**Reusable result:** ' + row['shortcut'], '',
                  'The remaining calculation follows this order:', '']
        lines += [f'{i}. {text}' for i, text in enumerate(row['remaining'], 1)]
        lines += ['', '[Owning source result](../../' + row['evidence'] + ').', '']
    lines += ['## Finite theorem premises', '',
              '| Existing theorem | Actual premise to supply | Conclusion and boundary |',
              '| --- | --- | --- |']
    for r in data['finite_theorem_reuse']:
        lines.append('| [' + r['result'] + '](../../' + r['path'] + ') | ' + r['premise']
                     + ' | ' + r['conclusion'] + ' ' + r['boundary'] + ' |')
    lines += ['', '## Why the two exact energy shortcuts do not give a rest mass', '',
        'For the protected spin factor, [h tensor I, I tensor A]=0. The spectral ',
        'projectors of H therefore commute with this spin-only probe. Its ideal ',
        'initial/final energy record has zero transfer under H evolution without ',
        'intervening repair. Other source ',
        'operators can change the energy-carrying multiplicity or number sector. ',
        'This does not imply a massless particle.', '',
        'For the eta isometry, applying H_6 V=V(H_4+U) to an H_4 eigenvector ',
        'with energy E proves the shifted eigenvalue E+U. This is an exact model ',
        'transition across two occupation sectors. Replacing H by H-(U/2)N makes ',
        'that shift zero and changes the generator on charge-changing records. ',
        'Across these sectors it is not a common additive energy zero. Within one ',
        'fixed-number sector the subtraction is a constant, so relative spectroscopy ',
        'is unchanged when the reference is shifted consistently. Neither convention ',
        'supplies a relativistic moving-particle ',
        'threshold or a physical mass by itself.', '']
    return '\n'.join(line.rstrip() for line in lines)
