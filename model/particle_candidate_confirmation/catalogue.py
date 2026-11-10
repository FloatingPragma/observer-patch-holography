"""Human labels preserve the domains and identifiers of the source catalogue."""
from copy import deepcopy


PROTECTED = {
    'exchange_cover': ('Exchange-cover memory', {
        'ordered_exchange_charge': ('Exchange-orientation bit',
            'A coherent orientation bit carries two protected charge sectors under the specified ordered-cover repair.'),
    }),
    'spin_and_memory': ('Spin, symmetry and memory sectors', {
        'selected_spin_charge': ('Spin-rotation charge',
            'A local spin operation creates two opposite conserved internal-charge channels. This charge is not identified with electric charge.'),
        'single_occupancy_spin_memory': ('One-particle-per-site spin memory',
            'Repair preserves the spin information in a fixed-occupation code. Identity repair supplies no winding invariant.'),
        'joint_positive_field_canonical_sector': ('Canonical spin fields',
            'A specified spin family admits two canonical field pairs with one positive-energy clock. Its physical particle dictionary is unassigned.'),
        'common_hidden_spin_factor': ('Shared spin subsystem',
            'A qutrit spin factor survives both the stated field dynamics and repair, with 105 multiplicity states.'),
        'cube_eta_range': ('Eta-pairing sector',
            'A rank-1820 symmetry sector is protected by the specified dynamics and repair. Its label is internal eta symmetry.'),
    }),
    'occupation_response': ('Occupation and current responses', {
        'one_hole_spin_readout': ('Mobile vacancy with a spin record',
            'One missing occupation can be located, while recorded transport carries the remaining spin information.'),
        'centered_full_CAR_charge_response': ('Conserved-number response family',
            'Particle addition, removal and spin probes share an interacting fermion algebra with conserved source occupation.'),
    }),
    'spatial_winding': ('Spatial winding textures', {
        'compact_opposite_texture_cores': ('Opposite winding cores',
            'A prescribed two-core texture has shell degrees +1 and -1 and four localized band modes. Spontaneous production from repair is not established.'),
    }),
    'momentum_winding': ('Momentum-space winding', {
        'kinetic_Weyl_nodes_and_free_fields': ('Oppositely wound Weyl nodes',
            'Eight kinetic nodes have four positive and four negative windings. A stated scaling family has a conditional free-field map; this does not identify eight observed massless species.'),
    }),
    'exclusion': ('Exclusion control', {
        'bare_native_source': ('Bare source has no nontrivial reducing charge',
            'The full branch commutant is scalar in the specified bare source-only domain. This entry is an exclusion, not a protected defect.'),
    }),
}

INTERNAL_LABELS = {
    'up_left': 'Up-type doublet component',
    'down_left': 'Down-type doublet component',
    'up_right_conjugate': 'Conjugate up-type singlet',
    'down_right_conjugate': 'Conjugate down-type singlet',
    'neutrino_left': 'Neutral lepton-doublet component',
    'charged_lepton_left': 'Charged lepton-doublet component',
    'charged_lepton_right_conjugate': 'Conjugate charged-lepton singlet',
}


def excitation_label(row):
    key = row['id']
    if key in INTERNAL_LABELS:
        return 'internal', 'Internal matter weights', INTERNAL_LABELS[key]
    if key.endswith(':regional_endpoint_irrep'):
        return 'charged', 'Charged endpoint and pair records', 'Charged regional endpoint'
    if key.endswith(':internal_pair_probe_irreps'):
        return 'charged', 'Charged endpoint and pair records', 'Compensated pair probe'
    if key.startswith('particle_primary_composite_spin:'):
        sign = 'conjugate' if any(n < 0 for n in row['exact']['native_occupation_increment'][:4]) else 'direct'
        if ':weak_pair_charge_six:' in key:
            return 'weak_composite', 'Neutral weak composites', 'Weak-pair trion (' + sign + ')'
        if ':color_triple:' in key:
            return 'color_composite', 'Neutral color composites', 'Color-singlet trion (' + sign + ')'
        return 'even_composite', 'Neutral even composites', 'Color determinant with singlet hole (' + sign + ')'
    if key == 'particle_primary_reduced_sector:isolated_A1u_cluster':
        return 'gauge', 'Gauge-sector response', 'Isolated odd-parity cube-scalar cluster'
    raise ValueError('unlabelled source excitation: ' + key)


def build(dictionary):
    originals = {row['id']: row for row in dictionary['protected_catalogue_rows']}
    covered = {key for _, rows in PROTECTED.values() for key in rows}
    if covered != set(originals):
        raise ValueError('protected catalogue coverage differs')
    protected = []
    for family, (label, rows) in PROTECTED.items():
        protected.append({'id': family, 'label': label, 'rows': [
            {'source_id': key, 'label': title, 'description': text,
             'source_record': deepcopy(originals[key]), 'evidence': originals[key]['evidence'],
             'physical_particle_id': None, 'first_physical_appearance': None}
            for key, (title, text) in rows.items()]})
    excitations = {}
    for row in dictionary['source_excitation_rows']:
        family, title, label = excitation_label(row)
        group = excitations.setdefault(family, {'id': family, 'label': title, 'rows': []})
        group['rows'].append({'source_id': row['id'], 'label': label,
            'description': row.get('best_guess', 'Exact internal weights; physical Lorentz spin and flavor are unassigned.'),
            'source_record': deepcopy(row), 'physical_particle_id': None,
            'first_physical_appearance': None})
    return {'protected_families': protected, 'excitation_families': list(excitations.values()),
            'physical_particle_map': [], 'first_appearance_order': None,
            'chronology_scope': 'Source stationarity is not a cosmological preparation or production history.'}


def render(value):
    data = value['labelled_catalogue']
    lines = ['# Protected objects and source excitation families', '',
        'The eleven protected-object entries retain their distinct construction domains. ',
        'They include conserved charges, protected memory, two kinds of winding and an ',
        'exclusion control. They are not eleven established topological particle species. ',
        'The sixteen excitation identifiers form a separate list. Their recorded ',
        '`protected_defect_id` fields are null: the defect-to-excitation attachment is unproved.', '',
        '## Protected-object families', '']
    for group in data['protected_families']:
        lines += ['### ' + group['label'], '']
        for row in group['rows']:
            source = row['source_record']
            lines += ['**' + row['label'] + '** (`' + row['source_id'] + '`). ' + row['description'], '',
                      'Domain: ' + source['source_domain'] + '.', '',
                      'Protected object: ' + source['protected_object'] + '.', '',
                      'Topology: ' + source['topology'] + '.', '',
                      '[Source proof](../../' + row['evidence'] + ').', '']
    lines += ['## Source excitation families', '',
              'Names below describe source components or conditional analogies. Repeated ',
              'electron, muon and tau comparison labels do not create three generations.', '']
    for group in data['excitation_families']:
        lines += ['### ' + group['label'], '']
        for row in group['rows']:
            lines += ['**' + row['label'] + '** (`' + row['source_id'] + '`). ' + row['description'], '']
    lines += ['## Early-universe comparison and source time', '',
        'The hot Big-Bang comparison depends on epoch and temperature. Around nucleosynthesis, ',
        'photons, electrons and positrons, neutrinos, protons and neutrons enter the thermal ',
        'history. At higher temperatures the Standard Model describes a quark-gluon plasma ',
        'with leptons and antiparticles; the populated degrees of freedom depend on temperature ',
        'and equilibration. This is not an observed inventory of the first instant. ',
        'See [PDG cosmology](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-bbang-cosmology.pdf) ',
        'and [CERN plasma experiments](https://home.cern/science/physics/heavy-ions-and-quark-gluon-plasma/).', '',
        'The primary source state is stationary under its declared Hamiltonian. For a fixed ',
        'source effect E, [rho,H]=0 gives Tr(rho exp(iHt) E exp(-iHt)) = Tr(rho E). ',
        'This statement concerns the source marginal. Retained-reference and conditional ',
        'multi-time records can change. A first detector click is not a proof of a first ',
        'particle-production event. The refinement index is not cosmological time.', '',
        'A formation history requires an initial preparation, a physical time/temperature ',
        'dictionary and complete production, annihilation and decay records, including ',
        'measurement backaction. Catalogue order, mass order and first-click order supply ',
        'none of these choices. See the [source-state proof](../particle_common_contract/PROOF.md) ',
        'and [candidate confirmation steps](CONFIRMATION_STEPS.md).', '']
    return '\n'.join(line.rstrip() for line in lines)
