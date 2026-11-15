"""Bind finite observer claims to executable sources and a recomputed assessment."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from tools.check_openwave_assessment_replay import compare

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT/'code/openwave_candidate'
OUTPUT = ROOT/'claims/openwave_joint_model.json'


def separately_validated_native_charge():
    path = ROOT/'claims/openwave_native_charge.json'
    data = json.loads(path.read_text())
    if data.get('full_criterion_met') is not True or data.get('row') != 1:
        raise ValueError('native charge criterion receipt absent')
    for relative, expected in data['source_sha256'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != expected:
            raise ValueError('stale native charge validation input: '+relative)
    return {'candidate_id': data['candidate_id'],
            'validated_rows': [1], 'ledger': str(path.relative_to(ROOT)),
            'ledger_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'scope': data['scope'],
            'transfer_to_this_canonical_source': False,
            'official_OpenWave_awards_added': 0}


def separately_validated_wave_rows():
    path = ROOT/'claims/openwave_wave_rows.json'
    data = json.loads(path.read_text())
    for relative, expected in data['source_sha256'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != expected:
            raise ValueError('stale wave validation input: '+relative)
    if data['same_source_transfer_between_rows'] is not False:
        raise ValueError('unproved wave source transfer')
    records = []
    for row in data['rows']:
        if row['full_criterion_met'] is not True:
            raise ValueError('incomplete wave criterion')
        records.append({'candidate_id': row['candidate_id'],
                        'validated_rows': [row['row']], 'ledger': str(path.relative_to(ROOT)),
                        'ledger_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                        'scope': row['scope'], 'transfer_to_this_canonical_source': False,
                        'official_OpenWave_awards_added': 0})
    return records


def separately_validated_common_charge():
    path = ROOT/'claims/openwave_common_charge.json'
    data = json.loads(path.read_text())
    if data.get('full_criterion_met') is not True or data.get('row') != 1:
        raise ValueError('common-source charge criterion receipt absent')
    for relative, expected in data['source_sha256'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != expected:
            raise ValueError('stale common-source charge input: '+relative)
    return {'candidate_id': data['candidate_id'], 'validated_rows': [1],
            'ledger': str(path.relative_to(ROOT)),
            'ledger_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'scope': data['scope'], 'transfer_to_this_canonical_source': False,
            'duplicates_existing_row_number': True,
            'official_OpenWave_awards_added': 0}


def observable_receipt(relative):
    path=BUNDLE/relative
    result=json.loads(path.read_text())
    no_award=(type(result.get('official_rows_promoted')) is int and result['official_rows_promoted']==0) or ('official_rows_promoted' not in result and result.get('official_rows_changed')==[] and result.get('full_criterion_met') is False)
    if not no_award:
        raise ValueError('observable experiment cannot award an official row')
    pins=result.get('source_sha256',{})
    if not pins:
        raise ValueError('observable source custody absent')
    for rel,expected in pins.items():
        source=(BUNDLE/rel).resolve()
        if not source.is_relative_to(BUNDLE.resolve()):
            raise ValueError('observable source outside bundle')
        if hashlib.sha256(source.read_bytes()).hexdigest()!=expected:
            raise ValueError('stale observable source: '+rel)
    return result


def c8_finite_claims():
    """Finite mechanisms instantiated on the declared active source; no awards."""
    return [
      {'id':'native_pinned_charge_spectrum','statement':'A position-only E2a pinned-charge energy is added to the same uniform-fiber H_m. Complete canonical1x1 spectra and selected source-overlap eigenstates retain negative branches; a complete modular triangle checks the charge-energy lift and every native Kraus map.','domain':'Supplied source, mass .23, coupling .29 and half-transit charge. Coherent eigenstate stationarity fails under the tested native monitoring. No physical atom, scattering threshold, hydrogen series or H0=L clock substitution.'},
      {'id':'native_position_current_readout','statement':'The same Hamiltonian supplies the exact regional continuity operator J=i[H,P]. Complete native outcomes account separately for committed regional-charge changes, and a declared complete current spectral PVM reproduces its expectation.','domain':'Whole port/transit position charge, distinct from endpoint-half energy charge. Supplied current apparatus; no C/E-only synthesis, physical electric/spin current label or reduced postmeasurement continuation inferred.'},
      {'id':'shared_many_load_native','statement':'A complete nonnegative Bose occupation sector with shared modular half-records implements signed port/transit hopping and literal occupation-balancing C/E instruments. Its one-load cyclic control equals the original Hamiltonian and every native Kraus map; three-unit tree controls retain complete outcomes and fixed Gauss residues.','domain':'Declared Bose statistics, fixed total, finite BFS sector and modulus. Modulus greater than twice an admitted absolute charge bound gives unique signed cut decoding; a net signed load alone does not. No signed sea, independent-unit commit substitution or many-load uniform-fiber theorem.'},
      {'id':'declared_uniform_fiber_preparation','statement':'Added dimension-N Fourier gates and reversible tree Gauss arithmetic prepare the uniform record fibers with ideal success one, including composite modulus. Complete tiny Gauss sectors and reversible work-register controls instantiate the prescription.','domain':'Additional preparation resources, not native-only C/E evolution; no hardware noise, efficient gate-count bound or general complexity lower bound asserted.'},
      {'id':'native_current_maxwell_bridge','statement':'The actual coherent native expectation current supplies a straight-edge Fourier source for the existing Maxwell mode executor. Exact instantaneous half-transit charge continuity and held-source Gauss/work identities pass independent endpoint, quadrature and augmented-exponential checks.', 'domain':'One-way supplied linear field and embedding, zero commit rate, finite source-hold step. The endpoint native charge mismatch is retained and compared under refinement. No quantum source-field state, reciprocal backreaction, outgoing radiation, derived dispersion or physical calibration.'},
      {'id':'continuing_record_apparatus','statement':'Individual half-record measurements retain every full poststate, including its code-complement coherence, for subsequent full native evolution. Whole-position closed-cut measurement branches intertwine with the reduced code.','domain':'Tiny complete modular apparatus controls. Code/complement weights are diagnostics rather than an extra dephasing measurement. Arbitrary terminal effects do not imply closed reduced continuation.'}]


def c9_finite_claims():
    return [
      {'id':'shared_car_native_instrument',
       'statement':'Direct full-occupation C/E/MIX Kraus maps on shared modular records form complete channels with CAR hopping and stated endpoint transport phases. Every one-load Hamiltonian and native Kraus map recovers the original model; complete work and external-reference controls are retained.',
       'domain':'Declared fixed-number Pauli/CAR statistical extension. This is not the non-trace-preserving exterior lift of each one-load Kraus operator. Local gauge conservation is modular; no physical spin, species or derived statistics.'},
      {'id':'car_all_character_transfer',
       'statement':'Normalized Gauss-fiber character isometries decompose the entire finite shared-record CAR source and intertwine its Hamiltonian and every C/E/MIX Kraus map. Uniform fibers are the zero-character member, not the entire source.',
       'domain':'Complete triangle sectors at the declared modulus and fixed number, with independent tree/cycle proofs and Fourier reconstructions. Character preparation is supplied. No unconditional global-ground inference from the zero-character Slater state.'},
      {'id':'car_global_reference',
       'statement':'The actual full fixed-number Hamiltonian has a positive excitation reference after subtracting its minimum and retaining the complete ground multiplet. Z3 has a ground doublet in nonzero characters; Z6 has a unique character-3 minimum on the tested triangle.',
       'domain':'Finite declared CAR source, mass .23, no potential. The lowest positive one-load band edge minimizes character 0, a different energy selection. Native channels preserve character and can deposit work; energetic preference does not implement cooling, character selection or a physical vacuum.'},
      {'id':'car_native_energy_response_support',
       'statement':'A complete ground/native-letter/final-energy experiment retains every energy eigenspace and all zero outcomes. Its supported excitation threshold is distinguished from the full spectral gap between invariant character sectors.',
       'domain':'Declared finite ground preparation, native instrument and energy PVM. Algebraic character exclusion and numerical within-sector support tolerance are separate. Neither spectral gap is a physical particle mass or GeV prediction.'}]


def c10_finite_claims():
    return [
      {'id':'one_interacting_all_number_source',
       'statement':'One complete 384-state shared-record CAR source carries the native instrument, periodic electric and cyclic face terms, integer boundary reference and charged reservoir. Its all-number spectra and unique N=3 ground reference are recomputed from that interacting Hamiltonian.',
       'domain':'Declared six-orbital triangle, modulus6, mass .23, electric coefficient .37 and face coefficient .61. Full states are retained because the electric term mixes the zero-field character sectors. No zero-field gap, physical vacuum, particle species or continuum field result is transferred.'},
      {'id':'same_source_complete_charged_response',
       'statement':'The boundary-dressed CAR reservoir swap retains both initial settings, every addition/removal/no-change outcome, full subnormalized source states and complete source-number/energy probabilities. All native and reservoir work uses the same interacting Hamiltonian and reference.',
       'domain':'Fair reservoir occupation, energy .23 times its occupation, zero extra boundary energy and ideal pulse are supplied apparatus choices. The pulse preserves unconditioned field energy; Hamiltonian exchange and monitoring work remain distinct. No autonomous battery, pair radiation or charge-conjugation symmetry of the fixed Gauss background is inferred.'},
      {'id':'continuing_degenerate_energy_instrument',
       'statement':'A complete grouped number/energy Lüders instrument retains within-eigenspace coherences as compact full coefficient matrices and every zero outcome. An actually populated degenerate branch continues under the same Hamiltonian and selected native instrument with its original subnormalization.',
       'domain':'Ideal spectral apparatus and explicit grouping tolerance. Degenerate basis rotation preserves the reconstructed state; rank-one measurement followed by forgetting is rejected by a coherent doublet control. The apparatus is not claimed synthesized by native C/E, and measured energy labels have no physical species or GeV identification.'}]


def c11_finite_claims():
    return [
      {'id':'all_number_boundary_charge_apparatus',
       'statement':'On the complete interacting carrier, known total occupation bounds each regional Pauli charge in an interval of width at most3. Modulus6 and the integer boundary record therefore reconstruct every regional charge; all complete boundary/cut projectors equal the corresponding number/charge projectors.',
       'domain':'All384 basis states,64 regions and2688 joint outcomes. Reading the boundary is part of the apparatus and disturbs cross-number coherence even at fixed fermion parity; it cannot be silently replaced by charge-only readout. No physical winding or exclusion of fractional species beyond the declared occupation domain follows.'},
      {'id':'interacting_external_phase_current',
       'statement':'An explicit external half-link phase probes the same interacting Hamiltonian. Its derivative current, regional continuity, static energy response and complete native charge/work ledgers satisfy the declared gauge covariance, including subnormalized continuation.',
       'domain':'Supplied external Peierls phase and gauge-family intervention on the finite source. Native record dressing is a scalar phase per classical branch, not a new branch channel. No physical magnetic field, spin, g-factor or internal continuous gauge algebra at cyclic wrap is inferred.'},
      {'id':'interacting_conditioned_modular_observer',
       'statement':'The actual selected interacting output has a faithful declared two-configuration corner with nontrivial modular correlations. A complete P/Q and sequential projective apparatus realizes the real correlation and retains every branch through subsequent native evolution.',
       'domain':'The corner is chosen after inspecting the known output and is an exploratory global configuration algebra, not a local spatial partial trace or a predeclared physical observer. Full-source and pure observed-branch faithfulness are rejected. The observable family is frozen from the premeasurement state and its modular parameter is injected; no autonomous H-clock, subjective time or physical time conversion is identified.'},
      {'id':'complete_charged_energy_unit_display',
       'statement':'Every already executed charged-response energy group and unconditional probability is retained under one common source-energy-to-GeV scale, with the actual interacting ground subtracted once and the complete work budget converted by that same scale.',
       'domain':'The most probable addition line is assigned the supplied0.000511GeV anchor. Its agreement is a calibration identity. Other lines remain unnamed source excitations, not predicted particle masses; source charge conjugation, rest-shell and species identification remain open.'}]


def generate():
    receipt = BUNDLE/'integration/assessment.json'
    assessment = json.loads(receipt.read_text())
    required={'integration/assess.py','platform_packet/packet.json'}
    required.update(str(p.relative_to(BUNDLE)) for p in (BUNDLE/'candidate/oph_openwave').glob('*.py'))
    if set(assessment['sources']) != required:
        raise ValueError('incomplete assessment source custody')
    for rel, expected in assessment['sources'].items():
        if hashlib.sha256((BUNDLE/rel).read_bytes()).hexdigest()!=expected:
            raise ValueError('stale assessment source: '+rel)
    with tempfile.TemporaryDirectory() as directory:
        replay = Path(directory)/'assessment.json'
        subprocess.run([sys.executable, str(BUNDLE/'integration/assess.py'), '--output', str(replay)], check=True, capture_output=True, text=True)
        compare(assessment, json.loads(replay.read_text()))
    rows = assessment['rows']
    if sorted(r['id'] for r in rows) != list(range(1,32)) or any(r['full_criterion_met'] for r in rows):
        raise ValueError('incomplete or promoted finite assessment')
    paths = [Path(__file__).resolve(), ROOT/'tools/check_openwave_assessment_replay.py', receipt, BUNDLE/'integration/assess.py', BUNDLE/'platform_packet/packet.json',
             BUNDLE/'platform_packet/criteria_frozen.json']
    paths += sorted((BUNDLE/'candidate/oph_openwave').glob('*.py'))
    experiment_root=BUNDLE/'validation/observables'
    modular=observable_receipt('validation/observables/modular_time/MODULAR_TIME.json')
    spatial=observable_receipt('validation/observables/spacetime/STACKED_WAVEPACKET.json')
    energy=observable_receipt('validation/observables/energy_gev/ENERGY_GEV.json')
    apparatus=observable_receipt('validation/observables/energy_gev/SPECTRAL_APPARATUS.json')
    if not all(v is True for result in (energy,apparatus) for v in result['checks'].values()):
        raise ValueError('native energy/apparatus check failed')
    if not all(v is True for v in modular['checks'].values()) or not all(v is True for v in modular['mutations'].values()):
        raise ValueError('modular observer check failed')
    archive=experiment_root/'modular_time'/modular['native_state_archive']['filename']
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=modular['native_state_archive']['sha256']:
        raise ValueError('native observer state archive changed')
    bridge_root=BUNDLE/'validation/bridges'
    bridge=json.loads((bridge_root/'NATIVE_MAXWELL_SOURCE.json').read_text())
    if bridge.get('full_physical_criterion_met') is not False or not bridge.get('source_sha256') or not bridge.get('checks') or not all(v is True for v in bridge['checks'].values()):
        raise ValueError('native current/field bridge scope or checks changed')
    for rel,expected in bridge['source_sha256'].items():
        source=(BUNDLE/rel).resolve()
        if not source.is_relative_to(BUNDLE.resolve()) or hashlib.sha256(source.read_bytes()).hexdigest()!=expected:
            raise ValueError('stale or external current/field bridge source: '+rel)
    paths += sorted(p for p in bridge_root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and '.pytest_cache' not in p.parts)
    paths += sorted(p for p in experiment_root.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                    and '.pytest_cache' not in p.parts)
    car_root=BUNDLE/'next/C9'
    car_receipts={name:observable_receipt('next/C9/'+name+'.json') for name in
                  ('FERMIONIC_SHARED_RECORDS','FERMIONIC_UNIFORM_REFERENCE','CHARACTER_REFERENCE','CHARACTER_MINIMA')}
    for name,result in car_receipts.items():
        if not result.get('checks') or not all(value is True for value in result['checks'].values()):
            raise ValueError('CAR source or reference check failed: '+name)
    paths += sorted(p for p in car_root.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                    and '.pytest_cache' not in p.parts)
    joint_root=BUNDLE/'next/C10'
    joint_receipts={name:observable_receipt('next/C10/'+name+'.json') for name in
                    ('JOINT_SOURCE','CANONICAL_CONSTRUCTION','CANONICAL_SOURCE_CONTRACT')}
    joint=joint_receipts['JOINT_SOURCE']; canonical=joint_receipts['CANONICAL_CONSTRUCTION']
    contract=joint_receipts['CANONICAL_SOURCE_CONTRACT']
    for name,result in (('joint',joint),('canonical',canonical)):
        if not result.get('checks') or not all(value is True for value in result['checks'].values()):
            raise ValueError('interacting source or continued apparatus check failed: '+name)
    if canonical['construction_id'] != contract['construction_id']:
        raise ValueError('canonical source identity differs from its contract')
    archive=joint_root/joint['state_archive']['filename']
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=joint['state_archive']['sha256']:
        raise ValueError('interacting source archive changed')
    paths += sorted(p for p in joint_root.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                    and '.pytest_cache' not in p.parts)
    continued_root=BUNDLE/'next/C11'
    continued_receipts={name:observable_receipt('next/C11/'+name+'.json') for name in
                        ('BOUNDARY_CHARGE','LINK_PHASE_PROBE','INTERACTING_MODULAR_OBSERVER','CHARGED_ENERGY_UNITS')}
    for name,result in continued_receipts.items():
        if not result.get('checks') or not all(value is True for value in result['checks'].values()):
            raise ValueError('continued common-source observable failed: '+name)
        if not all(value is True for value in result.get('mutations',{}).values()):
            raise ValueError('continued observer rejecting control failed: '+name)
        if result.get('construction_id',contract['construction_id'])!=contract['construction_id']:
            raise ValueError('continued observer source identity differs: '+name)
        if 'state_archive' in result:
            archive=continued_root/result['state_archive']['filename']
            if hashlib.sha256(archive.read_bytes()).hexdigest()!=result['state_archive']['sha256']:
                raise ValueError('continued observer archive changed: '+name)
    paths += sorted(p for p in continued_root.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                    and '.pytest_cache' not in p.parts)
    return {
      'schema':'oph.openwave_joint_model_claims.v1',
      'separately_validated_candidates':[separately_validated_native_charge(), separately_validated_common_charge(), *separately_validated_wave_rows()],
      'construction':'One declared all-number shared-record CAR source with periodic half-record electric energy, cyclic face term, recomputed interacting reference, complete native instrument and boundary-dressed charged reservoir. The canonical entry point binds one selected protocol and its continued grouped spectral apparatus. Earlier one-load, zero-field, bounded-integer, spatial and modular-time experiments retain their comparison domains.',
      'canonical_source':{'contract':'code/openwave_candidate/next/C10/CANONICAL_SOURCE_CONTRACT.json',
                          'receipt':'code/openwave_candidate/next/C10/CANONICAL_CONSTRUCTION.json',
                          'construction_id':contract['construction_id'],
                          'selected':canonical['selected'],
                          'scope':canonical['scope'],
                          'earlier_31_row_assessment_domain':'The earlier one-load assessment is preserved below. Its measured outcomes are not reassigned to this interacting all-number source.'},
      'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
      'assessment':str(receipt.relative_to(ROOT)),
      'canonical_continued_observables':[
        {'receipt':'code/openwave_candidate/next/C11/'+name+'.json',
         'checks':result['checks'],'scope':result['scope'],'physical_rows_promoted':0}
        for name,result in continued_receipts.items()],
      'finite_extension_experiments':[
        {'receipt':'code/openwave_candidate/next/C9/'+name+'.json',
         'parameters':result['parameters'],'checks':result['checks'],
         'scope':result['scope'],'physical_rows_promoted':0}
        for name,result in car_receipts.items()]+[
        {'receipt':'code/openwave_candidate/next/C10/JOINT_SOURCE.json',
         'parameters':joint['parameters'],'checks':joint['checks'],
         'scope':joint['scope'],'physical_rows_promoted':0},
        {'receipt':'code/openwave_candidate/next/C10/CANONICAL_CONSTRUCTION.json',
         'parameters':canonical['selected'],'checks':canonical['checks'],
         'scope':canonical['scope'],'physical_rows_promoted':0}],
      'observable_experiments':[
        {'id':'native_conditioned_modular_flow',
         'statement':'A declared two-port corner of the canonical native output is numerically faithful and has nontrivial modular observable correlations, checked with the standard finite GNS modular operator and a separate analytic reconstruction.',
         'receipt':'code/openwave_candidate/validation/observables/modular_time/MODULAR_TIME.json',
         'frequency_dimensionless':modular['modular_flow']['dimensionless_modular_frequency'],
         'inside_probability':modular['observer']['corner_probability'],
         'domain':'Supplied uniform Gauss-fiber preparation, native duration/rate, terminal P/Q instrument, conditioning and state-dependent observable. Full native density and both terminal branches are archived. Native numerical propagation has no certified enclosure. No physical or subjective time conversion, consciousness interpretation or raw C/E realization of the complex correlation is inferred.',
         'physical_rows_promoted':0},
        {'id':'explicit_stacked_ring_wavepacket',
         'statement':'An explicitly augmented ring graph has rank-three acoustic propagation. A prepared positive-band packet evolves under its full H_m at zero commit rate; all bins of the declared coarse spatial apparatus are retained and compared with a constant-velocity reference geodesic.',
         'receipt':'code/openwave_candidate/validation/observables/spacetime/STACKED_WAVEPACKET.json',
         'weighted_group_velocity':spatial['observables']['weighted_group_velocity'],
         'centroid_reference_max_distance':spatial['observables']['max_centroid_ray_distance'],
         'domain':'Supplied vertical edges, unit spacing/weights, periodic lattice, positive acoustic preparation, finite-grid phase convention, coarse position apparatus and flat anisotropic metric. This is an extension of the two-dimensional screen. The broad ensemble centroid is not an observed individual trajectory. No physical length/time calibration, dynamical gravity, modular/proper-time identification or transfer of the two-dimensional force kernel is asserted.',
         'physical_rows_promoted':0},
        {'id':'native_energy_and_declared_units','statement':'The complete canonical native spectrum and supplied interferometric projectors read positive/negative energy gaps. A separate complete sequential initial-energy/native-event/final-energy instrument retains every reference label and energy-conserving supported battery-shell column.',
         'receipts':['code/openwave_candidate/validation/observables/energy_gev/ENERGY_GEV.json','code/openwave_candidate/validation/observables/energy_gev/SPECTRAL_APPARATUS.json'],
         'domain':'Spectral preparation/readout and phase reference are declared resources. GeV values use a supplied511keV half-gap calibration; no electron mass prediction or vacuum/species identification. Sequential energy measurements change reference correlations and differ from the original unmeasured-reference spectral coupling. Supported battery isometry is not a reusable autonomous battery.',
         'physical_rows_promoted':0}],
      'finite_claims':[
        {'id':'additive_independent_loads','statement':'Independent prepared loads evolve under separate copies of the port/transit generator.','domain':'Product preparations; no interacting or exchange-statistics theorem.'},
        {'id':'integer_joint_gauss','statement':'Half-seam hopping and accepted literal native currents preserve the finite Gauss sector.','domain':'One load; bounded integer records; hard-wall coherent regulator.'},
        {'id':'monitored_native_instrument','statement':'The declared measured event instrument preserves trace and positivity with complete native, transit-noop and rejected-attempt outcomes.','domain':'Position and half-seam records measured; supplied discrete time and commit probability.'},
        {'id':'observer_complete_readout','statement':'The finite observer readout retains position probabilities, monitored-attempt outcomes and applied-current probabilities; sparse evolution agrees with finite density controls.', 'domain':'One-load finite sector. The carrier path-window error bound concerns coherent evolution of the specified initial basis preparation against the complete sector at the same cutoff; no repeated-channel or cutoff-removal bound.'},
        {'id':'bare_dark_mode_nontransfer','statement':'The oriented triangle bare incidence has rank 2 and transit kernel dimension 1; the cutoff-1 integer-record incidence in the initial-port Gauss sector has rank 6 and transit kernel dimension 0. The declared equal-amplitude BFS transit preparation leaks to ports.', 'domain':'Exact finite triangle rank control and the explicitly supplied BFS preparation at time .37 and mass .23. No transfer of bare dark-mode stability or particle-stability identification.'},
        {'id':'bounded_preparation_nonintertwiner','statement':'Uniform bounded integer fibers are isometric but fail to intertwine the declared bare and joint generators; a distinct complete modular control intertwines exactly.', 'domain':'Declared finite triangle controls. Modular and bounded integer state spaces are distinct constructions; no cutoff-limit or state-space equivalence follows.'},
        {'id':'integer_record_energy','statement':'The complete triangle unit-dipole integer family has quadratic minimum 1/2, above the real-flux minimum 1/3. First-order spectral energy removes cycle cost, while its massive variant charges cycles.', 'domain':'Supplied spectral quadratic forms and the declared effective half-record transit charge. No Coulomb force, physical confinement or selected interaction follows.'},
        {'id':'selected_interacting_instrument','statement':'Explicit selected record-energy diagonals alter the same finite observer Hamiltonian. Coherent evolution preserves its energy and complete monitored outcomes account for the resulting energy change.', 'domain':'Declared one-load first-order and massive variants. Monitoring generally changes system energy; no work/heat attribution, energy-conserving environment, or free-square identity is asserted.'},
        {'id':'static_flux_spin_nonidentification','statement':'Finite bare compact modes and antipodal-flux doublets pass direct residual controls. Rotation rephasing changes the cube without changing covariance, and adjacent flux removes the doublet.', 'domain':'Supplied bare static flux controls. No identified joint particle, physical spin, or exchange-statistics theorem.'},
        {'id':'regional_charge_instrument','statement':'The complete regional charge instrument measures sharp one-load outcomes 0 and 1, retains the complement and optional port/transit refinement, and obeys exact regional Gauss and coherent/committed continuity on the bounded integer construction.', 'domain':'Declared subdivided-vertex regions and one-load integer sector, including explicit interacting variants. Fractional expectation is not a sharp fractional outcome; no physical species or unrestricted fragmentation exclusion follows.'},
        {'id':'canonical_general_cut','statement':'The canonical ring half-incidence identity cancels every internal edge on an arbitrary region boundary. Integer circulation aliases modulo the register size but contributes zero to a complete cut.', 'domain':'Exact diagonal operator identity. Modular data alone do not distinguish integer Gauss sectors; integer endpoint-current interpretation requires separately justified source-sector ancestry.'},
        {'id':'local_coarse_modular_instrument','statement':'The separately named modular local-coarse C/E Kraus instrument preserves trace, positivity, its Gauss residues and an external quantum reference. Its per-seam generator retains the local rate when an unoccupied remote component is added.', 'domain':'Complete finite reachable one-load modular sectors and supplied coarse grouping, transit policy and rates. Global monitoring changes the measured triangle evolution. No integer-regulator equality, many-load independence theorem, native-clock derivation or physical row award follows.'},
        {'id':'certified_canonical_native_window','statement':'A sparse local-coarse modular instrument runs on the canonical ring with explicit rejected boundary events. Exact rational Taylor arithmetic and an analytic projection bound enclose every position and classical-quantum outcome probability relative to the complete modular sector.', 'domain':'One coherent interval and Bernoulli event from the specified initial zero-record basis state, at the same modulus. The executed 55-state radius-4 ring window has combined probability-effect error below .006. The polynomial norm defect is retained. No repeated-channel, continuous per-seam scheduling, infinite-integer, or physical row transfer follows.'},
        {'id':'uniform_fiber_native_reduction','statement':'On a connected modular Gauss sector, normalized uniform record fibers define an isometry intertwining the Hamiltonian and every local-coarse C/E Kraus operator. The reduced canonical ring therefore reproduces complete native outcomes, repeated channels, external references and the per-seam generator on this preparation subspace.', 'domain':'Explicit uniform-fiber preparation, including composite modulus6; connected graph and fixed source sector. Native sector-preserving operations do not deterministically prepare it from definite records; conditional success is bounded by the initial zero-holonomy weight. Additional preparation resources are allowed but must be declared. Arbitrary terminal effects compress to V*EV; continued reduced evolution needs instrument-level preservation. Individual record measurements generally leave the code. Numerical continuous propagation is distinct from the exact operator identity. No physical charge identification follows.'},
        {'id':'compact_native_persistence','statement':'A closed three-seam transit cycle in the reduced canonical graph lifts through the uniform-fiber isometry to an eigenstate with energy -m. Every native local endpoint event preserves its density, including supplied per-seam continuous scheduling. Canonical1x1 and2x2 controls and an explicit complete modular triangle verify persistence and reject broken-cycle, definite-record and full-monitor substitutions.', 'domain':'Spatially compact, flat nonpropagating mode with globally coherent modular record fibers and supplied transit-as-elsewhere rule. The -m eigenvalue supplies no vacuum/positive-energy electron dictionary, physical rest clock, magnetic moment or full particle-stability criterion.'}]+c8_finite_claims()+c9_finite_claims()+c10_finite_claims()+c11_finite_claims(),
      'row_assessment':[{'id':r['id'],'local_status':r['local_status'],'full_criterion_met':r['full_criterion_met'],'official_status':r['official_status'],'missing_prerequisites':r['missing_prerequisites']} for r in rows],
      'physical_premises_closed':[], 'official_OpenWave_awards_added':0,
      'comparison_domains':{
        'shared_register_native_programs':'The nine claims in claims/joint_observer.json retain their authenticated finite source domains. A shared event register does not implement independent additive loads. No transfer to the separate-copy generator is asserted.',
        'clock':'The port/transit probability oscillates at 2E; no physical electron rest clock or derived commit-rate calibration is identified.',
        'regulator':'Two finite cutoff comparisons provide no cutoff-removal theorem.',
        'records':'C and E outcome payloads are retained; global position and register monitoring is supplied. Native coarse local measurements are a different instrument. Separate seam-offset dynamics are not implemented.',
        'interaction':'The canonical all-number source selects periodic half-record energy and a cyclic face term. Earlier E2a/E2b, first-order and zero-field variants remain controls; their spectra and squared-generator identities do not transfer.',
        'local_native_candidate':'The same configuration-level native CAR instrument acts on the full interacting carrier. Its coherent ground is generally excited by native monitoring. Earlier uniform-fiber reductions require their invariant zero-field domain; complete-sector enumeration grows with graph cycle rank.',
        'continuum':'No physical particle dictionary, common relativistic continuum or complete 31-row validation follows.'}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args(); text=json.dumps(generate(),indent=2,sort_keys=True)+'\n'
    if args.write: OUTPUT.write_text(text)
    elif OUTPUT.read_text()!=text: raise SystemExit('stale finite joint claim ledger')
    print('PASS_SOURCE_BOUND_FINITE_JOINT_MODEL_CLAIMS')

if __name__=='__main__':main()
