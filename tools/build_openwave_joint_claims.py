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
    return {
      'schema':'oph.openwave_joint_model_claims.v1',
      'construction':'Active canonical port/transit modular construction with supplied mass, local-coarse C/E instrument and per-seam rate; definite and uniform-fiber preparations are explicit. Bounded integer fully monitored and record-energy variants retain separate comparison domains.',
      'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
      'assessment':str(receipt.relative_to(ROOT)),
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
        {'id':'uniform_fiber_native_reduction','statement':'On a connected modular Gauss sector, normalized uniform record fibers define an isometry intertwining the Hamiltonian and every local-coarse C/E Kraus operator. The reduced canonical ring therefore reproduces complete native outcomes, repeated channels, external references and the per-seam generator on this preparation subspace.', 'domain':'Explicit uniform-fiber preparation, including composite modulus6; connected graph and fixed source sector. Individual record-resolving measurements generally leave the subspace. No definite-record substitution, efficient preparation of the enormous joint state, or physical charge identification is asserted. Numerical continuous propagation is distinct from the exact operator identity.'},
        {'id':'compact_native_persistence','statement':'A closed three-seam transit cycle in the reduced canonical graph lifts through the uniform-fiber isometry to an eigenstate with energy -m. Every native local endpoint event preserves its density, including supplied per-seam continuous scheduling. Canonical1x1 and2x2 controls and an explicit complete modular triangle verify persistence and reject broken-cycle, definite-record and full-monitor substitutions.', 'domain':'Spatially compact, flat nonpropagating mode with globally coherent modular record fibers and supplied transit-as-elsewhere rule. The -m eigenvalue supplies no vacuum/positive-energy electron dictionary, physical rest clock, magnetic moment or full particle-stability criterion.'}],
      'row_assessment':[{'id':r['id'],'local_status':r['local_status'],'full_criterion_met':r['full_criterion_met'],'official_status':r['official_status'],'missing_prerequisites':r['missing_prerequisites']} for r in rows],
      'physical_premises_closed':[], 'official_OpenWave_awards_added':0,
      'comparison_domains':{
        'shared_register_native_programs':'The nine claims in claims/joint_observer.json retain their authenticated finite source domains. A shared event register does not implement independent additive loads. No transfer to the separate-copy generator is asserted.',
        'clock':'The port/transit probability oscillates at 2E; no physical electron rest clock or derived commit-rate calibration is identified.',
        'regulator':'Two finite cutoff comparisons provide no cutoff-removal theorem.',
        'records':'C and E outcome payloads are retained; global position and register monitoring is supplied. Native coarse local measurements are a different instrument. Separate seam-offset dynamics are not implemented.',
        'interaction':'No record-energy choice is selected by default. Explicit interacting variants do not inherit the free squared-generator identity.',
        'local_native_candidate':'The modular local-coarse instrument is the active proposed construction. Complete-sector enumeration grows with cycle rank; the definite-record sparse certificate covers one interval, while uniform-fiber reduction gives an exact different-preparation simulation of the whole instrument. Neither supplies the physical species dictionary.',
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
