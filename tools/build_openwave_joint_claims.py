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
    required={'integration/assess.py','candidate/oph_openwave/construction.py','candidate/oph_openwave/quantum_records.py','candidate/oph_openwave/geometry.py','platform_packet/packet.json','candidate/oph_openwave/observer_protocol.py'}
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
      'construction':'Portable port/transit incidence Hamiltonian with supplied mass; bounded integer half-seam Gauss sector and monitored native event instrument.',
      'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
      'assessment':str(receipt.relative_to(ROOT)),
      'finite_claims':[
        {'id':'additive_independent_loads','statement':'Independent prepared loads evolve under separate copies of the port/transit generator.','domain':'Product preparations; no interacting or exchange-statistics theorem.'},
        {'id':'integer_joint_gauss','statement':'Half-seam hopping and accepted literal native currents preserve the finite Gauss sector.','domain':'One load; bounded integer records; hard-wall coherent regulator.'},
        {'id':'monitored_native_instrument','statement':'The declared measured event instrument preserves trace and positivity with complete native, transit-noop and rejected-attempt outcomes.','domain':'Position and half-seam records measured; supplied discrete time and commit probability.'},
        {'id':'observer_complete_readout','statement':'The finite observer readout retains position probabilities, monitored-attempt outcomes and applied-current probabilities; sparse evolution agrees with finite density controls.', 'domain':'One-load finite sector. The carrier path-window error bound concerns coherent evolution of the specified initial basis preparation against the complete sector at the same cutoff; no repeated-channel or cutoff-removal bound.'},
        {'id':'bare_dark_mode_nontransfer','statement':'The oriented triangle bare incidence has rank 2 and transit kernel dimension 1; the cutoff-1 integer-record incidence in the initial-port Gauss sector has rank 6 and transit kernel dimension 0. The declared equal-amplitude BFS transit preparation leaks to ports.', 'domain':'Exact finite triangle rank control and the explicitly supplied BFS preparation at time .37 and mass .23. No transfer of bare dark-mode stability or particle-stability identification.'}],
      'row_assessment':[{'id':r['id'],'local_status':r['local_status'],'full_criterion_met':r['full_criterion_met'],'official_status':r['official_status'],'missing_prerequisites':r['missing_prerequisites']} for r in rows],
      'physical_premises_closed':[], 'official_OpenWave_awards_added':0,
      'comparison_domains':{
        'shared_register_native_programs':'The nine claims in claims/joint_observer.json retain their authenticated finite source domains. A shared event register does not implement independent additive loads. No transfer to the separate-copy generator is asserted.',
        'clock':'The port/transit probability oscillates at 2E; no physical electron rest clock or derived commit-rate calibration is identified.',
        'regulator':'Two finite cutoff comparisons provide no cutoff-removal theorem.',
        'records':'C and E outcome payloads are retained; separate seam-offset dynamics are not implemented.',
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
