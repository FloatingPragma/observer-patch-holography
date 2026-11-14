from pathlib import Path
from fractions import Fraction as F
import hashlib,json

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
read=lambda p:json.loads(p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    assert not (HERE/'MANIFEST.json').exists()
    native=read(HERE/'NATIVE_AUDIT.json');channel=read(HERE/'CHANNEL_AUDIT.json')
    assert len(native['experiments'])==36 and channel['all_checks_pass']
    assert native['raw_M6_pulses_replayed_including_pointer']==8000
    assert native['complete_fine_probabilities_recomputed']==54756
    assert F(native['complete_terminal_TV_upper_exact'])<F(32,10**16)
    pins={}
    for file in ('SOURCE_PINS.json','CHANNEL_SOURCE_PINS.json'):
        for name,digest in read(HERE/file).items():
            assert sha(ROOT/name)==digest,name
            assert name not in pins or pins[name]==digest,name
            pins[name]=digest
    review={'schema':'oph.charge_terminal_annex.independent_review.v1',
        'verdict':'PASS_SCOPED_TERMINAL_CHANNEL_AND_COMPLETE_NATIVE_INSTRUMENT',
        'producer_manifest':'openwave-tuning/next/CHAIN5/common_charge_native_review/MANIFEST.json',
        'producer_manifest_sha256':'7104fb54312af99ba577ed6e2a69c3bc6390954f2e4e329db3627285750a1a2b',
        'reviewed_parts':['PROOF.md::New annex: exact terminal-observable execution','TERMINAL_ANNEX_PROGRAMS.json','TERMINAL_ANNEX_FAILURES.json','TERMINAL_ANNEX_RECORDS.json','TERMINAL_ANNEX_RESULT.json'],
        'producer_independent':True,'producer_or_prior_reviewer_executor_imported':False,
        'first_independent_annex_review_read_before_decision':False,
        'exact_source_block_permutation_supported':True,
        'arbitrary_retained_reference_channel_supported':True,
        'all36_native_programs_complete_all_column_replay':True,
        'all_first_export_effects_and_fine_pointer_records_recomputed':True,
        'all78_reset_probabilities_and_fine_tapes_rechecked':True,
        'raw_native_pulses':native['raw_M6_pulses_replayed_including_pointer'],
        'fine_probabilities':native['complete_fine_probabilities_recomputed'],
        'complete_terminal_TV_upper_exact':native['complete_terminal_TV_upper_exact'],
        'complete_terminal_TV_display':native['complete_terminal_TV_display'],
        'output_domain':'L2 terminal regional charge/electric flux and source clock fine records with arbitrary quantum reference. Cell/internal outputs and private successful reset tapes are marginalized. Every first timeout/export remains a failure outcome.',
        'not_approved':['Equality of full source unitaries or retained field states','Equality of private native tapes','Intermediate measurements','L>2 deterministic cell reduction','Execution of original huge full-field prefix','A new full-row award'],
        'native_input_and_reset_scope':'All dirty logical input columns are checked for the word. Ideal blank preparation compares after marginalizing private successful reset records; first failures are never success-normalized.',
        'initial_failed_runs':'Both preserved PATH_FAILURE logs concern nested snapshot resolution only; no mathematical gate or tolerance changed.',
        'full_row_award':False,'official_rows_changed':[]}
    # Bind the actual immutable digest; the literal value is checked here.
    digest=sha(ROOT/review['producer_manifest'])
    assert digest=='7104fb54312af99ba577ed6e2a69c3bc6390954f2e4e329db3627285750a1a2b'
    review['producer_manifest_sha256']=digest
    (HERE/'REVIEW.json').write_text(json.dumps(review,indent=2)+'\n')
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted(HERE.iterdir()) if p.is_file()}
    result={'schema':'oph.charge_terminal_annex.review.freeze.v1','files':files,'source_pins':pins,
        'producer_manifest_sha256':digest,'verdict':review['verdict'],
        'producer_independent':True,'full_row_award':False,'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'manifest':sha(HERE/'MANIFEST.json'),'files':len(files),'pins':len(pins),
        'complete_terminal_TV':native['complete_terminal_TV_display']}),flush=True)


if __name__=='__main__':main()
