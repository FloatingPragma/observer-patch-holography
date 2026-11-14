"""Freeze this independent review after every required receipt passes."""
from review import HERE,ROOT,ROW,SIZE,EXPECTED,load,sha,pincheck
import json

def main():
    pins=pincheck()
    assert load(HERE/'CHECKS.json')['all_checks_pass']
    assert load(HERE/'SIZE3_CHECKS.json')['all_checks_pass']
    assert load(HERE/'GEOMETRY.json')['full_size3_integer_source_join']
    assert load(HERE/'GEOMETRY.json')['exact_all_direction_curl']
    assert load(HERE/'NATIVE.json')['all_settings']==32
    assert load(HERE/'SIZE3_NATIVE.json')['all_phase_targets']==1024
    assert load(HERE/'FRAME.json')['modes']==646
    assert len(load(HERE/'BOUNDS.json')['refinements'])==3
    assert load(HERE/'SOURCE.json')['symbolic_general_Hermitian_cap_Heisenberg_identity']
    for name,digest in load(HERE/'IMPORTS.json')['files'].items():
        assert sha(ROOT/name)==digest,name
        assert name not in pins or pins[name]==digest
        pins[name]=digest
    criterion=load(ROW/'CRITERION.json')
    assert criterion['literal_criterion']=='Maxwell recovered; transverse waves at c'
    upstream=ROOT/'openwave-tuning/next/CHAIN3/klein_gordon_completion/upstream/MODELS.md'
    assert sha(upstream)=='a69895e4af4f6bc180748c4e9a8ae11d745c75462650aea44c3a34bd687a139d'
    pins[str(upstream.relative_to(ROOT))]=sha(upstream)
    receipt_names=['GEOMETRY.json','NATIVE.json','BOUNDS.json','SOURCE.json','FRAME.json','SIZE3_NATIVE.json','CHECKS.json','SIZE3_CHECKS.json','IMPORTS.json']
    review={'schema':'oph.openwave.independent_full_criterion_review.v1','row':29,
      'reviewer':'/root/kg_complete','producer_independent':True,
      'candidate_id':criterion['candidate_id'],'criterion':criterion['literal_criterion'],
      'criterion_sha256':sha(ROW/'CRITERION.json'),'producer_manifest_sha256':EXPECTED[str((ROW/'MANIFEST.json').relative_to(ROOT))],
      'size3_manifest_sha256':EXPECTED[str((SIZE/'MANIFEST.json').relative_to(ROOT))],
      'full_criterion_met':True,'decision':'PASS: the literal complete Maxwell criterion is met on the declared free noncompact radiative construction, with finite native instrument errors.',
      'review_file':str((HERE/'REVIEW.md').relative_to(ROOT)),'review_sha256':sha(HERE/'REVIEW.md'),
      'scope':'Oriented vacuum Maxwell field dictionary and controlled radiative limit, two transverse polarizations, a=1/2 b=3 speed one in this source clock, spectral preparation, autonomous spectator reduction, complete terminal native field/clock instruments.',
      'source_clock_distinction':'This v4 field uses c=1 under its h-time convention. No equality with scalar Klein-Gordon c=1/sqrt48 dynamics, clock conversion or mass coin is inferred. Same routed Dirac/field-source KG transfer remains a separate connection.',
      'limitations':[
        'Coefficients selected to match admitted Dirac speed; no independent empirical c or alpha prediction.',
        'Free noncompact radiative source and first-moment POVMs; no compact-cosine, interacting Weyl, full continuous-field spectrum, or local old-edge encoder transfer.',
        'Huge L values are finite theorem bounds; only their reduced terminal instruments were literally executed.',
        'Native fixed precision leaves an explicit nonzero floor.',
        'Component arrays are before reset success weights; complete instruments compose actual rank-one resets and retain first timeouts.',
        'One source-specific full review does not establish inter-row unification or an official OpenWave award.'
      ],
      'resolved_concerns':[
        'Clock certificate includes its pointer; two source pointer terms yield three in total.',
        'Visible vacuum scalar retained; spectator scalars cancel only under proved autonomous reduction.',
        'All 646 fixed rational frame columns retained; no runtime eigensolver basis reselection.',
        'Initial reviewer return-arity error corrected; original attempt log retained and corrected full command passed.'
      ],
      'evidence_files':{str((HERE/name).relative_to(ROOT)):sha(HERE/name) for name in receipt_names},
      'official_rows_changed':[]}
    (HERE/'REVIEW.json').write_text(json.dumps(review,indent=2)+'\n')
    manifest={'schema':'oph.openwave.row29.independent_full_review.freeze.v1',
      'row':29,'reviewer':'/root/kg_complete','producer_independent':True,'full_criterion_met':True,
      'producer_manifest_sha256':review['producer_manifest_sha256'],'size3_manifest_sha256':review['size3_manifest_sha256'],
      'files':{str(p.relative_to(ROOT)):sha(p) for p in sorted(HERE.rglob('*')) if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts},
      'source_pins':pins,'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for group in ('files','source_pins'):
        for name,digest in manifest[group].items():assert sha(ROOT/name)==digest,name
    print(json.dumps({'full_criterion_met':True,'manifest':sha(HERE/'MANIFEST.json'),
                      'review':sha(HERE/'REVIEW.json'),'proof':sha(HERE/'REVIEW.md'),
                      'source_pins':len(pins),'review_files':len(manifest['files'])}),flush=True)
if __name__=='__main__':main()
