"""Refresh authenticated native projections while retaining exact legacy evidence.

This path is limited to the declared joint-observer extension. It does not
replace a scientific producer or the canonical full ledger builder. Every
pre-existing scientific/Lean/authority file must equal the supplied Git base.
Only explicitly reviewed presentation/projection files may differ. Both
affected ledger views are evaluated from their actual current producers.
"""
from copy import deepcopy
from contextlib import contextmanager
from functools import lru_cache
from hashlib import sha256,sha1
from pathlib import Path
import argparse
import importlib.util
import json
import os
import stat
import subprocess
from tools.joint_observer_results import NAMES

ROOT=Path(__file__).resolve().parents[1]
LEDGER='code/particles/runs/status/postdiction_ledger.json'
MARKDOWN='docs/POSTDICTION_LEDGER.md'
OUTPUTS={LEDGER,MARKDOWN}
PAPER='paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex'
PAPER_ROOTS=(
 'flagship/from_observer_consensus_to_standard_physics',
 'paper/deriving_the_particle_zoo_from_observer_consistency',
 'paper/executable_observer_construction',
 'paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency',
 'extra/yang_mills_gap_clay_problem',
 'cosmology/oph_cosmological_vacuum_and_structure_formation')
REVIEWED_EXISTING={
 'README.md','tools/native_physical_results.py','tools/test_native_physical_results.py',
 'tools/coherent_integration_results.py','tools/build_coherent_integration_claims.py',
 'claims/coherent_integration.json','model/publication_dependency_roles/paper_results.json',
 'paper/paper_release_manifest.json',
 *(p+ext for p in PAPER_ROOTS for ext in ('.tex','.pdf'))}
NEW_GLOBAL={
 'tools/joint_observer_results.py','tools/build_joint_observer_claims.py',
 'tools/test_joint_observer_results.py','tools/refresh_native_projection.py',
 'tools/test_refresh_native_projection.py','claims/joint_observer.json',
 'paper/tex_fragments/NATIVE_SHARED_SOURCE_FIELD_RECORDS.tex'}
PACKAGE_FILES={
 '__init__.py','__main__.py','source.py','README.md','PROOF.md','results.json','test_source.py'}
NEW_PACKAGE={f'model/{n}/{f}'for n in NAMES for f in PACKAGE_FILES}
NEW_PACKAGE|={
 'model/native_observer_program/contract.py','model/native_observer_program/reflection.py',
 'model/native_observer_program/test_program.py','model/native_observer_program/PAPER_RESULT.tex',
 'model/native_electric_gauss/PAPER_RESULT.tex','model/native_gauss_walk/PAPER_RESULT.tex'}
ENV=dict(os.environ,GIT_OPTIONAL_LOCKS='0')

def require(condition,message):
    if not condition:raise ValueError(message)

def git(root,*args):
    return subprocess.check_output(['git',*args],cwd=root,env=ENV,timeout=60)

def blob_oid(data):
    return sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()

def read_path(root,name,mode=None):
    p=Path(root)/name
    if p.is_symlink():
        require(mode=='120000','new symlink or regular-file type changed: '+name)
        return os.readlink(p).encode()
    require(p.is_file(),'missing retained source: '+name)
    require(mode!='120000','retained symlink replaced: '+name)
    actual='100755'if p.stat().st_mode&stat.S_IXUSR else'100644'
    require(mode is None or mode==actual,'retained executable mode changed: '+name)
    return p.read_bytes()

def base_tree(root,base):
    rows={}
    for record in git(root,'ls-tree','-rz','--full-tree',base).split(b'\0'):
        if not record:continue
        meta,path=record.split(b'\t',1);mode,kind,oid=meta.decode().split()
        require(kind=='blob','unsupported retained Git object: '+path.decode())
        rows[path.decode()]=(mode,oid)
    return rows

def current_names(root):
    return {p for p in git(root,'ls-files','--cached','--others','--exclude-standard','-z').decode().split('\0')if p}

def classify_changes(base_oids,current_oids):
    """Reject old operative drift, unknown additions and all deletions."""
    changes={}
    for name,old in base_oids.items():
        if name in OUTPUTS:continue
        require(name in current_oids,'deleted retained source: '+name)
        if current_oids[name]!=old:
            require(name in REVIEWED_EXISTING,'operative retained input changed: '+name)
            changes[name]=current_oids[name]
    for name,oid in current_oids.items():
        if name in base_oids or name in OUTPUTS:continue
        require(name in NEW_GLOBAL|NEW_PACKAGE,'unreviewed added input: '+name)
        changes[name]=oid
    return changes

def inventory(root):
    modes={p:mode for p,(mode,_)in base_tree(root,'HEAD').items()}
    return {p:sha256(read_path(root,p,modes.get(p))).hexdigest()
            for p in sorted(current_names(root)-OUTPUTS)}

def authenticate_boundary(root,base,review=None):
    require(git(root,'rev-parse','HEAD').decode().strip()==base,'checkout differs from authenticated base')
    tree=base_tree(root,base);names=current_names(root)
    data={p:read_path(root,p,tree.get(p,(None,None))[0])for p in sorted(names-OUTPUTS)}
    oids={p:blob_oid(b)for p,b in data.items()}
    changes=classify_changes({p:oid for p,(_,oid)in tree.items()},oids)
    hashes={p:sha256(data[p]).hexdigest()for p in sorted(changes)}
    report={'base_commit':base,'reviewed_changed_sha256':hashes,
            'retained_git_blob_count':len(tree),'all_other_base_blobs_exact':True,
            'canonical_full_builder_unchanged':True}
    if review is not None:require(report==review,'reviewed change manifest differs from actual source')
    return report

def retained_collector(old,fresh,receipt_sha):
    retained=deepcopy(fresh)
    require(set(retained['joint_observer_certificates']['source_packages'])==set(NAMES),
            'joint source package set differs')
    retained.pop('joint_observer_certificates')
    for n in NAMES:
        path=f'model/{n}/results.json'
        require(retained['receipts_sha256'].pop(path)==receipt_sha[path],
                'joint receipt binding differs: '+n)
    require(retained==old,'retained collector values differ')
    return True

def completion_projection(old,fresh,old_paper_sha,new_paper_sha):
    retained=deepcopy(fresh);renewed=[]
    for i,row in enumerate(retained['retained_foundations']):
        pins=row['source_sha256']
        if PAPER in pins:
            require(old['retained_foundations'][i]['source_sha256'].get(PAPER)==old_paper_sha,
                    'legacy paper provenance is not the authenticated base file')
            require(pins[PAPER]==new_paper_sha,'current paper provenance differs')
            pins[PAPER]=old_paper_sha;renewed.append(i)
    require(renewed==[4],'unexpected presentation-provenance scope')
    require(retained==old,'retained completion decision or scientific value differs')
    return {'path':PAPER,'old_sha256':old_paper_sha,'current_sha256':new_paper_sha,
            'foundation_indices':renewed}

def canonical_renderer(root):
    path=Path(root)/'code/particles/scripts/build_postdiction_ledger.py'
    spec=importlib.util.spec_from_file_location('native_projection_canonical_ledger',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module._render_md

@contextmanager
def source_ast_cache():
    """Cache pure parsing by complete text/arguments for this closure only.

    The original algorithm reads every source and resolves every import. Its
    AST traversal is read-only. Hashing and custody checks are unchanged, and
    edited source has a different cache key. The shared ast module is never
    modified; the original dependency is restored on every exit.
    """
    from model.particle_common_contract import provenance
    original=provenance.ast
    parse=lru_cache(maxsize=None)(original.parse)
    class ParserProxy:
        def __getattr__(self,name):
            return parse if name=='parse'else getattr(original,name)
    stats={}
    provenance.ast=ParserProxy()
    try:yield stats
    finally:
        provenance.ast=original
        info=parse.cache_info()
        stats.update(hits=info.hits,misses=info.misses,entries=info.currsize,
                     key='complete source text and parser arguments',
                     original_dependency_restored=provenance.ast is original,
                     global_ast_module_unchanged=True)
        parse.cache_clear()

def refresh(root,base,review,write=False):
    root=Path(root).resolve()
    require(root==ROOT.resolve(),'refresh root differs from the actual producer module root')
    boundary=authenticate_boundary(root,base,review)
    before=inventory(root)
    original=git(root,'show',base+':'+LEDGER);old=json.loads(original)
    if write:require((root/LEDGER).read_bytes()==original,'writer requires exact public ledger baseline')
    from tools.native_physical_results import collect
    from tools.build_completion_program_current import ledger_view
    from tools.build_joint_observer_claims import generate
    with source_ast_cache() as cache_stats:
        expected_claim=json.dumps(generate(),sort_keys=True,indent=2,ensure_ascii=False)+'\n'
        require((root/'claims/joint_observer.json').read_text()==expected_claim,'joint claim register differs')
        fresh=collect(root)
    new_receipts={f'model/{n}/results.json':sha256((root/f'model/{n}/results.json').read_bytes()).hexdigest()for n in NAMES}
    retained_collector(old['native_spatial_and_field_records'],fresh,new_receipts)
    completion=ledger_view(root)
    paper_renewal=completion_projection(old['executable_completion_program'],completion,
        sha256(git(root,'show',base+':'+PAPER)).hexdigest(),sha256((root/PAPER).read_bytes()).hexdigest())
    result=deepcopy(old);result['native_spatial_and_field_records']=fresh
    result['executable_completion_program']=completion
    text=json.dumps(result,sort_keys=True,indent=2)+'\n'
    markdown=canonical_renderer(root)(result)
    require(inventory(root)==before,'source changed during projection refresh')
    if write:
        (root/LEDGER).write_text(text);(root/MARKDOWN).write_text(markdown)
    else:
        require((root/LEDGER).read_text()==text,'ledger JSON differs from authenticated projection')
        require((root/MARKDOWN).read_text()==markdown,'canonical Markdown differs')
    require((root/LEDGER).read_text()==text and (root/MARKDOWN).read_text()==markdown,
            'stored output differs from evaluated canonical bytes')
    require(inventory(root)==before,'source changed during projection installation')
    return {'verdict':'PASS_AUTHENTICATED_NATIVE_PROJECTION_REFRESH',
        'mode':'write'if write else'check','base_commit':base,
        'retained_source_boundary':boundary,
        'input_inventory_sha256':sha256(json.dumps(before,sort_keys=True).encode()).hexdigest(),
        'input_count':len(before),'new_receipt_sha256':new_receipts,
        'all_retained_scientific_values_exact':True,
        'exact_same_file_presentation_renewal':paper_renewal,
        'collector_receipt_count':len(fresh['receipts_sha256']),
        'pure_source_parser_cache':cache_stats,
        'output_sha256':{p:sha256((root/p).read_bytes()).hexdigest()for p in sorted(OUTPUTS)},
        'original_full_ledger_builder_reexecuted':False,
        'original_scientific_owners_reexecuted':False,'CI_run':False,
        'acceptance_scope':'Actual current native collector and completion view; unchanged canonical renderer; exact authenticated legacy source and value reuse.'}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',required=True);p.add_argument('--review',type=Path,required=True)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--capture-review',action='store_true');g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    p.add_argument('--report',type=Path)
    a=p.parse_args()
    if a.capture_review:
        require(not a.review.exists(),'preserve existing review manifest')
        report=authenticate_boundary(ROOT,a.base)
        a.review.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    else:
        report=refresh(ROOT,a.base,json.loads(a.review.read_text()),write=a.write)
        if a.report:a.report.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps(report,sort_keys=True,indent=2))
if __name__=='__main__':main()
