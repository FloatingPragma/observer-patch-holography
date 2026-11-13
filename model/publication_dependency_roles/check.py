"""Authenticate three immutable presentation edges without changing their receipts."""
import ast
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
EXPECTED={
    ('model/local_repair_family/record_clock_results.json','docs/INTERPRETATION_AND_TERMS.md'),
    ('model/local_repair_family/record_clock_results.json','docs/research/COHERENT_SOURCE_CLOCKS.md'),
    ('model/common_source/chiral_candidate/results.json','docs/research/OPH_PARTICLE_MASS_RESPONSE.md')}
PROVIDERS={
    ('model/local_repair_family/record_clock_results.json','docs/INTERPRETATION_AND_TERMS.md'):'model/native_modular_reference/construction.py',
    ('model/local_repair_family/record_clock_results.json','docs/research/COHERENT_SOURCE_CLOCKS.md'):'model/compiled_native_interface/construction.py',
    ('model/common_source/chiral_candidate/results.json','docs/research/OPH_PARTICLE_MASS_RESPONSE.md'):'model/common_source/chiral_candidate/run.py'}


def hash_bytes(value):return hashlib.sha256(value).hexdigest()


def traced_uses(source, path, provider):
    tree=ast.parse(source)
    parents={child:node for node in ast.walk(tree) for child in ast.iter_child_nodes(node)}
    locations=[node for node in ast.walk(tree) if isinstance(node,ast.Constant) and node.value==path]
    assert len(locations)==1,'presentation path has an unclassified source use'
    current=locations[0]
    while not isinstance(current,ast.Assign):current=parents[current]
    target=current.targets[0]
    assert isinstance(target,ast.Name)
    if provider.endswith('chiral_candidate/run.py'):
        assert target.id=='NOTE'
        uses=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Name) and node.id=='NOTE' and isinstance(node.ctx,ast.Load):
                current=node
                while not isinstance(current,ast.FunctionDef):current=parents[current]
                uses.append(current.name)
        assert sorted(uses)==['_bindings','verify']
        return {'literal_assignment':'NOTE','only_consumers':uses,
            'read_domain':'_bindings hashes bytes; verify checks required custody membership. run obtains physical outputs independently from action, geometry, history and hierarchy controls.'}
    assert target.id=='SOURCE_PATHS'
    consumers=[]
    for node in ast.walk(tree):
        if isinstance(node,ast.Name) and node.id=='SOURCE_PATHS' and isinstance(node.ctx,ast.Load):
            current=node
            while not isinstance(current,(ast.FunctionDef,ast.Module)):current=parents[current]
            assert isinstance(current,ast.FunctionDef)
            # Every direct use is the returned source_sha256 comprehension.
            branch=node
            while not isinstance(branch,ast.DictComp):branch=parents[branch]
            assert isinstance(branch.value,ast.Call)
            assert ast.unparse(branch.value).endswith('.hexdigest()')
            consumers.append(current.name)
    assert len(consumers)==1
    return {'literal_assignment':'SOURCE_PATHS','only_consumers':consumers,
        'read_domain':'Returned SHA256 custody map only; downstream record_clock merges provider maps without reading presentation contents.'}


def verify(registry=None, *, override=None):
    registry=registry or json.loads((HERE/'roles.json').read_text())
    override=override or {}
    def read(name):
        p=(ROOT/name).resolve()
        assert p.is_relative_to(ROOT),'escaped source path'
        return override.get(name,p.read_bytes())
    edges=registry['edges']
    assert {(x['owner'],x['presentation_path']) for x in edges}==EXPECTED
    assert len(edges)==3
    reports=[]
    for row in edges:
        assert row['provider']==PROVIDERS[row['owner'],row['presentation_path']]
        assert row['role']=='retained_presentation_provenance' and row['scientific_input'] is False
        assert hash_bytes(read(row['owner']))==row['owner_sha256'],'immutable owner changed'
        assert hash_bytes(read(row['provider']))==row['provider_sha256'],'operative provider changed'
        assert hash_bytes(read(row['presentation_path']))==row['presentation_sha256'],'legacy presentation bytes changed'
        owner=json.loads(read(row['owner']))
        assert owner[row['legacy_map']][row['presentation_path']]==row['presentation_sha256']
        # Every original science pin remains checked, without a docs wildcard.
        for name,digest in owner[row['legacy_map']].items():
            assert hash_bytes(read(name))==digest,('original custody changed',name)
        use=traced_uses(read(row['provider']).decode(),row['presentation_path'],row['provider'])
        reports.append({**row,'read_trace':use,'legacy_bytes_authenticated':True})
    return {'schema':'oph.presentation_dependency_role_review.v1',
        'verdict':'PASS_THREE_EXACT_PRESENTATION_ROLES_WITH_IMMUTABLE_ANCESTRY',
        'edges':reports,'raw_legacy_documentary_policy_findings_retained':3,
        'raw_legacy_audit_exit_status_is_not_changed':True,
        'typed_scientific_evidence_uses_document_prose':False,
        'scientific_outputs_or_immutable_selection_goldens_changed':False,
        'catch_all_docs_exception':False,
        'scope':'Exact inherited edge roles only. The original raw no-docs audit retains its three findings; this checker authenticates all original bytes and does not award a physical claim.'}
