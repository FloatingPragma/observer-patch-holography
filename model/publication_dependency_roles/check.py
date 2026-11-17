"""Authenticate the exact current presentation edges and their scientific inputs."""
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
                if current.name == '_bindings':
                    assert isinstance(parents[node], ast.List), 'note contents used outside custody list'
                elif current.name == 'verify':
                    statement=node
                    while not isinstance(statement,ast.Expr):statement=parents[statement]
                    expected=ast.parse('required.add(NOTE.relative_to(ROOT).as_posix())').body[0]
                    assert ast.dump(statement)==ast.dump(expected), 'note contents used outside custody membership'
                else:
                    raise AssertionError('note has an unclassified consumer')
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
            assert len(branch.generators)==1 and not branch.generators[0].ifs
            generator=branch.generators[0]
            assert isinstance(generator.target,ast.Name)
            variable=generator.target.id
            expected=ast.parse(
                '{'+variable+': hashlib.sha256((ROOT/'+variable+').read_bytes()).hexdigest() '
                'for '+variable+' in SOURCE_PATHS}',mode='eval').body
            assert ast.dump(branch)==ast.dump(expected), 'source list used outside exact custody hashing'
            consumers.append(current.name)
    assert len(consumers)==1
    return {'literal_assignment':'SOURCE_PATHS','only_consumers':consumers,
        'read_domain':'Returned SHA256 custody map only; downstream record_clock merges provider maps without reading presentation contents.'}


def verify(registry=None, *, override=None):
    if registry is None:
        registry=json.loads((HERE/'roles.json').read_text())
    assert registry['schema']=='oph.exact_presentation_dependency_roles.v2'
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
        assert row['binding_map']==('source_bindings' if 'chiral_candidate' in row['owner'] else 'source_sha256')
        assert row['role']=='presentation_provenance' and row['scientific_input'] is False
        assert hash_bytes(read(row['owner']))==row['owner_sha256'],'owner binding changed'
        assert hash_bytes(read(row['provider']))==row['provider_sha256'],'operative provider changed'
        assert hash_bytes(read(row['presentation_path']))==row['presentation_sha256'],'presentation bytes changed'
        owner=json.loads(read(row['owner']))
        assert owner[row['binding_map']][row['presentation_path']]==row['presentation_sha256']
        # Every operative pin is checked; there is no document wildcard.
        for name,digest in owner[row['binding_map']].items():
            assert hash_bytes(read(name))==digest,('operative custody changed',name)
        use=traced_uses(read(row['provider']).decode(),row['presentation_path'],row['provider'])
        reports.append({**row,'read_trace':use,'current_bytes_authenticated':True})
    return {'schema':'oph.presentation_dependency_role_review.v2',
        'verdict':'PASS_THREE_EXACT_CURRENT_PRESENTATION_ROLES',
        'edges':reports,'exact_presentation_edge_count':3,
        'typed_scientific_evidence_uses_document_prose':False,
        'scientific_outputs_recomputed':False,
        'validation_mode':'Exact current bindings and AST role checks; no physical replay.',
        'catch_all_docs_exception':False,
        'scope':'Only the three enumerated presentation roles. Every operative input is authenticated; this check awards no physical claim.'}
