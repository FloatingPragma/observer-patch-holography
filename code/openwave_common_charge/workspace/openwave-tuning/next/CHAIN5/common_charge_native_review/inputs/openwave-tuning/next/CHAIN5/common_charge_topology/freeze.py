"""Freeze and verify the additive common-source charge connection."""
from bridge import *

def main():
    if '--verify' in sys.argv:
        manifest=read(HERE/'MANIFEST.json')
        for section in ('files','source_pins'):
            for path,digest in manifest[section].items():assert sha(ROOT/path)==digest,path
        print('FROZEN_COMMON_SOURCE_LOCAL_CHARGE_PASS',sha(HERE/'MANIFEST.json'));return
    pins=verify();embedding=read(HERE/'EMBEDDING.json');topology=read(HERE/'TOPOLOGY.json')
    dynamics=read(HERE/'DYNAMICS.json');meter=read(HERE/'NATIVE_RESULT.json')
    assert embedding['all_actual_edges_embedded']==288 and embedding['all_vertices']==96
    assert all(topology['checks'].values()) and dynamics['all_checks_pass'] and meter['all_checks_pass']
    assert topology['half_charge_counterexample']['declared_background_preserved']
    assert meter['native_2pi_integer_vs_half_control']['half_mutant_wrong_bit_probability']>1-2e-12
    assert not meter['source_prefix_native_whole_iterator_executed']
    # Import the executed entry points without running them, to freeze their
    # actual local module closure as well as the transitive source manifests.
    import embedding as embedding_code,topology as topology_code,dynamics as dynamics_code,native_flux
    import native_tomography
    runtime={}
    for module in tuple(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path:
            p=Path(path).resolve()
            if p.is_file() and p.is_relative_to(ROOT) and not p.is_relative_to(HERE):
                name=str(p.relative_to(ROOT));runtime[name]=sha(p)
                assert name not in pins or pins[name]==runtime[name]
                pins[name]=runtime[name]
    save('RUNTIME_IMPORTS.json',runtime)
    save('AUDIT_REQUEST.json',{'candidate_source':'native-local12-dirac-noncompact-quadratic-v4-proposed',
        'official_row':1,'required_producer_independent_full_audits':2,'producer_self_award':False,
        'criterion_sha256':'3af50af0955c64a0aa2019fa84af3cf0dd36e6c63a5abc76eb2c60c9b413cb6f',
        'review_questions':['Does the explicit actual-graph embedding and electric-curvature clutching meet the literal local topological-charge criterion?',
          'Does exact U Gauss law transfer through every g=1 signed route and leave no rounded-Q field substitution?',
          'Does the regional Fourier identity justify the native occupation instrument on the same source, state, reference and clock?',
          'Are every first export, pointer branch, reset timeout and full-source comparison bound retained?',
          'Are the global harmonic obstruction, one-particle domain and unexecuted huge native source iterator stated accurately?',
          'Do the fixed-background half-daughter, noninteger Q, wrong-reference and 2pi phase controls reject the relevant false claims?']})
    save('RESULT.json',{'schema':'oph.common_source.local_charge_connection.v1',
        'status':'PASS_CONSTRUCTIVE_LOCAL_CHARGE_CONNECTION_REQUIRES_INDEPENDENT_REVIEW',
        'official_validation_awarded':False,'all_packet_checks_pass':True,
        'candidate_source':'native-local12-dirac-noncompact-quadratic-v4-proposed','g':1,
        'source_manifest_sha256':'719036fe46769e68202731f765045dc3d75ebc8fe6e901e85319ce435776779f',
        'criterion':read(ROW1/'CRITERION.json'),
        'scope':'Local bounding-surface Chern charge on actual size2 common electric Gauss fiber; physical Q=n−G=−c1.',
        'same_source_preparation_clock':True,'extra_retained_history_register':False,'extra_dynamical_field':False,
        'global_torus_line_bundle_claim':False,'all_768_local_surface_checks':True,
        'signed_native_route_checks':dynamics['signed_native_routes'],
        'maximum_full_field_shift_diagnostic_residual':dynamics['maximum_full_field_current_shift_residual'],
        'native_meter_experiments':meter['complete_meter_experiments'],
        'native_fine_probabilities':meter['fine_probabilities_evaluated'],
        'maximum_complete_native_Gram_residual':meter['maximum_complete_Gram_residual'],
        'complete_source_meter_record_TV_upper':meter['maximum_complete_common_source_TV_bound'],
        'all_small_appended_native_columns_executed':True,'complete_large_native_source_iterator_executed':False,
        'review_gate':'Two independent complete criterion reviews; no automatic row award from producer checks.',
        'limitations':['One-particle native source; no two-particle creation implementation claimed.',
          'Local bundles only; continuous harmonic spectral support forbids a generic global compact bundle.',
          'Declared integer occupation and curvature normalization are architecture; uniqueness is not proved.',
          'Full source native word is constructive and bounded, not literally executed.',
          'No interacting spatial continuum, empirical coupling or charge-mass theorem.']})
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted(HERE.rglob('*'))
           if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts}
    manifest={'schema':'oph.frozen.common_source_charge.v1','files':files,'source_pins':dict(sorted(pins.items())),
              'local_file_count':len(files),'source_pin_count':len(pins),'runtime_import_count':len(runtime),
              'source_manifest_sha256':'719036fe46769e68202731f765045dc3d75ebc8fe6e901e85319ce435776779f',
              'official_validation_awarded':False}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('MANIFEST',sha(HERE/'MANIFEST.json'));print('PROOF',sha(HERE/'PROOF.md'));print('RESULT',sha(HERE/'RESULT.json'))
if __name__=='__main__':main()
