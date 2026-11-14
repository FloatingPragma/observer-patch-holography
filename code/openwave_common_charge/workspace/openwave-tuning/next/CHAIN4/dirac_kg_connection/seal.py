"""Pin the constructive free-source connection without a row self-award."""
from connection import *
import sympy,scipy

def main():
    pins=verify()
    for name in ('GEOMETRY_CHECKS.json','NATIVE_RESULT.json','SPATIAL_CHECKS.json','CONTROLS.json'):
        assert read(HERE/name)['all_checks_pass'],name
    assert read(HERE/'FREE_FACTOR.json')['g1_original_word_delegation']
    assert len(read(HERE/'DISPERSION.json')['controlled_cases'])==120
    assert read(HERE/'SPATIAL_NATIVE.json')['all_spin_flavor_columns']==16
    assert read(HERE/'NATIVE_RESULT.json')['fine_probabilities_evaluated']==908544
    assert len(read(HERE/'CONTROLS.json')['finite_size2_small_coupling_connection'])==9
    # Import the actual pointer and static-field consumer closure; no external
    # service or producer mutation participates in this freeze.
    import native_tomography,spatial_check,geometry_checks,controls,transfer
    native_tomography.pointer()
    runtime={}
    for module in tuple(sys.modules.values()):
        filename=getattr(module,'__file__',None)
        if not filename:continue
        p=Path(filename).resolve()
        if p.is_file() and p.is_relative_to(ROOT) and not p.is_relative_to(HERE):
            rel=str(p.relative_to(ROOT));digest=sha(p)
            assert rel not in pins or pins[rel]==digest,rel
            pins[rel]=digest;runtime[rel]=digest
    for p in (ROOT/'reverse-engineering-reality/AGENTS.md',ROOT/'reverse-engineering-reality/docs/policies/STYLE_GUIDE.md'):
        pins[str(p.relative_to(ROOT))]=sha(p)
    (HERE/'RUNTIME_IMPORTS.json').write_text(json.dumps({'files':runtime,'python':sys.version,'numpy':np.__version__,
        'scipy':scipy.__version__,'sympy':sympy.__version__,'numerical_threads':1},indent=2)+'\n')
    native_result=read(HERE/'NATIVE_RESULT.json');refinement=read(HERE/'REFINEMENT_CONTRACT.json')
    result={'schema':'oph.dirac_field_same_source_kg_connection.v1',
      'status':'PASS_SCOPED_SAME_ROUTED_SOURCE_KG_FREE_INTERVENTION',
      'candidate_id':'native-local12-dirac-noncompact-flavor-coupling-v1',
      'source_parameter':'g scales every actual signed route current and physical longitudinal charge; g^2 scales Coulomb energy.',
      'source_equalities':{'size2_g1':'Literal frozen field25 iterator, with task29 flavor mass-slot substitution.',
          'g0':'Exact free routed Dirac factor tensor the same autonomous noncompact field and native clock.',
          'larger_spatial_family':'Same LOCAL12 root-port path rule; dyadic native signed endpoint compiler. Coupled larger-L native field implementation not claimed.'},
      'field_clock_speed_join':{'lattice_field_time_per_macro':1,'physical_time_per_macro':'h','native_clock_phase_per_macro':'h/8',
          'field_physical_generator':'H_lattice/h','Maxwell_speed_squared':'a b rho = (1/2)*3*(2/3)=1','Dirac_speed_squared':'1'},
      'mass_coefficients_exact':list(map(str,MASSES)),'mass_ancestry':'Executed finite register/analyzer effect with retrospective architecture provenance; exact frozen rational midpoints, no GeV calibration in this connection.',
      'retained_degrees':{'spin':4,'physical_flavors':3,'flavor_padding':1,'source_size2_coexact_modes':190,'source_size2_harmonic_modes':3,'Brillouin_corners':8},
      'executed':{'route_sizes':[2,3,4],'native_spatial_member_L':4,'native_spatial_source_pulses':read(HERE/'SPATIAL_NATIVE.json')['source_native_pulses'],
          'native_spatial_spin_flavor_columns':16,'complete_native_source_cases':4,'native_controlled_columns_per_case':32,
          'spin_preparations':28,'complete_tomography_component_records':5376,'evaluated_fine_probabilities':908544,'controlled_dispersion_cases':120},
      'native_maximum_tomography_residual':native_result['maximum_tomographic_operator_error'],
      'finite_unexecuted_refinement':{'L':refinement['L'],'macros':refinement['macro_count'],'native_pulses_upper':refinement['native_pulses_upper'],
          'precision_bits':refinement['precision_bits'],'complete_observed_TV_bound':refinement['complete_observed_record_TV_bound']},
      'finite_size2_small_coupling_transfer':read(HERE/'CONTROLS.json')['finite_size2_small_coupling_connection'],
      'all_failures_retained':True,'success_normalized':False,
      'limitations':['g=0 is a neutral/free source intervention, not interacting charged-particle KG or a g=1 mass pole.',
        'The scalar c=1/sqrt48 source and full90 Wilson Hamiltonian are not identified with the routed Dirac source.',
        'L=2^18 native word and full tensor-field holding iterators are finite unexecuted constructions.',
        'Complete records refer to the declared fresh terminal apparatus, not a transferred full-field or per-macro failure tape.',
        'Size-two field error/cost is not transferred to a larger-L interacting field register.',
        'Coefficient and clock choices retain their input ancestry; no parameter-free empirical spectrum or unique construction is asserted.'],
      'physical_rows_validated':[],'independent_reviews':[],'official_rows_changed':[],
      'source_pins':pins}
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    validation={'schema':'oph.dirac_field_kg.review_request.v1','producer_status':result['status'],
      'literal_criterion':'emergent omega^2 = k^2 + m^2 for matter waves','full_criterion_met':False,
      'reason':'Producer evidence frozen for two independent source-specific full-criterion reviews.',
      'requested_review':'Assess actual same-source free intervention, spatial/native transfer, complete clock/readout and controlled limits. Do not infer interacting g=1 or scalar/Wilson equality.',
      'independent_reviews':[],'official_rows_changed':[]}
    (HERE/'VALIDATION.json').write_text(json.dumps(validation,indent=2)+'\n')
    manifest={'schema':'oph.dirac_field_kg.connection_freeze.v1','candidate_id':result['candidate_id'],
      'files':{str(p.relative_to(ROOT)):sha(p) for p in sorted(HERE.rglob('*')) if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts},
      'source_pins':pins,'physical_rows_validated':[],'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for section in ('files','source_pins'):
        for name,digest in manifest[section].items():assert sha(ROOT/name)==digest,name
    print(json.dumps({'manifest':sha(HERE/'MANIFEST.json'),'proof':sha(HERE/'PROOF.md'),'result':sha(HERE/'RESULT.json'),
                      'files':len(manifest['files']),'source_pins':len(pins),'runtime_imports':len(runtime)}),flush=True)
if __name__=='__main__':main()
