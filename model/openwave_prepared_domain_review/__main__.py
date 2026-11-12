"""Review receipts without invoking their scientific producers."""
import argparse,hashlib,json
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PACKAGES=('openwave_calibrated_hadrons','openwave_nuclear_exchange','openwave_nuclear_refinement','openwave_reachable_averaging','openwave_meson_ring','openwave_dirac_neutral_current','openwave_geometry_refinement','openwave_relational_clock')

def car_identity():
    checks=0
    for b in range(5):
        total=np.zeros((64,64),dtype=np.int64)
        for word in combinations(range(6),b):
            A=np.zeros((64,64),dtype=np.int64)
            for mask in range(64):
                state=mask;coefficient=1
                for mode in reversed(word):
                    if not(state>>mode)&1:coefficient=0;break
                    coefficient*=(-1)**((state&((1<<mode)-1)).bit_count());state^=1<<mode
                if coefficient:A[state,mask]=coefficient
            total+=A.T@A
        expected=np.diag([comb(mask.bit_count(),b)if mask.bit_count()>=b else 0 for mask in range(64)])
        assert np.array_equal(total,expected);checks+=1
    return {'independent_exact_full_Fock_modes':6,'integer_operator_identities_checked':checks,'identity':'sum_I a_I^dagger a_I = binomial(N,b), b=0..4; every matrix entry checked'}

def run():
    data={};pins={};count=0
    for name in PACKAGES:
        base=ROOT/'model'/name;d=json.loads((base/'results.json').read_text());data[name]=d
        for rel,want in d['source_sha256'].items():
            actual=hashlib.sha256((ROOT/rel).read_bytes()).hexdigest();assert actual==want,(name,rel);count+=1
        for path in base.iterdir():
            if path.is_file()and path.suffix in('.py','.md','.tex','.json'):pins[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    h=data['openwave_calibrated_hadrons']
    assert len(h['fits'])==3
    for fit in h['fits']:
        lo,hi=map(F,fit['computed_energy_interval_GeV']);assert lo<=F(fit['target_GeV'])<=hi
        assert hi-lo<F(1,10**8)
    assert all(not row['physical_assignment_supported']for row in h['held_out_comparisons'])
    assert F(h['exact_beta_available_energy_GeV'])>0
    w=data['openwave_reachable_averaging'];assert not w['source_truncation_or_added_penalty']
    assert F(w['exact_resonant_block_centred_norm'])==F(5,4)
    coefficient=F(w['complete_primitive_operator_bound'])+F(44,35)*(F(w['complete_current_times_primitive_operator_bound'])+F(5,4)*F(w['complete_primitive_operator_bound']))
    assert coefficient==F(w['uniform_first_peak_error_coefficient'])
    assert coefficient/F(2**30)<F(1,318)
    ring=data['openwave_meson_ring'];assert F(ring['initial_quadratic_emission_coefficient_exact'])==F(3,512)
    assert ring['exact_memory_kernel_parent_dimension']==6
    assert not ring['physical_strange_decay_lifetime_or_observed_meson_identity_claimed']
    neutral=data['openwave_dirac_neutral_current'];assert len(neutral['exact_native_spin_contractions'])==32
    for member in neutral['complete_massive_members']:
        assert member['full_free_source']['full_CAR_modes']==288
        gram=np.array(member['native_current_intertwiner_Gram_real']);assert np.max(abs(gram-np.eye(2)))<1e-12
    for row in neutral['complete_massive_members'][0]['complete_joint_flight_records']['complete_records']:
        assert row['production_filter_failure_probability']>0
        assert abs(row['original_record_probability_sum']-1)<1e-12
    control=neutral['autonomous_unfiltered_short_time_control']
    lo,hi=map(F,control['total_three_body_probability_interval']);assert 0<lo<hi
    n=data['openwave_nuclear_refinement']
    for row in n['same_source_graphs']:
        L=row['vertices'];assert row['nuclear_total_low_dimension_including_magnetic_components']==50*L+400*comb(L,2)
        assert any(x['spin']==0 and x['isospin']==0 and x['strictly_bound']for x in row['binding_channels'])
    g=data['openwave_geometry_refinement'];assert g['full_source_dimension']==64*25*50
    assert 0<g['certified_compression_comparisons'][-1]['full_arrival_instrument_reference_stable_TV_upper']<.000358
    assert not g['infinite_geometry_cap_or_spatial_continuum_limit_proved']
    c=data['openwave_relational_clock'];assert all(abs(x-1)<1e-10 for x in c['complete_probability_totals'])
    assert c['delay_in_calibrated_native_clock_cycles_interval'][0]>0
    assert c['actual_aggregate_arrival_clock_phase_response_cycles_interval'][0]>0
    assert not c['physical_seconds_Einstein_proper_time_or_fundamental_global_clock_assumed']
    for path in HERE.iterdir():
        if path.suffix in('.py','.md'):pins[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    return {'schema':'oph.independent_prepared_domain_review.v1','verdict':'PASS_REVIEWED_FINITE_ALGEBRA_AND_SOURCE_DOMAINS','source_sha256':pins,
      'direct_producer_source_bindings_checked':count,'independent_CAR_number_lemma_check':car_identity(),
      'reviewed_packages':list(PACKAGES),'review_scope':'Independent proof/source review plus exact CAR operator check, rational fit/averaging arithmetic, complete representation census, record normalization and source-custody checks. Numerical producers and their focused tests are separately replayed; this review does not add a Lean proof or an official physical criterion award.',
      'domain_findings':['Hadron targets are calibration inputs; all held-out mass failures remain.',
      'Complete local nuclear census retains an extra bound S0,I0 channel and finite graph dependence.',
      'Number-sector norms bound the weak primitive and its current image, not the entire evolving-state particle number.',
      'Meson volume bounds compare common embedded instruments, not distinct finite-ring and infinite-line momentum spectral measurements. The exact injection coefficient is3/512.',
      'The complete Dirac source retains all144 channels, tiny positive production failures, active-port work and its declared scale. Its controlled flight is distinct from the short-time full-contact theorem.',
      'Geometry cap comparison is fixed-master and admits native00/10, right0. The right-cell clock experiment instead uses full-master00/01/10/11 blocks.',
      'Clock setting outcomes are alternative X/Y measurements, with full failure and independent phase-calibration runs.']}

def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True);group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');a=parser.parse_args();r=run();p=HERE/'results.json'
    if a.write:p.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
    else:assert json.loads(p.read_text())==r
    print(r['verdict'],r['direct_producer_source_bindings_checked'],'source bindings')
if __name__=='__main__':main()
