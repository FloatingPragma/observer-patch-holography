"""Source-bound exact finite transfer-energy histogram."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from pathlib import Path

from tools import strict_json
from ..native_exchange_clock import packet as clock
from ..native_vacuum_response import packet as ground
from ..native_fair_local_preparation import packet as preparing
from ..native_correlation_sectors import records
from . import energy as e

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple('model/native_energy_identification/'+name for name in (
    '__init__.py', 'energy.py', 'packet.py', 'test_energy.py', 'PROOF.md', 'README.md'))


def source_bindings():
    paths = set(ground.source_bindings()) | set(preparing.source_bindings()) | set(SOURCE_FILES)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def local_programs(first, second, lags, blocks=64):
    """Same fair H66 source, native preparation, complete four-E scores."""
    pair = [{'op':'repair', 'selector':{'kind':'fixed','seam':list(edge)},
             'law':'E', 'repeat':1} for edge in (first,second)]
    source = preparing.p.source
    programs = []
    for lag in lags[1:]:
        body = [{'op':'repair', 'selector':{'kind':'uniform'}, 'law':'mixture', 'repeat':lag}]
        program = source.program(biased=False, commands=preparing.p.commands(blocks)+pair+body+pair,
                                 capacity=1, overflow='export_oldest')
        programs.append({'lag_attempts':lag, 'program_digest':source.digest(program),
                         'resources':program['declared_resources']})
    return programs


def produce():
    before = source_bindings()
    family, first, second, probe = ground.local_case()
    local = e.exact_case(family, probe, 2, 16)
    variance = e.response.preparation.variance(12,2)
    error = F(1,2**63)/variance
    local.update({'occupation':2, 'probe_seams':(first,second),
        'working_reference':'I/66 on the fair total-two carrier',
        'preparation':preparing.p.budget(64),
        'stationary_probe_variance':variance,
        'normalized_prepared_score_error_upper':error,
        'complete_setting_programs':local_programs(first,second,local['detector']['lag_attempts']),
        'estimation':e.estimation_budget(local['detector'],variance,error,F(1,100),F(1,100)),
        'resolution':e.resolution_budget(30,5,2,16,F(1,3))})
    value = {'schema':'oph.native_energy_identification.v1', 'source_bindings':before,
        'native_complete_record_attachment':records.pair_certificate(),
        'native_generator_attachment':clock.check_source_generator(),
        'local_prepared_histogram':local,
        'growing_source_detector_budgets':[
            {'level':level, 'ports':240*4**level,
             'detector':e.plan(630*4**level,6,2,32),
             'resolution':e.resolution_budget(630*4**level,6,2,32,F(1,4)),
             'global_reference_preparation_supplied_by_this_budget':False}
            for level in (0,2,20)],
        'theorem':{
            'effects':'B_k(Q)=binomial(n,k)*Q^k*(I-Q)^(n-k), k=0,...,n',
            'complete_positive_partition':'B_k>=0; sum_k B_k=I',
            'counted_energy':'H=-2*M*log(P)=-(2*M/b)*log(Q)',
            'native_lag_identity':'<u,Q^j u>=C(b*j)/v_q',
            'full_vector_amplification':'sum_k sum_(j>=1)|c_kj|=3^n-1',
            'finite_sample_boundary':'signed estimate, without clipping or renormalization',
            'scale_conversion':'physical energy requires one common hbar/time conversion',
            'degree_domain':'primitive degree r; the four-E attachment uses r=2 and2<=q<=V-2'},
        'scope':{
            'complete_transfer_space_energy_effects':True,
            'exact_native_record_statistics':True,
            'local_preparation_from_arbitrary_input_with_untouched_reference':True,
            'all_source_raw_records_and_resource_stops_retained':True,
            'single_run_native_coherent_energy_measurement_executed':False,
            'all_displayed_independent_repetitions_executed':False,
            'all_global_canonical_references_prepared':False,
            'energy_bins_identify_particle_species':False,
            'physical_energy_unit_or_mass_ratio_identified':False,
            'interacting_relativistic_continuum_constructed':False,
            'new_canonical_physical_premise_discharged':False}}
    if before != source_bindings():
        raise ArithmeticError('source bytes changed during the energy replay')
    return clock.encode(value)


def verify(path=OUTPUT):
    if clock.canonical(strict_json.load(Path(path))) != clock.canonical(produce()):
        raise ValueError('native energy histogram receipt differs from exact replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(clock.canonical(produce()))
        print('Native transfer-energy histogram written.')
    else:
        verify()
        print('Native transfer-energy histogram verified.')
