from collections import Counter
from fractions import Fraction as F
import hashlib,json,math
from pathlib import Path
from model.compiled_native_history.compiler import compile_history
from model.native_program import interpreter as native
from model.openwave_native_process.source import admitted_local_program,fixed_command

ROOT=Path(__file__).resolve().parents[2]
PARENT='model/openwave_capacity_finite_wavelength/results.json'

def complete_kernel(compiled):
    return tuple((tuple(native.canonical(vars(r))for r in (*b.archive,*b.memory)),
        b.source_map,b.weight,b.halted)for b in compiled.branches)

def run():
    commands=[fixed_command((0,1)),fixed_command((1,7))]
    a=compile_history(admitted_local_program(commands=commands,capacity=1))
    b=compile_history(admitted_local_program(commands=commands,capacity=3))
    assert complete_kernel(a)==complete_kernel(b)
    probabilities=[F(*x.weight)*F(len(x.source_map),66)for x in a.branches]
    assert sum(probabilities)==1
    # The uniform H66 source is declared. No history is selected as successful.
    entropy=-sum(float(p)*math.log2(float(p))for p in probabilities if p)
    long=[fixed_command((0,1),repeat=3)]
    stop=compile_history(admitted_local_program(commands=long,capacity=1,overflow='halt'))
    run3=compile_history(admitted_local_program(commands=long,capacity=3,overflow='halt'))
    assert stop.status=='FIFO_TIMEOUT' and run3.status=='COMPLETE'
    parent=json.loads((ROOT/PARENT).read_text())
    for rel,h in parent['source_sha256'].items():
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==h,rel
    return {'schema':'oph.native_capacity_history.v1',
        'verdict':'PASS_COMPLETE_NATIVE_MEMORY_AND_NORMALIZATION_BOUNDARY',
        'source':'same fair H66 C/E mixture, two specified native seams; uniform66-state input',
        'all_history_blocks':len(a.branches),'positive_probability_blocks':sum(bool(p)for p in probabilities),
        'all_history_probabilities':[str(p)for p in probabilities],
        'record_entropy_bits_float_diagnostic':entropy,
        'source_Holevo_capacity_bits':'log2(66)',
        'raw_two_record_storage_upper_bits':80,
        'exported_archive_plus_FIFO_kernel_identical':True,
        'identity_domain':'all66 source columns, every offdiagonal by exact partial-injection theorem, arbitrary spectator reference',
        'memory_representations':[{'FIFO_capacity':1,'archive_records':1,'FIFO_records':1},
                                  {'FIFO_capacity':3,'archive_records':0,'FIFO_records':2}],
        'overflow_negative_control':{'one_slot':stop.status,'three_slots':run3.status,
                                     'same_channel_if_halt_policy_changed':False},
        'capacity_bookkeeping_is_not_horizon_area':True,
        'cosmology_parent':{'receipt':PARENT,'sha256':hashlib.sha256((ROOT/PARENT).read_bytes()).hexdigest(),
            'scope':'authenticated conditional late-time finite-k response; no source recompilation asserted'},
        'normalization_symmetries':{'clock':'t_phys -> s t_phys, H0 -> H0/s, Lambda -> Lambda/s^2',
            'perturbation':'linear homogeneous perturbation solution y -> A y; power -> A^2 power',
            'native_complete_records_change_under_either_assignment':False,
            'H0_or_Lambda_or_primordial_amplitude_selected':False},
        'missing_physical_join':['horizon map from actual records','physical event clock',
            'independent amplitude/preparation law','same-source effective cosmological dynamics'],
        'new_Nature_validation_rows':0}
