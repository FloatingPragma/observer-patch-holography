"""Reuse the existing source-prepared count-derived transfer reference."""
from fractions import Fraction as F
from pathlib import Path

from model.native_phase_experiment import independent, verify
from tools import strict_json
from . import spectral

ROOT=Path(__file__).resolve().parents[2]


def certificate():
    """Exact state/clock identity plus the owning numerical admission check.

    The source phase experiment's numerical preparation is not promoted to
    an outward-rounded complete-channel certificate.
    """
    oracle=independent.compute(ROOT,compact=False)
    verified=verify.verify(strict_json.load(ROOT/'model/native_phase_experiment/results.json'),oracle)
    counts=spectral.source_counts()
    states=tuple(tuple(row) for row in oracle['carrier']['states'])
    if states!=tuple(spectral.synthesis.tomo.compiled.quantum.STATES):
        raise ArithmeticError('the prepared and measured source configuration charts differ')
    transition={(i,j):F(x) for i,j,x in oracle['preparation']['P']['entries']}
    rho={(i,j):F(x) for i,j,x in oracle['preparation']['rho']['entries']}
    trace=sum(F(counts[i][i],60) for i in range(66))
    if trace!=56 or F(oracle['preparation']['trace_P'])!=trace:
        raise ArithmeticError('the count-derived preparation has a different transfer trace')
    for i in range(66):
        for j in range(66):
            if transition.get((i,j),F())!=F(counts[i][j],60) or rho.get((i,j),F())!=F(counts[i][j],60*56):
                raise ArithmeticError('the existing prepared reference is not the measured fair P/56')
    observable={(i,j):F(x) for i,j,x in oracle['preparation']['observable']['entries']}
    if (any(i==j or value!=1 or observable.get((j,i))!=1 for (i,j),value in observable.items())
            or len({i for i,_ in observable})!=len(observable) or not observable):
        raise ArithmeticError('the source exchange observable lost its norm-one partial-swap form')
    contrast=sum((value*rho.get((j,i),F()) for (i,j),value in observable.items()),F())
    if contrast!=F(1,168):
        raise ArithmeticError('the two prepared-reference distinction changed')
    return {'same_source_basis_states':66,'same_transfer_entries_checked':66**2,
            'reference':'rho=P/56','trace_P':56,'faithful_reference_floor':'1/84',
            'exact_modular_identity':'K=-log(rho)=H/60+log(56)*I; H=-60*log(P)',
            'exact_modular_action':'rho^(i*s) A rho^(-i*s)=U_H(s/60) A U_H(s/60)^*, U_H(t)=exp(-i*t*H)',
            'controlled_transfer_reference_stationary':True,
            'unread_native_repair_reference_stationary':False,
            'source_preparation_owner':'model/native_phase_experiment',
            'preparation_ideal_timeout_probability':str(F(23,128)**16),
            'preparation_numerical_frobenius_admission':'1/5000000000',
            'owning_independent_numerical_verifier_passed':verified['all_pass'],
            'diagonal_tree_reference_separation':{
                'norm_one_source_exchange_observable':'ports0,1 occupation exchange',
                'count_derived_rho_expectation':'1/168','diagonal_tree_reference_expectation':'0',
                'half_trace_distance_lower':'1/336',
                'regional_tree_reference_modular_properties_transferred_to_rho':False},
            'preparation_outward_rounded_channel_error_supplied':False,
            'proper_code_load_and_export_interfaces_supplied':True,
            'physical_pulse_to_attempt_time_conversion_derived':False}
