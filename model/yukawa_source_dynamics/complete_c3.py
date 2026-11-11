"""All-outcome family channel from the original exclusive native records.

The source routing R is the actual product of its two half-turns. The
conditional family routing is inherited from the coherent source. Every
one of its 128 original outcomes is retained. Tracing those records produces
an exact dephasing channel; an analyzer outcome is not a mass vertex.
"""
from functools import lru_cache
import numpy as np
import sympy as s
from model.standard_model_source_map import koide_source as source
from model.charged_lepton_response import phase_selection


def clean(a):
    return s.ImmutableMatrix(a.applyfunc(lambda x:s.simplify(s.expand(x))))


@lru_cache(maxsize=1)
def exact():
    records=source.native_partition()
    geometry=phase_selection.native_geometry()
    R=geometry['R']
    branches=[]; channel=s.zeros(9); gram=s.zeros(3)
    for row in records:
        group=row['group']; exponent=1 if group==1 else -1 if group==2 else 0
        weight=s.Rational(row['probability']); rotation=R**exponent
        branches.append((tuple(row['raw']),exponent,weight,rotation))
        gram+=weight*rotation.H*rotation
        channel+=weight*s.kronecker_product(rotation.conjugate(),rotation)
    gram=clean(gram); channel=clean(channel)
    if gram!=s.eye(3): raise ArithmeticError('complete source isometry is not normalized')
    twirl=clean(sum((s.kronecker_product((R**i).conjugate(),R**i) for i in range(3)),s.zeros(9))/3)
    r=s.Rational(19,24)
    if clean(twirl*twirl)!=twirl or s.trace(twirl)!=3:
        raise ArithmeticError('complete cyclic fixed algebra changed')
    if clean(channel-r*s.eye(9)-(1-r)*twirl)!=s.zeros(9):
        raise ArithmeticError('full native channel is not the derived dephasing law')
    # Pure dephasing has no antisymmetric Hilbert--Schmidt component. Its
    # principal real logarithm is log(r)(I-twirl), with no Hamiltonian term.
    if channel!=channel.T or clean(channel*twirl)!=twirl:
        raise ArithmeticError('source channel acquired a coherent generator')
    # Native preparation has sqrt(w) amplitudes. Summing K^dag K keeps all
    # modes, so every singular value of the full Stinespring map is exactly1.
    # Orthogonal record Fourier transforms cannot change this Gram.
    analyzer=tuple(source.analyzer_block(i,0)[:3,:3] for i in range(8))
    effects=tuple(clean(a.H*a) for a in analyzer)
    if clean(sum(effects,s.zeros(3)))!=s.eye(3):
        raise ArithmeticError('an analyzer complement was discarded')
    return dict(R=R,records=tuple(branches),channel=channel,twirl=twirl,
                gram=gram,contraction=r,analyzer_effects=effects)


def report():
    x=exact(); R=x['R']; r=x['contraction']
    a=s.sqrt(source.WEIGHTS[0]); b=s.sqrt(source.WEIGHTS[1])
    singlet=s.simplify((a+2*b)**2/3)
    doublet=s.simplify((a-b)**2/3)
    return dict(native_raw_outcomes=len(x['records']),total_probability_exact=str(sum(t[2] for t in x['records'])),
                actual_family_rotation=[[str(v) for v in row] for row in R.tolist()],
                raw_shift_probabilities_exact={'zero':'31/36','plus':'5/72','minus':'5/72'},
                native_grouping_probabilities_exact=list(map(str,source.WEIGHTS)),
                canonical_full_compression_gram='I3',full_stinespring_singular_values_exact=['1','1','1'],
                channel_formula='E=(19/24) Id+(5/24) Twirl_C3',fixed_operator_dimension=3,
                channel_spectrum_exact={'1':3,'19/24':6},
                principal_generator='log(19/24)*(Id-Twirl_C3)',coherent_Hamiltonian_exact='0 modulo scalar',
                generator_time='one declared record-routing invocation; no physical clock',
                all_positive_integer_powers='E^n=Twirl_C3+(19/24)^n*(Id-Twirl_C3)',
                all_continuous_embedding_times='E_t=Twirl_C3+(19/24)^t*(Id-Twirl_C3), t>=0',
                comparison_phase_or_N_selected=False,
                analyzer_outcomes=8,analyzer_completeness_exact=True,
                default_analyzer_effect_eigenvalues_exact=[str(doublet),str(doublet),str(singlet)],
                default_analyzer_hierarchy_ratio=float(singlet/doublet),
                default_phase='z=1; inherited source default',
                all_outcome_channel_depends_on_scalar_branch_phase=False,
                equal_plus_minus_weight_implies_pair_degeneracy=True,
                raw_record_probabilities_native=True,controlled_family_routing_primitive_source_forced=False,
                original_record_classical_channel_derives_coherent_phase=False,
                three_charged_sectors_distinguished=False,physical_mass_or_Yukawa_derived=False)
