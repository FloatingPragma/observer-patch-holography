import numpy as np
from fractions import Fraction as F
from .source import chain,block,effective,complement_bounds
from .potential import constrained_potential
from .certificate import one_chain


def test_full_block_inverse_and_current_are_reconstructed():
    d=chain(3);b=block(d,1)
    E=float(np.linalg.eigvalsh(b['H'])[0]-.1)
    r=effective(b,E)
    complete=np.linalg.inv(E*np.eye(len(b['H']))-b['H'])
    coarse=np.linalg.inv(E*np.eye(len(b['P']))-r['H_effective'])
    assert np.max(abs(complete[np.ix_(b['P'],b['P'])]-coarse))<2e-8
    assert complement_bounds(b,F.from_float(E)+F(1,4))['gap']>0


def test_exact_complete_internal_string_shift_and_failures():
    d=constrained_potential()
    assert d['complete_internal_singlet_dimension']==40
    assert d['exact_slope_GeV_per_fm']=='1'
    assert F(d['singlet_success_on_uniform_pinned_input'])+F(d['coloured_failure_on_uniform_pinned_input'])==1


def test_native_full_current_and_complete_time_records():
    d=one_chain(3)
    assert abs(d['native_axial_with_exact_vertex']-d['native_gA'])<2e-8
    diagonal={int(i):F(v)for i,v in d['complete_physical_instrument']['complete_colour_projector_diagonal_nonzero']}
    size=d['complete_physical_instrument']['input_dimension']
    assert sum(diagonal.values())==size
    for row in d['complete_finite_time_records']:
        assert abs(row['coarse_probability']+row['retained_complement_probability']-1)<2e-11
        assert abs(sum(x[1] for x in row['nonzero_success_native_occupations'])-float(F(row['preparation_success'])))<2e-12

        success=dict(row['nonzero_success_native_occupations'])
        assert all(float(v)/size-success.get(i,0)>=-2e-12 for i,v in diagonal.items())
        assert abs(sum(float(v)/size-success.get(i,0)for i,v in diagonal.items())+sum(success.values())-1)<2e-11
