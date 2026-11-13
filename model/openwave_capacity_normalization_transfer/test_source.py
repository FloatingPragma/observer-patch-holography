from fractions import Fraction as F
import pytest
from .source import exact_mode,normalization,validated

def test_complete_regular_mode_and_normalization_symmetries():
    d=exact_mode();assert d['all_seven_evolution_and_Hamiltonian_constraint_residuals']==['0']*8
    n=normalization();assert n['scale_residuals_exact']==['0']*3
    assert n['examples'][2]['density_multiplier_exact']=='4'
    assert not n['absolute_dark_energy_or_primordial_amplitude_fixed_by_dimensionless_flow']

def test_radiation_and_dust_power_law_controls():
    # H=C a^-p makes the complete improper integral exact.
    for p,expected in((F(2),F(1)),(F(3,2),F(9,10))):
        T=F(3,2)/(p+1);phi=F(3,2)-T
        assert phi==expected and phi-p*T==0

def test_validated_refinement_enclosures_overlap_and_narrow():
    coarse=validated(F(1,64));fine=validated()
    for a,b in zip(coarse['records'],fine['records'],strict=True):
        for kind in('capacity','lambda'):
            lo,hi=map(F,a[kind]['metric_transfer_interval']);ll,hh=map(F,b[kind]['metric_transfer_interval'])
            assert max(lo,ll)<min(hi,hh) and hh-ll<(hi-lo)/2
    present=fine['records'][0]
    assert F(present['capacity']['metric_transfer_interval'][0])>F(present['lambda']['metric_transfer_interval'][1])
    assert F(fine['initial_integral_tail_upper_exact'])<F(1,10**13)
    with pytest.raises(ValueError):validated(F(1))
    with pytest.raises(ValueError):validated(tail_t=F(1))
