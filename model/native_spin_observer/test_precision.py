from fractions import Fraction as F
import mpmath
import pytest
from . import precision as p
from ..certified_reference import certificate as cert


def test_radical_and_machin_enclosures_independently():
    ctx=mpmath.mp.clone();ctx.dps=100
    def mp(q):return ctx.mpf(q.numerator)/q.denominator
    for bits in (8,80,240):
        for n in (2,5,9):
            lo,hi=p.sqrt_interval(n,bits);assert lo*lo<=n<=hi*hi
        lo,hi=p.pi_interval(bits)
        assert mp(lo)<=ctx.pi<=mp(hi)
        assert hi-lo<F(1,1<<bits)


@pytest.mark.parametrize('kind',['X','Y','Z','I','H','J'])
def test_dyadic_commands_converge_without_binary64_floor(kind):
    exact=p.compiler._source()[2];ctx=mpmath.mp.clone();ctx.prec=350
    def mp(q):return ctx.mpf(q.numerator)/q.denominator
    coeff=exact[kind] if kind!='H' else [(a[0]+b[0],a[1]+b[1]) for a,b in zip(exact['X'],exact['Z'],strict=True)]
    target=[(mp(a)+mp(b)*ctx.sqrt(5))*ctx.pi/7/(ctx.sqrt(2) if kind=='H' else 1) for a,b in coeff]
    for bits in (24,80,160):
        row=p.dyadic_command(kind,cert._pi(1,7),bits)
        for s,e,x in zip(row['integrated_port_coefficients'],row['port_error_bounds'],target,strict=True):
            q=F(s);assert q.denominator&(q.denominator-1)==0
            assert abs(mp(q)-x)<=mp(F(e))<=mp(F(1,1<<bits))
        assert F(row['generator_and_complete_output_error'])<=29*F(1,1<<bits)


def test_zero_and_invalid_precision():
    assert set(p.dyadic_command('J',F(),80)['integrated_port_coefficients'])=={'0'}
    for value in (True,0,-2,1.5):
        with pytest.raises(ValueError):p.dyadic_command('X',F(1),value)
    with pytest.raises(ValueError):p.dyadic_command('X',.1,80)
    with pytest.raises(ValueError):p.angle_interval(cert.PiMultiple(False),80)
    assert p.bits_for(1000,F(1,100))>10
