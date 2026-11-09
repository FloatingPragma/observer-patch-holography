from fractions import Fraction as F
import numpy as np
import pytest
from model.particle_ground_spectral_records import records
from . import cutoff as c


def test_unchanged_native208_matrix_and_all_sectors():
    assert c.rotor_blocks(6)==records.rational_blocks()
    for k in (6,12,18):
        blocks=c.rotor_blocks(k)
        assert sorted(i for ids,_,_ in blocks for i in ids)==list(range(16*(2*k+1)))
        assert {q for _,q,_ in blocks}==set(range(5))


def test_exact_retained_energy_monotonicity_and_cutoff_effect():
    results=[c.sector_intervals(k) for k in (6,12,18)]
    assert results[1][2][1]<results[0][2][0]
    for a,b in zip(results,results[1:]):
        for q in a:assert b[q][0]<=a[q][1]
    for r in results:assert r[1]==r[3]


def test_parent_enclosures_include_larger_retained_sources():
    for k in (6,12):
        report=c.rotor_limit_control(k);larger=c.sector_intervals(k+6)
        for q,(lo,hi) in report['uncut_sector_intervals'].items():
            assert lo<larger[q][0]<=larger[q][1]<=hi+F(1,1<<24)
    row=c.rotor_limit_control(18)
    assert row['uncut_adjacent_charge_threshold_intervals'][1][0]>1


def test_two_block_bound_against_coupled_hermitian_spectra():
    for a,b,z in [(-2,10,3),(0,1,4),(3,9,0),(-4,8,2)]:
        h=np.array([[a,z],[z,b]],dtype=float)
        assert float(c.block_floor(a,b,z))<=np.linalg.eigvalsh(h)[0]+1e-14


@pytest.mark.parametrize('args',[(0,0,1),(1,0,1),(0,1,-1),(0.,1,1)])
def test_invalid_block_hypotheses_rejected(args):
    with pytest.raises(ValueError):c.block_floor(*args)


def test_exact_receipt_keeps_limit_domains_separate():
    r=c.certificate()
    assert r['fixed_volume_ground_energy_convergence']
    assert not r['joint_volume_spacing_cutoff_limit_proved']
    assert not r['physical_particle_or_Einstein_identification']
