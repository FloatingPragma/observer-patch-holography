import json
from fractions import Fraction as F

import numpy as np
import pytest

from ..defect_field_bridge import dynamics_packet
from . import packet


@pytest.fixture(scope='module')
def retained():
    report=json.loads(packet.FIELD_REPORT.read_text())
    arrays=dynamics_packet.load_arrays(report,packet.FIELD_ARRAYS)
    return report,arrays


def test_same_field_clock_join_retains_all_records_and_quantum_outputs(retained):
    report,arrays=retained
    rows,clocks=packet.append_to_field(report,arrays)
    assert len(rows)==len(clocks)==24
    assert sum(row['new_clock_pointer_shape'][0]*row['new_clock_pointer_shape'][1] for row in rows)==96
    assert any(row['native_record']['abort'] for row in rows)
    seen=0
    for row in rows:
        matrix=clocks[row['array_prefix']+'_clock_work_and_purifiers']
        assert matrix.shape==(2,2,4,4)
        assert abs(np.linalg.norm(matrix)**2-1)<1e-12
        assert abs(sum(row['joint_clock_marginals'].values())-row['prior_probability_mass'])<1e-12
        assert row['source_clock_resources']['fresh_pointer_blanks']==2
        inverse=row['record_only_inverse']
        if 'energy_gap_GeV_interval' in inverse:
            seen+=1
            assert inverse['record_probability_grid']=='1/1000000000'
            assert not inverse['source_energies_or_program_duration_used_by_inverse']
            lo,hi=map(F,inverse['energy_gap_GeV_interval'])
            assert 0<lo<hi
    assert seen>0


def test_missing_old_record_is_not_silently_reconstructed(retained):
    report,arrays=retained
    bad=dict(arrays);bad.pop(next(iter(bad)))
    with pytest.raises(ValueError,match='missing or extra'): packet.append_to_field(report,bad)


def test_wrong_prior_quantum_probability_is_rejected(retained):
    report,arrays=retained
    bad=dict(arrays)
    name=next(key for key in bad if key.endswith('_all_pointer_probabilities') and bad[key].sum()>.1)
    bad[name]=np.zeros_like(bad[name])
    with pytest.raises(ValueError,match='disagree'): packet.append_to_field(report,bad)
