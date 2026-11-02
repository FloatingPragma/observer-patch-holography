import numpy as np
import pytest
from . import source,packet


def test_full_fock_current_identity_and_signed_group():
    r=source.source_certificate()
    assert r['full_Fock_continuity_dimension']==65536
    assert r['density_and_twice_spin_continuity_checked'] is True


def test_every_neutral_state_is_counted_once(result):
    r=result['spectrum']['global']
    assert len(r['blocks'])==60
    assert r['coverage']['balanced_dimension']==4900
    assert r['coverage']['full_N8_dimension']==12870
    assert r['lowest_excitation_multiplicity']==3


@pytest.mark.parametrize('s,e,k',[(True,0,0),(0,False,0),(0,0,True),(-1,0,0),(3,2,0),(0,0,8)])
def test_wrong_domains_rejected_even_with_warmed_cache(s,e,k):
    source.block(0,0,0)
    with pytest.raises(ValueError):source.block(s,e,k)


def test_non_tuple_and_duplicate_columns_are_rejected():
    w=packet.load_witness();cols=tuple(map(int,w['global_0_0_0_columns']))
    source.block(0,0,0,columns=cols)
    with pytest.raises(ValueError):source.block(0,0,0,columns=list(cols))
    with pytest.raises(ValueError):source.block(0,0,0,columns=(cols[0],)*len(cols))


def test_conserved_channels_are_removed_by_exact_operators():
    from ..centered_particle_spectrum import packet as old
    w=old.load_witness();g=source.parent.reduced('ground',tuple(map(int,w['ground_columns'])))
    for s,e,forbidden in ((1,0,0),(0,1,7)):
        b=source.block(s,e,kind='probe')
        image=source.probe_operator(s)@g['full_basis']
        assert not np.any(source.translation_sum(b['up'],b['down'],forbidden)@image)
        assert source.transition(s,g,b).shape==(136,32)
