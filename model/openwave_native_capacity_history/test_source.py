from .source import run

def test_complete_history_and_halt_control():
    result=run()
    assert result['all_history_blocks']==625
    assert result['exported_archive_plus_FIFO_kernel_identical']
    assert result['overflow_negative_control']['one_slot']=='FIFO_TIMEOUT'
    assert not result['normalization_symmetries']['H0_or_Lambda_or_primordial_amplitude_selected']
