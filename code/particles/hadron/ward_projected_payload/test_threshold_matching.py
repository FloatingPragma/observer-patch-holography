"""The independent Ward payload retains the same conditional QCD convention."""
import importlib.util
from pathlib import Path

import pytest

import spectral_modules as spectral


def _transmutation():
    path=Path(__file__).resolve().parents[2]/'qcd/derive_lambda_qcd_source_transmutation.py'
    spec=importlib.util.spec_from_file_location('checked_transmutation',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('alpha,nl',[(.1,3),(.2,4),(.3,5)])
def test_down_and_up_matching_retain_the_same_direction(alpha,nl):
    down=spectral._decouple_down(alpha,nl)
    assert down>alpha
    assert down==_transmutation().decouple_down(alpha,nl)
    assert abs(spectral._decouple_up(down,nl)-alpha)<1e-12


def test_lambda_grid_is_the_transmutation_corner_range():
    result=_transmutation().build()
    source=result['source_input']['mz_run_gev_display']
    assert source==spectral.LAMBDA_LANE_MZ_DISPLAY_GEV
    expected=[*result['lambda3_interval_gev'][:1],result['central']['lambda3_gev'],
              result['lambda3_interval_gev'][1]]
    for key,value in zip(('lane_lo','lane_central','lane_hi'),expected,strict=True):
        assert spectral.LAMBDA3_GRID[key]==value/source
