import hashlib
import json
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np

from . import full_clock, measure


def test_complete_head_clock_and_independent_native_effects():
    record=full_clock.run()
    path = Path(full_clock.__file__).with_suffix(".json")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == "9c98105b9902f60b01671efdacf695aad80a8567e5177bb99ef9fd8bffd2129a"
    saved = json.loads(path.read_text())
    improved = record["certificate"]
    archived = saved["certificate"]
    assert improved.keys() == archived.keys()
    assert improved["inverse_operator_norm_upper"] <= archived["inverse_operator_norm_upper"]
    assert improved["inverse_residual_operator_norm_upper"] <= archived["inverse_residual_operator_norm_upper"]
    assert improved["smallest_singular_value_lower"] >= archived["smallest_singular_value_lower"]
    assert {k: v for k, v in record.items() if k != "certificate"} == {
        k: v for k, v in saved.items() if k != "certificate"
    }
    assert record["rank"]==record["free_real_state_parameters"]==31
    assert len(set(tuple(w) for w in record["kept_clock_bins"]))==31
    assert all(0<t<12 for _,_,t in record["kept_clock_bins"])
    archive=record["affine_clock"]
    denom=int(archive["denominator"])
    a=np.array(archive["numerators"],dtype=float)/denom
    offsets=np.array(archive["offset_numerators"],dtype=float)/denom
    s=full_clock.station
    board=s.iso.pairing.formation.federation(1)
    ports=s.iso.pairing.formation.region(board)
    effects=np.zeros((31,2,4,4),complex)
    for q in (-1,1):
        for parity in range(2):
            h=s.iso.ranging.generator(board,ports,q,8,bool(parity))
            k=s.iso.ranging.kernels(h,ports,.5,12)
            hits={route:s.iso.ranging.first_hits(k,*route) for route in measure.ROUTES}
            for row,(i,j,t) in enumerate(record["kept_clock_bins"]):
                for step in range(1,t):
                    amplitude=hits[j,i][t-step]@hits[i,j][step]
                    effects[row,parity] += amplitude.conj().T@amplitude/2
    # Pauli/2 columns coincide with the original survey basis. The final column
    # parametrizes odd weight directly, and the affine offset is at I_8/8.
    dirs=s.directions()
    dirs[-1]=np.stack((-np.eye(4),np.eye(4)))/4
    expected=np.einsum("rbxy,dbyx->rd",effects,dirs).real
    np.testing.assert_allclose(expected,a,rtol=0,atol=1e-13)
    np.testing.assert_allclose(effects.trace(axis1=2,axis2=3).sum(axis=1)/8,
                               offsets,rtol=0,atol=1e-13)
    true_bound=record["certificate"]["smallest_singular_value_lower"]
    assert 5e-8<true_bound<np.linalg.svd(a,compute_uv=False)[-1]


def test_simultaneous_two_state_entropy_width():
    """Check the uniform finite-shot guarantee against the archived calibration."""
    record = json.loads(Path(full_clock.__file__).with_suffix(".json").read_text(),
                        parse_float=Decimal)
    with localcontext() as ctx:
        ctx.prec = 60
        bins = Decimal(len(record["kept_clock_bins"]))
        starts_per_route = 2**72
        failure = Decimal("0.01")
        bin_error = ((4*bins/failure).ln()/(2*starts_per_route)).sqrt()
        gain = Decimal("1.576e7")
        entry_error = Decimal("1.629e-15")
        offset_error = Decimal("8.145e-16")
        assert gain > record["certificate"]["inverse_operator_norm_upper"]
        assert entry_error > record["affine_clock"]["entry_error_upper"]
        assert offset_error > record["affine_clock"]["offset_error_upper"]
        hs_error = gain*(bins.sqrt()*(bin_error+offset_error) + bins*entry_error)
        t = Decimal(2).sqrt()*hs_error
        entropy_error = t*(Decimal(7)/t).ln() - (1-t)*(1-t).ln()
        assert bin_error < Decimal("3.160e-11")
        assert t < Decimal("0.003922")
        radius = Decimal("0.03328")
        assert entropy_error < radius
        assert Decimal(2).ln() - 4*radius > Decimal("0.5600")
