import numpy as np
from scipy import linalg as la
import pytest

from . import overlap as o


def test_complete_kernel_and_random_cross_site_leakage():
    v=o.build()
    rng=np.random.default_rng(5701)
    f=rng.normal(size=len(v["raw"]["points"]))
    row=o.profile_response(f)
    literal=v["W"].conj().T @ (np.repeat(f,2)[:,None]*v["c"]) @ v["W"].conj()
    assert la.norm(literal-row["vertex"]) < 1e-11
    for x,y in ((0,0),(2,173),(160,160),(169,171),(342,342)):
        px=np.zeros(343);py=np.zeros(343);px[x]=1;py[y]=1
        lx=v["Q"].conj().T @ (np.repeat(px,2)[:,None]*(v["c"]@v["W"].conj()))
        ly=v["Q"].conj().T @ (np.repeat(py,2)[:,None]*(v["c"]@v["W"].conj()))
        bx,by=v["density"][x],v["density"][y]
        expected=v["pair"].conj().T @ ((bx if x==y else 0)-bx@by) @ v["pair"]
        assert la.norm(lx.conj().T@ly-expected) < 1e-11


def test_arbitrary_complex_band_frame_changes_no_strength():
    v=o.build();rng=np.random.default_rng(771)
    u,_=la.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    profiles=o.diagnostic_profiles()
    for f in profiles.values():
        row=o.profile_response(f)
        transformed=u.conj().T@row["vertex"]@u.conj()
        w=v["W"]@u
        direct=w.conj().T @ (np.repeat(f,2)[:,None]*v["c"]) @ w.conj()
        assert la.norm(transformed-direct) < 1e-11
        assert np.allclose(la.svdvals(transformed),row["singulars"],atol=1e-11)


def test_rotated_complete_archive_is_accepted_but_changed_sites_are_rejected():
    v=o.build();rng=np.random.default_rng(98)
    u,_=la.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    packet={"points":np.asarray(v["raw"]["points"]),
            "band_energies":v["energies"][v["mask"]], "W":v["W"]@u,
            "pair":u.conj().T@v["pair"]@u.conj(),
            "density":u.conj().T@v["density"]@u,
            "kernel":u.conj().T@v["kernel"]@u.conj()}
    o.validate_archive(packet)
    packet["points"]=packet["points"].copy()
    packet["points"][0,0]+=1
    with pytest.raises(ArithmeticError,match="site alphabet"):
        o.validate_archive(packet)


def test_uniform_no_hierarchy_and_local_patch_two_strengths():
    profiles=o.diagnostic_profiles()
    uniform=o.profile_response(profiles["uniform"])
    assert np.allclose(uniform["singulars"],1,atol=1e-11)
    assert uniform["leakage"] < 1e-11
    core=o.profile_response(profiles["positive_core"])
    assert np.count_nonzero(core["singulars"]>1e-11)==2
    patch=o.profile_response(profiles["positive_core_patch"])
    assert 16000 < patch["singulars"][0]/patch["singulars"][-1] < 17000
    assert np.allclose(patch["singulars"][::2],patch["singulars"][1::2],atol=1e-11)
    scaled=o.profile_response(profiles["positive_core_patch"]/np.sqrt(27))
    assert np.allclose(scaled["singulars"]*np.sqrt(27),patch["singulars"],atol=1e-11)
    assert patch["leakage"] > 0.007


def test_all_residues_and_complete_full_intervention():
    rows=o.core_residue()
    assert len(rows)==4
    assert all(row["raw_core_norm_squared"]<1 for row in rows)
    assert all(abs(row["residue_normalized_norm_squared"]-1)<1e-11 for row in rows)
    f=o.diagnostic_profiles()["positive_core_patch"]/np.sqrt(27)
    row=o.full_intervention(f)
    assert row["full_Nambu_rows"]==1372 and row["retained_low_rows"]==8
    assert row["maximum_spectral_difference"] < row["Schur_self_energy_norm_bound_diagnostic"]+1e-11
    assert row["residue_metric_distance_from_identity_bound_diagnostic"] < 1e-10


def test_geometric_first_hop_scale_without_family_parameters():
    rows=o.texture_controls()
    assert [row["texture_amplitude"] for row in rows]==[32,64,128]
    coefficients=[row["weak_strength_times_amplitude_squared"] for row in rows]
    assert all(0.249 < x < 0.251 for x in coefficients)
    assert abs(coefficients[-1]-.25)<abs(coefficients[0]-.25)


@pytest.mark.parametrize("value", [np.ones(342),np.full(343,np.nan),np.ones(343,dtype=complex)])
def test_incomplete_or_invalid_scalar_intervention_rejected(value):
    with pytest.raises(ValueError):
        o.profile_response(value)
    with pytest.raises(ValueError):
        o.full_intervention(value)
