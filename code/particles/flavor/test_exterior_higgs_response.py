"""Mutation tests for the exact full-doublet response and its claim boundary."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import pytest


HERE=Path(__file__).resolve().parent


def load(name):
    spec=importlib.util.spec_from_file_location(name,HERE/(name+".py"))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


producer=load("exterior_higgs_response_certificate")
verifier=load("verify_exterior_higgs_response")


@pytest.fixture(scope="module")
def receipt():
    return json.loads((producer.ROOT/producer.OUTPUT).read_text())


def test_independent_car_verifier_accepts(receipt):
    assert verifier.verify(receipt)["verified"] is True


def test_fresh_producer_matches(receipt):
    assert producer.build()==receipt


@pytest.mark.parametrize("sector,index,row,column",[
    ("up_higgs_component_matrices",0,1,0),
    ("up_higgs_component_matrices",1,0,0),
    ("up_higgs_component_matrices",0,0,2),
    ("up_higgs_component_matrices",1,4,2),
    ("down_conjugate_higgs_component_matrices",0,0,0),
    ("down_conjugate_higgs_component_matrices",1,1,0),
    ("down_conjugate_higgs_component_matrices",0,4,2),
    ("down_conjugate_higgs_component_matrices",1,0,1),
])
def test_charged_neutral_and_color_tensor_corruption(receipt,sector,index,row,column):
    changed=deepcopy(receipt)
    changed[sector][index][row][column]+=1
    with pytest.raises(verifier.VerificationError,match="TENSOR"):
        verifier.verify(changed)


@pytest.mark.parametrize("claim",[
    "source_selected_record_state","source_selected_Higgs_dependent_quark_kernel",
    "physical_kinetic_normalization_emitted","g_u_g_d_selected",
    "family_hierarchy_derived","quark_masses_emitted","mass_scheme_emitted",
    "empirical_prediction",
])
def test_physical_promotion_rejected(receipt,claim):
    changed=deepcopy(receipt)
    changed["claims"][claim]=True
    with pytest.raises(verifier.VerificationError,match="SPEC_PROJECTION:claims"):
        verifier.verify(changed)


@pytest.mark.parametrize("key,value",[
    ("within_sector_singular_value_multiplicities",[1,1,1]),
    ("canonical_up","g_u I3"),
    ("canonical_down","g_d I3 / sqrt(z_Q z_u Z_H)"),
    ("CKM_identifiable_on_triply_degenerate_branch",True),
    ("up","diag(g_u,g_c,g_t)"),
    ("normalization_convention","all-left Weyl coefficient"),
])
def test_family_and_normalization_corruption(receipt,key,value):
    changed=deepcopy(receipt)
    changed["family_response"][key]=value
    with pytest.raises(verifier.VerificationError,match="FAMILY_RESPONSE"):
        verifier.verify(changed)


@pytest.mark.parametrize("key,value",[
    ("color_commutant_constraint_rank",6),
    ("color_commutant_complex_dimension",3),
    ("shared_metric_generalized_root_multiplicity",1),
    ("kinetic_congruence_matrix",[["1","0","0"],["0","2","0"],["0","0","3"]]),
])
def test_forged_exact_checks(receipt,key,value):
    changed=deepcopy(receipt)
    changed["exact_checks"][key]=value
    with pytest.raises(verifier.VerificationError,match="EXACT_CHECKS"):
        verifier.verify(changed)


def test_higgs_conjugation_is_not_optional(receipt):
    changed=deepcopy(receipt)
    changed["equations"]["full_down"]="I_color tensor column(-h0,-h1)"
    with pytest.raises(verifier.VerificationError,match="SPEC_PROJECTION:equations"):
        verifier.verify(changed)


def test_cancellation_is_conditional_not_selected(receipt):
    changed=deepcopy(receipt)
    changed["hypotheses"][-1]="The source selects Y=g R."
    with pytest.raises(verifier.VerificationError,match="SPEC_PROJECTION:hypotheses"):
        verifier.verify(changed)


def test_no_empirical_numeric_input(receipt):
    changed=deepcopy(receipt)
    changed["numeric_input_paths"].append("code/particles/data/particle_reference_values.json")
    with pytest.raises(verifier.VerificationError,match="NUMERIC_INPUT_PATHS"):
        verifier.verify(changed)


def test_no_unlisted_ancestry(receipt):
    changed=deepcopy(receipt)
    changed["source_pins"]["target-masses.json"]="0"*64
    with pytest.raises(verifier.VerificationError,match="PIN_PATHS"):
        verifier.verify(changed)


def test_no_unreported_output(receipt):
    changed=deepcopy(receipt)
    changed["masses"]=[1,2,3,4,5,6]
    with pytest.raises(verifier.VerificationError,match="RECEIPT_KEYS"):
        verifier.verify(changed)


def test_source_pin_required(receipt):
    changed=deepcopy(receipt)
    changed["source_pins"][verifier.INPUTS[0]]="0"*64
    with pytest.raises(verifier.VerificationError,match="SOURCE_PIN"):
        verifier.verify(changed)


@pytest.mark.parametrize("kind",["scalar","charges","family"])
def test_rehashed_source_semantic_mutation_rejected(receipt,tmp_path,kind):
    for source in verifier.PIN_PATHS:
        dest=tmp_path/source
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(producer.ROOT/source,dest)
    changed=deepcopy(receipt)
    path=verifier.INPUTS[1 if kind=="family" else 0]
    data=json.loads((tmp_path/path).read_text())
    if kind=="scalar":
        data["exterior_matter_contract"]["one_scalar"]="color_block"
    elif kind=="charges":
        data["exterior_matter_contract"]["block_trace_charges"]["color_block"]="-1/6"
    else:
        data["attachment"]["family_dimension"]=4
    (tmp_path/path).write_text(json.dumps(data))
    changed["source_pins"][path]=hashlib.sha256((tmp_path/path).read_bytes()).hexdigest()
    with pytest.raises(verifier.VerificationError,match="SOURCE_(SCALAR|CHARGES|FAMILY)"):
        verifier.verify(changed,root=tmp_path)


def test_closed_numeric_read_boundary(monkeypatch):
    observed=set()
    read_text,read_bytes=Path.read_text,Path.read_bytes
    def text_spy(path,*args,**kwargs):
        observed.add(path.resolve())
        return read_text(path,*args,**kwargs)
    def bytes_spy(path,*args,**kwargs):
        observed.add(path.resolve())
        return read_bytes(path,*args,**kwargs)
    monkeypatch.setattr(Path,"read_text",text_spy)
    monkeypatch.setattr(Path,"read_bytes",bytes_spy)
    producer.build()
    assert observed=={(producer.ROOT/p).resolve() for p in producer.SOURCES+producer.IMPLEMENTATION}
