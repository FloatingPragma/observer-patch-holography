import copy
import importlib.util
import json
import shutil
import hashlib
from pathlib import Path
import pytest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
RECEIPT=HERE/"runtime"/"k7_gauge_kinetic_cross_control.json"
spec=importlib.util.spec_from_file_location("cross_verifier",HERE/"verify_k7_gauge_kinetic_cross_control_independent.py")
ver=importlib.util.module_from_spec(spec);spec.loader.exec_module(ver)
pspec=importlib.util.spec_from_file_location("cross_producer",HERE/"k7_gauge_kinetic_cross_control.py")
prod=importlib.util.module_from_spec(pspec);pspec.loader.exec_module(prod)

def receipt(): return json.loads(RECEIPT.read_text())
def rejected(mutator):
 d=copy.deepcopy(receipt());mutator(d)
 try: ver.validate(d,ROOT)
 except (ValueError,KeyError,TypeError): return
 raise AssertionError("hostile mutation passed")

def test_clean_packet_passes(): assert ver.validate(receipt(),ROOT)
def test_exact_quadratic_field_operations():
 a=prod.quadratic(1,1); sq=prod.q2mul(a,a)
 assert sq==prod.quadratic(3,2)
 assert prod.q2sub(sq,prod.quadratic(3,2))==prod.quadratic()
 assert prod.q2zero(prod.quadratic()) and not prod.q2zero(sq)
 assert prod.q2scale(a, __import__("fractions").Fraction(3,5))==prod.quadratic(__import__("fractions").Fraction(3,5),__import__("fractions").Fraction(3,5))
def test_mixed_scale_is_control_only():
 d=receipt(); c=d["controls"]["mixed_scale_probe"]
 assert c["classification"]=="NONPHYSICAL_MIXED_SCALE_CONTROL" and not c["eligible_for_framework_cross_verdict"]
def test_b_test_is_control_only():
 d=receipt(); c=d["controls"]["b_test_branch_mismatch_probe"]
 assert c["classification"]=="RGE_BRANCH_MISMATCH_CONTROL" and not c["eligible_for_framework_cross_verdict"]
def test_b_test_scale_conflict_is_preserved():
 assert receipt()["k7_candidates"]["b_test_holonomy_bundle"]["scale"].startswith("CONFLICTED:")
def test_type_three_output_is_recorded_but_not_promoted():
 r=receipt()["k7_candidates"]["rge_bundle"]
 assert r["found"] and r["displayed_outputs"]["alpha_em_inverse_MZ"]=="131.19"
 assert r["exact_vector"].startswith("NOT AVAILABLE") and not r["same_object_contract"]
def test_normalization_covariance_recorded(): assert receipt()["controls"]["normalization_covariance"]["zero_status_preserved"]
def test_no_admissible_vector_is_not_an_off_plane_claim():
 d=receipt(); assert d["comparison"]["admissible_candidate"] is None and not d["comparison"]["same_object_contract"]

# Hostile controls alter pinned premises, exact values, types and trust
# boundaries. These are selected semantic attacks, not a prose-completeness proof.
def test_hostile_alpha_zero_relabelled_mz(): rejected(lambda d:d["k7_candidates"]["type1_structural_bundle"]["alpha_em_inverse"].update(scale="M_Z"))
def test_hostile_mixed_scale_promotion(): rejected(lambda d:d["controls"]["mixed_scale_probe"].update(eligible_for_framework_cross_verdict=True))
def test_hostile_untransformed_gut_u1(): rejected(lambda d:d["oph_contract"].update(u1_normalization="GUT inserted into Y row"))
def test_hostile_transform_x_only(): rejected(lambda d:d["controls"]["normalization_covariance"].update(all_x_k_b_first_entries_transformed=False))
def test_hostile_transform_k_only(): rejected(lambda d:d["controls"]["normalization_covariance"].update(hostile_partial_row_transform_rejected=False))
def test_hostile_swap_y_and_weak(): rejected(lambda d:d["oph_contract"].update(vector_role="swapped alpha_2 and alpha_Y"))
def test_hostile_alpha_s_instead_inverse(): rejected(lambda d:d["controls"]["mixed_scale_probe"].update(xY=["a","b","sqrt(2)/12"]))
def test_hostile_float_sqrt2(): rejected(lambda d:d["controls"]["mixed_scale_probe"]["exact_D"].update(sqrt2_part="1849.??"))
def test_hostile_beta_mutation(): rejected(lambda d:d["oph_contract"]["beta_column"].__setitem__(0,"7"))
def test_hostile_kinetic_mutation(): rejected(lambda d:d["oph_contract"]["kinetic_column"].__setitem__(0,"4"))
def test_hostile_ng_nh_mutation(): rejected(lambda d:d["oph_contract"].update(nG=2))
def test_hostile_mssm_beta_substitution(): rejected(lambda d:d["oph_contract"].update(beta_column=["33/5","1","-3"]))
def test_hostile_btest_branch_promotion(): rejected(lambda d:d["controls"]["b_test_branch_mismatch_probe"].update(eligible_for_framework_cross_verdict=True))
def test_hostile_public_measurement_input(): rejected(lambda d:d["trust_boundary"].update(public_measurement_consumed_as_cross_input=True))
def test_hostile_sealed_column_read(): rejected(lambda d:d["trust_boundary"].update(oph_sealed_comparison_opened=True))
def test_hostile_mixed_verdict_promotion(): rejected(lambda d:d.update(primary_verdict="K7_GAUGE_VECTOR_EXACTLY_OFF_OPH_MATTER_TRACE_RG_PLANE"))
def test_hostile_exact_near_zero_claim(): rejected(lambda d:d["controls"]["mixed_scale_probe"]["exact_D"].update(rational_part="0",sqrt2_part="0"))
def test_hostile_external_bytes_mutation_and_rehashed_json(tmp_path):
 mini=tmp_path/"repo"
 for rel in ver.OPH_PINS:
  dest=mini/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,dest)
 d=receipt()
 for s in d["pins"]["k7_sources"]:
  dest=mini/s["vendored_path"];dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/s["vendored_path"],dest)
 target=mini/d["pins"]["k7_sources"][-1]["vendored_path"]
 target.write_bytes(target.read_bytes()+b"\\n-- mutation")
 d["pins"]["k7_sources"][-1]["sha256"]=hashlib.sha256(target.read_bytes()).hexdigest()
 try: ver.validate(d,mini)
 except ValueError: return
 raise AssertionError("mutated source with updated JSON digest passed pinned blob check")
def test_hostile_external_sha_rehash_does_not_change_commit_pin(): rejected(lambda d:d["pins"].update(k7_commit="f"*40))
def test_hostile_erratum_omitted(): rejected(lambda d:d["pins"]["k7_sources"].pop(0))
def test_hostile_lean_arithmetic_does_not_fix_scale(): rejected(lambda d:d["k7_candidates"]["type1_structural_bundle"]["alpha_em_inverse"].update(scale="M_Z"))
def test_hostile_port_branch_identified_with_k7(): rejected(lambda d:d["trust_boundary"].update(oph_port_response_identified_with_k7=True))
def test_hostile_matter_branch_selected(): rejected(lambda d:d["oph_contract"].update(matter_branch_selected_physical=True))
def test_hostile_off_plane_global_falsification(): rejected(lambda d:d["trust_boundary"].update(off_plane_does_not_falsify_oph_globally=False))

def replace_field(doc, path, value):
 for key in path[:-1]: doc=doc[key]
 doc[path[-1]]=value

@pytest.mark.parametrize("path,value", [
 (("controls","mixed_scale_probe","zero_status"),"EXACT_ZERO"),
 (("controls","b_test_branch_mismatch_probe","zero_status"),"EXACT_ZERO"),
 (("controls","b_test_branch_mismatch_probe","xY"),["1","2","3"]),
 (("comparison","exact_determinant"),{"rational_part":"1","sqrt2_part":"0"}),
 (("comparison","exact_zero_status"),"EXACT_NONZERO"),
 (("oph_contract","integer_zero_locus"),[69.5,-333.5,218.5]),
 (("oph_contract","integer_zero_locus"),["69","-333","218"]),
 (("oph_contract","nH"),True),
 (("oph_contract","nG"),3.0),
 (("oph_contract","scheme"),"MS-bar"),
 (("oph_contract","thresholds"),"physical threshold derived"),
 (("controls","normalization_covariance","raw_determinant_row_factor"),"5/3"),
 (("controls","normalization_covariance","Y_to_GUT_replayed"),False),
 (("controls","normalization_covariance","zero_test_input","xY"),["0","0","1"]),
 (("k7_candidates","b_test_holonomy_bundle","ratio_GUT_normalized"),"14:7:3"),
 (("k7_candidates","b_test_holonomy_bundle","normalization"),"unconverted GUT"),
 (("k7_candidates","b_test_holonomy_bundle","field_content"),"MSSM-free SM branch"),
 (("k7_candidates","b_test_holonomy_bundle","scale"),"CONFLICTED: actually all at M_Z"),
 (("k7_candidates","b_test_holonomy_bundle","rge_branch"),"OPH one-loop SM"),
 (("k7_candidates","b_test_holonomy_bundle","physical_role"),"measured comparison column"),
 (("k7_candidates","type1_structural_bundle","sin2_theta_W","scale"),"alpha_em^-1(0)"),
 (("k7_candidates","type1_structural_bundle","alpha_s","scheme"),"MS-bar"),
 (("k7_candidates","type1_structural_bundle","alpha_s","status"),"target-blind prediction"),
 (("k7_candidates","rge_bundle","displayed_outputs","alpha_s_all_MSSM"),"0.1224"),
 (("k7_candidates","rge_bundle","exact_vector"),"NOT AVAILABLE except this exact vector"),
 (("k7_candidates","rge_bundle","public_data_used_in_construction"),False),
 (("k7_candidates","rge_bundle","direction_of_running"),"M_Z to M_GUT"),
 (("k7_candidates","rge_bundle","field_content"),"OPH SM"),
])
def test_hostile_contradictory_control_and_source_typing(path,value):
 rejected(lambda d:replace_field(d,path,value))

@pytest.mark.parametrize("field", ["alpha_1_GUT_inverse","alpha_2_inverse","alpha_3_inverse","alpha_Y_inverse","alpha_2_inverse_Y","alpha_3_inverse_Y"])
def test_hostile_each_duplicate_btest_component(field):
 rejected(lambda d:d["k7_candidates"]["b_test_holonomy_bundle"]["values"].update({field:"22*sqrt(2)"}))

@pytest.mark.parametrize("path", [
 ("k7_candidates","type1_structural_bundle",field,"public_data_used_in_construction")
 for field in ("alpha_em_inverse","sin2_theta_W","alpha_s")
]+[("k7_candidates","b_test_holonomy_bundle","public_data_used_in_construction")])
@pytest.mark.parametrize("unsupported_claim", [False,True,"target-blind"])
def test_historical_input_status_cannot_be_overclaimed(path,unsupported_claim):
 rejected(lambda d:replace_field(d,path,unsupported_claim))

def test_historical_unknown_is_separate_from_cross_input_boundary():
 d=receipt()
 assert all(d["k7_candidates"]["type1_structural_bundle"][field]["public_data_used_in_construction"].startswith("UNKNOWN;") for field in ("alpha_em_inverse","sin2_theta_W","alpha_s"))
 assert d["k7_candidates"]["b_test_holonomy_bundle"]["public_data_used_in_construction"].startswith("UNKNOWN;")
 assert d["trust_boundary"]["public_measurement_consumed_as_cross_input"] is False
 assert d["trust_boundary"]["public_measurement_values_parsed_by_cross_control"] is False

@pytest.mark.parametrize("field", ["oph_repository","k7_repository","k7_lean_repository","oph_commit","k7_commit","k7_fetched_main_head","k7_lean_commit","k7_lean_fetched_main_head"])
def test_hostile_each_top_level_custody_field(field):
 rejected(lambda d:d["pins"].update({field:"invented/provenance"}))

@pytest.mark.parametrize("field,value", [("repository","invented/repository"),("commit","0"*40),("path","elsewhere.md"),("git_blob","0"*40),("vendored_path","elsewhere.md"),("sha256","0"*64),("bytes",9173.0)])
def test_hostile_each_external_custody_field(field,value):
 rejected(lambda d:d["pins"]["k7_sources"][0].update({field:value}))

@pytest.mark.parametrize("path", [("pins",),("pins","k7_sources",0),("comparison",)])
def test_hostile_extra_custody_or_comparison_claim(path):
 rejected(lambda d:replace_field(d,path+("unverified_claim",),"promoted"))

@pytest.mark.parametrize("text", [
 '{"comparison":{"same_object_contract":true,"same_object_contract":false}}',
 '{"zero_status":"EXACT_ZERO","zero_status":"EXACT_NONZERO"}',
 '{"value":NaN}', '{"value":Infinity}', '{"value":-Infinity}',
])
def test_ambiguous_or_nonfinite_json_rejected(text):
 with pytest.raises(ValueError): ver.load_receipt(text)

def test_strict_json_load_accepts_committed_receipt():
 assert ver.validate(ver.load_receipt(RECEIPT.read_text()),ROOT)

def leaf_paths(value,path=()):
 if isinstance(value,dict):
  for key,item in value.items(): yield from leaf_paths(item,path+(key,))
 elif isinstance(value,list):
  for index,item in enumerate(value): yield from leaf_paths(item,path+(index,))
 else: yield path

@pytest.mark.parametrize("path", list(leaf_paths(receipt())),ids=lambda p:"/".join(map(str,p)))
def test_each_published_leaf_rejects_contradictory_replacement(path):
 rejected(lambda d:replace_field(d,path,"HOSTILE_CONTRADICTORY_CLAIM"))

def object_paths(value,path=()):
 if isinstance(value,dict):
  yield path
  for key,item in value.items(): yield from object_paths(item,path+(key,))
 elif isinstance(value,list):
  for index,item in enumerate(value): yield from object_paths(item,path+(index,))

@pytest.mark.parametrize("path",list(object_paths(receipt())),ids=lambda p:"/".join(map(str,p)) or "root")
def test_each_object_rejects_added_unreviewed_claim(path):
 rejected(lambda d:replace_field(d,path+("unreviewed_claim",),"physical agreement established"))

@pytest.mark.parametrize("field,value",[("kinetic_column",["10/3",2,2]),("beta_column",["41/6","-19/6",-7.0]),("cofactors",["-23/3",37.0,"-218/9"])])
def test_exact_arithmetic_fields_reject_float_or_integer_coercion(field,value):
 rejected(lambda d:d["oph_contract"].update({field:value}))

@pytest.mark.parametrize("field",["kinetic_column","beta_column","cofactors","integer_zero_locus"])
def test_exact_vectors_cannot_be_replaced_by_mapping_keys(field):
 rejected(lambda d:d["oph_contract"].update({field:{value:None for value in d["oph_contract"][field]}}))

def test_arithmetic_is_checked_beyond_semantic_custody():
 d=receipt(); d["controls"]["mixed_scale_probe"]["exact_D"]["sqrt2_part"]="0"
 assert ver.semantic_digest(d)==ver.REVIEWED_SEMANTICS_SHA256
 with pytest.raises(ValueError,match="exact value changed"): ver.validate(d,ROOT)
