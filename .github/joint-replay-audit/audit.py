from pathlib import Path
import hashlib,json,re,subprocess,sys,math
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUT=ROOT/"audit-artifacts";OUT.mkdir(exist_ok=True)
manifest=json.loads((HERE/"manifest.json").read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for name in ("input_sha256","snapshot_sha256"):
    for path,digest in manifest[name].items():
        assert sha(ROOT/path)==digest,(name,path,"source bytes changed")
report={"phase":sys.argv[1],"git_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),"input_files":len(manifest["input_sha256"]),"snapshot_files":len(manifest["snapshot_sha256"]),"all_inputs_exact":True}
if sys.argv[1] != "before":
    changes=[];float_changes=[]
    def same(a,b,path=""):
        assert type(a) is type(b),("type",path)
        if isinstance(a,dict):
            assert a.keys()==b.keys(),("keys",path)
            for k in a:same(a[k],b[k],path+"/"+k)
        elif isinstance(a,list):
            assert len(a)==len(b),("length",path)
            for i,(x,y) in enumerate(zip(a,b)):same(x,y,path+"/"+str(i))
        elif a!=b:
            if isinstance(a,float):
                assert math.isfinite(a) and math.isfinite(b),("nonfinite",path)
                float_changes.append({"path":path,"before":a,"after":b,"absolute_delta":abs(a-b),"relative_delta":abs(a-b)/max(abs(a),abs(b),sys.float_info.min)})
                return
            assert isinstance(a,str) and re.fullmatch("[0-9a-f]{64}",a) and re.fullmatch("[0-9a-f]{64}",b),("nonhash scientific difference",path,a,b)
            changes.append({"path":path,"before":a,"after":b})
    fresh=json.loads((ROOT/manifest["target_output"]).read_text())
    same(json.loads((HERE/"before.json").read_text()),fresh)
    allowed={"/input_sha256/extra/yang_mills_gap_clay_problem.tex","/input_sha256/model/particle_joint_reconstruction/candidate_results.json"}
    assert {x["path"] for x in changes}==allowed,("unexpected exact changed leaves",changes)
    n=[0]
    def pins(x):
        if isinstance(x,dict):
            for k,v in x.items():
                if k in ("source_sha256","input_sha256","source_bindings") and isinstance(v,dict):
                    for rel,digest in v.items():
                        if isinstance(digest,str) and re.fullmatch("[0-9a-f]{64}",digest):
                            assert sha(ROOT/rel)==digest,("stale produced source pin",rel)
                            n[0]+=1
                pins(v)
        elif isinstance(x,list):
            for v in x:pins(v)
    pins(fresh)
    report.update({"changed_hash_leaves":changes,"all_other_keys_types_values_exact":not float_changes,"all_nonfloat_keys_types_values_exact":True,"float_differences":float_changes,"supplemental_float_acceptance":"PENDING_REVIEW" if float_changes else "EXACT","live_produced_bindings":n[0],"result_sha256":sha(ROOT/manifest["target_output"])})
(OUT/(sys.argv[1]+"-source-report.json")).write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
print(json.dumps(report,sort_keys=True))
