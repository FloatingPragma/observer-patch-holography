"""Current replay routing with one explicitly enumerated successor.

The original verifier, artifact/source guards, group partition and every other
replay declaration retain their ordering. The protected particle readout is
checked by its current numerical successor.
"""
import argparse,gzip,hashlib,json,os,subprocess,sys
from pathlib import Path
from . import verify as retained
from tools import strict_json
try:
    from tools.pin_scope import bound_pins, exempt, scoped
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1]))
    from tools.pin_scope import bound_pins, exempt, scoped

BASE_REPLAYS_SHA256='87c475d86192b569fe969485ee0bab67bf16d7bfe3251713b9ffb62a3ee1fdfa'
LEGACY=('model.protected_particle_readout',None,'model/protected_particle_readout/results.json')
CURRENT=('model.protected_readout_numerics',None,'model/protected_readout_numerics/results.json')
ROUTES={LEGACY:CURRENT}


def digest(value):return hashlib.sha256(value).hexdigest()

def replay_digest(rows):return digest(json.dumps(tuple(rows),separators=(',',':')).encode())

def current_replays(rows=None):
    rows=retained.REPLAYS if rows is None else tuple(rows)
    if (replay_digest(rows)!=BASE_REPLAYS_SHA256
            or any(rows.count(old)!=1 or new in rows for old,new in ROUTES.items())):
        raise ValueError('unreviewed replay inventory/order or duplicate successor')
    return tuple(ROUTES.get(row,row) for row in rows)


def verify_successor_source_bindings():
    from .protected_readout_numerics import model
    from .protected_particle_readout import packet
    model.specification()
    value=strict_json.load(retained.ROOT/CURRENT[2]);pins=dict(packet.source_bindings())
    owned={p.name:model.digest(p)for p in sorted(model.HERE.iterdir())if p.suffix in('.py','.md')}
    pins.update({str((model.HERE/name).relative_to(retained.ROOT)):sha for name,sha in owned.items()})
    for name in('contract.json','retained.json'):pins[str((model.HERE/name).relative_to(retained.ROOT))]=model.digest(model.HERE/name)
    spec,_,_=model.specification();pins[spec['retained_array']]=spec['retained_array_sha256']
    # Prose and test entries (tools/pin_scope.py) are neither pinned nor compared.
    pins,owned=bound_pins(pins),bound_pins(owned)
    if bound_pins(value.get('input_sha256'))!=scoped(pins) or bound_pins(value.get('owned_source_sha256'))!=scoped(owned):
        raise ValueError('numerical successor source scope or custody changed')


# The retained verifier's artifact and source-map checks, without prose and
# test pins (tools/pin_scope.py), code or JSON receipts in package manifests.
def verify_receipt_source_bindings(digests=None):
    """Check live replay source maps without recalculating their science.

    Recognized maps use their producers' repository or package-relative
    path convention. Local-only maps reject missing or escaping files.
    The mixed ``sources`` maps also contain external custody metadata and
    aliases; those and frozen captures outside REPLAYS retain their owning
    verifiers. This does not replace schema validation or scientific replay.
    Prose and test entries (tools/pin_scope.py) are never compared.
    """
    digests = {} if digests is None else digests
    root = retained.ROOT.resolve()
    repository_fields = {"sources", "source_sha256", "source_bindings",
                         "input_sha256", "source_binding_sha256",
                         "source_control_bindings", "source_pins"}
    package_levels = {"implementation_sha256": 0, "source_interface_hashes": 1}
    # The endpoint audit's generic field name is scoped to its actual location.
    endpoint_pins = "root.retrospective_endpoint_audit.source.pins"

    def walk(value, receipt, location="root"):
        if isinstance(value, dict):
            for key, child in value.items():
                here = f"{location}.{key}"
                recognized = key in repository_fields or key in package_levels or here == endpoint_pins
                local_only = recognized and key != "sources"
                if local_only and (not isinstance(child, dict) or not child):
                    raise ValueError(f"Empty or malformed local source map: {receipt} at {here}")
                if recognized and isinstance(child, dict):
                    base = ((root/receipt).parents[package_levels[key]]
                            if key in package_levels else root)
                    for relative, expected in child.items():
                        if isinstance(relative, str) and exempt(relative):
                            continue
                        if not isinstance(relative, str) or "://" in relative:
                            if local_only:
                                raise ValueError(f"Invalid local source path: {receipt} at {here}")
                            continue
                        path = Path(relative)
                        if path.is_absolute():
                            if local_only:
                                raise ValueError(f"Absolute local source path: {receipt} at {here}[{relative}]")
                            continue
                        try:
                            target = (base/path).resolve()
                            local_file = target.is_relative_to(root) and target.is_file()
                        except (OSError, ValueError):
                            local_file = False
                        if not local_file:
                            if local_only:
                                raise ValueError(f"Missing or escaping local source path: {receipt} at {here}[{relative}]")
                            continue
                        if target not in digests:
                            digests[target] = hashlib.sha256(target.read_bytes()).hexdigest()
                        if digests[target] != expected:
                            raise ValueError(
                                f"Receipt source digest mismatch: {receipt} at {here}[{relative}]")
                walk(child, receipt, here)
        elif isinstance(value, list):
            for index, child in enumerate(value):
                walk(child, receipt, f"{location}[{index}]")

    for _, _, receipt in retained.REPLAYS:
        path = retained.ROOT/receipt
        raw = path.read_bytes()
        if path.suffix == ".gz":
            raw = gzip.decompress(raw)
        walk(strict_json.loads(raw.decode("utf-8")), receipt)


# Non-artifact files that legitimately sit inside a packet: the manifest that
# binds it, repository/VCS metadata, and interpreter caches. Every other file
# must be bound by the packet manifest, whatever its suffix, except prose and
# test modules (tools/pin_scope.py), code, and JSON receipts. Code is not
# byte-pinned, and every receipt carries its own source map and replay.
UNPINNED_SUFFIXES = (".py", ".c", ".json", ".json.gz")
PACKET_METADATA = frozenset({"manifest.json", ".gitattributes", ".gitignore"})
PACKET_CACHE_DIRS = frozenset({"__pycache__", ".pytest_cache"})


def verify_artifacts():
    covered = set()
    digests = {}
    for name in retained.PACKETS:
        folder = retained.ROOT / name
        manifest = json.loads((folder / "manifest.json").read_text())
        files = manifest["files"]
        if not files:
            raise ValueError(f"Empty artifact manifest: {name}")
        for relative, expected in files.items():
            if exempt(f"{name}/{relative}") or relative.endswith(UNPINNED_SUFFIXES):
                continue
            target = (folder / relative).resolve()
            if not target.is_relative_to(folder.resolve()):
                raise ValueError(f"Artifact outside its package: {name}/{relative}")
            actual = hashlib.sha256(target.read_bytes()).hexdigest()
            if actual != expected:
                raise ValueError(f"Artifact digest mismatch: {name}/{relative}")
            covered.add(target)
            digests[target] = actual
    for name in retained.PACKETS:
        for target in (retained.ROOT / name).rglob("*"):
            if not target.is_file() or target.name in PACKET_METADATA:
                continue
            if exempt(target.relative_to(retained.ROOT).as_posix()) or target.name.endswith(UNPINNED_SUFFIXES):
                continue
            if PACKET_CACHE_DIRS.intersection(target.relative_to(retained.ROOT / name).parts):
                continue
            if target.resolve() not in covered:
                raise ValueError(f"Unbound artifact: {target.relative_to(retained.ROOT)}")
    verify_receipt_source_bindings(digests)
    return len(covered)



def verify_all_hashes():
    current_replays()
    count=verify_artifacts()
    verify_successor_source_bindings()
    return count


def selected_groups(index=None,count=None):
    return retained.selected_replay_groups(index,count,replays=current_replays())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hashes-only',action='store_true')
    p.add_argument('--shard-index',type=int);p.add_argument('--shard-count',type=int)
    a=p.parse_args()
    try:groups=selected_groups(a.shard_index,a.shard_count)
    except ValueError as error:p.error(str(error))
    count=verify_all_hashes();print(f'Verified {count} retained artifacts and the numerical successor source',flush=True)
    if a.hashes_only:return
    environment=dict(os.environ)
    for name in('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):environment.setdefault(name,'1')
    for group in groups:
        if group==retained.EMPIRICAL_REPLAYS:
            subprocess.run([sys.executable,'-m','model.replay_empirical'],cwd=retained.ROOT,env=environment,check=True)
            continue
        (module,command,receipt),=group
        argv=[sys.executable,'-m',module]
        if command is not None:argv.append(command)
        argv.extend(['--check',receipt]);print(f'Replaying {receipt}',flush=True)
        subprocess.run(argv,cwd=retained.ROOT,env=environment,check=True)
    print(f'All {sum(map(len,groups))} selected current receipts reproduced',flush=True)


if __name__=='__main__':main()
