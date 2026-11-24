"""Check a prospective OW-30 registration and execute only its first public release."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
PUBLIC = 'https://github.com/openwave-labs/openwave.git'
BASE = 'd4fc3e4046c472f4ad38468df8dddfa1f97bac97'
PIN = '4fe10f9ac3e81423155fb3ac7baaf0fdb4c0b11c'
TASK = 'openwave/xperiments/m10_oph/research/tasks/m10_5_task_details.md'
RER_PUBLIC = 'https://github.com/FloatingPragma/observer-patch-holography.git'
SOURCE_POINTER = TASK.removesuffix('.md')+'.source'
CANON = 'openwave/xperiments/m10_oph/research/m10_theory_canonical.md'
PREFIX = b'OPH-OW30-collective-confirmation-v1:'
PILOTS = {312141,312142,312143,927016}


def git(repo, *args):
    return subprocess.run(['git','-C',str(repo),*args],check=True,capture_output=True).stdout


def file_hash(path):
    with path.open('rb') as stream:
        digest=sha256()
        for chunk in iter(lambda:stream.read(1024*1024),b''):
            digest.update(chunk)
        return digest.hexdigest()


def frozen_sources():
    manifest=HERE/'files.json'
    task=(HERE/'TASK.md').read_bytes()
    declared=re.findall(rb'Frozen source manifest SHA-256 `([0-9a-f]{64})`',task)
    if declared!=[file_hash(manifest).encode()]:
        raise ValueError('Task does not bind this exact source manifest')
    sections=json.loads(manifest.read_text())
    for path, expected in sections['working_inputs'].items():
        if file_hash(ROOT/path)!=expected:
            raise ValueError('Changed local source: '+path)
    for path, expected in sections['adopted_native_inputs'].items():
        content=git(ROOT,'show',PIN+':'+path)
        if sha256(content).hexdigest()!=expected or file_hash(ROOT/path)!=expected:
            raise ValueError('Native generator/board does not match adopted bytes: '+path)
    return task


def first_exact(repo, base, head, path, terms):
    """Identify the first first-parent landing, including replacement or restoration."""
    lineage=git(repo,'rev-list','--first-parent',head).decode().splitlines()
    if base not in lineage:
        raise ValueError('Complete first-parent chain through the recorded baseline is required')
    def exact(at):
        return bool(git(repo,'ls-tree',at,'--',path)) and git(repo,'show',at+':'+path)==terms
    if exact(base):
        raise ValueError('Task was already present at the baseline')
    changes=git(repo,'log','--first-parent','--reverse','--format=%H',base+'..'+head,'--',path).decode().splitlines()
    for rev in changes:
        if exact(rev):
            if exact(rev+'^1'):
                raise ValueError('Registration predecessor already carried these bytes')
            return rev
    return None


def official_anchor(checkout):
    terms=frozen_sources()
    ref='refs/remotes/origin/main'
    # Read-only network check pins the local history to the official default branch.
    public=git(checkout,'ls-remote',PUBLIC,'refs/heads/main').decode().split()[0]
    if public!=git(checkout,'rev-parse',ref).decode().strip():
        raise ValueError('Fetch official OpenWave main into origin/main, then retry')
    if git(checkout,'rev-parse','--is-shallow-repository').strip()==b'true':
        raise ValueError('A complete first-parent history is required')
    anchor=first_exact(checkout,BASE,public,TASK,terms)
    if anchor:
        changes=git(checkout,'log','--first-parent','--format=%H',anchor+'..'+public,'--',TASK,SOURCE_POINTER)
        if changes.strip():
            raise ValueError('Task or source pointer changed after registration; prior release remains historical')
        if not git(checkout,'ls-tree',public,'--',TASK) or git(checkout,'show',public+':'+TASK)!=terms:
            raise ValueError('Registration has been withdrawn or amended')
        public_sources(checkout,anchor,public)
        for rev in (anchor,public):
            canonical=git(checkout,'show',rev+':'+CANON).decode()
            pinned=next((line for line in canonical.splitlines() if line.startswith('| Pinned candidate')), '')
            if PIN not in pinned:
                raise ValueError('Upstream candidate pin changed; stop for candidate adjudication')
    return anchor


def public_sources(checkout,anchor,public):
    """The upstream merge binds the public RER release in a separate Git blob."""
    at_landing=git(checkout,'show',anchor+':'+SOURCE_POINTER)
    current=git(checkout,'show',public+':'+SOURCE_POINTER)
    if at_landing!=current or not re.fullmatch(rb'[0-9a-f]{40}\n',at_landing):
        raise ValueError('Missing, invalid or changed public source pointer')
    revision=at_landing.decode().strip()
    remote_main=git(ROOT,'ls-remote',RER_PUBLIC,'refs/heads/main').decode().split()[0]
    if remote_main!=git(ROOT,'rev-parse','refs/remotes/origin/main').decode().strip():
        raise ValueError('Fetch official RER main before verifying upstream source custody')
    git(ROOT,'merge-base','--is-ancestor',revision,remote_main)
    manifest=HERE/'files.json'
    for path in (manifest,HERE/'TASK.md'):
        released=git(ROOT,'show',revision+':'+path.relative_to(ROOT).as_posix())
        if released!=path.read_bytes():
            raise ValueError('Public RER registration bytes differ: '+path.name)
    sections=json.loads(manifest.read_text())
    for group in ('working_inputs','adopted_native_inputs'):
        for path,expected in sections[group].items():
            if sha256(git(ROOT,'show',revision+':'+path)).hexdigest()!=expected:
                raise ValueError('Public release source differs: '+path)
    return revision


def settings(anchor):
    if not re.fullmatch('[0-9a-f]{40}',anchor):
        raise ValueError('Invalid registration commit')
    seeds=[]
    for lev in range(3):
        digest=sha256(PREFIX+anchor.encode('ascii')+b':'+str(lev).encode('ascii')).digest()
        seed=(1<<63)|int.from_bytes(digest[:8],'big')
        if seed in seeds or seed in PILOTS:
            raise ValueError('Selector collision: retain failure, no replacement seeds')
        seeds.append(seed)
    return seeds


def clean_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def run_bank(anchor, dest, source_commit=None):
    sys.path.insert(0,str(ROOT))
    import numpy as np
    from importlib.metadata import version
    from unittest.mock import patch
    from model.candidate_laws.interval_echo.openwave.row_30_collective import experiment as ex
    from model.candidate_laws.interval_echo.openwave.row_30_collective import formation, holdout, coupling_meter
    from model.executable_model import federation

    requirements=(ROOT/'model/requirements.txt').read_text().splitlines()
    for name in ('numpy','scipy'):
        expected=next(line.split('==')[1] for line in requirements if line.startswith(name+'=='))
        if version(name)!=expected:
            raise ValueError('Registered acquisition requires '+name+'=='+expected)
    if dest.exists():
        raise ValueError('Output directory already exists; use replay for an acquired release')
    dest.mkdir(parents=True)
    seeds=settings(anchor)
    clean_json(dest/'draw.json',dict(registration=anchor,source_commit=source_commit,adopted_pin=PIN,protocol_sha256=file_hash(HERE/'TASK.md'),seeds=seeds))
    meter=json.loads((ex.HERE/'coupling_meter.json').read_text())
    values,_=coupling_meter.measure(np.array(meter['all_counts']),
                                   np.array(meter['hand_calibration'])/meter['calibration_trials'])
    h,j,m,c2=map(float,values[:4])
    meter_ok=bool(h>0 and j>h and abs(values[6])<.005 and abs(values[7])<.005)
    if ex.K!=8 or ex.candidate.native_law.__name__!='model.candidate_laws.contrast_memory_spherical_v1.law':
        raise ValueError('Unexpected effective generator binding or pinned coupling')
    checks=ex.native_checks()
    checks['commissioning_column_residual']=formation.native()
    if checks['commissioning_column_residual']>1e-11:
        raise ValueError('Native commissioning columns changed')
    dense=[ex.exact_check(n) for n in (5,6)]
    boards=[]
    for lev,seed in enumerate(seeds):
        folder=dest/str(lev)
        folder.mkdir()
        try:
            board=federation(lev)
            setup=formation.prepare(lev,seed)
            route=max((p for closed,p in ex.grow(board,seed) if not closed),key=len)
            def keep(stream,key,clicks):
                data=np.asarray(clicks,dtype='<u8',order='C')
                np.save(folder/(key+'.npy'),data,allow_pickle=False)
                stream.hashes[key]=sha256(data.tobytes()).hexdigest()
            with patch.object(ex.CountStream,'keep',keep):
                obs=ex.spectroscopy(len(route),np.random.default_rng([seed,lev,47]))
            native,off=[ex.transport(route,coupling,seed) for coupling in (ex.K,0)]
            # Preserve both positive outcome counts behind the published count differences.
            for tag,receipt,coupling in [('native',native,ex.K),('control',off,0)]:
                rng=np.random.default_rng(seed)
                n=len(route)
                for flight in receipt['flights']:
                    kernel=ex.walls(n,coupling,flight['t'])
                    arms=[]
                    for site in (receipt['source'],receipt['reference']):
                        initial=np.zeros(n+1)
                        initial[-1]=1-n%2
                        initial[site]=1
                        p=.5+kernel@(initial-.5)
                        arms.append(rng.binomial(receipt['trials_per_arm'],np.clip(p,0,1)))
                    np.testing.assert_array_equal(arms[0]-arms[1],flight['differences'])
                    np.save(folder/(tag+'_flight_'+str(flight['t'])+'.npy'),np.asarray(arms,dtype='<u8'))
            score=holdout.score(obs['modes'],obs,native,off,m,c2)
            ordinary_modes=[r for r in obs['modes'] if 0<r['k']<=.50]
            ordinary=[abs(r['omega']-(m*m+c2*r['k']**2)**.5) for r in ordinary_modes]
            continuum_errors=[r['omega']-(m*m+c2*r['k']**2)**.5 for r in obs['modes']]
            cells=len(set(np.ravel(route)//12))
            controls=obs['same_spatial_source_changed_contrast'].values()
            min_changed_at_two=min(x['changed_at_two'] for x in controls)
            gates=dict(changed_reader=min_changed_at_two>.40,meter=meter_ok,commissioned=bool(setup['dimers_enrolled']==len(route) and cells>1),
                       lattice=all(r['omega']>0 for r in obs['modes']) and score['lattice_max']<.0048,closure=score['sum_closure']<.0048,
                       live=score['min_native_amplitude']>.40,changed=score['max_changed_amplitude']<.02,
                       transport=score['far_count_difference']>.02,quiet=score['far_changed_count_difference']<.005,
                       shape=score['shape_rms_lattice']<.5*score['shape_rms_parabola'],
                       continuum=len(ordinary)>=2 and np.ptp([r['k'] for r in ordinary_modes])>=.02 and max(ordinary)<.0048)
            boards.append(dict(level=lev,seed=seed,route=route,setup=setup,observer=obs,transport=native,
                               changed_transport=off,score=score,low_k_errors=ordinary,continuum_residuals=continuum_errors,min_changed_at_two=min_changed_at_two,gates={k:bool(v) for k,v in gates.items()}))
        except Exception as err:
            failed=dict(level=lev,seed=seed,error=type(err).__name__,detail=str(err),gates={'records_complete':False})
            clean_json(folder/'failure.json',failed)
            boards.append(failed)
    passed=all(all(row['gates'].values()) for row in boards)
    clean_json(dest/'result.json',dict(registration=anchor,m=m,c_squared=c2,native=checks,
                                      finite_spectrum_checks=dense,boards=boards,passed=passed))
    return passed


def compare_release(a,b):
    import numpy as np
    files=lambda root:{p.relative_to(root) for p in root.rglob('*') if p.is_file()}
    if files(a)!=files(b):
        raise ValueError('Release file inventory differs')
    def compare(x,y):
        if isinstance(x,dict):
            if not isinstance(y,dict) or x.keys()!=y.keys():
                raise ValueError('Missing result keys')
            for k in x: compare(x[k],y[k])
        elif isinstance(x,list):
            if not isinstance(y,list) or len(x)!=len(y): raise ValueError('Result length differs')
            for i,j in zip(x,y): compare(i,j)
        elif isinstance(x,float):
            if not np.isfinite(y) or abs(x-y)>1e-8: raise ValueError('Numerical replay differs')
        elif x!=y:
            raise ValueError('Discrete result differs')
    for file in files(b):
        if file.suffix=='.npy':
            np.testing.assert_array_equal(np.load(a/file,allow_pickle=False),np.load(b/file,allow_pickle=False))
        else:
            compare(json.loads((b/file).read_text()),json.loads((a/file).read_text()))


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('command',choices=['check-sources','status','acquire','replay'])
    cli.add_argument('--openwave',type=Path)
    cli.add_argument('--release',type=Path)
    args=cli.parse_args()
    if args.command=='check-sources':
        frozen_sources()
        print('Verified reviewed inputs and native files at the adopted M10.3 pin.')
        return
    if not args.openwave: cli.error('--openwave requires an official Git checkout')
    anchor=official_anchor(args.openwave)
    if args.command=='status':
        print('First exact public registration: '+(anchor or 'absent. No settings derived.'))
        return
    if not anchor: raise SystemExit('No exact public registration. No settings derived.')
    if not args.release: cli.error('--release is required')
    source=git(args.openwave,'show',anchor+':'+SOURCE_POINTER).decode().strip()
    if args.command=='acquire':
        passed=run_bank(anchor,args.release,source)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            other=Path(tmp)/'run'
            passed=run_bank(anchor,other,source)
            compare_release(args.release,other)
    print('All prospective gates passed.' if passed else 'Prospective release rejected.')
    if not passed: raise SystemExit(2)


if __name__=='__main__':
    main()
