#!/usr/bin/env python3
"""Postrun dependency custody. This does not retroactively extend PREREG."""
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import sys
import numpy
import scipy
import headline as h

if not (h.HERE/'source_results.json').exists():
    raise SystemExit('Collect final custody only after the complete source replay.')
runtime=set()
for module in tuple(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name:
        p=Path(name).resolve()
        if p.is_file() and p.is_relative_to(h.RER):runtime.add(p)
runtime|={h.ROOT/p for p in h.dependencies() if p.startswith('reverse-engineering-reality/')}

base=h.RER/'code/a5_closure'
sys.path.insert(0,str(base))
import sl2f5_port_spin_bridge_certificate
cited={h.RER/p for p in (
    'model/observer_particles/spin_matter_join.py',
    'model/defect_field_bridge/REPRESENTATION.md',
    'model/defect_field_bridge/COHERENT_EXTENSION.md',
    'model/native_spin_observer/README.md',
    'model/native_observer_program/README.md',
    'model/native_observer_program/PROOF.md',
    'model/native_observer_program/source.py',
    'model/native_electric_gauss/PROOF.md',
    'model/native_electric_gauss/source.py',
    'model/observer/README.md',
    'model/observer/native_loop_readout.py',
    'model/observer/holonomy_reconstruction.py',
    'model/native_circulation/README.md',
    'model/particle_common_contract/PROOF.md',
    'Lean/EventAlgebra/PauliMinimalCoupling.lean',
    'Lean/EventAlgebra/Robertson.lean',
    'Lean/Screen/SL2F5PortCover.lean',
    'Lean/Screen/PortFrameGram.lean',
    'Lean/ObserverPatchHolography/CoreAxioms.lean',
    'Lean/Screen/A5PortAction.lean')}
for module in tuple(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name:
        p=Path(name).resolve()
        if p.is_file() and p.is_relative_to(base):cited.add(p)
pending=[base/'manifests/coset_carrier_reference.json',
         base/'manifests/super_tannakian_matter_reference.json']
seen=set()
while pending:
    p=pending.pop()
    if p in seen:continue
    seen.add(p);cited.add(p)
    if p.suffix!='.json':continue
    try:value=json.loads(p.read_text())
    except (ValueError,UnicodeError):continue
    def walk(x):
        if isinstance(x,dict):
            for y in x.values():walk(y)
        elif isinstance(x,list):
            for y in x:walk(y)
        elif (isinstance(x,str) and len(x)<240 and '\n' not in x
              and ('/' in x or x.endswith(('.json','.py','.lean','.md')))):
            for root in (base,h.RER,p.parent):
                target=(root/x).resolve()
                if target.is_relative_to(h.RER) and target.is_file():
                    pending.append(target);break
    walk(value)

result={'schema':'oph.native_field_spin_postrun_custody.v1',
        'collected_utc':datetime.now(timezone.utc).isoformat(),
        'timing':'Collected after evidence. Transitive and cited files below were not all preregistered; this is custody, not retroactive preregistration.',
        'preregistration_sha256':h.sha(h.HERE/'PREREG.json'),
        'runtime_dependency_sha256':{str(p.relative_to(h.ROOT)):h.sha(p) for p in sorted(runtime)},
        'cited_proof_dependency_sha256':{str(p.relative_to(h.ROOT)):h.sha(p) for p in sorted(cited)},
        'libraries':{'python':platform.python_version(),'numpy':numpy.__version__,
                     'scipy':scipy.__version__}}
(h.HERE/'SOURCE_PINS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({'runtime_files':len(runtime),'cited_files':len(cited),
                  'SOURCE_PINS_sha256':h.sha(h.HERE/'SOURCE_PINS.json')}))
