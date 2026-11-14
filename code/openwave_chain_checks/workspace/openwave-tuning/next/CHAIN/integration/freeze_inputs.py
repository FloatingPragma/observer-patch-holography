"""Maintainer utility to freeze reviewed inputs; not run by integrate.py."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
CHAIN=HERE.parent
ROOT=HERE.parents[3]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    frozen={}
    def add(p,expected=None):
        p=p.resolve()
        digest=sha(p)
        if expected is not None and digest!=expected:
            raise ValueError('input pin mismatch: '+str(p))
        name=str(p.relative_to(ROOT))
        if name in frozen:
            assert frozen[name]==digest
            return
        frozen[name]=digest
        if p.suffix=='.json':
            data=json.loads(p.read_text())
            # These maps use repository-root paths in the admitted packages.
            for key in ('source_sha256','source_pins','source_dependencies'):
                for relative,value in data.get(key,{}).items():
                    if isinstance(value,str) and len(value)==64:
                        candidates=[ROOT/relative,ROOT/'openwave-tuning'/relative,p.parent/relative]
                        candidates=[c for c in candidates if c.exists()]
                        if not candidates:raise ValueError('unresolved declared pin: '+relative)
                        matching=[c for c in candidates if sha(c)==value]
                        if not matching:raise ValueError('declared source drift: '+relative)
                        add(matching[0],value)
            for relative,value in data.get('files',{}).items():
                if isinstance(value,str) and len(value)==64 and (p.parent/relative).exists():
                    add(p.parent/relative,value)
    # Require all three opposing reviews before freezing.
    for name in ('audit_clock_force','audit_field','audit_particle'):
        if not (CHAIN/name/'AUDIT_CONTRACT.json').exists():
            raise ValueError('audit missing: '+name)
    for name in ('clock','field_limit','force','dynamical_field','particle','audit_clock_force','audit_field','audit_particle','coordination'):
        for p in sorted((CHAIN/name).glob('*')):
            if p.is_file() and p.suffix in ('.py','.json','.npz','.md'):
                add(p)
    op=ROOT/'openwave-tuning/next/ROW1_NATIVE/routes/quantum_operational'
    for name in ('OPERATIONAL_QUANTUM_REVIEW.json','operational_readout.py','operational_states.npz'):
        add(op/name)
    # This ledger is historical evidence with explicitly classified source drift.
    # Bind its bytes without treating every historical source pin as a live claim.
    ledger=ROOT/'openwave-tuning/next/ROW1_NATIVE/ALL31_READY.json'
    frozen[str(ledger.relative_to(ROOT))]=sha(ledger)
    for name in ('ROW1_SCOPE.json','SOURCE_CLOCK_MAP.json','EXTERNAL_SOURCE_MAP.json'):
        add(HERE/name)
    for relative,digest in json.loads((HERE/'ROW1_SCOPE.json').read_text())['inputs'].items():
        add(ROOT/relative,digest)
    data={'schema':'oph.chain.integration_frozen_inputs.v1','files':dict(sorted(frozen.items())),
          'policy':'All admitted producer and audit bytes are immutable inputs. External LOCAL12 comparison metadata is a separately hashed observation, not a live admitted dependency.'}
    (HERE/'INPUT_PINS.json').write_text(json.dumps(data,indent=2)+'\n')
    print(len(frozen),'frozen inputs')


if __name__=='__main__':main()
