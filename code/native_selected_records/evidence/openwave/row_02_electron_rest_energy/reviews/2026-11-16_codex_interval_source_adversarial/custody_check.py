"""Read-only source/certificate custody and finite dyadic endpoint verification."""
import argparse
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
CERTIFICATE=HERE.parents[1]/'work/2026-11-16_codex_interval_source_certificate/CERTIFICATE.json'
SOURCE=HERE.parents[2]/'row_01_charge_quantization/work/2026-11-16_codex_common_reference_menu'


def require(ok):
    if not ok:raise AssertionError('Source custody or endpoint encoding failed')


def walk(value):
    if isinstance(value,dict):
        if 'lo_dyadic' in value:
            for key in ('lo_dyadic','hi_dyadic'):
                sign,mantissa,exponent=value[key]
                require(sign in (0,1) and int(mantissa)>=0 and type(exponent) is int)
                require(int(mantissa)!=0 or exponent==0)
            return 1
        return sum(walk(x) for x in value.values())
    if isinstance(value,list):return sum(walk(x) for x in value)
    return 0


def run():
    manifest=SOURCE/'MANIFEST.sha256';digest=hashlib.sha256(manifest.read_bytes()).hexdigest()
    require(digest=='a9ced490f0bc6fee5ed08aacb6ef9c8613758ee1af243ce3e98aa153928ebf2c')
    for line in manifest.read_text().splitlines():
        expected,name=line.split('  ',1);require(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==expected)
    certificate=json.loads(CERTIFICATE.read_text())
    for name,expected in certificate['source_files'].items():
        require(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==expected)
    return {'passed':True,'producer_native_source_manifest_sha256':digest,
            'manifest_leaves':len(manifest.read_text().splitlines()),'bound_source_receipts':len(certificate['source_files']),
            'certificate_file_sha256':hashlib.sha256(CERTIFICATE.read_bytes()).hexdigest(),
            'finite_dyadic_interval_records':walk(certificate),'zero_endpoint_special_values':False,
            'scope':'Read-only source custody and endpoint shape/finite encoding audit'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--verify',type=Path,default=HERE/'CUSTODY.json');mode.add_argument('--record-to',type=Path)
    args=parser.parse_args();result=run()
    if args.record_to:
        with args.record_to.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    else:require(result==json.loads(args.verify.read_text()))
    print(json.dumps(result))
