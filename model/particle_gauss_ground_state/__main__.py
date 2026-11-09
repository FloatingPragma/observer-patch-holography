import json
from .gauss import certificate
print(json.dumps(certificate(include_records=True),indent=2,allow_nan=False))
