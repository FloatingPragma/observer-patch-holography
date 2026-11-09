import json
from .cone import certificate

print(json.dumps(certificate(), indent=2, sort_keys=True, allow_nan=False))
