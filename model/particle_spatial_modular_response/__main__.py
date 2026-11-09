import json
from .regional import certificate

print(json.dumps(certificate(), indent=2, sort_keys=True, allow_nan=False))
