import json

from .ground import certificate

print(json.dumps(certificate(), indent=2, sort_keys=True))
