import json
from .history import certificate

print(json.dumps(certificate(), indent=2))
