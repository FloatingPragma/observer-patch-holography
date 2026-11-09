import json
from .state import certificate

print(json.dumps(certificate(), indent=2))
