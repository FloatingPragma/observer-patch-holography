import json
from .smooth import certificate

print(json.dumps(certificate(), indent=2))
