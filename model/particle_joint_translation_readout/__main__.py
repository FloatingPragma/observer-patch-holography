import json
from .readout import certificate

print(json.dumps(certificate(), indent=2, sort_keys=True))
