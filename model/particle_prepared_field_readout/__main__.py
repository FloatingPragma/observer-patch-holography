import json

from .experiment import certificate

print(json.dumps(certificate(), indent=2))
