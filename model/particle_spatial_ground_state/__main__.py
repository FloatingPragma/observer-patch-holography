import json
from fractions import Fraction
from .ground import certificate
print(json.dumps(certificate(),indent=2,default=lambda value: str(value) if isinstance(value,Fraction) else value))
