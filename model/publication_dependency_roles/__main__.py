from pathlib import Path
import json
from model.openwave_pure_su2_source.custody import pin_package,execute
from .check import verify,HERE


def pins():
    registry=json.loads((HERE/'roles.json').read_text())
    return pin_package(HERE,['model/publication_dependency_roles/roles.json']+
        [x['owner'] for x in registry['edges']]+[x['provider'] for x in registry['edges']])


if __name__=='__main__':execute(HERE,pins,verify)
