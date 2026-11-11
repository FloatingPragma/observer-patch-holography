"""Public replay with immutable research comparison and closed source custody."""
from pathlib import Path
from . import model
from ..certificate import public_certificate,run_cli

HERE=Path(__file__).resolve().parent

def certificate():
    return public_certificate('native_orthogonal_wilson_spatial_join',model.certificate,HERE)

def main():
    run_cli(certificate,HERE)

if __name__=='__main__':main()
