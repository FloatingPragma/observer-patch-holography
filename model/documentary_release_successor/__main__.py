import argparse
import json
from .model import validate_current_science


def main():
    parser = argparse.ArgumentParser(description='Check current science with retained documentary identity and actual current PDFs.')
    parser.add_argument('--check', action='store_true', required=True)
    parser.parse_args()
    print(json.dumps(validate_current_science(), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
