import argparse
from .overlap import check, write

parser = argparse.ArgumentParser()
parser.add_argument("--write", action="store_true")
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
result = write() if args.write else check()
row = result["profiles"]["positive_core_patch"]
print("Complete site kernels:", result["sites"], "; band rank:", result["complete_band_rank"])
print("Patch singular values:", row["singular_values"])
print("Patch strength ratio:", row["largest_to_smallest_ratio"])
