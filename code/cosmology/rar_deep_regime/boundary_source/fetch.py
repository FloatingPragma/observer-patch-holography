"""Reacquire the author-deposited Leroy et al. source columns from VizieR."""

import argparse
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
CATALOG = "J/AJ/136/2782"
COLUMNS = {
    "table7": "Name,r,r.n,SigmaHI,e_SigmaHI,SigmaH2,e_SigmaH2,Sigma*,e_Sigma*",
    "sample": "Name,Dist,i,Res,r25,l*,l_logMH2,logMH2",
}


def canonical(table, raw):
    lines = [COLUMNS[table].replace(",", "\t")]
    for line in raw.splitlines():
        if line.startswith(("DDO ", "Ho ", "IC ", "NGC ")):
            lines.append("\t".join(field.strip() for field in line.split("\t")))
    if len(lines) != {"table7": 688, "sample": 24}[table]:
        raise ValueError(f"Unexpected row count for {table}: {len(lines) - 1}")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cached-prefix", type=Path)
    args = parser.parse_args()
    for table, columns in COLUMNS.items():
        if args.cached_prefix:
            raw = Path(f"{args.cached_prefix}_{table}.txt").read_text()
        else:
            query = urlencode(
                {
                    "-source": f"{CATALOG}/{table}",
                    "-out": columns,
                    "-out.max": "unlimited",
                }
            )
            with urlopen(
                f"https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?{query}", timeout=60
            ) as reply:
                raw = reply.read().decode("utf-8")
        (HERE / f"leroy_{table}.dat").write_text(canonical(table, raw))


if __name__ == "__main__":
    main()
