#!/usr/bin/env python3
"""Render the credited prior work of each paper from the claims ledger.

``claims/prior_work.json`` holds the verified bibliographic entries for the prior
work that the registered claims build on. ``claims/claim_registry.yaml`` records,
in ``claim_context[<claim>].prior_work.citations``, each established result a claim
uses as ``{"work": <name>, "keys": [<entry keys>]}``; a result of the project
itself carries ``"keys": []`` and ``"internal": true``. This tool

* validates the bibliography (author, title, year, venue and a resolvable
  identifier per entry; unique keys) and every citation key the registry uses;
* renders, for every paper that owns registered claims, a section listing the
  credited prior work with citations (``paper/tex_fragments/PRIOR_WORK_*.tex``);
* renders the bibliography entries that paper does not already carry
  (``paper/tex_fragments/PRIOR_WORK_*_BIBITEMS.tex`` for ``thebibliography``
  papers, ``extra/*_prior_work.bib`` for BibTeX papers). An entry the paper
  already lists, matched by DOI, arXiv identifier or title, is cited under the
  paper's own key instead of being listed twice.

``--check`` fails when a rendered file differs from its committed bytes.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "claims" / "claim_registry.yaml"
BIBLIOGRAPHY = ROOT / "claims" / "prior_work.json"
FRAGMENTS = ROOT / "paper" / "tex_fragments"

REQUIRED_FIELDS = ("key", "type", "author", "title", "year")
IDENTIFIER_FIELDS = ("doi", "eprint", "isbn", "url")
VENUE_FIELDS = ("journal", "booktitle", "publisher", "howpublished", "institution", "school")
ENTRY_TYPES = {"article", "book", "inproceedings", "incollection", "misc", "phdthesis", "techreport"}


@dataclass(frozen=True)
class Target:
    paper: str          # owner paper path used by claim rows
    stem: str           # fragment stem
    bib_mode: str       # "bibitems" or "bibtex"
    bib_sources: tuple  # files whose existing references are matched
    extra_owners: tuple = ()  # further owner paths whose claims this paper states


def _target(paper: str, bib_mode: str, bib_sources: tuple, extra_owners: tuple = ()) -> Target:
    stem = Path(paper).stem.upper()
    return Target(paper, stem, bib_mode, bib_sources, extra_owners)


SHARED_REFERENCES = "paper/tex_fragments/OBSERVERS_REFERENCES.tex"
TARGETS = (
    _target("flagship/from_observer_consensus_to_standard_physics.tex", "bibitems",
            ("flagship/from_observer_consensus_to_standard_physics.tex",)),
    _target("paper/observers_are_all_you_need.tex", "bibitems",
            ("paper/observers_are_all_you_need.tex", SHARED_REFERENCES)),
    _target("paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex", "bibitems",
            ("paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex", SHARED_REFERENCES)),
    _target("paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex", "bibitems",
            ("paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex", SHARED_REFERENCES)),
    _target("paper/reality_as_consensus_protocol.tex", "bibitems",
            ("paper/reality_as_consensus_protocol.tex",)),
    _target("paper/deriving_the_particle_zoo_from_observer_consistency.tex", "bibitems",
            ("paper/deriving_the_particle_zoo_from_observer_consistency.tex",)),
    _target("paper/screen_microphysics_and_observer_synchronization.tex", "bibitems",
            ("paper/screen_microphysics_and_observer_synchronization.tex",),
            extra_owners=("paper/tex_fragments/SCREEN_SPECTRUM_THEOREMS.tex",)),
    _target("extra/machine_checked_finite_event_algebras.tex", "bibtex",
            ("extra/machine_checked_finite_event_algebras.bib",)),
    _target("extra/observable_normal_forms.tex", "bibtex",
            ("extra/observable_normal_forms.bib",)),
    _target("extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex", "bibitems",
            ("extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex",)),
    _target("extra/yang_mills_gap_clay_problem.tex", "bibitems",
            ("extra/yang_mills_gap_clay_problem.tex",)),
    _target("extra/fine_structure_constant_derivation.tex", "bibitems",
            ("extra/fine_structure_constant_derivation.tex",)),
    _target("cosmology/oph_dark_matter_paper.tex", "bibitems",
            ("cosmology/oph_dark_matter_paper.tex",)),
    _target("cosmology/oph_cosmology_finite_source_cmb_program.tex", "bibitems",
            ("cosmology/oph_cosmology_finite_source_cmb_program.tex",)),
    _target("cosmology/oph_black_hole_information_ledger.tex", "bibitems",
            ("cosmology/oph_black_hole_information_ledger.tex",)),
)
# Papers without a bibliography environment (three cosmology working papers)
# keep their claims' citations in the ledger only.


class PriorWorkError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PriorWorkError(message)


def load_bibliography(path: Path = BIBLIOGRAPHY) -> dict[str, dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    require(data.get("schema") == "oph.prior_work.v1", f"{path.name}: schema must be oph.prior_work.v1")
    entries = data.get("entries")
    require(isinstance(entries, list) and entries, f"{path.name}: entries must be a nonempty list")
    by_key: dict[str, dict] = {}
    for entry in entries:
        key = entry.get("key", "<missing>")
        for field in REQUIRED_FIELDS:
            require(isinstance(entry.get(field), str) and entry[field].strip(), f"{path.name}: {key}: missing {field}")
        require(re.fullmatch(r"[a-z][a-z0-9]*", key) is not None, f"{path.name}: {key}: keys are lowercase ASCII letters and digits")
        require(key not in by_key, f"{path.name}: duplicate key {key}")
        require(entry["type"] in ENTRY_TYPES, f"{path.name}: {key}: unknown type {entry['type']}")
        require(re.fullmatch(r"\d{4}", entry["year"]) is not None, f"{path.name}: {key}: year must have four digits")
        require(any(entry.get(field) for field in IDENTIFIER_FIELDS), f"{path.name}: {key}: needs a doi, eprint, isbn or url")
        require(any(entry.get(field) for field in VENUE_FIELDS) or entry.get("eprint"), f"{path.name}: {key}: needs a venue or an eprint")
        require(entry.get("verified") is True, f"{path.name}: {key}: every entry must be verified against a primary source")
        by_key[key] = entry
    require([entry["key"] for entry in entries] == sorted(by_key), f"{path.name}: entries must be sorted by key")
    return by_key


def claims_and_context(registry: Path = REGISTRY) -> tuple[list[dict], dict]:
    data = json.loads(registry.read_text(encoding="utf-8"))
    return data["claims"], data["claim_context"]


def check_citations(claims: list[dict], context: dict, bibliography: dict) -> None:
    used: set[str] = set()
    for claim in claims:
        claim_id = claim["claim_id"]
        records = context[claim_id]["prior_work"]["citations"]
        require(isinstance(records, list), f"{claim_id}: prior_work.citations must be a list")
        for record in records:
            require(
                isinstance(record, dict) and set(record) <= {"work", "keys", "internal"}
                and isinstance(record.get("work"), str) and record["work"].strip()
                and isinstance(record.get("keys"), list),
                f"{claim_id}: each prior-work citation is a work name with a key list",
            )
            keys = record["keys"]
            require(len(keys) == len(set(keys)), f"{claim_id}: duplicate key for {record['work']!r}")
            unknown = [key for key in keys if key not in bibliography]
            require(not unknown, f"{claim_id}: {record['work']!r} cites keys without an entry: {unknown}")
            if keys:
                require("internal" not in record, f"{claim_id}: {record['work']!r} is cited and marked internal")
            else:
                require(record.get("internal") is True,
                        f"{claim_id}: {record['work']!r} has no citation; cite it or mark it internal")
            used.update(keys)
    unused = sorted(set(bibliography) - used)
    require(not unused, f"prior_work.json entries that no claim cites: {unused[:8]}")


# ---------------------------------------------------------------- formatting

_TEX_ESCAPES = {"&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_", "^": r"\^{}", "~": r"\textasciitilde{}"}

# The papers set Latin Modern through 8-bit T1 fonts, which carry no glyph above
# U+00FF; every non-ASCII character is therefore written as a TeX command.
_TEX_CHARACTERS = {
    "\u2013": "--", "\u2014": "---", "\u2018": "`", "\u2019": "'", "\u201c": "``", "\u201d": "''",
    "\u2026": r"\ldots{}", "\u00a0": "~", "\u00df": r"{\ss}", "\u00f8": r"{\o}", "\u00d8": r"{\O}",
    "\u00e6": r"{\ae}", "\u00c6": r"{\AE}", "\u0153": r"{\oe}", "\u0152": r"{\OE}", "\u0142": r"{\l}",
    "\u0141": r"{\L}", "\u0131": r"{\i}", "\u00e5": r"{\aa}", "\u00c5": r"{\AA}", "\u00d7": r"$\times$",
    "\u2032": "'", "\u00b7": r"$\cdot$", "\u2212": "-", "\u03b1": r"$\alpha$", "\u03b2": r"$\beta$",
    "\u03b3": r"$\gamma$", "\u03bc": r"$\mu$", "\u03c0": r"$\pi$", "\u00b1": r"$\pm$",
}
_COMBINING = {
    "\u0301": "'", "\u0300": "`", "\u0302": "^", "\u0303": "~", "\u0308": '"', "\u030c": "v",
    "\u030b": "H", "\u0327": "c", "\u0328": "k", "\u0304": "=", "\u0306": "u", "\u030a": "r",
    "\u0307": ".",
}


def _unicode_to_tex(char: str) -> str:
    if ord(char) < 128:
        return char
    if char in _TEX_CHARACTERS:
        return _TEX_CHARACTERS[char]
    import unicodedata

    base, *marks = unicodedata.normalize("NFD", char)
    if marks and ord(base) < 128 and all(mark in _COMBINING for mark in marks):
        out = base
        for mark in marks:
            command = _COMBINING[mark]
            out = (rf"\{command}{{{out}}}" if command.isalpha() else rf"\{command}{{{out}}}")
        return "{" + out + "}"
    raise PriorWorkError(f"no TeX form for character {char!r} (U+{ord(char):04X})")


def tex(text: str) -> str:
    out = "".join(_TEX_ESCAPES.get(char, char) for char in text)
    if out.count("$") % 2:
        out = out.replace("$", r"\$")
    return "".join(_unicode_to_tex(char) for char in out)


def authors_short(author: str) -> str:
    if author.startswith("{") and author.endswith("}"):
        return author[1:-1]
    names = [name.strip() for name in author.split(" and ") if name.strip()]
    rendered = []
    for name in names:
        if name.lower() == "others":
            rendered.append("others")
            continue
        if "," in name:
            last, first = [part.strip() for part in name.split(",", 1)]
        else:
            parts = name.split()
            last, first = parts[-1], " ".join(parts[:-1])
        initials = " ".join(f"{part[0]}." for part in re.split(r"[\s\-]+", first) if part and part[0].isalpha())
        rendered.append(f"{last}, {initials}" if initials else last)
    if len(rendered) > 6:
        return ", ".join(rendered[:6]) + ", et al."
    if len(rendered) == 1:
        return rendered[0]
    return ", ".join(rendered[:-1]) + ", and " + rendered[-1]


def identifier(entry: dict) -> str:
    if entry.get("doi"):
        return rf"\url{{https://doi.org/{entry['doi']}}}"
    if entry.get("eprint"):
        return rf"arXiv:{entry['eprint']}. \url{{https://arxiv.org/abs/{entry['eprint']}}}"
    if entry.get("isbn"):
        return f"ISBN {entry['isbn']}."
    return rf"\url{{{entry['url']}}}"


def bibitem(cite_key: str, entry: dict) -> str:
    head = f"{tex(authors_short(entry['author']))} ({entry['year']})."
    title = entry.get("title_tex") or tex(entry["title"])
    title = "".join(_unicode_to_tex(char) for char in title)
    if entry["type"] in {"book", "phdthesis", "techreport"}:
        body = rf"\emph{{{title}}}."
        venue = entry.get("publisher") or entry.get("institution") or entry.get("school")
        if venue:
            body += f" {tex(venue)}."
    else:
        body = f"``{title}.''"
        if entry.get("journal"):
            body += rf" \emph{{{tex(entry['journal'])}}}"
            if entry.get("volume"):
                body += rf" \textbf{{{tex(entry['volume'])}}}"
            if entry.get("pages"):
                body += f", {tex(entry['pages'])}"
            body += "."
        elif entry.get("booktitle"):
            body += rf" In \emph{{{tex(entry['booktitle'])}}}"
            if entry.get("pages"):
                body += f", {tex(entry['pages'])}"
            body += "."
            if entry.get("publisher"):
                body += f" {tex(entry['publisher'])}."
        elif entry.get("howpublished"):
            body += f" {tex(entry['howpublished'])}."
    return f"\\bibitem{{{cite_key}}}\n{head} {body} {identifier(entry)}\n"


def bibtex_entry(cite_key: str, entry: dict) -> str:
    fields = []
    for field in ("author", "title", "journal", "booktitle", "publisher", "institution", "school",
                  "howpublished", "volume", "number", "pages", "year", "doi", "eprint", "isbn", "url"):
        if entry.get(field):
            if field in {"doi", "url", "eprint", "isbn"}:
                value = entry[field]
            elif field == "title" and entry.get("title_tex"):
                value = entry["title_tex"]
            elif field == "author" and entry[field].startswith("{"):
                value = entry[field][1:-1]
            else:
                value = tex(entry[field])
            if field == "author" and entry[field].startswith("{"):
                value = "{" + value + "}"
            fields.append(f"  {field} = {{{value}}}")
    return f"@{entry['type']}{{{cite_key},\n" + ",\n".join(fields) + "\n}\n"


# ---------------------------------------------------------------- matching existing references

def _normal(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def existing_references(target: Target, root: Path = ROOT) -> list[tuple[str, str]]:
    """(key, block text) for every reference the paper already carries."""
    blocks: list[tuple[str, str]] = []
    for source in target.bib_sources:
        text = (root / source).read_text(encoding="utf-8")
        if source.endswith(".bib"):
            for match in re.finditer(r"@\w+\s*\{\s*([^,\s]+)\s*,(.*?)\n\}", text, re.S):
                blocks.append((match.group(1), match.group(2)))
        else:
            parts = re.split(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", text)
            for index in range(1, len(parts) - 1, 2):
                blocks.append((parts[index], parts[index + 1].split("\\end{thebibliography}")[0]))
    return blocks


def _block_titles(block: str) -> set[str]:
    """Normalized title fields of one reference: BibTeX, quoted or emphasized."""
    titles = re.findall(r"\btitle\s*=\s*[{\"](.+?)[}\"]\s*,?\s*$", block, re.M | re.I)
    titles += re.findall(r"``(.+?)''", block, re.S)
    titles += re.findall(r"\\emph\{([^{}]+)\}", block)
    return {_normal(title) for title in titles}


def _first_surname(author: str) -> str:
    if author.startswith("{") and author.endswith("}"):
        return _normal(author[1:-1])
    first = author.split(" and ")[0].strip()
    return _normal(first.split(",")[0] if "," in first else first.split()[-1])


def match_existing(entry: dict, blocks: list[tuple[str, str]]) -> str | None:
    doi = (entry.get("doi") or "").lower()
    eprint = (entry.get("eprint") or "").lower()
    title = _normal(entry["title"])
    for key, block in blocks:
        lower = block.lower()
        if doi and doi in lower:
            return key
        if eprint and eprint in lower:
            return key
        if len(title) >= 24 and title in _normal(block):
            return key
    # A short title is too common to match inside a reference. It matches when
    # it is that reference's whole title and the reference names the same
    # first author: Kemeny and Snell's "Finite Markov Chains", for instance.
    surname = _first_surname(entry["author"])
    for key, block in blocks:
        if title in _block_titles(block) and surname and surname in _normal(block):
            return key
    return None


# ---------------------------------------------------------------- rendering

def owned(target: Target, claims: list[dict]) -> list[dict]:
    owners = (target.paper, *target.extra_owners)
    return [claim for claim in claims
            if claim["owner_paper"] == target.paper
            or any(owner.endswith("/") and claim["owner_paper"].startswith(owner) for owner in target.extra_owners)
            or claim["owner_paper"] in owners]


def render_target(target: Target, claims: list[dict], context: dict, bibliography: dict,
                  root: Path = ROOT) -> dict[Path, str]:
    """Rendered files for one paper; empty when its claims cite no prior work."""
    mine = owned(target, claims)
    works: dict[str, list[str]] = {}
    keys: list[str] = []
    for claim in mine:
        for record in context[claim["claim_id"]]["prior_work"]["citations"]:
            if not record["keys"]:
                continue
            listed = works.setdefault(record["work"], [])
            for key in record["keys"]:
                if key not in listed:
                    listed.append(key)
                if key not in keys:
                    keys.append(key)
    if not keys:
        return {}
    blocks = existing_references(target, root)
    cite_key: dict[str, str] = {}
    new_entries: list[tuple[str, dict]] = []
    for key in sorted(keys):
        entry = bibliography[key]
        found = match_existing(entry, blocks)
        if found:
            cite_key[key] = found
        else:
            cite_key[key] = f"pw{key}"
            new_entries.append((cite_key[key], entry))
    lines = [
        "% Generated by tools/build_prior_work.py from the claims ledger; do not edit.",
        r"\section*{Credited prior work}",
        "",
        "The results of this paper build on the following established work. Each",
        "entry names an established result as it is used here and cites its sources.",
        "",
        r"\begin{itemize}\setlength{\itemsep}{0pt}",
    ]
    for name in sorted(works, key=lambda w: (w.lower(), w)):
        cites = ",".join(cite_key[key] for key in works[name])
        lines.append(rf"\item {tex(name)}~\cite{{{cites}}}.")
    lines += [r"\end{itemize}", ""]
    files = {root / "paper" / "tex_fragments" / f"PRIOR_WORK_{target.stem}.tex": "\n".join(lines)}
    if target.bib_mode == "bibtex":
        bib_path = root / target.bib_sources[0].replace(".bib", "_prior_work.bib")
        files[bib_path] = "% Generated by tools/build_prior_work.py from the claims ledger; do not edit.\n" + "\n".join(
            bibtex_entry(k, e) for k, e in new_entries)
    else:
        items = "".join(bibitem(k, e) + "\n" for k, e in new_entries)
        files[root / "paper" / "tex_fragments" / f"PRIOR_WORK_{target.stem}_BIBITEMS.tex"] = (
            "% Generated by tools/build_prior_work.py from the claims ledger; do not edit.\n" + items)
    return files


def renderings(root: Path = ROOT) -> dict[Path, str]:
    bibliography = load_bibliography(root / "claims" / "prior_work.json")
    claims, context = claims_and_context(root / "claims" / "claim_registry.yaml")
    check_citations(claims, context, bibliography)
    files: dict[Path, str] = {}
    for target in TARGETS:
        files.update(render_target(target, claims, context, bibliography, root))
    return files


def stale(root: Path = ROOT) -> list[str]:
    rendered = renderings(root)
    out = [str(path.relative_to(root)) for path, text in rendered.items()
           if not path.is_file() or path.read_text(encoding="utf-8") != text]
    managed = set(rendered)
    for path in sorted((root / "paper" / "tex_fragments").glob("PRIOR_WORK_*.tex")):
        if path not in managed:
            out.append(f"{path.relative_to(root)} (no longer rendered)")
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="compare rendered files; write nothing")
    args = parser.parse_args(argv)
    try:
        if args.check:
            problems = stale()
            if problems:
                print("prior-work renderings differ: " + ", ".join(problems), file=sys.stderr)
                return 1
            print("prior work OK: bibliography verified, citations resolve, paper sections current")
            return 0
        for path, text in renderings().items():
            path.write_text(text, encoding="utf-8")
            print(f"wrote {path.relative_to(ROOT)}")
    except PriorWorkError as error:
        print(f"prior work: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
