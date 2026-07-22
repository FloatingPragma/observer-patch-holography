# OPH Book Source

This directory contains the canonical Markdown source for *Reverse Engineering
Reality*, the OPH book:

- `prologue.md`
- `chapter-01-*.md` through `chapter-20-*.md`
- `appendix-*.md`
- `epilogue.md`

## Contents

- [Prologue](prologue.md)

1. [Consistency](chapter-01-consistency.md)
2. [Lineage](chapter-02-lineage.md)
3. [The Screen](chapter-03-screen.md)
4. [Entropy](chapter-04-entropy.md)
5. [Algebra](chapter-05-algebra.md)
6. [Overlap](chapter-06-overlap.md)
7. [Recovery](chapter-07-recovery.md)
8. [Holography](chapter-08-holography.md)
9. [Entanglement](chapter-09-entanglement.md)
10. [Error Correction](chapter-10-error-correction.md)
11. [Maximum Entropy](chapter-11-maxent.md)
12. [Symmetry](chapter-12-symmetry.md)
13. [de Sitter Space](chapter-13-desitter.md)
14. [The Standard Model](chapter-14-standard-model.md)
15. [Relativity](chapter-15-relativity.md)
16. [Matter](chapter-16-matter.md)
17. [Darwinian Selection](chapter-17-darwin.md)
18. [The Strange Loop](chapter-18-strangeloop.md)
19. [Synthesis](chapter-19-synthesis.md)
20. [Metaphysics](chapter-20-metaphysics.md)

- [Appendix: Concept Glossary](appendix-concept-glossary.md)
- [Appendix: Extended Interludes](appendix-extended-interludes.md)
- [Epilogue](epilogue.md)

The list above follows the build order: prologue, numbered chapters,
appendices, epilogue. Narration exports for the audiobook are generated from
these files and live in
[`../../book/audiobook-chapters/`](../../book/audiobook-chapters/).
The technical
[chapter ledger](../docs/BOOK_CHAPTER_LEDGER.md) lives with the project docs and
is not included in the reader-facing book build.

## Print PDF Build

The polished print-style PDF is built with:

```bash
cd reverse-engineering-reality
python3 tools/build_book_pdf.py
```

This writes the finished PDF inside this repo as:

- `book/reverse-engineering-reality-book.pdf`

The build also writes temporary intermediate files under:

- `../temp/book_pdf_build/`

## Requirements

The builder expects these tools to be available in `PATH`:

- `pandoc`
- `tectonic`
- `rsvg-convert`

## Layout And Styling

The book PDF builder is implemented in:

- `tools/build_book_pdf.py`

The LaTeX header and print styling live in:

- `tools/book_pdf_header.tex`

Current output style:

- trade-book trim (`6in x 9in`)
- two-sided layout
- full-page cover using `assets/book-cover.svg`
- curated inline diagrams from `assets/book_diagrams/`
- chapter-level table of contents
- widow/orphan control
- running headers
- tightened long-table handling for summary tables

## Notes

- The builder automatically combines the prologue, numbered chapters, appendices, and epilogue in order.
- Unnumbered front/back matter gets corrected running heads (`Prologue`, `Epilogue`).
- SVG assets referenced by the book are converted to PDF during the build.
- A few long inline scientific expressions are normalized into TeX-safe forms before compilation.
- `book/reverse-engineering-reality-book.pdf` is a release artifact. Rebuild it before publishing
  whenever book chapters or public release materials change.

## License

The OPH book source and generated book artifact are part of the OPH public
repository. See the main [LICENSE](../LICENSE).
