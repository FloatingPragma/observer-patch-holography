# Reverse Engineering Reality

*An Observer-First Theory of Everything.*

Forty-one chapters, a prologue, an epilogue and four appendices, in 48 files.

## Reading order

The two-digit filename prefix is the reading order and the only ordering rule. Sorting
`*.md` by name gives front matter, prologue, chapters 1 to 41, epilogue, appendices A to
D. `tools/build_book_pdf.py` does exactly that and rejects a file with no prefix or a
duplicated one.

The ten parts are a property of the chapter ranges rather than of the directory:

| Part | Chapters | Title |
|---|---|---|
| One | 1 to 6 | The Only Place to Stand |
| Two | 7 to 12 | The Protocol |
| Three | 13 to 16 | The Machine |
| Four | 17 to 22 | What an Observer Finds |
| Five | 23 to 25 | Why Gravity Is Different |
| Six | 26 to 28 | The Parts List |
| Seven | 29 to 31 | The Numbers and the One Line |
| Eight | 32 to 35 | Coming Back Down |
| Nine | 36 to 37 | The Loop Closes |
| Ten | 38 to 41 | Being One of These |

## Building

    python3 tools/build_book_pdf.py

Produces `book/reverse-engineering-reality-book.pdf` on a 6x9 trade trim through pandoc
and tectonic.

## Conventions

Paragraphs are single lines in the source, separated by blank lines, and wrapping is left
to the renderer. Straight quotes, American spelling. Bold marks a key term at first use
and is not used for emphasis. There are no em-dashes anywhere, by rule.


## Interpretation and source scope

The book is written under one premise, stated in its front matter: Observer Patch Holography is a consistent and complete theory of everything, and the book explains it from the inside, filling each gap with the theory's best-justified mechanism. The papers, the [OPH interpretation and working terms](../docs/INTERPRETATION_AND_TERMS.md) and the model documents carry the proofs, their assumptions and the computations that are work in progress. Statements about experience and ethics rest on further steps that the last chapters name.

The [current construction summary](../model/completion/CONTINUUM_FIELD_AND_ENTROPY.md) records what has been proved: fixed-amplitude small-ball entropy control on the free source, conditional spin-packet readout and exact current spectral bounds. Sharp spin-particle mass shells, measured-species identification, a common area response and fixed-volume generalized-entropy stationarity remain open. The [mass and energy table](../model/completion/PARTICLE_MASS_TABLE.md) distinguishes certified model energies from calibrated physical masses. These current results do not establish the book's assumed completeness premise.

[Appendix D](47-appendix-d.md#the-continuing-ce-campaign) includes the focused C/E campaign replay and navigation to the current [completion plan](../docs/COMPLETION_PROGRAM.md).
