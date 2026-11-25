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

The book follows the proposed observer history described in its front matter. The main candidate is [Interval Echo](../model/candidate_laws/interval_echo/README.md). Small observers occupy the ports of twelve-port cells, and the laws under test act on their registers. Chapter 28 describes the [signed material excitation](../model/candidate_laws/interval_echo/openwave/electron_field/README.md) and the experiment that distinguishes its orbital state from the electron’s spin. Its ports, seams and cells are the hardware of a computation and no space, and a test compares what its observers record with what we record. The papers distinguish its theorems, its computed records and the parts that its makers supply. Independent repair and clock examples retain their own source assumptions. Statements about experience and ethics rest on further steps named in the last chapters.

The official [scorecard](../model/evaluation/SCORECARD.md) states what the correct law must do and shows how each candidate meets it, and the [observables register](../docs/registers/OBSERVABLES.json) lists what it reproduces. Physical spectra, a common continuum and the regional entropy and stress estimates required for Einstein dynamics need their own source attachment.

[Appendix D](47-appendix-d.md) explains reproduction and the scope of the checks.
