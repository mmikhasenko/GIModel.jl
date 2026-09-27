# Agentic Godfrey-Isgur

A learning track from elementary quantum mechanics to the Godfrey-Isgur model.
`part1/` contains nine pen-and-paper theory sheets; `part2/` begins the
computational track. Each sheet includes problems, worked solutions, and short
post-solution concept checks.

Problem statements specify the inputs, conventions, and requested result.
Necessary assumptions belong in the statement; short `Hint:` paragraphs supply
an identity or a starting method. Students should be able to begin from the
statement and its cited earlier results without reading the worked solution.
Conceptual design questions specify the expected output, such as a table,
controlled comparison, or brief pseudocode.

## Compile

A LaTeX installation with `latexmk` is required. From this directory, run:

```sh
make
```

The generated PDFs are written to the gitignored `pdf/` directory. To remove
LaTeX intermediate files while keeping the PDFs, run:

```sh
make clean
```
