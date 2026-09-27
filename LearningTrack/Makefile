SOURCES := $(sort $(wildcard part1/*.tex) $(wildcard part2/*.tex))
PDFS := $(patsubst %.tex,pdf/%.pdf,$(notdir $(SOURCES)))

.PHONY: all clean
all: $(PDFS)

pdf/%.pdf: part1/%.tex course-common.tex | pdf
	latexmk -cd -pdf -interaction=nonstopmode -halt-on-error -outdir=../pdf $<

pdf/%.pdf: part2/%.tex course-common.tex | pdf
	latexmk -cd -pdf -interaction=nonstopmode -halt-on-error -outdir=../pdf $<

pdf:
	mkdir -p pdf

clean:
	latexmk -cd -c -outdir=../pdf $(SOURCES)
