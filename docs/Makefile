# Render the Quarto tutorials in docs/.
#
# Requires: quarto (https://quarto.org) and Julia with the GIPaper project
# instantiated (the tutorial cells `Pkg.activate ../GIPaper`).
#
# Quick start:
#   make            # render + open in browser (default)
#   make preview    # live-reloading server, auto-opens, Ctrl-C to stop
#   make render     # just build the HTML
#   make clean      # remove rendered artifacts
#
# Target a different tutorial:  make QMD=other_tutorial.qmd

QMD  ?= strong_decays_tutorial.qmd
HTML := $(QMD:.qmd=.html)

# Browser opener: `open` on macOS, `xdg-open` on Linux.
OPEN := $(shell command -v open >/dev/null 2>&1 && echo open || echo xdg-open)

.DEFAULT_GOAL := open
.PHONY: open preview render clean help

## open:    render the tutorial and open it in the default browser
open: render
	$(OPEN) $(HTML)

## preview: live-reloading server (auto-opens browser; Ctrl-C to stop)
preview:
	quarto preview $(QMD)

## render:  build the standalone HTML (only if the source changed)
render: $(HTML)

$(HTML): $(QMD)
	quarto render $(QMD) --to html

## clean:   remove the rendered HTML and its _files/ directory
clean:
	rm -rf $(HTML) $(basename $(QMD))_files

## help:    list targets
help:
	@grep -E '^## ' $(MAKEFILE_LIST) | sed 's/^## /  /'
