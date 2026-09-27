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
.PHONY: open preview render clean help dashboard dashboard-preview check

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

## check:   validate the paper manifest (link-check + status rollup)
check:
	julia ../scripts/check_paper_manifest.jl

## dashboard: validate the manifest, then render the reproduction dashboard + open
dashboard: check
	quarto render paper_dashboard.qmd --to html
	$(OPEN) paper_dashboard.html

## dashboard-preview: live-reloading dashboard (Ctrl-C to stop)
dashboard-preview:
	quarto preview paper_dashboard.qmd

## clean:   remove the rendered HTML and its _files/ directory
clean:
	rm -rf $(HTML) $(basename $(QMD))_files paper_dashboard.html paper_dashboard_files

## help:    list targets
help:
	@grep -E '^## ' $(MAKEFILE_LIST) | sed 's/^## /  /'
