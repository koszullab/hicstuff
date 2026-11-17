---
title: hicstuff
hide:
  - navigation
  - toc
  - footer
---

<div class="cx-landing" markdown>

<h1 class="cx-masthead"><img src="assets/mark.svg" alt="" width="64" height="64"><span>hicstuff</span></h1>

<p class="cx-tagline">Hi-C contact maps, from reads to matrices.</p>

[hicstuff](https://github.com/koszullab/hicstuff) aligns Hi-C reads, pairs
them, filters out the uninformative 3C events and builds the contact map, in
one command, for any genome. It also views, rebins, subsamples and converts
contact maps, computes the distance law, and exposes every step as a Python
library. Its recommended mapping mode, `5spm`, aligns the reads as the 4DN
pipeline does, with bwa running in-process: `pip install hicstuff` is all it
needs.
{ .cx-lede }

<div class="cx-install" markdown>

```bash
pip install hicstuff
```

</div>

<div class="cx-actions" markdown>

[Why hicstuff?](why.md){ .cx-button }
[Quickstart :lucide-arrow-right:](quickstart.md){ .cx-button .cx-button--filled }

</div>

<div class="cx-badges" markdown>

[![PyPI version](https://img.shields.io/pypi/v/hicstuff.svg)](https://pypi.org/project/hicstuff/)
[![Python versions](https://img.shields.io/pypi/pyversions/hicstuff.svg)](https://pypi.org/project/hicstuff/)
[![bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/recipes/hicstuff/README.html)
[![CI](https://github.com/koszullab/hicstuff/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/koszullab/hicstuff/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/koszullab/hicstuff/branch/master/graph/badge.svg)](https://codecov.io/gh/koszullab/hicstuff)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2620608.svg)](https://doi.org/10.5281/zenodo.2620608)

</div>

<p class="cx-kicker">Start here</p>

<div class="cx-tiles" markdown>

[<span class="cx-glyph">:lucide-rocket:</span> <strong>Quickstart</strong> <span class="cx-tile-text">From two FASTQ files to a contact map, with nothing but pip.</span>](quickstart.md){ .cx-tile }
[<span class="cx-glyph">:lucide-lightbulb:</span> <strong>Concepts</strong> <span class="cx-tile-text">Fragments, bins, pairs and 3C events.</span>](concepts.md){ .cx-tile }
[<span class="cx-glyph">:lucide-square-terminal:</span> <strong>CLI reference</strong> <span class="cx-tile-text">Every command and option.</span>](cli.md){ .cx-tile }
[<span class="cx-glyph">:lucide-braces:</span> <strong>Python API</strong> <span class="cx-tile-text">The functions behind the commands.</span>](api/hicstuff.md){ .cx-tile }

</div>

</div>
