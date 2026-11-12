# Codon accessibility atlas: figure code

Code and data for the figures in *Synonymous codons constrain which insecticide resistance mutations evolve across arthropods*.

## Reproducing the figures

    python3 src/fig99_build_and_verify.py      # all nine figures, about 40 s
    python3 src/fig96_export_tiff.py           # the production TIFFs

`fig99` rebuilds every figure from `data/` and `results/`, then checks each
written PDF against the journal specification (column width, page height,
font sizes, embedded TrueType, nothing outside the page) and exits non-zero if
one fails. `fig96` writes the 600 dpi RGB TIFF beside each vector PDF.

Fonts are in `assets/fonts/`. The figures are laid out for Arial; Liberation
Sans is metrically compatible and is loaded automatically, so the published
metrics reproduce on a machine without Arial. Without it matplotlib falls back
to DejaVu Sans, which is wider, and the line breaks differ.

## Which script makes which figure

| Figure | Script | What it shows |
|---|---|---|
| Fig. 1 | `src/fig01_dataset_and_context.py` | the atlas, and the codon each species carries at VGSC L1014 |
| Fig. 2 | `src/fig02_forest.py` | per-site forest and the specification ladder |
| Fig. 3 | `src/step8_04_figure2.py` | accessibility against tolerability, and the odds ratios |
| Fig. 4 | `src/step8_12_figure4_supply.py` | the mutational-supply exponent in context, and what more data settles |
| Fig. 5 | `src/step8_11_figure5_enumeration.py` | which sites the genetic code permits the question to be asked at |
| Supplementary Fig. 1 | `src/step8_06_figureS1.py` | composition of the atlas, and independent-assembly replication |
| Supplementary Fig. 2 | `src/step8_08_figureS2.py` | the estimate under every stratification |
| Supplementary Fig. 3 | `src/step8_10_figureS3_matched_nulls.py` | rapid adaptation against deep time, under matched nulls |
| Supplementary Fig. 4 | `src/step8_05_figure3.py` | the transition spectrum, and the three nulls |

`src/fig00_style.py` holds the shared style and the specification checks.
`src/atlas_common.py` holds the paths and `src/codon_engine.py` the
genetic-code operations.

## Contents

| | |
|---|---|
| `src/` | the 12 figure and driver scripts, the 26 scripts that compute a value a figure plots, and the 6 modules they import |
| `data/` | the codon atlas and the species tree, which the deposit starts from |
| `results/` | the analysis tables every plotted number is read from, including `canonical_quantities.csv` |
| `figures/` | each figure as published, in vector PDF, 600 dpi RGB TIFF and 600 dpi PNG, with its source data |
| `MANIFEST.csv` | every file with its SHA-256, and for each script which of the three it is |
| `VERSIONS.txt` | the environment that drew the figures, and the pinned one that computed the tables |

Not included: the raw genome assemblies and proteomes the atlas was built
from, which are public downloads named in the Methods with their accessions,
and the scripts that build the atlas from them. Everything downstream of
`data/processed/codon_atlas.csv` is here and runs.

Scripts name the upstream step that writes each of their inputs, and a few of
those steps are among the ones not included. Their outputs are shipped in
`results/`, so nothing here has to be regenerated and no such message is
reached on a clean run.

## Software

Python 3.11 with NumPy, SciPy, pandas, Matplotlib, statsmodels, Biopython and
openpyxl. Versions in `VERSIONS.txt`.

## Licence

MIT. See LICENSE.
