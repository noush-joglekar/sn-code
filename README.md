# SnISOr-Seq

We are interested in cell-type specific isoform expression in the brain,
however, studies of splicing patterns are difficult to conduct in human because
human brain tissue is often frozen, and isolating spliced RNA from frozen tissue
is prohibitive on multiple levels.
Using a combination of microfluidics, PCR-based artifact removal,
target enrichment, and long-read sequencing, we developed 
**S**ingle-**n**uclei **Iso**form **R**NA **seq**uencing (*SnISOr-Seq*)
and applied it to the analysis of human adult frontal cortex samples.
This is all original code generated for processing of data described
in our [manuscript](https://www.biorxiv.org/content/10.1101/2021.12.29.474385v1).

Most of our processing workflow has already been
[published](https://github.com/noush-joglekar/scisorseqr) but
this repository contains code for initial data processing with scripts for:

- Extending barcode detection to allow mismatches
- Error correction to prevent 'molecule inflation' in UMI mismatches
- Coordination analysis for exon pairs
- Coordination analysis for exon-endsites

The structure is as follows:

```bash

lr_scripts
├── bashScripts
│   ├── parallelizeMolInflation.sh
│   ├── run_11merMatching.sh
│   └── v2.0a_bc_finder.awk
├── pythonScripts
│   ├── barcodeUMI_finder_11mer.py
│   ├── minEditDistance_post11merMatching.py
│   └── moleculeInflation_10x_forParallel.py
└── rScripts
    ├── endSite_exon_chiSq_countMatrix.R
    ├── endSite_exon_countMatrix_Generate.R
    ├── exon_coord_chiSq_countMatrix.R
    └── get10xMoleculeInfo.R

```

In addition, to ensure maximal reproducibility,
we have provided the input as well as code for producing panels
the main figures in our text. Code to generate isoform plots can be
found in the form of an R-package [ScisorWiz](https://github.com/ans4013/ScisorWiz)

```bash

plotGeneration
├── figure4.R
├── figure5.R
└── figure6.R

```

If you want to explore this dataset and query expression of disease-associated exons
in our human frontal cortex samples, we have made an interactive portal which allows
you to upload lists of exons and see their cell-type specific expression. Please visit
[isoformAtlas.com](https://isoformatlas.com/) --> Access Data --> 
[Exons](https://noush-joglekar.shinyapps.io/snisorDisease/) to view.

Feel free to [email](mailto:anj2026@med.cornell.edu) with questions
about the data or approach
