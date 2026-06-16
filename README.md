# Kataegis Detection Suite

A species-agnostic R-based toolkit for detecting and visualizing **kataegis** - localized clusters of hypermutation - from somatic variant data in Mutation Annotation Format (MAF) files. The suite includes a scriptable back-end library and an interactive Shiny web application for exploratory and publication-ready analysis.

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Running the Shiny App](#running-the-shiny-app)
- [Detection Module (`kataegis_detect.R`)](#detection-module-kataegis_detectr)
- [Visualization Module (`plots.R`)](#visualization-module-plotsr)
- [MAF File Format](#maf-file-format)
- [Typical Scripted Workflow](#typical-scripted-workflow)
- [Output Files](#output-files)
- [Algorithm Notes](#algorithm-notes)

---

## Overview

Kataegis is a mutational phenomenon characterized by clusters of closely spaced, strand-coordinated mutations - typically C>T or C>G transitions at TpC dinucleotides - often attributed to APOBEC cytidine deaminase activity. This suite provides:

- **Flexible detection** using either a fixed intermutation distance (IMD) threshold or a chromosome-adaptive dynamic cutoff
- **Interactive analysis** via a Shiny app supporting per-sample and cohort-aggregated detection
- **Publication-quality rainfall plots** using `ggplot2` and Base R rendering engines
- **Downstream utilities** for exporting results as TSV, BED, and PNG files

---

## Project Structure

```
kataegis_detection_suite
├── app.R                # Shiny web application (UI + server)
├── kataegis_detect.R    # Detection library (core algorithm + utilities)
├── plots.R              # Rainfall plot rendering functions (ggplot2 + Base R)
├── packages.txt         # Required R package list
└──  Readme.md            # This file
```

---

## Requirements

**R version:** ≥ 4.1 recommended

**Required packages** (listed in `packages.txt`):

| Package           | Purpose                                      |
|-------------------|----------------------------------------------|
| `shiny`           | Web application framework                    |
| `shinycssloaders` | Spinner feedback during computation          |
| `DT`              | Interactive tables in the Shiny UI           |
| `data.table`      | Fast MAF parsing and in-memory manipulation  |
| `dplyr`           | Data frame manipulation utilities            |
| `ggplot2`         | Publication-quality rainfall plots           |
| `scales`          | Axis formatting for log-scale plots          |

---

## Installation

Install all dependencies from the R console:

```r
install.packages(c("shiny", "shinycssloaders", "DT",
                   "data.table", "dplyr", "ggplot2", "scales"))
```

Or use the provided `packages.txt`:

```r
pkgs <- readLines("packages.txt")
install.packages(pkgs)
```

---

## Running the Shiny App

From the project directory in R:

```r
shiny::runApp("app.R")
```

Or from the terminal:

```bash
Rscript -e "shiny::runApp('app.R')"
```

### App Walkthrough

1. **Upload MAF File** - accepts `.maf`, `.tsv`, or `.txt` (case-insensitive). Required columns are listed in the [MAF File Format](#maf-file-format) section.
2. **Upload Chromosome Lengths** *(optional)* - accepts `.txt`, `.tsv`, or `.fai`. If `.fai`, only the first two columns (chromosome name, size) are used.
3. **Select Samples** - choose one, multiple, or all samples from the MAF. Leave blank to include all.
4. **Set Detection Mode:**
   - `Per-sample` - runs detection independently on each selected sample
   - `Aggregated` - pools all selected samples into a single pseudo-sample for cohort-level detection (analogous to `katdetectr` cohort mode)
5. **Set Analysis Parameters:**
   - `Use Dynamic Cutoff` - enables chromosome-adaptive IMD thresholds (see [Algorithm Notes](#algorithm-notes))
   - `Minimum Mutations per Cluster` - minimum SNPs to form a kataegis region (default: 6)
   - `Maximum Avg Intermutation Distance` - IMD ceiling in bp when using fixed cutoff (default: 1000 bp)
6. **Click ▶ Run Analysis**
7. **Adjust Plot Settings** - log-scale Y axis, facet by chromosome - without rerunning detection
8. **Export Results:**
   - `Download Plot (PNG)` - saves the rendered rainfall plot
   - `Download Kataegis TSV` - saves the detected regions table

Note - The shiny app has the advantage of being able to pool/aggregate the samples in any order you want. Feel free to explore both the app and the attached scripts. The shiny app sources both the detection and plotting script.
---

## Detection Module (`kataegis_detect.R`)

**Author:** Harish Kothandaraman  
**Inspired by:** `maftools::detect_kataegis` (changepoint approach) and `katdetectr` (IMD function)

Source this file to use the detection functions in a script:

```r
source("kataegis_detect.R")
```

### Functions

#### `detect_kataegis()`

Core detection engine using a sliding deque (double-ended queue) algorithm.

```r
detect_kataegis(
  maf_file,
  sample_id        = NULL,     # NULL processes all samples in the MAF
  min_mutations    = 6,        # Minimum cluster size
  max_avg_distance = 1000,     # Fixed IMD cutoff in bp (used if use_dynamic_cutoff = FALSE)
  use_dynamic_cutoff = FALSE,  # Enable chromosome-adaptive thresholds
  chrom_sizes      = NULL      # data.table with columns: Chromosome, Size
)
```

Returns a list of region objects, each containing sample ID, chromosome, start/end positions, mutation count, span (bp), average IMD, IMD cutoff used, raw mutations, and mutation spectrum counts.

#### `kataegis_to_dataframe()`

Converts the detection list into a flat `data.frame` suitable for downstream analysis.

```r
kataegis_to_dataframe(
  kataegis_list,
  include_spectrum = TRUE,  # Add per-type mutation counts (C>A, C>G, C>T, T>A, T>C, T>G)
  group_other      = TRUE   # Aggregate non-canonical types into an "Other" column
)
```

When `include_spectrum = TRUE`, the output also includes `pct_C_to_T` and `pct_C_to_G` proportion columns.

#### `get_kataegis_mutations()`

Prints a summary and returns the raw mutations for a specific region ID.

```r
get_kataegis_mutations(kataegis_list, region_id = 1)
```

#### `export_kataegis_bed()`

Exports detected regions as a 0-indexed BED file, compatible with the UCSC Genome Browser.

```r
export_kataegis_bed(kataegis_list, output_file = "kataegis_regions.bed")
```

#### `subset_kataegis_by_bed()`

Filters detection results to regions overlapping a provided BED file. Useful for restricting analysis to a genomic region of interest before plotting.

```r
subset_kataegis_by_bed(
  kataegis_list,
  bed_file,
  overlap_type = "any"   # Options: "any", "complete", "start"
)
```

---

## Visualization Module (`plots.R`)

Source this file to access rainfall plot functions:

```r
source("plots.R")
```

### Functions

#### `plot_rainfall_genome_ggplot2()` *(recommended)*

High-quality ggplot2-based renderer. Best for publication figures.

```r
plot_rainfall_genome_ggplot2(
  maf_file,
  sample_id,
  log_scale               = TRUE,
  highlight_kataegis      = NULL,    # Pass detect_kataegis() output to overlay arrows
  show_legend             = TRUE,
  chrom_sizes             = NULL,
  include_other_mutations = TRUE,
  facet_by_chromosome     = FALSE    # TRUE = per-chromosome panels; FALSE = genome-wide
)
```

**Features:**
- Genome-wide view with cumulative chromosome coordinates and vertical dividers
- Per-chromosome faceted view (`facet_wrap`, 4 columns, free X scales)
- Dashed horizontal reference line at IMD = 1000 bp
- Kataegis region arrows positioned at geometric mean IMD of cluster mutations
- Mutation type color palette: C>A (#4DBBD5), C>G (#E64B35), C>T (#00A087), T>A (#3C5488), T>C (#F39B7F), T>G (#8491B4)

#### `plot_rainfall_genome_base()`

Base R renderer. Faster for quick exploratory plots; no external graphics dependencies.

```r
plot_rainfall_genome_base(
  maf_file,
  sample_id,
  chromosomes             = NULL,    # NULL = genome-wide; character vector = multi-panel
  log_scale               = TRUE,
  highlight_kataegis      = NULL,
  ncol                    = 3,
  show_legend             = TRUE,
  chrom_sizes             = NULL,
  include_other_mutations = TRUE
)
```

---

## MAF File Format

The MAF file must be tab-delimited and include the following columns:

| Column                  | Description                                     |
|-------------------------|-------------------------------------------------|
| `Chromosome`            | Chromosome name (e.g., `chr1`, `1`)             |
| `Start_Position`        | 1-based mutation start coordinate               |
| `End_Position`          | 1-based mutation end coordinate                 |
| `Reference_Allele`      | Reference base (e.g., `C`)                      |
| `Tumor_Seq_Allele2`     | Alternate allele (e.g., `T`)                    |
| `Tumor_Sample_Barcode`  | Sample identifier                               |
| `Variant_Type`          | Must include `SNP` rows (indels are ignored)    |
| `Variant_Classification`| Mutation classification (e.g., `Missense_Mutation`) |

> TCGA-format MAF files are directly compatible. Only `SNP` rows are used for kataegis detection and rainfall plotting.
> All the SNPs are considered.

---

## Typical Scripted Workflow

```r
source("kataegis_detect.R")
source("plots.R")

# 1. Detect kataegis with chromosome-adaptive cutoff
results <- detect_kataegis(
  maf_file           = "cohort.maf",
  sample_id          = "TCGA-AB-1234",
  use_dynamic_cutoff = TRUE,
  min_mutations      = 6
)

# 2. Inspect as a data frame
df <- kataegis_to_dataframe(results, include_spectrum = TRUE)
head(df)

# 3. Export BED for UCSC Genome Browser
export_kataegis_bed(results, output_file = "kataegis.bed")

# 4. Rainfall plot for a single sample (ggplot2)
p <- plot_rainfall_genome_ggplot2(
  maf_file           = "cohort.maf",
  sample_id          = "TCGA-AB-1234",
  highlight_kataegis = results,
  log_scale          = TRUE
)
ggsave("rainfall_TCGA-AB-1234.png", p, width = 14, height = 6, dpi = 300)

# 5. Subset by a region of interest before plotting
roi_results <- subset_kataegis_by_bed(results, "regions_of_interest.bed", overlap_type = "any")
```

---

## Output Files

| File                     | Format | Contents                                           |
|--------------------------|--------|----------------------------------------------------|
| `kataegis_plot_<date>.png` | PNG    | Rendered rainfall plot                             |
| `kataegis_results.tsv`   | TSV    | Detected regions with mutation spectrum columns    |
| `kataegis.bed`           | BED    | 0-indexed regions for UCSC Genome Browser          |

---

## Algorithm Notes

### Deque Algorithm

Detection uses a sliding double-ended queue (deque) over position-sorted SNPs per chromosome. A window of `min_mutations` is initialized and extended forward as long as the average intermutation distance (IMD) stays ≤ the cutoff. When the threshold is exceeded, the window slides forward. Each qualifying window is saved as a kataegis region, and the queue reseeds from the end of the last region to allow adjacent clusters.

### Dynamic (Chromosome-Adaptive) Cutoff

When `use_dynamic_cutoff = TRUE`, the IMD threshold for each chromosome is calculated as:

$$\text{IMDcutoff} = \frac{-\ln\!\left(1 - \left(\frac{0.01}{W}\right)^{1/(N-1)}\right)}{\lambda}$$

where:
- \(W\) = chromosome width (true length if provided, else observed span)
- \(N\) = number of SNPs on the chromosome
- \(\lambda = \ln(2) / \text{median(IMDs)}\) - the modeled sample mutation rate

The resulting cutoff is capped at 1000 bp to prevent overly permissive thresholds on low-coverage chromosomes. This formula originates from the **Pan-Cancer Analysis of Whole Genomes (PCAWG) Consortium** and is demonstrated as a custom IMD cutoff in the `katdetectr` vignette. It is more sensitive than a fixed 1000 bp threshold on densely mutated chromosomes.

---

## Acknowledgements & License

### Attribution

**Deque detection algorithm (`kataegis_detect.R`)**

The core sliding deque algorithm was originally devised by **Moritz Goretzky (WWU Münster)** and underlies the detection logic in [`maftools::detect_kataegis`](https://github.com/PoisonAlien/maftools) (Mayakonda et al.). The detection scaffolding in this project is structurally inspired by that implementation. `maftools` is distributed under the **MIT License**.

> Detection logic adapted from `maftools` (MIT License) - Anand Mayakonda / PoisonAlien.  
> Original deque algorithm: Moritz Goretzky, WWU Münster.

**Dynamic IMD cutoff (`IMDcutoffFun`)**

The `IMDcutoffFun`, `modelSampleRate`, and `nthroot` helper functions in `kataegis_detect.R` are directly adapted from the custom IMD cutoff example provided in the [`katdetectr` vignette](https://bioconductor.org/packages/release/bioc/vignettes/katdetectr/inst/doc/General_overview.html) (Rens et al., *GigaScience* 2023). The vignette attributes the original formula to the **Pan-Cancer Analysis of Whole Genomes (PCAWG) Consortium**. `katdetectr` is distributed under the **GPL-3 License**.

The cutoff formula is:

$$\text{IMDcutoff} = \frac{-\ln\!\left(1 - \left(\frac{0.01}{W}\right)^{1/(N-1)}\right)}{\lambda}$$

where \(W\) is the segment/chromosome width, \(N\) is the number of variants, and \(\lambda = \ln(2) / \text{median(IMDs)}\) is the modeled sample mutation rate. The cutoff is capped at 1000 bp.

> Dynamic IMD cutoff adapted from `katdetectr` vignette (GPL-3 License) - Rens JA et al. / ErasmusMC-CCBC.  
> Original formula: Pan-Cancer Analysis of Whole Genomes (PCAWG) Consortium.

### Key Extensions Beyond Source Packages

This suite adds the following capabilities not present in either maftools or katdetectr:

- Species-agnostic design - no hardcoded reference genome (hg19/hg38)
- Multi-sample MAF processing in a single run
- Chromosome-adaptive IMD cutoffs applied per-chromosome within a sample
- Interactive Shiny application with per-sample and cohort-aggregated detection modes
- ggplot2-based rainfall renderer with faceting, cumulative genome coordinates, and publication-ready output
- BED export and BED-based region subsetting utilities

### License

Given the structural similarity of `kataegis_detect.R` to `maftools` (MIT) and the adaptation of the IMD cutoff from `katdetectr` (GPL-3), this project is released under the **GPL-3 License** to ensure compatibility with the more restrictive of the two upstream licenses.

See [https://www.gnu.org/licenses/gpl-3.0.html](https://www.gnu.org/licenses/gpl-3.0.html) for full license terms.

### References

- Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. (2018). Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Research*. doi:10.1101/gr.239244.118
- Rens JA, et al. (2023). Katdetectr: an R/Bioconductor package utilizing unsupervised changepoint analysis for robust kataegis detection. *GigaScience*. doi:10.1093/gigascience/giad081
- ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium. (2020). Pan-cancer analysis of whole genomes. *Nature*, 578, 82–93. doi:10.1038/s41586-020-1969-6
