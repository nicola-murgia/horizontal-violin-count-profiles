# Horizontal Count-Scaled Violin Profiles

R script generating **horizontal count-scaled violin profiles** of tyrosine hydroxylase-positive (TH+) neurons in the substantia nigra pars compacta (SNpc).

**Paper:** Li et al. (2025) *"Selective abundance of the stemness-promoting cluster miR-290-295 within the adult substantia nigra dopamine neurons is neuroprotective via preservation of protein synthesis"* — bioRxiv  
📄 <https://www.biorxiv.org/content/10.1101/2025.08.05.668587v1>

---

## Description

This script produces publication-ready vector graphics (SVG) displaying TH+ neuron counts across 11 serial coronal sections of the mouse SNpc. Each profile is a **count-scaled violin** where the full vertical thickness at every section position equals the exact neuronal count (1 neuron = 1 SVG pixel). The shape is drawn with smooth cubic Bézier curves using shape-preserving (PCHIP) interpolation. Individual counts are indicated by points and numerical labels. Figures are exported as editable SVG files. No normalization, data transformation, group averaging, or inferential statistical analysis is performed.

The script is fully self-contained and requires **no additional R packages** beyond base R 4.1.1+.

---

## Data

TH-positive neurons in the substantia nigra pars compacta (SNpc) were quantified across 11 sequential coronal sections in wild-type (WT; _n_ = 6), homozygous (Homo; _n_ = 5), and heterozygous (Het; _n_ = 10) animals. Each animal was considered an independent biological replicate, and section order was preserved.

---

## Output

The script generates an editable SVG file (`TH_SNc.svg`) with:
- Horizontal count-scaled violin profiles for each animal
- Individual neuronal-count dots at each section
- Numerical count labels above each dot
- Grid lines for all 11 coronal sections
- Legend and axis labels

Output is saved to an `output/` directory in the current working directory by default.

---

## Usage

### Basic run
```bash
Rscript horizontal_violin_count_profiles.R
```

### Custom output path
```bash
Rscript horizontal_violin_count_profiles.R --output /path/to/figures
```

### In RStudio
Open the script and click **Source**, or run:
```r
source("horizontal_violin_count_profiles.R")
```

---

## Parameters

Key parameters can be edited at the top of the script:

| Parameter | Default | Description |
|---|---|---|
| `ANIMAL_CENTER_GAP_PX` | 150 | Vertical spacing between animal profiles |
| `VIOLIN_FILL_OPACITY` | 0.36 | Fill opacity of the violin shapes |
| `POINT_RADIUS_PX` | 2.4 | Radius of individual count dots |
| `SHOW_COUNT_LABELS` | TRUE | Toggle numeric count labels |
| `PIXELS_PER_NEURON` | 1 | Scale factor (1 = 1 neuron per SVG pixel) |

---

## Citation

If you use this code in your work, please cite the associated preprint:

> Li, Z., Xu, Y., Murgia, N., Kovzel, N., Ng, D. S., Liu, Y., Kang, X., Jiang, H., Domanskyi, A., & Vinnikov, I. A. (2025). Selective abundance of the stemness-promoting cluster miR-290-295 within the adult substantia nigra dopamine neurons is neuroprotective via preservation of protein synthesis. *bioRxiv*.  
> <https://doi.org/10.1101/2025.08.05.668587>

---

## License

[MIT](LICENSE)