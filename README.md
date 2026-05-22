# MS2 Spectrum Visualization (DDA + PRM)

This directory provides two scripts with matched figure style:

- `dda_vis.py`: DDA visualization from `pepXML + mzXML`
- `prm_vis.py`: PRM visualization from `target list + mzXML`

Both scripts export one annotated PDF spectrum per selected scan and keep the
same visual layout (top peptide cleavage map + lower MS2 peak panel).

## Environment

```bash
cd /home/Share/Codex_CRC/MS/20260520_vis
rtk mamba env create -f environment.yml
rtk mamba activate dda-vis
```

Or run directly:

```bash
rtk mamba run -n dda-vis python dda_vis.py --help
rtk mamba run -n dda-vis python prm_vis.py --help
```

## DDA Workflow (`dda_vis.py`)

### Input

- Raw spectra: `mzXML`
- Search result: `pepXML` (MS-GF+ style)
- Target peptide sequence: exact match on unmodified peptide

### Example

```bash
rtk mamba run -n dda-vis python dda_vis.py \
  --pepxml 20240401_SW620_kit_DDA_1.pepXML \
  --mzxml 20240401_SW620_kit_DDA_1.mzXML \
  --peptide IHFISPNIYCCGAGTAADTDMTTQLISSNLELHSLSTGR \
  --outdir outputs
```

### DDA Key options

- `--max-qvalue 0.01` (default): keep `QValue <= 0.01`; disable via `--max-qvalue None`
- `--frag-tol-ppm 20` (default): fragment matching tolerance
- `--max-labels-per-series 20` (default): cap b/y labels in lower panel
- `--topmap-min-rel-int 0.02` (default): minimum relative intensity for top-map b/y labels
- `--intensity-scale {absolute,relative}` (default `absolute`)
- `--sample-id <text>`: override sample ID in output filename
- `--annotate-neutral-losses` / `--no-annotate-neutral-losses` (default enabled)

## PRM Workflow (`prm_vis.py`)

### Core behavior

- Match MS2 scans to target precursor m/z + charge
- Select one representative scan per target:
  - scan with highest `basePeakIntensity`
- Render with the same style as `dda_vis.py`

### `--prm-list` input modes

`--prm-list` supports both:

1. **File path** (TSV)
2. **Direct peptide text**

#### Mode 1: TSV file path

Example:

```bash
rtk mamba run -n dda-vis python prm_vis.py \
  --prm-list /home/Share/Codex_CRC/MS/20260521_PRM/PRM_list.txt \
  --mzxml /home/Share/Codex_CRC/MS/20260521_PRM/HCT116_PRM_1.mzXML \
  --outdir prm_outputs \
  --sample-id HCT116_PRM_1
```

Supported TSV header patterns:

- Standard: `gene`, `peptide_sequence`, `charge`
- Flexible:
  - peptide column can be `peptide_sequence` / `peptide` / `sequence`
  - charge column can be `charge` / `z`
  - ID column can be `gene` / `peptide_id` / first column

If your first column is peptide ID, it is supported.

#### Mode 2: direct peptide text

With explicit charge:

```bash
rtk mamba run -n dda-vis python prm_vis.py \
  --prm-list "ELTYPQQQLRDDDVGELGR/2,QVEWGAQLWVLYAGVERPVSR/3" \
  --mzxml /home/Share/Codex_CRC/MS/20260521_PRM/HCT116_PRM_1.mzXML \
  --outdir prm_outputs
```

Without explicit charge:

```bash
rtk mamba run -n dda-vis python prm_vis.py \
  --prm-list "ELTYPQQQLRDDDVGELGR QVEWGAQLWVLYAGVERPVSR" \
  --direct-charges 2,3,4 \
  --mzxml /home/Share/Codex_CRC/MS/20260521_PRM/HCT116_PRM_1.mzXML \
  --outdir prm_outputs
```

### PRM Key options

- `--precursor-tol-ppm 10` (default): precursor matching tolerance
- `--isotope-errors 0,1` (default): allow monoisotopic and +1 isotope targeting
- `--direct-charges 2,3,4` (default): used only for direct peptide mode without `/z`
- `--frag-tol-ppm 20` (default): fragment matching tolerance
- `--max-labels-per-series 20` (default)
- `--topmap-min-rel-int 0.02` (default)
- `--intensity-scale {absolute,relative}` (default `absolute`)
- `--annotate-neutral-losses` / `--no-annotate-neutral-losses` (default enabled)

## Figure/Output Details

- Output file pattern:
  - `protein_or_id-sample_id-z{charge}-peptide-scan_id.pdf`
- Figure content:
  - lower panel: MS2 peaks and matched b/y annotations
  - upper panel: peptide cleavage map with b/y labels
  - neutral loss labels:
    - `o` for `-H2O`
    - `*` for `-NH3`
- Lower panel does not repeat peptide sequence title.
- Peak-to-label connector lines are removed.

## Illustrator Compatibility

- PDF export uses:
  - `matplotlib.rcParams["pdf.fonttype"] = 42`
  - `matplotlib.rcParams["ps.fonttype"] = 42`
- This keeps text as editable text objects in Illustrator (instead of paths) in normal cases.
