# DDA Peptide Spectrum Visualization

Generate one annotated MS2 spectrum PDF per scan for a target peptide using:
- Raw spectra: `mzXML`
- Identification: `pepXML` (MS-GF+)

## Environment (mamba)

```bash
cd /home/Share/Codex_CRC/MS/20260520_vis
rtk mamba env create -f environment.yml
rtk mamba activate dda-vis
```

Or run without activation:

```bash
rtk mamba run -n dda-vis python dda_vis.py --help
```

## Usage

```bash
rtk mamba run -n dda-vis python dda_vis.py \
  --pepxml 20240401_SW620_kit_DDA_1.pepXML \
  --mzxml 20240401_SW620_kit_DDA_1.mzXML \
  --peptide IHFISPNIYCCGAGTAADTDMTTQLISSNLELHSLSTGR \
  --outdir outputs
```

## Output

- One PDF per scan ID.
- Filename pattern:
  - `protein_id-sample_id-peptide_sequence-scan_id.pdf`
- Figure contains:
  - `scan_id`
  - precursor `m/z`
  - precursor charge `z`
  - matched b/y ion highlights and labels in spectrum panel (default 20 ppm)
  - no connector line from peak to text label
  - top peptide cleavage map:
    - bold peptide sequence
    - `y` cleavage lines/labels above (red), `b` cleavage lines/labels below (blue)
    - charge states shown as superscript in labels (example: `y7²⁺`)
  - lower panel does not repeat peptide sequence title

## Key options

- `--max-qvalue 0.01` (default): keep `QValue <= 0.01`
  - Disable filter: `--max-qvalue None`
- `--frag-tol-ppm 20` (default)
- `--max-labels-per-series 20` (default)
- `--topmap-min-rel-int 0.02` (default): min relative intensity for top backbone cleavage labels
- `--sample-id <text>`: override sample ID in filenames
- `--annotate-neutral-losses` / `--no-annotate-neutral-losses`:
  - annotate neutral-loss ions in peak labels (`o` for H2O loss, `*` for NH3 loss)
  - default is enabled

## Illustrator Compatibility

- PDF text is exported with `pdf.fonttype=42`, so labels are kept as editable text in Illustrator.
