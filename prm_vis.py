#!/usr/bin/env python3
"""PRM spectrum visualization from mzXML + PRM target list.

Output one PDF per PRM target (gene + peptide + charge), using the MS2 scan
with the highest basePeakIntensity among scans matched by precursor m/z and charge.
Visualization style is shared with dda_vis.py.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from pyteomics import mzxml

from dda_vis import AA_MASS, H2O, PROTON, PSM, ppm_error, render_spectrum

NEUTRON = 1.0033548378


@dataclass(frozen=True)
class PRMTarget:
    gene: str
    peptide: str
    charge: int
    precursor_mz: float

    @property
    def key(self) -> Tuple[str, str, int]:
        return (self.gene, self.peptide, self.charge)


@dataclass
class SelectedScan:
    target: PRMTarget
    scan_id: int
    precursor_mz_obs: float
    base_peak_intensity: float
    mz_array: np.ndarray
    int_array: np.ndarray


def peptide_neutral_mass(peptide: str) -> float:
    total = H2O
    for aa in peptide:
        if aa not in AA_MASS:
            raise ValueError(f"Unsupported amino acid '{aa}' in peptide {peptide}")
        total += AA_MASS[aa]
    return total


def precursor_mz_from_peptide(peptide: str, charge: int) -> float:
    neutral_mass = peptide_neutral_mass(peptide)
    return (neutral_mass + charge * PROTON) / charge


def parse_int_list(text: str, name: str) -> List[int]:
    values = [int(x.strip()) for x in text.split(",") if x.strip() != ""]
    if not values:
        raise ValueError(f"No valid {name} provided.")
    return values


def add_target(targets: List[PRMTarget], seen: set, gene: str, peptide: str, charge: int) -> None:
    gene_clean = gene.strip() if gene else peptide
    peptide_clean = peptide.strip().upper()
    if not gene_clean or not peptide_clean:
        return
    key = (gene_clean, peptide_clean, charge)
    if key in seen:
        return
    seen.add(key)
    targets.append(
        PRMTarget(
            gene=gene_clean,
            peptide=peptide_clean,
            charge=charge,
            precursor_mz=precursor_mz_from_peptide(peptide_clean, charge),
        )
    )


def parse_prm_list_file(prm_list_path: Path) -> List[PRMTarget]:
    targets: List[PRMTarget] = []
    seen = set()
    with prm_list_path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = [x.strip() for x in (reader.fieldnames or [])]
        lower = {x.lower(): x for x in fieldnames}

        pep_col = lower.get("peptide_sequence") or lower.get("peptide") or lower.get("sequence")
        charge_col = lower.get("charge") or lower.get("z")
        id_col = lower.get("gene") or lower.get("peptide_id") or (fieldnames[0] if fieldnames else None)

        if pep_col and charge_col and id_col:
            for row in reader:
                peptide = (row.get(pep_col) or "").strip()
                charge_text = (row.get(charge_col) or "").strip()
                item_id = (row.get(id_col) or "").strip()
                if not peptide or not charge_text:
                    continue
                charge = int(charge_text)
                add_target(targets, seen, item_id or peptide, peptide, charge)
            return targets

    # Fallback: parse as plain TSV (supports first column as peptide ID).
    with prm_list_path.open("r", newline="") as fh:
        raw = csv.reader(fh, delimiter="\t")
        for i, row in enumerate(raw, start=1):
            if len(row) < 3:
                continue
            item_id = row[0].strip()
            peptide = row[1].strip()
            charge_text = row[2].strip()
            if not item_id or not peptide or not charge_text:
                continue
            try:
                charge = int(charge_text)
            except ValueError:
                # likely header row in non-standard format
                if i == 1:
                    continue
                raise
            add_target(targets, seen, item_id, peptide, charge)
    return targets


def parse_prm_text_targets(prm_text: str, default_charges: List[int]) -> List[PRMTarget]:
    targets: List[PRMTarget] = []
    seen = set()
    chunks = [x.strip() for x in re.split(r"[,\s;]+", prm_text.strip()) if x.strip()]
    if not chunks:
        return targets

    for token in chunks:
        m = re.fullmatch(r"([A-Za-z]+)(?:[/:\-]z?(\d+))?$", token)
        if not m:
            raise ValueError(
                f"Invalid peptide token '{token}'. Use formats like PEPTIDE or PEPTIDE/2."
            )
        peptide = m.group(1).upper()
        charge_text = m.group(2)
        charges = [int(charge_text)] if charge_text else default_charges
        for z in charges:
            add_target(targets, seen, peptide, peptide, z)
    return targets


def parse_prm_input(prm_list_or_peptide: str, default_charges: List[int]) -> List[PRMTarget]:
    as_path = Path(prm_list_or_peptide)
    if as_path.exists() and as_path.is_file():
        return parse_prm_list_file(as_path)
    return parse_prm_text_targets(prm_list_or_peptide, default_charges)


def extract_precursor(scan: Dict) -> Tuple[Optional[float], Optional[int]]:
    prec = scan.get("precursorMz")
    if not isinstance(prec, list) or not prec:
        return None, None
    first = prec[0]
    try:
        mz = float(first.get("precursorMz"))
    except Exception:
        return None, None
    charge_raw = first.get("precursorCharge")
    charge = int(charge_raw) if charge_raw is not None else None
    return mz, charge


def best_target_for_scan(
    precursor_mz_obs: float,
    precursor_charge_obs: Optional[int],
    targets_by_charge: Dict[int, List[PRMTarget]],
    precursor_tol_ppm: float,
    isotope_errors: List[int],
) -> Optional[PRMTarget]:
    if precursor_charge_obs is not None and precursor_charge_obs in targets_by_charge:
        candidates = targets_by_charge[precursor_charge_obs]
    else:
        candidates = [t for rows in targets_by_charge.values() for t in rows]

    best: Optional[PRMTarget] = None
    best_abs_ppm = float("inf")
    for target in candidates:
        for iso in isotope_errors:
            theo = target.precursor_mz + (iso * NEUTRON / target.charge)
            err = ppm_error(precursor_mz_obs, theo)
            abs_err = abs(err)
            if abs_err > precursor_tol_ppm:
                continue
            if abs_err < best_abs_ppm:
                best_abs_ppm = abs_err
                best = target
    return best


def select_best_scans(
    mzxml_path: Path,
    targets: List[PRMTarget],
    precursor_tol_ppm: float,
    isotope_errors: List[int],
) -> Dict[Tuple[str, str, int], SelectedScan]:
    targets_by_charge: Dict[int, List[PRMTarget]] = {}
    for t in targets:
        targets_by_charge.setdefault(t.charge, []).append(t)

    selected: Dict[Tuple[str, str, int], SelectedScan] = {}
    with mzxml.read(str(mzxml_path)) as reader:
        for scan in reader:
            if int(scan.get("msLevel", 0)) != 2:
                continue

            precursor_mz_obs, precursor_charge_obs = extract_precursor(scan)
            if precursor_mz_obs is None:
                continue

            target = best_target_for_scan(
                precursor_mz_obs=precursor_mz_obs,
                precursor_charge_obs=precursor_charge_obs,
                targets_by_charge=targets_by_charge,
                precursor_tol_ppm=precursor_tol_ppm,
                isotope_errors=isotope_errors,
            )
            if target is None:
                continue

            mz_array = np.asarray(scan.get("m/z array", []), dtype=float)
            int_array = np.asarray(scan.get("intensity array", []), dtype=float)
            if mz_array.size == 0 or int_array.size == 0:
                continue
            order = np.argsort(mz_array)
            mz_array = mz_array[order]
            int_array = int_array[order]

            scan_id = int(scan.get("num", -1))
            base_peak_intensity = float(scan.get("basePeakIntensity", 0.0) or 0.0)
            old = selected.get(target.key)
            if old is None or base_peak_intensity > old.base_peak_intensity:
                selected[target.key] = SelectedScan(
                    target=target,
                    scan_id=scan_id,
                    precursor_mz_obs=precursor_mz_obs,
                    base_peak_intensity=base_peak_intensity,
                    mz_array=mz_array,
                    int_array=int_array,
                )
    return selected


def render_selected_scans(
    selected: Dict[Tuple[str, str, int], SelectedScan],
    sample_id: str,
    outdir: Path,
    frag_tol_ppm: float,
    max_labels_per_series: int,
    topmap_min_rel_int: float,
    intensity_scale: str,
    annotate_neutral_losses: bool,
) -> List[Path]:
    outputs: List[Path] = []
    for key in sorted(selected.keys()):
        row = selected[key]
        target = row.target
        psm = PSM(
            scan_id=row.scan_id,
            charge=target.charge,
            precursor_neutral_mass=peptide_neutral_mass(target.peptide),
            precursor_mz=row.precursor_mz_obs,
            peptide=target.peptide,
            modified_peptide=target.peptide,
            protein_id=target.gene,
            qvalue=None,
            modifications=[],
        )
        out_pdf = render_spectrum(
            psm=psm,
            mz_array=row.mz_array,
            int_array=row.int_array,
            sample_id=f"{sample_id}-z{target.charge}",
            outdir=outdir,
            frag_tol_ppm=frag_tol_ppm,
            max_labels_per_series=max_labels_per_series,
            topmap_min_rel_int=topmap_min_rel_int,
            intensity_scale=intensity_scale,
            annotate_neutral_losses=annotate_neutral_losses,
        )
        outputs.append(out_pdf)
    return outputs


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Visualize PRM spectra from target list and mzXML.")
    p.add_argument(
        "--prm-list",
        required=True,
        type=str,
        help=(
            "PRM input: (1) TSV file path, or (2) direct peptide string list, "
            "e.g. 'PEPTIDE/2,ANOTHERPEP/3' or 'PEPTIDE' (uses --direct-charges)."
        ),
    )
    p.add_argument(
        "--direct-charges",
        default="2,3,4",
        help="Default charges for direct peptide input when charge is omitted. Default: 2,3,4",
    )
    p.add_argument("--mzxml", required=True, type=Path, help="Path to PRM mzXML file")
    p.add_argument("--outdir", type=Path, default=Path("prm_vis_output"), help="Output directory for PDFs")
    p.add_argument(
        "--sample-id",
        default=None,
        help="Sample ID in output filename. Default: mzXML basename without extension.",
    )
    p.add_argument(
        "--precursor-tol-ppm",
        type=float,
        default=10.0,
        help="Precursor m/z tolerance for matching PRM targets (ppm). Default: 10",
    )
    p.add_argument(
        "--isotope-errors",
        default="0,1",
        help="Comma-separated isotope errors allowed for precursor matching. Default: 0,1",
    )
    p.add_argument(
        "--frag-tol-ppm",
        type=float,
        default=20.0,
        help="Fragment ion match tolerance in ppm. Default: 20",
    )
    p.add_argument(
        "--max-labels-per-series",
        type=int,
        default=20,
        help="Maximum number of text labels per ion series (b and y). Default: 20",
    )
    p.add_argument(
        "--topmap-min-rel-int",
        type=float,
        default=0.02,
        help="Minimum relative intensity for top backbone labels. Default: 0.02",
    )
    p.add_argument(
        "--intensity-scale",
        choices=("absolute", "relative"),
        default="absolute",
        help="Spectrum intensity scale: absolute or relative. Default: absolute",
    )
    nl_group = p.add_mutually_exclusive_group()
    nl_group.add_argument(
        "--annotate-neutral-losses",
        dest="annotate_neutral_losses",
        action="store_true",
        help="Annotate neutral-loss ions: H2O as 'o', NH3 as '*'. Default: enabled.",
    )
    nl_group.add_argument(
        "--no-annotate-neutral-losses",
        dest="annotate_neutral_losses",
        action="store_false",
        help="Disable annotation of neutral-loss ions.",
    )
    p.set_defaults(annotate_neutral_losses=True)
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    sample_id = args.sample_id if args.sample_id else args.mzxml.stem
    args.outdir.mkdir(parents=True, exist_ok=True)

    default_charges = parse_int_list(args.direct_charges, "direct charges")
    targets = parse_prm_input(args.prm_list, default_charges)
    if not targets:
        print("No valid PRM targets found in list.")
        return

    isotope_errors = parse_int_list(args.isotope_errors, "isotope errors")

    selected = select_best_scans(
        mzxml_path=args.mzxml,
        targets=targets,
        precursor_tol_ppm=args.precursor_tol_ppm,
        isotope_errors=isotope_errors,
    )
    outputs = render_selected_scans(
        selected=selected,
        sample_id=sample_id,
        outdir=args.outdir,
        frag_tol_ppm=args.frag_tol_ppm,
        max_labels_per_series=args.max_labels_per_series,
        topmap_min_rel_int=args.topmap_min_rel_int,
        intensity_scale=args.intensity_scale,
        annotate_neutral_losses=args.annotate_neutral_losses,
    )

    print(f"PRM targets in list: {len(targets)}")
    print(f"Rendered PDFs (one per matched target): {len(outputs)}")
    missing = [t for t in targets if t.key not in selected]
    if missing:
        print(f"Targets with no matched MS2 scan: {len(missing)}")
        preview = ", ".join(f"{x.gene}:{x.peptide}/z{x.charge}" for x in missing[:10])
        tail = " ..." if len(missing) > 10 else ""
        print(f"Missing preview: {preview}{tail}")
    print(f"Output directory: {args.outdir.resolve()}")
    if outputs:
        print("Example output file:")
        print(f"  {outputs[0]}")


if __name__ == "__main__":
    main()
