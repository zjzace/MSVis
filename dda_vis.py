#!/usr/bin/env python3
"""DDA spectrum visualization for a target peptide from mzXML + pepXML.

Output one PDF per scan, with precursor metadata and b/y ion annotations.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["text.usetex"] = False
import matplotlib.pyplot as plt
import numpy as np
from pyteomics import mzxml, pepxml

PROTON = 1.00727646688
H2O = 18.0105646837
NH3 = 17.026549101

AA_MASS = {
    "A": 71.037113805,
    "R": 156.10111105,
    "N": 114.04292747,
    "D": 115.026943065,
    "C": 103.009184505,
    "E": 129.042593135,
    "Q": 128.05857754,
    "G": 57.021463735,
    "H": 137.058911875,
    "I": 113.084064015,
    "L": 113.084064015,
    "K": 128.09496305,
    "M": 131.040484645,
    "F": 147.068413945,
    "P": 97.052763875,
    "S": 87.032028435,
    "T": 101.047678505,
    "W": 186.07931298,
    "Y": 163.063328575,
    "V": 99.068413945,
}

B_COLOR = "#1C75BC"
Y_COLOR = "#BE1E2D"
BY_COLOR = "#7F3F98"


@dataclass
class PSM:
    scan_id: int
    hit_rank: int
    charge: int
    precursor_neutral_mass: float
    precursor_mz: float
    peptide: str
    modified_peptide: str
    protein_id: str
    qvalue: Optional[float]
    pepqvalue: Optional[float]
    modifications: List[Dict[str, float]]


def parse_optional_float(value: str) -> Optional[float]:
    text = value.strip().lower()
    if text in {"none", "null", "na", "nan"}:
        return None
    return float(value)


def safe_text(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("._-") or "NA"


def get_qvalue(hit: Dict) -> Optional[float]:
    scores = hit.get("search_score", {}) or {}
    qv = scores.get("QValue")
    return float(qv) if qv is not None else None


def get_pepqvalue(hit: Dict) -> Optional[float]:
    scores = hit.get("search_score", {}) or {}
    pqv = scores.get("PepQValue")
    return float(pqv) if pqv is not None else None


def normalized_sample_key(path: Path) -> str:
    name = path.name.lower()
    return re.sub(r"\.(pepxml|mzxml)$", "", name)


def pair_sample_files(pepxml_paths: List[Path], mzxml_paths: List[Path]) -> List[Tuple[str, Path, Path]]:
    if len(pepxml_paths) != len(mzxml_paths):
        raise ValueError(
            f"--pepxml and --mzxml must contain the same number of files "
            f"(got {len(pepxml_paths)} and {len(mzxml_paths)})."
        )

    pep_map: Dict[str, Path] = {}
    mz_map: Dict[str, Path] = {}

    for p in pepxml_paths:
        key = normalized_sample_key(p)
        if key in pep_map:
            raise ValueError(f"Duplicate pepXML basename after normalization: {p} and {pep_map[key]}")
        pep_map[key] = p

    for p in mzxml_paths:
        key = normalized_sample_key(p)
        if key in mz_map:
            raise ValueError(f"Duplicate mzXML basename after normalization: {p} and {mz_map[key]}")
        mz_map[key] = p

    pep_keys = set(pep_map.keys())
    mz_keys = set(mz_map.keys())
    if pep_keys != mz_keys:
        only_pep = sorted(pep_keys - mz_keys)
        only_mz = sorted(mz_keys - pep_keys)
        raise ValueError(
            "basename matching failed between pepXML and mzXML inputs. "
            f"Only in pepXML: {only_pep or '[]'}; only in mzXML: {only_mz or '[]'}"
        )

    return [(k, pep_map[k], mz_map[k]) for k in sorted(pep_keys)]


def score_passes(
    hit: Dict,
    score_type: str,
    qvalue_threshold: float,
    pepqvalue_threshold: float,
) -> bool:
    if score_type == "none":
        return True
    if score_type in {"qvalue", "both"}:
        qv = get_qvalue(hit)
        if qv is None or qv > qvalue_threshold:
            return False
    if score_type in {"pepqvalue", "both"}:
        pqv = get_pepqvalue(hit)
        if pqv is None or pqv > pepqvalue_threshold:
            return False
    return True


def select_psms(
    pepxml_path: Path,
    target_peptide: str,
    score_type: str,
    qvalue_threshold: float,
    pepqvalue_threshold: float,
    max_hit_rank: Optional[int],
) -> List[PSM]:
    selected: List[PSM] = []

    with pepxml.read(str(pepxml_path)) as reader:
        for sq in reader:
            scan_id = int(sq.get("start_scan"))
            charge = int(sq.get("assumed_charge"))
            neutral_mass = float(sq.get("precursor_neutral_mass"))
            precursor_mz = (neutral_mass + charge * PROTON) / charge

            for hit in sq.get("search_hit", []):
                hit_rank = int(hit.get("hit_rank", 9999))
                if max_hit_rank is not None and hit_rank > max_hit_rank:
                    continue
                if hit.get("peptide") != target_peptide:
                    continue

                if not score_passes(hit, score_type, qvalue_threshold, pepqvalue_threshold):
                    continue

                proteins = hit.get("proteins") or []
                protein_id = proteins[0].get("protein", "NA") if proteins else "NA"
                modifications = hit.get("modifications") or []
                qv = get_qvalue(hit)
                pqv = get_pepqvalue(hit)

                selected.append(
                    PSM(
                    scan_id=scan_id,
                    hit_rank=hit_rank,
                    charge=charge,
                    precursor_neutral_mass=neutral_mass,
                    precursor_mz=precursor_mz,
                    peptide=hit.get("peptide", target_peptide),
                    modified_peptide=hit.get("modified_peptide", target_peptide),
                    protein_id=protein_id,
                    qvalue=qv,
                    pepqvalue=pqv,
                    modifications=modifications,
                )
                )

    selected.sort(key=lambda x: (x.scan_id, x.hit_rank))
    return selected


def residue_masses_with_mods(peptide: str, modifications: List[Dict[str, float]]) -> Tuple[List[float], float, float]:
    masses = []
    for aa in peptide:
        if aa not in AA_MASS:
            raise ValueError(f"Unsupported amino acid '{aa}' in peptide {peptide}")
        masses.append(AA_MASS[aa])

    nterm_delta = 0.0
    cterm_delta = 0.0

    for mod in modifications:
        pos = int(mod["position"])
        mod_mass = float(mod["mass"])
        if 1 <= pos <= len(peptide):
            masses[pos - 1] = mod_mass
        elif pos == 0:
            # pepXML stores absolute modified N-terminus mass.
            nterm_delta = mod_mass - 1.00782503223
        elif pos == len(peptide) + 1:
            # pepXML stores absolute modified C-terminus mass.
            cterm_delta = mod_mass - 17.00273965163

    return masses, nterm_delta, cterm_delta


def theoretical_by_ions(
    peptide: str,
    modifications: List[Dict[str, float]],
    max_frag_charge: int = 2,
    annotate_neutral_losses: bool = True,
) -> List[Tuple[str, float, str]]:
    if len(peptide) < 2:
        return []

    masses, nterm_delta, cterm_delta = residue_masses_with_mods(peptide, modifications)
    prefix = np.cumsum(masses)
    total = prefix[-1]

    ions: List[Tuple[str, float, str]] = []
    n = len(peptide)
    for i in range(1, n):
        b_neutral = prefix[i - 1] + nterm_delta
        y_len = n - i
        y_neutral = (total - prefix[i - 1]) + H2O + cterm_delta

        for z in range(1, max_frag_charge + 1):
            b_mz = (b_neutral + z * PROTON) / z
            y_mz = (y_neutral + z * PROTON) / z
            ions.append((f"b{i}^{z}+", b_mz, "b"))
            ions.append((f"y{y_len}^{z}+", y_mz, "y"))

            # Keep neutral-loss annotation limited to 1+ to reduce label crowding.
            if annotate_neutral_losses and z == 1:
                b_h2o_mz = ((b_neutral - H2O) + PROTON)
                b_nh3_mz = ((b_neutral - NH3) + PROTON)
                y_h2o_mz = ((y_neutral - H2O) + PROTON)
                y_nh3_mz = ((y_neutral - NH3) + PROTON)

                if b_neutral > H2O:
                    ions.append((f"b{i}o^{z}+", b_h2o_mz, "b"))
                if b_neutral > NH3:
                    ions.append((f"b{i}*^{z}+", b_nh3_mz, "b"))
                if y_neutral > H2O:
                    ions.append((f"y{y_len}o^{z}+", y_h2o_mz, "y"))
                if y_neutral > NH3:
                    ions.append((f"y{y_len}*^{z}+", y_nh3_mz, "y"))

    return ions


def ppm_error(obs: float, theo: float) -> float:
    return (obs - theo) / theo * 1e6


def match_theoretical_ions(
    mz_array: np.ndarray,
    int_array: np.ndarray,
    theoretical: Iterable[Tuple[str, float, str]],
    frag_tol: float,
    frag_tol_unit: str,
) -> List[Dict]:
    matches: List[Dict] = []
    if mz_array.size == 0:
        return matches

    for ion_name, ion_mz, ion_series in theoretical:
        if frag_tol_unit == "ppm":
            tol = ion_mz * frag_tol * 1e-6
        else:
            tol = frag_tol
        left = np.searchsorted(mz_array, ion_mz - tol, side="left")
        right = np.searchsorted(mz_array, ion_mz + tol, side="right")
        if left >= right:
            continue

        local_idx = np.argmax(int_array[left:right])
        peak_idx = int(left + local_idx)
        obs_mz = float(mz_array[peak_idx])
        err_ppm = ppm_error(obs_mz, ion_mz)
        matches.append(
            {
                "ion": ion_name,
                "series": ion_series,
                "theoretical_mz": float(ion_mz),
                "peak_idx": peak_idx,
                "obs_mz": obs_mz,
                "intensity": float(int_array[peak_idx]),
                "ppm_error": float(err_ppm),
            }
        )

    return matches


def parse_ion_name(ion_name: str) -> Tuple[str, int, int, str]:
    m = re.fullmatch(r"([by])(\d+)([o*]?)\^(\d+)\+", ion_name)
    if not m:
        raise ValueError(f"Unexpected ion name format: {ion_name}")
    nl = m.group(3) if m.group(3) else ""
    return m.group(1), int(m.group(2)), int(m.group(4)), nl


def ion_label(ion_name: str) -> str:
    series, idx, charge, nl = parse_ion_name(ion_name)
    return f"{series}{idx}{nl}{charge}+"


def ion_label_parts(ion_name: str) -> Tuple[str, str]:
    series, idx, charge, nl = parse_ion_name(ion_name)
    return f"{series}{idx}{nl}", f"{charge}+"


def draw_charge_label(
    ax: plt.Axes,
    x: float,
    y: float,
    ion_name: str,
    color: str,
    ha: str,
    base_size: float = 8.0,
    sup_size: float = 6.2,
    sup_dy: float = 0.048,
    sup_dx_scale: float = 0.115,
    sup_dx_offset: float = 0.0,
) -> None:
    base, sup = ion_label_parts(ion_name)
    ax.text(
        x,
        y,
        base,
        fontsize=base_size,
        ha=ha,
        va="baseline",
        color=color,
        fontfamily="DejaVu Sans",
    )
    if ha == "left":
        sup_x = x + sup_dx_scale * len(base) + sup_dx_offset
    elif ha == "right":
        sup_x = x + 0.006 + sup_dx_offset
    else:
        sup_x = x + 0.08 + sup_dx_offset
    ax.text(
        sup_x,
        y + sup_dy,
        sup,
        fontsize=sup_size,
        ha="left",
        va="baseline",
        color=color,
        fontfamily="DejaVu Sans",
    )


def draw_charge_label_centered(
    ax: plt.Axes,
    x: float,
    y: float,
    ion_name: str,
    color: str,
    base_size: float = 9.0,
    sup_size: float = 7.0,
    sup_raise_points: float = 5.0,
    sup_right_points: float = 4.0,
) -> None:
    # Draw a compact "base + superscript charge" label centered over (x, y).
    base, sup = ion_label_parts(ion_name)
    base_w = base_size * 0.58 * len(base)
    sup_w = sup_size * 0.58 * len(sup)
    total_w = base_w + sup_w
    start_x = -0.5 * total_w

    ax.annotate(
        base,
        xy=(x, y),
        xytext=(start_x, 0),
        textcoords="offset points",
        ha="left",
        va="bottom",
        fontsize=base_size,
        color=color,
        fontfamily="DejaVu Sans",
    )
    ax.annotate(
        sup,
        xy=(x, y),
        xytext=(start_x + base_w + sup_right_points, sup_raise_points),
        textcoords="offset points",
        ha="left",
        va="bottom",
        fontsize=sup_size,
        color=color,
        fontfamily="DejaVu Sans",
    )


def pick_peak_labels(
    matches: List[Dict],
    mz_array: np.ndarray,
    max_labels_per_series: int,
    min_mz_gap: float,
) -> List[Dict]:
    selected: List[Dict] = []
    for series in ("b", "y"):
        candidates = [m for m in matches if m["series"] == series]
        candidates.sort(key=lambda x: x["intensity"], reverse=True)
        picked: List[Dict] = []
        used_mz: List[float] = []
        for m in candidates:
            mz = float(mz_array[m["peak_idx"]])
            if any(abs(mz - x) < min_mz_gap for x in used_mz):
                continue
            picked.append(m)
            used_mz.append(mz)
            if len(picked) >= max_labels_per_series:
                break
        selected.extend(picked)

    # Merge labels that land on the same peak.
    by_peak: Dict[int, List[Dict]] = {}
    for m in selected:
        by_peak.setdefault(m["peak_idx"], []).append(m)

    merged: List[Dict] = []
    for peak_idx, rows in by_peak.items():
        rows.sort(key=lambda x: (x["series"], -x["intensity"]))
        merged.append({"peak_idx": peak_idx, "items": rows})
    merged.sort(key=lambda x: x["peak_idx"])
    return merged


def select_cleavage_top_labels(
    matches: List[Dict],
    peptide_len: int,
    base_peak_intensity: float,
    min_rel_intensity: float,
) -> Dict[int, Dict[str, Dict]]:
    # cleavage index i means break between residues i and i+1 (1-based), 1..len-1
    cleavage_map: Dict[int, Dict[str, Dict]] = {i: {} for i in range(1, peptide_len)}
    if peptide_len < 2 or base_peak_intensity <= 0:
        return cleavage_map

    abs_min_int = base_peak_intensity * min_rel_intensity
    for m in matches:
        if m["intensity"] < abs_min_int:
            continue
        series, idx, _charge, nl = parse_ion_name(m["ion"])
        if nl:
            # Top backbone map should only reflect primary b/y ions.
            continue
        cleavage = idx if series == "b" else peptide_len - idx
        if cleavage < 1 or cleavage >= peptide_len:
            continue
        current = cleavage_map[cleavage].get(series)
        if current is None:
            cleavage_map[cleavage][series] = m
            continue
        cur_err = abs(current["ppm_error"])
        new_err = abs(m["ppm_error"])
        if new_err < cur_err or (new_err == cur_err and m["intensity"] > current["intensity"]):
            cleavage_map[cleavage][series] = m
    return cleavage_map


def draw_top_backbone(ax: plt.Axes, peptide: str, cleavage_labels: Dict[int, Dict[str, Dict]]) -> None:
    n = len(peptide)
    if n < 2:
        return

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    # Sequence letters
    for i, aa in enumerate(peptide):
        ax.text(
            i,
            0.50,
            aa,
            ha="center",
            va="center",
            fontsize=18,
            color="black",
            fontfamily="DejaVu Sans",
            fontweight="bold",
        )

    # Cleavage markers and labels in the legacy style:
    # b ions below, y ions above the peptide sequence.
    for i in range(1, n):
        x = i - 0.5
        cell = cleavage_labels[i]

        if "y" in cell:
            y_line = 0.78
            h_len = 0.68
            ax.plot([x, x], [0.58, y_line], color=Y_COLOR, linewidth=0.8)
            ax.plot([x, x + h_len], [y_line, y_line], color=Y_COLOR, linewidth=0.8)
            draw_charge_label(
                ax=ax,
                x=x + 0.018,
                y=y_line - 0.086,
                ion_name=cell["y"]["ion"],
                color=Y_COLOR,
                ha="left",
                base_size=8.0,
                sup_size=6.0,
                sup_dy=0.024,
                sup_dx_scale=0.134,
                sup_dx_offset=0.058,
            )

        if "b" in cell:
            b_line = 0.22
            h_len = 0.68
            ax.plot([x, x], [0.42, b_line], color=B_COLOR, linewidth=0.8)
            ax.plot([x - h_len, x], [b_line, b_line], color=B_COLOR, linewidth=0.8)
            draw_charge_label(
                ax=ax,
                x=x - 0.300,
                y=b_line + 0.034,
                ion_name=cell["b"]["ion"],
                color=B_COLOR,
                ha="right",
                base_size=8.0,
                sup_size=6.0,
                sup_dy=0.030,
                sup_dx_scale=0.106,
                sup_dx_offset=0.000,
            )


def render_spectrum(
    psm: PSM,
    mz_array: np.ndarray,
    int_array: np.ndarray,
    sample_id: str,
    outdir: Path,
    frag_tol: float,
    frag_tol_unit: str,
    max_labels_per_series: int,
    topmap_min_rel_int: float,
    intensity_scale: str,
    annotate_neutral_losses: bool,
) -> Path:
    theoretical = theoretical_by_ions(
        psm.peptide,
        psm.modifications,
        max_frag_charge=2,
        annotate_neutral_losses=annotate_neutral_losses,
    )
    matches = match_theoretical_ions(mz_array, int_array, theoretical, frag_tol, frag_tol_unit)

    peak_series: Dict[int, str] = {}
    for m in matches:
        idx = m["peak_idx"]
        if idx not in peak_series:
            peak_series[idx] = m["series"]
        elif peak_series[idx] != m["series"]:
            peak_series[idx] = "by"

    ymax_raw = float(int_array.max()) if int_array.size else 1.0

    if intensity_scale == "relative":
        plot_int = (int_array / ymax_raw) * 100.0 if ymax_raw > 0 else int_array
        y_label = "Relative ion abundance (%)"
        y_max = 106.0
    else:
        plot_int = int_array
        y_label = "Ion intensity"
        y_max = ymax_raw * 1.06 if ymax_raw > 0 else 1.0

    labels = pick_peak_labels(matches, mz_array, max_labels_per_series=max_labels_per_series, min_mz_gap=12.0)
    cleavage_labels = select_cleavage_top_labels(
        matches=matches,
        peptide_len=len(psm.peptide),
        base_peak_intensity=ymax_raw,
        min_rel_intensity=topmap_min_rel_int,
    )
    fig = plt.figure(figsize=(14, 7.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.5, 5], hspace=0.02)
    ax_top = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1])
    draw_top_backbone(ax_top, psm.peptide, cleavage_labels)

    ax.vlines(mz_array, 0, plot_int, color="#8F8F8F", linewidth=0.5, alpha=0.5)
    if matches:
        b_idx = [m["peak_idx"] for m in matches if peak_series[m["peak_idx"]] == "b"]
        y_idx = [m["peak_idx"] for m in matches if peak_series[m["peak_idx"]] == "y"]
        by_idx = [m["peak_idx"] for m in matches if peak_series[m["peak_idx"]] == "by"]

        if b_idx:
            ax.vlines(mz_array[b_idx], 0, plot_int[b_idx], color=B_COLOR, linewidth=0.9, alpha=0.95, label="b ions")
        if y_idx:
            ax.vlines(mz_array[y_idx], 0, plot_int[y_idx], color=Y_COLOR, linewidth=0.9, alpha=0.95, label="y ions")
        if by_idx:
            ax.vlines(mz_array[by_idx], 0, plot_int[by_idx], color=BY_COLOR, linewidth=0.9, alpha=0.95, label="b/y overlap")

    for j, row in enumerate(labels):
        idx = row["peak_idx"]
        mz = float(mz_array[idx])
        inten = float(plot_int[idx])
        label_items = row["items"]
        label_text = "/".join(ion_label(item["ion"]) for item in label_items)
        series_set = {item["series"] for item in row["items"]}
        if series_set == {"b"}:
            color = B_COLOR
        elif series_set == {"y"}:
            color = Y_COLOR
        else:
            color = BY_COLOR
        yoff = 2.5 + (j % 3) * 2.0
        label_y = min(inten + yoff, y_max * 0.98)

        if len(label_items) == 1:
            # Keep charge as superscript style in intensity labels, consistent with top map.
            draw_charge_label_centered(
                ax=ax,
                x=mz,
                y=label_y,
                ion_name=label_items[0]["ion"],
                color=color,
                base_size=9.0,
                sup_size=7.0,
                sup_raise_points=5.0,
                sup_right_points=1.0,
            )
        else:
            ax.text(
                mz,
                label_y,
                label_text,
                fontsize=9,
                ha="center",
                va="bottom",
                color=color,
                fontfamily="DejaVu Sans",
            )

    ax.set_xlabel("m/z")
    ax.set_ylabel(y_label)
    ax.set_ylim(0, y_max)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(width=0.6, labelsize=9)

    meta = f"scan_id={psm.scan_id} | precursor_m/z={psm.precursor_mz:.6f} | z={psm.charge}"
    ax.text(0.01, 0.98, meta, transform=ax.transAxes, va="top", ha="left", fontsize=10)

    if psm.qvalue is not None:
        ax.text(0.01, 0.92, f"QValue={psm.qvalue:.4g}", transform=ax.transAxes, va="top", ha="left", fontsize=9)
    if psm.pepqvalue is not None:
        pepq_y = 0.86 if psm.qvalue is not None else 0.92
        ax.text(0.01, pepq_y, f"PepQValue={psm.pepqvalue:.4g}", transform=ax.transAxes, va="top", ha="left", fontsize=9)

    if matches:
        ax.legend(loc="upper right", frameon=False, fontsize=9)

    fig.subplots_adjust(top=0.93, bottom=0.08, left=0.07, right=0.99)

    fname = (
        f"{safe_text(psm.protein_id)}-"
        f"{safe_text(sample_id)}-"
        f"{safe_text(psm.peptide)}-"
        f"{psm.scan_id}-rank{psm.hit_rank}.pdf"
    )
    outpath = outdir / fname
    fig.savefig(outpath)
    plt.close(fig)
    return outpath


def parse_mzxml_and_render(
    mzxml_path: Path,
    psms: List[PSM],
    sample_id: str,
    outdir: Path,
    frag_tol: float,
    frag_tol_unit: str,
    max_labels_per_series: int,
    topmap_min_rel_int: float,
    intensity_scale: str,
    annotate_neutral_losses: bool,
) -> Tuple[int, int, List[int], List[Path]]:
    psms_by_scan: Dict[int, List[PSM]] = {}
    for psm in psms:
        psms_by_scan.setdefault(psm.scan_id, []).append(psm)
    wanted = set(psms_by_scan.keys())
    found_scans: set[int] = set()
    outputs: List[Path] = []

    with mzxml.read(str(mzxml_path)) as reader:
        for scan in reader:
            scan_id = int(scan.get("num", -1))
            if scan_id not in wanted:
                continue

            ms_level = int(scan.get("msLevel", 0))
            if ms_level != 2:
                continue

            mz_array = np.asarray(scan.get("m/z array", []), dtype=float)
            int_array = np.asarray(scan.get("intensity array", []), dtype=float)
            if mz_array.size == 0 or int_array.size == 0:
                continue

            order = np.argsort(mz_array)
            mz_array = mz_array[order]
            int_array = int_array[order]

            for psm in psms_by_scan[scan_id]:
                out_pdf = render_spectrum(
                    psm=psm,
                    mz_array=mz_array,
                    int_array=int_array,
                    sample_id=sample_id,
                    outdir=outdir,
                    frag_tol=frag_tol,
                    frag_tol_unit=frag_tol_unit,
                    max_labels_per_series=max_labels_per_series,
                    topmap_min_rel_int=topmap_min_rel_int,
                    intensity_scale=intensity_scale,
                    annotate_neutral_losses=annotate_neutral_losses,
                )
                outputs.append(out_pdf)
            found_scans.add(scan_id)

            if found_scans == wanted:
                break

    missing = sorted(wanted - found_scans)
    return len(wanted), len(found_scans), missing, outputs


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Visualize DDA spectra for a target peptide from pepXML and mzXML.")
    p.add_argument("--pepxml", required=True, nargs="+", type=Path, help="One or more MS-GF+ pepXML files")
    p.add_argument("--mzxml", required=True, nargs="+", type=Path, help="One or more mzXML files")
    p.add_argument("--peptide", required=True, help="Target unmodified peptide sequence for exact matching")
    p.add_argument("--outdir", type=Path, default=Path("dda_vis_output"), help="Output directory for PDFs")
    p.add_argument(
        "--sample-id",
        default=None,
        help="Override sample ID in output filename for single-sample runs only.",
    )
    p.add_argument(
        "--score-type",
        choices=("qvalue", "pepqvalue", "both", "none"),
        default="qvalue",
        help="Score filter field(s). Default: qvalue",
    )
    p.add_argument(
        "--qvalue-threshold",
        type=float,
        default=0.1,
        help="QValue threshold. Used when --score-type is qvalue or both. Default: 0.1",
    )
    p.add_argument(
        "--pepqvalue-threshold",
        type=float,
        default=0.1,
        help="PepQValue threshold. Used when --score-type is pepqvalue or both. Default: 0.1",
    )
    p.add_argument(
        "--max-hit-rank",
        type=int,
        default=None,
        help="Maximum hit rank to keep. Default: None (keep all ranks).",
    )
    p.add_argument("--frag-tol", type=float, default=20.0, help="Fragment ion match tolerance value. Default: 20")
    p.add_argument(
        "--frag-tol-unit",
        choices=("ppm", "da"),
        default="ppm",
        help="Fragment ion tolerance unit. Default: ppm. Use 'da' for TPP-like matching.",
    )
    p.add_argument(
        "--frag-tol-ppm",
        type=float,
        default=None,
        help="Deprecated alias of --frag-tol when --frag-tol-unit ppm.",
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
    if args.frag_tol_ppm is not None:
        args.frag_tol = args.frag_tol_ppm
        args.frag_tol_unit = "ppm"
    args.outdir.mkdir(parents=True, exist_ok=True)
    sample_pairs = pair_sample_files(args.pepxml, args.mzxml)

    print(f"Target peptide: {args.peptide}")
    print(f"Sample pairs by basename: {len(sample_pairs)}")
    grand_outputs = 0

    for sample_key, pepxml_path, mzxml_path in sample_pairs:
        sample_id = args.sample_id if (args.sample_id and len(sample_pairs) == 1) else mzxml_path.stem
        sample_outdir = args.outdir

        psms = select_psms(
            pepxml_path=pepxml_path,
            target_peptide=args.peptide,
            score_type=args.score_type,
            qvalue_threshold=args.qvalue_threshold,
            pepqvalue_threshold=args.pepqvalue_threshold,
            max_hit_rank=args.max_hit_rank,
        )
        if not psms:
            print(f"[{sample_key}] No matching PSMs found with current filters.")
            continue

        unique_scans = len({x.scan_id for x in psms})
        total_hits = len(psms)
        found, outputs, missing = 0, [], []
        total, found, missing, outputs = parse_mzxml_and_render(
            mzxml_path=mzxml_path,
            psms=psms,
            sample_id=sample_id,
            outdir=sample_outdir,
            frag_tol=args.frag_tol,
            frag_tol_unit=args.frag_tol_unit,
            max_labels_per_series=args.max_labels_per_series,
            topmap_min_rel_int=args.topmap_min_rel_int,
            intensity_scale=args.intensity_scale,
            annotate_neutral_losses=args.annotate_neutral_losses,
        )
        grand_outputs += len(outputs)

        print(f"[{sample_key}] pepXML={pepxml_path}")
        print(f"[{sample_key}] mzXML={mzxml_path}")
        print(f"[{sample_key}] Matched hits in pepXML: {total_hits} across {unique_scans} scan IDs")
        print(f"[{sample_key}] Scan IDs rendered from mzXML: {found}/{total}")
        print(f"[{sample_key}] Output PDFs: {len(outputs)}")
        if missing:
            preview = ", ".join(map(str, missing[:20]))
            tail = " ..." if len(missing) > 20 else ""
            print(f"[{sample_key}] Missing MS2 scans in mzXML: {len(missing)} [{preview}{tail}]")
        if outputs:
            print(f"[{sample_key}] Example output: {outputs[0]}")

    print(f"Output directory root: {args.outdir.resolve()}")
    print(f"Total PDFs generated: {grand_outputs}")


if __name__ == "__main__":
    main()
