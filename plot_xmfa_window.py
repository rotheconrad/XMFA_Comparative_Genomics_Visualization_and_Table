#!/usr/bin/env python3
import argparse
import gzip
import re
from collections import OrderedDict
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap, BoundaryNorm

# Illustrator-friendly output
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["savefig.transparent"] = False


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


# --------------------------------------------------
# XMFA sequence name normalization
# --------------------------------------------------

def normalize_sid(label: str) -> str:
    return label.split(".")[0]


# --------------------------------------------------
# GFF3 parsing
# --------------------------------------------------

def parse_gff3_features(
    gff_path: str,
    seqid_wanted: str,
    start: int,
    end: int,
    feature_types: Tuple[str, ...]
) -> List[dict]:
    feats = []
    with open_text(gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, ftype, fstart, fend, score, strand, phase, attrs = parts
            if seqid != seqid_wanted:
                continue
            if ftype not in feature_types:
                continue

            s = int(fstart)
            e = int(fend)
            if e < start or s > end:
                continue

            ad = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    ad[k] = v

            feats.append({
                "seqid": seqid,
                "type": ftype,
                "start": s,
                "end": e,
                "strand": strand,
                "attrs": ad
            })
    return feats


# --------------------------------------------------
# XMFA parsing
# --------------------------------------------------

XMFA_HEADER_RE = re.compile(
    r"^>\s*(?P<idx>\S+):(?P<start>\d+)-(?P<end>\d+)\s+(?P<strand>[+-])\s+(?P<label>\S+)$"
)


def read_xmfa_blocks(xmfa_path: str):
    block = {}
    cur_sid = None

    def flush():
        nonlocal block
        if not block:
            return None
        out = {}
        for sid, rec in block.items():
            out[sid] = {
                "start": rec["start"],
                "end": rec["end"],
                "strand": rec["strand"],
                "label": rec["label"],
                "seq": "".join(rec["seq_chunks"])
            }
        block = {}
        return out

    with open_text(xmfa_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                continue

            if line == "=":
                b = flush()
                if b is not None:
                    yield b
                cur_sid = None
                continue

            if line.startswith(">"):
                m = XMFA_HEADER_RE.match(line)
                if not m:
                    raise ValueError(f"Unrecognized XMFA header: {line}")

                label = m.group("label")
                sid = normalize_sid(label)

                block[sid] = {
                    "start": int(m.group("start")),
                    "end": int(m.group("end")),
                    "strand": m.group("strand"),
                    "label": label,
                    "seq_chunks": []
                }
                cur_sid = sid
                continue

            if cur_sid is None:
                continue

            block[cur_sid]["seq_chunks"].append(line.strip())

    b = flush()
    if b is not None:
        yield b


def ref_span_for_block(block: dict, ref_sid: str) -> Optional[Tuple[int, int]]:
    if ref_sid not in block:
        return None
    s = block[ref_sid]["start"]
    e = block[ref_sid]["end"]
    return (s, e) if s <= e else (e, s)


def slice_alignment_by_ref_window(
    block: dict,
    ref_sid: str,
    win_start: int,
    win_end: int,
    sids_order: List[str]
):
    if ref_sid not in block:
        return None, None

    ref_rec = block[ref_sid]
    ref_s = ref_rec["start"]
    ref_e = ref_rec["end"]
    ref_aln = ref_rec["seq"]

    ref_lo, ref_hi = (ref_s, ref_e) if ref_s <= ref_e else (ref_e, ref_s)
    if ref_hi < win_start or ref_lo > win_end:
        return None, None

    refpos_to_fullcol = {}
    ref_pos = ref_lo
    selected_cols = []

    for j, c in enumerate(ref_aln):
        if c not in "-.":
            refpos_to_fullcol[ref_pos] = j
            if win_start <= ref_pos <= win_end:
                selected_cols.append(j)
            ref_pos += 1

    if not selected_cols:
        return None, None

    j0, j1 = min(selected_cols), max(selected_cols)

    sliced = OrderedDict()
    for sid in sids_order:
        if sid in block:
            sliced[sid] = block[sid]["seq"][j0:j1 + 1]

    refpos_to_slicecol = {}
    for pos, fullcol in refpos_to_fullcol.items():
        if j0 <= fullcol <= j1:
            refpos_to_slicecol[pos] = fullcol - j0

    return sliced, refpos_to_slicecol


# --------------------------------------------------
# Alignment helpers
# --------------------------------------------------

BASE_MISMATCH_CODE = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5}


def find_snp_columns(aln_seqs: Dict[str, str], ref_sid: str) -> np.ndarray:
    ref = aln_seqs[ref_sid].upper()
    L = len(ref)
    snp = np.zeros(L, dtype=bool)

    for j in range(L):
        r = ref[j]
        if r in "-.":
            continue
        for sid, s in aln_seqs.items():
            if sid == ref_sid:
                continue
            c = s[j].upper()
            if c in "-." or c != r:
                snp[j] = True
                break
    return snp


def build_difference_matrix(aln: Dict[str, str], ref_sid: str, ordered: List[str]):
    ref = aln[ref_sid].upper()
    L = len(ref)

    M = np.zeros((len(ordered), L), dtype=np.int8)
    informative = np.zeros(L, dtype=bool)

    for i, sid in enumerate(ordered):
        seq = aln[sid].upper()
        for j, (r, c) in enumerate(zip(ref, seq)):
            if c in "-.":
                M[i, j] = 6
                if sid != ref_sid:
                    informative[j] = True
            elif sid == ref_sid:
                M[i, j] = 0
            elif r == c:
                M[i, j] = 0
            else:
                M[i, j] = BASE_MISMATCH_CODE.get(c, 5)
                if r not in "-.":
                    informative[j] = True

    return M, informative


def make_difference_colormap():
    colors = [
        "#d9d9d9",  # ref match
        "#4daf4a",  # A
        "#377eb8",  # C
        "#ff7f00",  # G
        "#e41a1c",  # T
        "#984ea3",  # N/other
        "#ffffff",  # gap
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5, 1), cmap.N)
    return cmap, norm


# --------------------------------------------------
# Feature helpers
# --------------------------------------------------

def feature_to_cols(feature_start: int, feature_end: int, refpos_to_col: dict):
    coords = [refpos_to_col[p] for p in range(feature_start, feature_end + 1) if p in refpos_to_col]
    if not coords:
        return None
    return min(coords), max(coords)


def get_repeat_color(feat: dict) -> str:
    text = ";".join(f"{k}={v}" for k, v in feat["attrs"].items()).lower()
    if "simple_repeat" in text or "motif:" in text:
        return "#bdbdbd"
    if "cacta" in text:
        return "#e41a1c"
    if "gypsy" in text:
        return "#377eb8"
    if "copia" in text:
        return "#4daf4a"
    if "dna" in text:
        return "#984ea3"
    return "#ff7f00"


def assign_gene_rows(gene_feats: List[dict], max_rows: int = 4):
    rows_end = [-1] * max_rows
    gene_row = {}

    genes = [g for g in gene_feats if g["type"] == "gene"]
    genes = sorted(genes, key=lambda x: (x["start"], x["end"]))

    for feat in genes:
        gid = feat["attrs"].get("ID", feat["attrs"].get("Name", "gene"))
        placed = False
        for i in range(max_rows):
            if feat["start"] > rows_end[i]:
                gene_row[gid] = i
                rows_end[i] = feat["end"]
                placed = True
                break
        if not placed:
            gene_row[gid] = max_rows - 1

    return gene_row


def pick_exon_parent(exon_feat: dict) -> str:
    parent = exon_feat["attrs"].get("Parent", "")
    return parent.split(",")[0] if parent else parent


def parent_matches_gene(parent: str, gene_id: str) -> bool:
    return parent == gene_id or parent.startswith(gene_id + ".") or parent.startswith(gene_id + "-")


# --------------------------------------------------
# Quantitative summaries
# --------------------------------------------------

def build_sample_summary(
    aln: Dict[str, str],
    ref_sid: str,
    ordered: List[str],
    diffM: np.ndarray,
    informative_mask: np.ndarray,
    snp_cols: np.ndarray,
    plotted_fullcols: np.ndarray,
    informative_only: bool,
    win_start: int,
    win_end: int
) -> pd.DataFrame:
    rows = []
    ref = aln[ref_sid].upper()

    fullcols_used = plotted_fullcols if informative_only else np.arange(diffM.shape[1])

    for sid in ordered:
        seq = aln[sid].upper()
        seq_used = "".join(seq[j] for j in fullcols_used)
        ref_used = "".join(ref[j] for j in fullcols_used)

        total_cols = len(seq_used)
        gap_cols = sum(c in "-." for c in seq_used)
        match_cols = 0
        mismatch_cols = 0
        A_mismatch = 0
        C_mismatch = 0
        G_mismatch = 0
        T_mismatch = 0
        N_other = 0

        for r, c in zip(ref_used, seq_used):
            if c in "-.":
                continue
            if sid == ref_sid or c == r:
                match_cols += 1
            else:
                mismatch_cols += 1
                if c == "A":
                    A_mismatch += 1
                elif c == "C":
                    C_mismatch += 1
                elif c == "G":
                    G_mismatch += 1
                elif c == "T":
                    T_mismatch += 1
                else:
                    N_other += 1

        rows.append({
            "sample": sid,
            "window_start": win_start,
            "window_end": win_end,
            "informative_only": informative_only,
            "displayed_columns": total_cols,
            "displayed_snp_columns": int(total_cols if informative_only else snp_cols.sum()),
            "displayed_informative_columns": int(total_cols if informative_only else informative_mask.sum()),
            "match_columns": match_cols,
            "mismatch_columns": mismatch_cols,
            "gap_columns": gap_cols,
            "A_mismatch": A_mismatch,
            "C_mismatch": C_mismatch,
            "G_mismatch": G_mismatch,
            "T_mismatch": T_mismatch,
            "N_or_other_mismatch": N_other
        })

    return pd.DataFrame(rows)


def build_column_summary(
    aln: Dict[str, str],
    ref_sid: str,
    ordered: List[str],
    plotted_fullcols: np.ndarray,
    informative_only: bool,
    refpos_to_col: dict
) -> pd.DataFrame:
    ref = aln[ref_sid].upper()
    col_to_refpos = {col: pos for pos, col in refpos_to_col.items()}

    rows = []
    for plot_idx, fullcol in enumerate(plotted_fullcols):
        ref_base = ref[fullcol].upper()
        row = {
            "plot_column_index": int(plot_idx),
            "full_alignment_column": int(fullcol),
            "reference_position": int(col_to_refpos[fullcol]) if fullcol in col_to_refpos else np.nan,
            "reference_base": ref_base,
            "informative_only_mode": informative_only
        }

        distinct = set()
        n_gap = 0
        n_match = 0
        n_mismatch = 0

        for sid in ordered:
            base = aln[sid][fullcol].upper()
            row[sid] = base
            if base in "-.":
                n_gap += 1
            elif base == ref_base:
                n_match += 1
                distinct.add(base)
            else:
                n_mismatch += 1
                distinct.add(base)

        row["n_gap"] = n_gap
        row["n_match_to_ref"] = n_match
        row["n_nonref"] = n_mismatch
        row["n_distinct_non_gap_states"] = len(distinct)
        rows.append(row)

    return pd.DataFrame(rows)


# --------------------------------------------------
# Legend
# --------------------------------------------------

def draw_legend_panel(ax_leg):
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(0, 1)
    ax_leg.axis("off")

    y = 0.97
    ax_leg.text(0.02, y, "Legend", fontsize=13, fontweight="bold", va="top")
    y -= 0.08

    ax_leg.text(0.02, y, "Gene track", fontsize=11, fontweight="bold", va="top")
    y -= 0.05
    ax_leg.plot([0.05, 0.42], [y, y], color="black", lw=1.2)
    ax_leg.add_patch(Rectangle((0.16, y - 0.018), 0.08, 0.036, color="black"))
    ax_leg.annotate(
        "",
        xy=(0.49, y),
        xytext=(0.42, y),
        arrowprops=dict(arrowstyle="-|>", lw=1.2, color="black", mutation_scale=11)
    )
    ax_leg.text(0.55, y, "gene backbone + exon box + strand arrow", fontsize=9, va="center")
    y -= 0.08

    ax_leg.text(0.02, y, "Repeat track", fontsize=11, fontweight="bold", va="top")
    y -= 0.05
    repeat_items = [
        ("simple repeat", "#bdbdbd"),
        ("CACTA", "#e41a1c"),
        ("Gypsy", "#377eb8"),
        ("Copia", "#4daf4a"),
        ("DNA transposon", "#984ea3"),
        ("other repeat", "#ff7f00"),
    ]
    for label, color in repeat_items:
        ax_leg.add_patch(Rectangle((0.05, y - 0.015), 0.06, 0.03, color=color))
        ax_leg.text(0.15, y, label, fontsize=9, va="center")
        y -= 0.045

    y -= 0.02
    ax_leg.text(0.02, y, "SNP track", fontsize=11, fontweight="bold", va="top")
    y -= 0.05
    ax_leg.vlines(0.08, y - 0.025, y + 0.025, color="black", linewidth=1.2)
    ax_leg.text(0.15, y, "polymorphic / informative column marker", fontsize=9, va="center")
    y -= 0.08

    ax_leg.text(0.02, y, "Alignment tiles", fontsize=11, fontweight="bold", va="top")
    y -= 0.05
    base_items = [
        ("reference match", "#d9d9d9"),
        ("A mismatch", "#4daf4a"),
        ("C mismatch", "#377eb8"),
        ("G mismatch", "#ff7f00"),
        ("T mismatch", "#e41a1c"),
        ("N / other", "#984ea3"),
        ("gap", "#ffffff"),
    ]
    for label, color in base_items:
        edge = "black" if color == "#ffffff" else "none"
        ax_leg.add_patch(Rectangle((0.05, y - 0.015), 0.06, 0.03, color=color, ec=edge, lw=0.5))
        ax_leg.text(0.15, y, label, fontsize=9, va="center")
        y -= 0.045


# --------------------------------------------------
# Plotting helpers
# --------------------------------------------------

def build_col_to_refpos(refpos_to_col: dict):
    return {col: pos for pos, col in refpos_to_col.items()}


def nearest_refpos_for_plotted_cols(plotted_fullcols: np.ndarray, col_to_refpos: dict):
    vals = [col_to_refpos.get(int(c), np.nan) for c in plotted_fullcols]
    vals = np.array(vals, dtype=float)

    if np.all(np.isnan(vals)):
        return vals

    valid_idx = np.where(~np.isnan(vals))[0]
    for i in range(len(vals)):
        if np.isnan(vals[i]):
            nearest = valid_idx[np.argmin(np.abs(valid_idx - i))]
            vals[i] = vals[nearest]

    return vals.astype(int)


def format_mb_ticks(vals):
    return [f"{v / 1_000_000:.3f}" for v in vals]


# --------------------------------------------------
# Plotting
# --------------------------------------------------

def plot_window(
    aln: Dict[str, str],
    ref_sid: str,
    sids_order: List[str],
    seqid: str,
    win_start: int,
    win_end: int,
    gene_feats: List[dict],
    repeat_feats: List[dict],
    refpos_to_col: dict,
    out_path: str,
    title: str,
    informative_only: bool = False,
    svg_out: Optional[str] = None,
    add_legend: bool = True,
    stats_prefix: Optional[str] = None
):
    ordered = [sid for sid in sids_order if sid in aln]
    if ref_sid not in ordered:
        raise RuntimeError(f"Reference sid '{ref_sid}' not present in alignment slice.")

    diffM_full, informative_mask = build_difference_matrix(aln, ref_sid, ordered)
    snp_cols = find_snp_columns(aln, ref_sid)

    if informative_only:
        plotted_fullcols = np.where(informative_mask)[0]
        if len(plotted_fullcols) == 0:
            raise RuntimeError("No informative columns found in requested window.")
        diffM = diffM_full[:, plotted_fullcols]
        snp_cols_plot = np.ones(len(plotted_fullcols), dtype=bool)
    else:
        plotted_fullcols = np.arange(diffM_full.shape[1])
        diffM = diffM_full
        snp_cols_plot = snp_cols

    nrows, ncols = diffM.shape
    cmap, norm = make_difference_colormap()

    if stats_prefix:
        build_sample_summary(
            aln=aln,
            ref_sid=ref_sid,
            ordered=ordered,
            diffM=diffM,
            informative_mask=informative_mask,
            snp_cols=snp_cols,
            plotted_fullcols=plotted_fullcols,
            informative_only=informative_only,
            win_start=win_start,
            win_end=win_end
        ).to_csv(f"{stats_prefix}.sample_summary.tsv", sep="\t", index=False)

        build_column_summary(
            aln=aln,
            ref_sid=ref_sid,
            ordered=ordered,
            plotted_fullcols=plotted_fullcols,
            informative_only=informative_only,
            refpos_to_col=refpos_to_col
        ).to_csv(f"{stats_prefix}.column_summary.tsv", sep="\t", index=False)

    if add_legend:
        fig = plt.figure(figsize=(19, 9))
        gs = fig.add_gridspec(
            nrows=4,
            ncols=2,
            width_ratios=[5.8, 1.9],
            height_ratios=[0.85, 0.65, 0.45, 5.20],
            wspace=0.06,
            hspace=0.08
        )
        ax_gene = fig.add_subplot(gs[0, 0])
        ax_rep = fig.add_subplot(gs[1, 0], sharex=ax_gene)
        ax_snp = fig.add_subplot(gs[2, 0], sharex=ax_gene)
        ax_aln = fig.add_subplot(gs[3, 0], sharex=ax_gene)
        ax_leg = fig.add_subplot(gs[:, 1])
    else:
        fig = plt.figure(figsize=(16, 9))
        gs = fig.add_gridspec(
            nrows=4,
            ncols=1,
            height_ratios=[0.85, 0.65, 0.45, 5.20],
            hspace=0.08
        )
        ax_gene = fig.add_subplot(gs[0, 0])
        ax_rep = fig.add_subplot(gs[1, 0], sharex=ax_gene)
        ax_snp = fig.add_subplot(gs[2, 0], sharex=ax_gene)
        ax_aln = fig.add_subplot(gs[3, 0], sharex=ax_gene)
        ax_leg = None

    # ----------------------------
    # Gene track
    # ----------------------------
    ax_gene.set_ylim(0, 1)
    ax_gene.set_yticks([])
    ax_gene.set_ylabel("Genes", fontsize=10)

    gene_list = [g for g in gene_feats if g["type"] == "gene"]
    gene_map = {g["attrs"].get("ID", g["attrs"].get("Name", "gene")): g for g in gene_list}
    gene_rows = assign_gene_rows(gene_feats, max_rows=4)
    row_y = [0.08, 0.30, 0.52, 0.74]

    for gene_id, gene_feat in sorted(gene_map.items(), key=lambda kv: (kv[1]["start"], kv[1]["end"])):
        coords = feature_to_cols(gene_feat["start"], gene_feat["end"], refpos_to_col)
        if coords is None:
            continue
        x0_full, x1_full = coords

        if informative_only:
            keep_set = set(plotted_fullcols)
            gene_cols = [c for c in range(x0_full, x1_full + 1) if c in keep_set]
            if not gene_cols:
                continue
            x0 = np.searchsorted(plotted_fullcols, min(gene_cols))
            x1 = np.searchsorted(plotted_fullcols, max(gene_cols))
        else:
            x0, x1 = x0_full, x1_full

        row = gene_rows.get(gene_id, 0)
        y = row_y[row]
        ax_gene.plot([x0, x1], [y + 0.06, y + 0.06], color="black", lw=1.0)

        label = gene_feat["attrs"].get("Name", gene_id)
        ax_gene.text((x0 + x1) / 2, y + 0.12, label, ha="center", va="bottom", fontsize=8)

        gene_width = max(1.0, x1 - x0)
        arrow_len = min(max(3.0, gene_width * 0.18), max(3.0, gene_width * 0.8))
        if gene_feat["strand"] == "+":
            arrow_start = max(x0, x1 - arrow_len)
            arrow_end = x1
        else:
            arrow_start = min(x1, x0 + arrow_len)
            arrow_end = x0

        ax_gene.annotate(
            "",
            xy=(arrow_end, y + 0.06),
            xytext=(arrow_start, y + 0.06),
            arrowprops=dict(
                arrowstyle="-|>",
                lw=1.0,
                color="black",
                shrinkA=0,
                shrinkB=0,
                mutation_scale=11
            ),
            annotation_clip=True
        )

    exons = [f for f in gene_feats if f["type"] == "exon"]
    for exon in exons:
        parent = pick_exon_parent(exon)
        coords = feature_to_cols(exon["start"], exon["end"], refpos_to_col)
        if coords is None:
            continue
        x0_full, x1_full = coords

        if informative_only:
            keep_set = set(plotted_fullcols)
            exon_cols = [c for c in range(x0_full, x1_full + 1) if c in keep_set]
            if not exon_cols:
                continue
            x0 = np.searchsorted(plotted_fullcols, min(exon_cols))
            x1 = np.searchsorted(plotted_fullcols, max(exon_cols))
        else:
            x0, x1 = x0_full, x1_full

        gene_row = 0
        for gid in gene_rows:
            if parent_matches_gene(parent, gid):
                gene_row = gene_rows[gid]
                break

        y = row_y[gene_row]
        ax_gene.add_patch(Rectangle((x0, y), max(1, x1 - x0 + 1), 0.12, color="black"))

    # ----------------------------
    # Repeat track
    # ----------------------------
    ax_rep.set_ylim(0, 1)
    ax_rep.set_yticks([])
    ax_rep.set_ylabel("Repeats", fontsize=10)

    for feat in repeat_feats:
        coords = feature_to_cols(feat["start"], feat["end"], refpos_to_col)
        if coords is None:
            continue
        x0_full, x1_full = coords

        if informative_only:
            keep_set = set(plotted_fullcols)
            rep_cols = [c for c in range(x0_full, x1_full + 1) if c in keep_set]
            if not rep_cols:
                continue
            x0 = np.searchsorted(plotted_fullcols, min(rep_cols))
            x1 = np.searchsorted(plotted_fullcols, max(rep_cols))
        else:
            x0, x1 = x0_full, x1_full

        ax_rep.add_patch(
            Rectangle(
                (x0, 0.2),
                max(1, x1 - x0 + 1),
                0.6,
                color=get_repeat_color(feat),
                alpha=0.9
            )
        )

    # ----------------------------
    # SNP track
    # ----------------------------
    ax_snp.set_ylim(0, 1)
    ax_snp.set_yticks([])
    ax_snp.set_ylabel("SNPs", fontsize=10)

    snp_x = np.where(snp_cols_plot)[0]
    ax_snp.vlines(snp_x + 0.5, 0.18, 0.82, color="black", linewidth=0.25, alpha=0.7)

    # ----------------------------
    # Alignment matrix
    # ----------------------------
    X, Y = np.meshgrid(np.arange(ncols + 1), np.arange(nrows + 1))
    ax_aln.pcolormesh(
        X,
        Y,
        diffM,
        cmap=cmap,
        norm=norm,
        shading="flat",
        linewidth=0.0,
        edgecolors="none",
        rasterized=False
    )

    ax_aln.set_ylim(nrows, 0)
    ax_aln.set_yticks(np.arange(nrows) + 0.5)
    ax_aln.set_yticklabels(ordered, fontsize=10)

    # bottom x-axis = Wm82a6 reference position (Mb)
    col_to_refpos = build_col_to_refpos(refpos_to_col)
    plotted_refpos = nearest_refpos_for_plotted_cols(plotted_fullcols, col_to_refpos)

    n_ticks = min(8, max(2, ncols))
    tick_positions = np.linspace(0, max(0, ncols - 1), n_ticks, dtype=int)

    if len(plotted_refpos) > 0:
        tick_refpos = [int(plotted_refpos[min(i, len(plotted_refpos) - 1)]) for i in tick_positions]
        tick_labels = format_mb_ticks(tick_refpos)
    else:
        tick_labels = [""] * len(tick_positions)

    ax_aln.set_xticks(tick_positions + 0.5)
    ax_aln.set_xticklabels(tick_labels, rotation=0, fontsize=9)
    ax_aln.tick_params(axis="x", length=4, width=0.8, pad=2)
    ax_aln.set_xlabel("Wm82a6 reference position (Mb)", fontsize=10, labelpad=6)

    for x in tick_positions + 0.5:
        ax_aln.axvline(x, color="black", lw=0.25, alpha=0.15, zorder=0)
        ax_gene.axvline(x, color="black", lw=0.25, alpha=0.10, zorder=0)
        ax_rep.axvline(x, color="black", lw=0.25, alpha=0.10, zorder=0)
        ax_snp.axvline(x, color="black", lw=0.25, alpha=0.10, zorder=0)

    for ax in [ax_gene, ax_rep, ax_snp]:
        ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
        for spine in ax.spines.values():
            spine.set_visible(False)

    ax_aln.spines["top"].set_visible(False)
    ax_aln.spines["right"].set_visible(False)

    if ax_leg is not None:
        draw_legend_panel(ax_leg)

    fig.suptitle(title, fontsize=16, y=0.955)

    if add_legend:
        fig.subplots_adjust(top=0.92, bottom=0.10, left=0.07, right=0.98, hspace=0.04, wspace=0.05)
    else:
        fig.subplots_adjust(top=0.90, bottom=0.10, left=0.10, right=0.985, hspace=0.08)

    plt.savefig(out_path, bbox_inches="tight")
    if svg_out:
        plt.savefig(svg_out, bbox_inches="tight")
    elif out_path.lower().endswith(".pdf"):
        plt.savefig(out_path[:-4] + ".svg", bbox_inches="tight")

    plt.close(fig)


# --------------------------------------------------
# MAIN
# --------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Plot an XMFA window with gene/repeat overlays, legend, and quantitative TSV outputs."
    )
    ap.add_argument("--xmfa", required=True)
    ap.add_argument("--ref-sid", required=True)
    ap.add_argument("--seqid", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--genes-gff", required=True)
    ap.add_argument("--repeats-gff", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--svg-out", default=None)
    ap.add_argument("--title", default="XMFA window")
    ap.add_argument("--order", nargs="+", required=True)
    ap.add_argument("--informative-only", action="store_true")
    ap.add_argument(
        "--repeat-types",
        nargs="+",
        default=["repeat_region", "transposable_element", "similarity", "match"]
    )
    ap.add_argument("--no-legend", action="store_true")
    ap.add_argument("--stats-prefix", default=None)

    args = ap.parse_args()

    if args.start > args.end:
        raise ValueError("--start must be <= --end")

    gene_feats = parse_gff3_features(
        args.genes_gff,
        args.seqid,
        args.start,
        args.end,
        feature_types=("gene", "exon")
    )

    repeat_feats = parse_gff3_features(
        args.repeats_gff,
        args.seqid,
        args.start,
        args.end,
        feature_types=tuple(args.repeat_types)
    )

    merged = OrderedDict((sid, "") for sid in args.order)
    refpos_to_col_merged = {}
    col_offset = 0
    any_hit = False

    for block in read_xmfa_blocks(args.xmfa):
        rs = ref_span_for_block(block, args.ref_sid)
        if rs is None:
            continue

        rlo, rhi = rs
        if rhi < args.start or rlo > args.end:
            continue

        sliced, refpos_to_col = slice_alignment_by_ref_window(
            block, args.ref_sid, args.start, args.end, args.order
        )
        if sliced is None:
            continue

        any_hit = True
        current_len = len(next(iter(sliced.values())))

        for sid in merged:
            if sid in sliced:
                merged[sid] += sliced[sid]
            else:
                merged[sid] += "-" * current_len

        for pos, col in refpos_to_col.items():
            refpos_to_col_merged[pos] = col + col_offset

        col_offset += current_len

    if not any_hit:
        raise RuntimeError(
            "No XMFA blocks overlapped the requested window. "
            "Check --ref-sid, coordinates, and XMFA sequence labels."
        )

    merged = OrderedDict((sid, seq) for sid, seq in merged.items() if seq)

    if args.ref_sid not in merged:
        raise RuntimeError(f"Reference sid '{args.ref_sid}' was not found after slicing.")

    lengths = {sid: len(seq) for sid, seq in merged.items()}
    if len(set(lengths.values())) != 1:
        raise RuntimeError(f"Sliced alignment sequences have unequal lengths: {lengths}")

    plot_title = args.title
    if args.informative_only and "informative" not in plot_title.lower():
        plot_title = f"{plot_title}; informative columns"

    plot_window(
        aln=merged,
        ref_sid=args.ref_sid,
        sids_order=list(merged.keys()),
        seqid=args.seqid,
        win_start=args.start,
        win_end=args.end,
        gene_feats=gene_feats,
        repeat_feats=repeat_feats,
        refpos_to_col=refpos_to_col_merged,
        out_path=args.out,
        title=plot_title,
        informative_only=args.informative_only,
        svg_out=args.svg_out,
        add_legend=(not args.no_legend),
        stats_prefix=args.stats_prefix
    )


if __name__ == "__main__":
    main()