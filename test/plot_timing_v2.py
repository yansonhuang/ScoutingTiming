#!/usr/bin/env python3
"""Plot 1D timing distributions and timing vs cell/jet properties.
 
Cell level (EB and EBCleaned separately):
  - Cell timing 1D distribution
  - Cell timing vs flags (integer value, log colour scale)
  - Cell timing vs cell energy (2D)
 
Jet level (EB timing, matched jets only):
  - Jet timing vs pT (2D)
  - Jet timing vs energy (2D)
  - Jet timing vs nCell (2D, log x)
"""
 
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import uproot
import awkward as ak
import mplhep as hep
 
plt.style.use(hep.style.CMS)
 
def save(fig, path):
    fig.savefig(path, dpi=180, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")
 
# ── Reference lines ───────────────────────────────────────────────────────────
# Each entry: (t_value, color, linestyle, linewidth, label, key)
# key is the string used to select which lines to draw — must be unique.
REF_LINES = [
    ( 0.0,  "gray",   ":",  1.0, "t = 0 ns",    "0"),
    (-4.5,  "orange", "--", 1.5, "t = −4.5 ns", "-4.5"),
    (-5.5,  "red",    "--", 1.5, "t = −5.5 ns", "-5.5"),
    ( 2.0,  "green",  "--", 1.5, "t = +2 ns",   "+2"),
]
 
def add_hlines(ax, which=("0", "-4.5", "-5.5", "+2")):
    """Draw horizontal reference lines on ax. which = tuple of keys."""
    drawn = []
    for t, color, ls, lw, label, key in REF_LINES:
        if key in which:
            ax.axhline(t, color=color, linestyle=ls, linewidth=lw, label=label)
            drawn.append(label)
    return drawn
 
def add_vlines(ax, which=("0", "-4.5", "-5.5", "+2")):
    """Draw vertical reference lines on ax. which = tuple of keys."""
    drawn = []
    for t, color, ls, lw, label, key in REF_LINES:
        if key in which:
            ax.axvline(t, color=color, linestyle=ls, linewidth=lw, label=label)
            drawn.append(label)
    return drawn
 
# ── Plot helpers ──────────────────────────────────────────────────────────────
 
def plot_1d(h, bins, xlabel, label, color, path, ref_lines=("0", "-4.5", "-5.5", "+2")):
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.stairs(h, bins, linewidth=2, color=color, label=label)
    add_vlines(ax, ref_lines)
    ax.set_yscale("log")
    ax.set(xlabel=xlabel, ylabel="Hits / bin")
    ax.legend(frameon=False, fontsize=14)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, path)
 
def plot_2d(h2d, binsx, binsy, xlabel, ylabel, title, path,
            cmap="Blues", hlines=("0", "-4.5", "-5.5"), vlines=None, lognorm=True):
    fig, ax = plt.subplots(figsize=(12, 8))
    norm = mcolors.LogNorm(vmin=1) if lognorm else None
    im = ax.pcolormesh(binsx, binsy, h2d.T, cmap=cmap, norm=norm)
    plt.colorbar(im, ax=ax, label="Entries")
    drawn = add_hlines(ax, hlines)
    if vlines:
        for xv, c, ls in vlines:
            ax.axvline(xv, color=c, linestyle=ls, linewidth=1.5)
    if drawn:
        ax.legend(frameon=False, fontsize=13)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, path)
 
# ── Main ─────────────────────────────────────────────────────────────────────
 
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path,
                        default=Path("test_6files_v2.root"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs_timing"))
    parser.add_argument("--tree",       default="testScoutingTiming/timeTree")
    parser.add_argument("--chunk-size", type=int, default=100_000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
 
    # ── Bins ──────────────────────────────────────────────────────────────────
    bins_t_cell  = np.linspace(-25,  25, 101)   # cell timing 1D
    bins_t_jet   = np.linspace(-25,  25, 101)   # jet timing 1D / 2D y-axis
    bins_t_2d    = np.linspace(-25,  25, 101)   # 2D y-axis for cell plots
    bins_flags   = np.arange(0, 300_001, 2000)  # flag integer (0..300000 in steps)
    bins_e_cell  = np.linspace(0.5,  50,  100)  # cell energy
    bins_pt_jet  = np.linspace(0,   500,  51)   # jet pT
    bins_e_jet   = np.linspace(0,   600,  61)   # jet energy
    bins_ncell   = np.concatenate([              # log-spaced nCell bins
        [0], np.logspace(0, 2, 51)
    ])
 
    # ── EB accumulators ───────────────────────────────────────────────────────
    h_eb_t1d        = np.zeros(len(bins_t_cell) - 1)
    h_eb_t_flags    = np.zeros((len(bins_flags) - 1, len(bins_t_2d) - 1))
    h_eb_t_e        = np.zeros((len(bins_e_cell) - 1, len(bins_t_2d) - 1))
    # Hits with any flag beyond kGood-only (flags != 1)
    h_eb_t1d_nongood = np.zeros(len(bins_t_cell) - 1)
    # kGood-only hits (flags == 1)
    h_eb_t1d_kgood   = np.zeros(len(bins_t_cell) - 1)
 
    # ── EBCleaned accumulators ────────────────────────────────────────────────
    h_cl_t1d     = np.zeros(len(bins_t_cell) - 1)
    h_cl_t_flags = np.zeros((len(bins_flags) - 1, len(bins_t_2d) - 1))
    h_cl_t_e     = np.zeros((len(bins_e_cell) - 1, len(bins_t_2d) - 1))
 
    # ── Jet accumulators ──────────────────────────────────────────────────────
    h_jet_t1d    = np.zeros(len(bins_t_jet) - 1)
    h_jet_t_pt   = np.zeros((len(bins_pt_jet) - 1, len(bins_t_jet) - 1))
    h_jet_t_e    = np.zeros((len(bins_e_jet) - 1,  len(bins_t_jet) - 1))
    h_jet_t_nc   = np.zeros((len(bins_ncell) - 1,  len(bins_t_jet) - 1))
 
    branches = [
        "caloCell_ecalTime",        "caloCell_flags",        "caloCell_e",
        "caloCellCleaned_ecalTime", "caloCellCleaned_flags", "caloCellCleaned_e",
        "pfJet_weightedTime",       "pfJet_pt",              "pfJet_e",
        "pfJet_nCell",
    ]
 
    tree = uproot.open(f"{args.input}:{args.tree}")
    n_entries = tree.num_entries
    print(f"Total entries: {n_entries:,}")
 
    n_jet_matched = 0
    n_jet_delayed = 0
 
    for start in range(0, n_entries, args.chunk_size):
        stop = min(start + args.chunk_size, n_entries)
        print(f"  {stop:,} / {n_entries:,} ({100*stop/n_entries:.1f}%)", flush=True)
        chunk = tree.arrays(branches, entry_start=start, entry_stop=stop, library="ak")
 
        # ── EB cells ──
        c_t = ak.to_numpy(ak.flatten(chunk["caloCell_ecalTime"]))
        c_f = ak.to_numpy(ak.flatten(chunk["caloCell_flags"])).astype(np.uint32)
        c_e = ak.to_numpy(ak.flatten(chunk["caloCell_e"]))
 
        h_eb_t1d        += np.histogram(c_t,              bins=bins_t_cell)[0]
        h_eb_t_flags    += np.histogram2d(c_f.astype(float), c_t,
                                           bins=[bins_flags, bins_t_2d])[0]
        h_eb_t_e        += np.histogram2d(c_e, c_t,       bins=[bins_e_cell, bins_t_2d])[0]
        # Non-kGood-only: any hit where flags != 0x1 (kGood alone)
        nongood = c_f != 1
        kgood   = c_f == 1
        h_eb_t1d_nongood += np.histogram(c_t[nongood], bins=bins_t_cell)[0]
        h_eb_t1d_kgood   += np.histogram(c_t[kgood],   bins=bins_t_cell)[0]
 
        # ── EBCleaned cells ──
        cl_t = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_ecalTime"]))
        cl_f = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_flags"])).astype(float)
        cl_e = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_e"]))
 
        h_cl_t1d     += np.histogram(cl_t,          bins=bins_t_cell)[0]
        h_cl_t_flags += np.histogram2d(cl_f, cl_t,  bins=[bins_flags,  bins_t_2d])[0]
        h_cl_t_e     += np.histogram2d(cl_e, cl_t,  bins=[bins_e_cell, bins_t_2d])[0]
 
        # ── Jets (matched, sentinel filtered) ──
        j_t  = ak.to_numpy(ak.flatten(chunk["pfJet_weightedTime"]))
        j_pt = ak.to_numpy(ak.flatten(chunk["pfJet_pt"]))
        j_e  = ak.to_numpy(ak.flatten(chunk["pfJet_e"]))
        j_nc = ak.to_numpy(ak.flatten(chunk["pfJet_nCell"])).astype(int)
 
        matched = (j_nc > 0) & (j_t > -100)
        n_jet_matched += matched.sum()
        n_jet_delayed += (matched & (j_t > 2.0)).sum()
 
        h_jet_t1d  += np.histogram(j_t[matched], bins=bins_t_jet)[0]
        h_jet_t_pt += np.histogram2d(j_pt[matched], j_t[matched],
                                      bins=[bins_pt_jet, bins_t_jet])[0]
        h_jet_t_e  += np.histogram2d(j_e[matched],  j_t[matched],
                                      bins=[bins_e_jet,  bins_t_jet])[0]
        h_jet_t_nc += np.histogram2d(j_nc[matched].astype(float), j_t[matched],
                                      bins=[bins_ncell,   bins_t_jet])[0]
 
        del chunk, c_t, c_f, c_e, cl_t, cl_f, cl_e
        del j_t, j_pt, j_e, j_nc, matched
 
    print(f"\nMatched jets (nCell>0): {n_jet_matched:,}")
    print(f"Delayed jets (t>2 ns):  {n_jet_delayed:,}  "
          f"({100*n_jet_delayed/max(n_jet_matched,1):.3f}%)")
 
    # ── EB cell plots ─────────────────────────────────────────────────────────
    plot_1d(h_eb_t1d, bins_t_cell,
            xlabel="Crystal ECAL time [ns]",
            label="EB rec hits (E > 0.5 GeV)",
            color="steelblue",
            ref_lines=("0", "-4.5", "-5.5", "+2"),
            path=args.output_dir / "eb_cell_time_1d.png")
 
    # Overlay: all EB hits vs kGood-only vs non-kGood-only
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.stairs(h_eb_t1d,         bins_t_cell, linewidth=2, color="steelblue",
              label="All EB hits (E > 0.5 GeV)")
    ax.stairs(h_eb_t1d_kgood,   bins_t_cell, linewidth=2, color="mediumseagreen",
              label="kGood-only (flags = 0x1)")
    ax.stairs(h_eb_t1d_nongood, bins_t_cell, linewidth=2, color="crimson",
              label="Non-kGood-only (flags ≠ 0x1)")
    add_vlines(ax, ("0", "-4.5", "-5.5", "+2"))
    ax.set_yscale("log")
    ax.set(xlabel="Crystal ECAL time [ns]", ylabel="Hits / bin")
    ax.legend(frameon=False, fontsize=14)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "eb_cell_time_1d_flagcheck.png")
 
    plot_2d(h_eb_t_flags, bins_flags, bins_t_2d,
            xlabel="caloCell_flags (integer)",
            ylabel="Crystal ECAL time [ns]",
            title="EB cell time vs flags (integer)",
            cmap="Blues", hlines=("0", "-4.5", "-5.5"),
            path=args.output_dir / "eb_cell_time_vs_flags.png")
 
    plot_2d(h_eb_t_e, bins_e_cell, bins_t_2d,
            xlabel="Crystal energy [GeV]",
            ylabel="Crystal ECAL time [ns]",
            title="EB cell time vs energy",
            cmap="Blues", hlines=("0", "-4.5", "-5.5"),
            path=args.output_dir / "eb_cell_time_vs_energy.png")
 
    # ── EBCleaned cell plots ───────────────────────────────────────────────────
    plot_1d(h_cl_t1d, bins_t_cell,
            xlabel="Crystal ECAL time [ns]",
            label="EBCleaned rec hits (E > 0.5 GeV)",
            color="darkorange",
            ref_lines=("0", "-4.5", "-5.5"),
            path=args.output_dir / "ebcleaned_cell_time_1d.png")
 
    plot_2d(h_cl_t_flags, bins_flags, bins_t_2d,
            xlabel="caloCellCleaned_flags (integer)",
            ylabel="Crystal ECAL time [ns]",
            title="EBCleaned cell time vs flags (integer)",
            cmap="Oranges", hlines=("0", "-4.5", "-5.5"),
            path=args.output_dir / "ebcleaned_cell_time_vs_flags.png")
 
    plot_2d(h_cl_t_e, bins_e_cell, bins_t_2d,
            xlabel="Crystal energy [GeV]",
            ylabel="Crystal ECAL time [ns]",
            title="EBCleaned cell time vs energy",
            cmap="Oranges", hlines=("0", "-4.5", "-5.5"),
            path=args.output_dir / "ebcleaned_cell_time_vs_energy.png")
 
    # ── Jet plots ─────────────────────────────────────────────────────────────
    plot_1d(h_jet_t1d, bins_t_jet,
            xlabel="ET-weighted jet ECAL time [ns]",
            label=f"Matched jets (nCell > 0): {n_jet_matched:,}",
            color="mediumseagreen",
            ref_lines=("0", "-4.5", "-5.5", "+2"),
            path=args.output_dir / "jet_time_1d.png")
 
    plot_2d(h_jet_t_pt, bins_pt_jet, bins_t_jet,
            xlabel="Jet $p_T$ [GeV]",
            ylabel="ET-weighted jet ECAL time [ns]",
            title="Jet time vs $p_T$ (matched jets)",
            cmap="Blues", hlines=("0", "-4.5", "-5.5", "+2"),
            path=args.output_dir / "jet_time_vs_pt.png")
 
    plot_2d(h_jet_t_e, bins_e_jet, bins_t_jet,
            xlabel="Jet energy [GeV]",
            ylabel="ET-weighted jet ECAL time [ns]",
            title="Jet time vs energy (matched jets)",
            cmap="Blues", hlines=("0", "-4.5", "-5.5", "+2"),
            path=args.output_dir / "jet_time_vs_energy.png")
 
    # nCell: log x axis — built manually to set log scale
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.pcolormesh(bins_ncell, bins_t_jet, h_jet_t_nc.T,
                       cmap="Blues", norm=mcolors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label="Jets")
    add_hlines(ax, ("0", "-4.5", "-5.5", "+2"))
    ax.set_xscale("log")
    ax.set(xlabel="Matched EB crystals per jet (nCell)",
           ylabel="ET-weighted jet ECAL time [ns]",
           title="Jet time vs nCell (matched jets)")
    ax.legend(frameon=False, fontsize=13)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "jet_time_vs_ncell.png")
 
    print(f"\nAll plots written to {args.output_dir}/")
 
if __name__ == "__main__":
    main()