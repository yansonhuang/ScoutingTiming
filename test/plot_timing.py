#!/usr/bin/env python3
"""Plot ECAL timing distributions - memory efficient chunked version."""

import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
import awkward as ak
import mplhep as hep

plt.style.use(hep.style.CMS)

def save(fig, path):
    fig.savefig(path, dpi=180)
    fig.savefig(path.with_suffix(".pdf"), dpi=180)
    plt.close(fig)
    print(f"Saved {path}")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path,
                        default=Path("Validation_60files.root"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs_validation"))
    parser.add_argument("--tree", default="testScoutingTiming/timeTree")
    parser.add_argument("--chunk-size", type=int, default=100_000)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Opening file...")
    tree = uproot.open(f"{args.input}:{args.tree}")
    n_entries = tree.num_entries
    print(f"Total entries: {n_entries:,}")

    # ── Define histogram bins once ─────────────────────────────────────────
    bins_time      = np.linspace(-25, 25, 101)
    bins_time_zoom = np.linspace(-10, 10,  81)
    bins_cell_time = np.linspace(-25, 25, 101)
    bins_nc        = np.arange(0, 102, 1)
    bins_pt        = np.linspace(0, 500,  51)
    bins_t2d       = np.linspace(-10, 10,  81)
    bins_pt_spec   = np.linspace(0, 600, 121)


    branches = ["pfJet_weightedTime", "pfJet_pt", "pfJet_nCell",
                "caloCell_ecalTime"]

    # ── accumulators ──
    h_jet_time_all      = np.zeros(len(bins_time) - 1)
    h_jet_time_zoom_all = np.zeros(len(bins_time_zoom) - 1)
    h_cell_time         = np.zeros(len(bins_cell_time) - 1)
    h_ncell             = np.zeros(len(bins_nc) - 1)
    h_pt_spec           = np.zeros(len(bins_pt_spec) - 1)
    h_2d                = np.zeros((len(bins_pt) - 1, len(bins_t2d) - 1))

    n_jet_all     = 0
    n_jet_matched = 0
    n_jet_delayed = 0

    for start in range(0, n_entries, args.chunk_size):
        stop = min(start + args.chunk_size, n_entries)
        pct  = 100 * stop / n_entries
        print(f"  Processing {stop:,} / {n_entries:,} ({pct:.1f}%)", flush=True)

        chunk = tree.arrays(branches, entry_start=start, entry_stop=stop, library="ak")

        jet_time  = ak.to_numpy(ak.flatten(chunk["pfJet_weightedTime"]))
        jet_pt    = ak.to_numpy(ak.flatten(chunk["pfJet_pt"]))
        jet_ncell = ak.to_numpy(ak.flatten(chunk["pfJet_nCell"])).astype(int)

        matched = jet_ncell > 0
        prompt  = matched & (jet_time > -100)

        n_jet_all     += len(jet_time)
        n_jet_matched += matched.sum()
        n_jet_delayed += (matched & (jet_time > 2.0)).sum()

        h_jet_time_all      += np.histogram(jet_time[prompt], bins=bins_time)[0]
        h_jet_time_zoom_all += np.histogram(jet_time[prompt], bins=bins_time_zoom)[0]
        h_ncell             += np.histogram(jet_ncell,        bins=bins_nc)[0]
        h_pt_spec           += np.histogram(jet_pt,           bins=bins_pt_spec)[0]
        h_2d                += np.histogram2d(jet_pt[prompt], jet_time[prompt],
                                              bins=[bins_pt, bins_t2d])[0]

        cell_time = ak.to_numpy(ak.flatten(chunk["caloCell_ecalTime"]))
        h_cell_time += np.histogram(cell_time, bins=bins_cell_time)[0]

        del chunk, jet_time, jet_pt, jet_ncell, cell_time, matched, prompt

    print(f"Total jets:            {n_jet_all:,}")
    print(f"Jets with nCell > 0:   {n_jet_matched:,}  ({100*n_jet_matched/n_jet_all:.1f}%)")
    print(f"Jets with time > 2 ns: {n_jet_delayed:,}  ({100*n_jet_delayed/n_jet_all:.3f}%)")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 1: full range
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.stairs(h_jet_time_all, bins_time, linewidth=2, color="steelblue",
              label=f"Jets with nCell > 0: {n_jet_matched:,}")
    ax.axvline(2.0, color="gray", linestyle="--", linewidth=1.5, label="t = 2 ns")
    ax.set(xlabel="ET-weighted jet ECAL time [ns]", ylabel="Jets / 0.5 ns", xlim=(-25, 25))
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=16)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "jet_weighted_time.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 2: zoomed
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.stairs(h_jet_time_zoom_all, bins_time_zoom, linewidth=2, color="steelblue",
              label=f"Jets with nCell > 0: {n_jet_matched:,}")
    ax.axvline(2.0, color="gray", linestyle="--", linewidth=1.5, label="t = 2 ns")
    ax.set(xlabel="ET-weighted jet ECAL time [ns]", ylabel="Jets / 0.25 ns", xlim=(-10, 10))
    ax.legend(frameon=False, fontsize=16)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "jet_weighted_time_zoom.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 3: Crystal ECAL time
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.stairs(h_cell_time, bins_cell_time, linewidth=2, color="darkorange",
              label=f"EB rec hits (E > 0.5 GeV)")
    ax.set(xlabel="Crystal ECAL time [ns]", ylabel="Crystals / 0.5 ns", xlim=(-25, 25))
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=16)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "crystal_ecal_time.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 4: nCell per jet
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    frac_zero = h_ncell[0] / h_ncell.sum()
    ax.stairs(h_ncell, bins_nc, linewidth=2, color="mediumseagreen",
              label=f"All jets: {n_jet_all:,}\nFraction nCell=0: {frac_zero:.2f}")
    ax.set(xlabel="Matched EB crystals per jet (dR<0.4, E>0.5 GeV)",
           ylabel="Jets", xlim=(0, 100))
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=16)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "jet_ncell.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 5: Jet time vs pT 2D
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.pcolormesh(bins_pt, bins_t2d, h_2d.T, cmap="Blues",
                       norm=matplotlib.colors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label="Jets")
    ax.axhline(2.0, color="crimson", linestyle="--", linewidth=1.5, label="t = 2 ns")
    ax.axhline(0.0, color="gray",    linestyle=":",  linewidth=1.0)
    ax.set(xlabel="Jet $p_T$ [GeV]", ylabel="ET-weighted jet ECAL time [ns]",
           xlim=(0, 500), ylim=(-10, 10))
    ax.legend(frameon=False, fontsize=16)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "jet_time_vs_pt_2d.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 6: Jet pT spectrum
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.stairs(h_pt_spec, bins_pt_spec, linewidth=2, color="steelblue",
              label=f"All jets: {n_jet_all:,}")
    ax.set(xlabel="Jet $p_T$ [GeV]", ylabel="Jets / 5 GeV", xlim=(0, 600))
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=16)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "jet_pt_spectrum.png")

    print(f"\nAll plots written to {args.output_dir}/")

if __name__ == "__main__":
    main()