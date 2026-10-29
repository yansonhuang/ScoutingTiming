#!/usr/bin/env python3
"""Plot timing vs angular position (no cuts) for EB and EBCleaned collections."""
 
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
 
def plot_time_vs_eta(h2d, bins_eta, bins_time, title, path):
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.pcolormesh(bins_eta, bins_time, h2d.T,
                       cmap="Blues", norm=mcolors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label="Hits")
    ax.axhline(0,    color="gray",   linestyle=":",  linewidth=1)
    ax.axhline(-3.0, color="orange", linestyle="--", linewidth=1.5, label="t = −3 ns")
    ax.axhline(-8.0, color="red",    linestyle="--", linewidth=1.5, label="t = −8 ns")
    ax.set(xlabel=r"$\eta$", ylabel="ECAL time [ns]", title=title)
    ax.legend(frameon=False, fontsize=14)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, path)
 
def plot_time_vs_phi(h2d, bins_phi, bins_time, title, path):
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.pcolormesh(bins_phi, bins_time, h2d.T,
                       cmap="Blues", norm=mcolors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label="Hits")
    ax.axhline(0,    color="gray",   linestyle=":",  linewidth=1)
    ax.axhline(-3.0, color="orange", linestyle="--", linewidth=1.5, label="t = −3 ns")
    ax.axhline(-8.0, color="red",    linestyle="--", linewidth=1.5, label="t = −8 ns")
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
                  [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.set(xlabel=r"$\phi$ [rad]", ylabel="ECAL time [ns]", title=title)
    ax.legend(frameon=False, fontsize=14)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, path)
 
def plot_occupancy(h2d, bins_eta, bins_phi, title, path, cmap="Blues"):
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.pcolormesh(bins_eta, bins_phi, h2d.T,
                       cmap=cmap, norm=mcolors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label="Entries")
    ax.set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
                  [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.set(xlabel=r"$\eta$", ylabel=r"$\phi$ [rad]", title=title)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, path)
 
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path,
                        default=Path("test_6files_v2.root"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs_timing_angular"))
    parser.add_argument("--tree", default="testScoutingTiming/timeTree")
    parser.add_argument("--chunk-size", type=int, default=100_000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
 
    # ── Bins ──────────────────────────────────────────────────────────────
    bins_cell_eta  = np.linspace(-1.5,  1.5,  91)
    bins_jet_eta   = np.linspace(-5.0,  5.0, 101)
    bins_phi       = np.linspace(-np.pi, np.pi, 73)
    bins_time_cell = np.linspace(-25,   25,  101)
    bins_time_jet  = np.linspace(-25,   25,  101)
 
    # ── Accumulators ── EB ────────────────────────────────────────────────
    h_cell_t_eta   = np.zeros((len(bins_cell_eta)-1, len(bins_time_cell)-1))
    h_cell_t_phi   = np.zeros((len(bins_phi)-1,      len(bins_time_cell)-1))
    h_cell_occ     = np.zeros((len(bins_cell_eta)-1, len(bins_phi)-1))
    h_jet_t_eta    = np.zeros((len(bins_jet_eta)-1,  len(bins_time_jet)-1))
    h_jet_t_phi    = np.zeros((len(bins_phi)-1,      len(bins_time_jet)-1))
    h_jet_occ      = np.zeros((len(bins_jet_eta)-1,  len(bins_phi)-1))
 
    # ── Accumulators ── EBCleaned ─────────────────────────────────────────
    h_cl_t_eta     = np.zeros((len(bins_cell_eta)-1, len(bins_time_cell)-1))
    h_cl_t_phi     = np.zeros((len(bins_phi)-1,      len(bins_time_cell)-1))
    h_cl_occ       = np.zeros((len(bins_cell_eta)-1, len(bins_phi)-1))
 
    branches = [
        "caloCell_eta",            "caloCell_phi",            "caloCell_ecalTime",
        "caloCellCleaned_eta",     "caloCellCleaned_phi",     "caloCellCleaned_ecalTime",
        "pfJet_eta",               "pfJet_phi",               "pfJet_weightedTime",
        "pfJet_nCell",
    ]
 
    tree = uproot.open(f"{args.input}:{args.tree}")
    n_entries = tree.num_entries
    print(f"Total entries: {n_entries:,}")
 
    for start in range(0, n_entries, args.chunk_size):
        stop = min(start + args.chunk_size, n_entries)
        print(f"  {stop:,} / {n_entries:,} ({100*stop/n_entries:.1f}%)", flush=True)
        chunk = tree.arrays(branches, entry_start=start, entry_stop=stop, library="ak")
 
        # ── EB cells ──
        c_eta  = ak.to_numpy(ak.flatten(chunk["caloCell_eta"]))
        c_phi  = ak.to_numpy(ak.flatten(chunk["caloCell_phi"]))
        c_time = ak.to_numpy(ak.flatten(chunk["caloCell_ecalTime"]))
 
        h_cell_t_eta += np.histogram2d(c_eta,  c_time, bins=[bins_cell_eta, bins_time_cell])[0]
        h_cell_t_phi += np.histogram2d(c_phi,  c_time, bins=[bins_phi,      bins_time_cell])[0]
        h_cell_occ   += np.histogram2d(c_eta,  c_phi,  bins=[bins_cell_eta, bins_phi])[0]
 
        # ── EBCleaned cells ──
        cl_eta  = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_eta"]))
        cl_phi  = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_phi"]))
        cl_time = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_ecalTime"]))
 
        h_cl_t_eta += np.histogram2d(cl_eta, cl_time, bins=[bins_cell_eta, bins_time_cell])[0]
        h_cl_t_phi += np.histogram2d(cl_phi, cl_time, bins=[bins_phi,      bins_time_cell])[0]
        h_cl_occ   += np.histogram2d(cl_eta, cl_phi,  bins=[bins_cell_eta, bins_phi])[0]
 
        # ── Jets (matched only: nCell > 0, sentinel filtered) ──
        j_eta   = ak.to_numpy(ak.flatten(chunk["pfJet_eta"]))
        j_phi   = ak.to_numpy(ak.flatten(chunk["pfJet_phi"]))
        j_time  = ak.to_numpy(ak.flatten(chunk["pfJet_weightedTime"]))
        j_ncell = ak.to_numpy(ak.flatten(chunk["pfJet_nCell"])).astype(int)
        matched = (j_ncell > 0) & (j_time > -100)
 
        h_jet_t_eta += np.histogram2d(j_eta[matched],  j_time[matched],
                                       bins=[bins_jet_eta, bins_time_jet])[0]
        h_jet_t_phi += np.histogram2d(j_phi[matched],  j_time[matched],
                                       bins=[bins_phi, bins_time_jet])[0]
        h_jet_occ   += np.histogram2d(j_eta[matched],  j_phi[matched],
                                       bins=[bins_jet_eta, bins_phi])[0]
 
        del chunk, c_eta, c_phi, c_time, cl_eta, cl_phi, cl_time
        del j_eta, j_phi, j_time, j_ncell, matched
 
    # ── EB cell plots ─────────────────────────────────────────────────────
    plot_time_vs_eta(h_cell_t_eta, bins_cell_eta, bins_time_cell,
                     "EB cell ECAL time vs η (no cuts)",
                     args.output_dir / "eb_cell_time_vs_eta.png")
 
    plot_time_vs_phi(h_cell_t_phi, bins_phi, bins_time_cell,
                     "EB cell ECAL time vs φ (no cuts)",
                     args.output_dir / "eb_cell_time_vs_phi.png")
 
    plot_occupancy(h_cell_occ, bins_cell_eta, bins_phi,
                   "EB cell η-φ occupancy (no cuts)",
                   args.output_dir / "eb_cell_occupancy.png", cmap="Blues")
 
    # ── EBCleaned cell plots ───────────────────────────────────────────────
    plot_time_vs_eta(h_cl_t_eta, bins_cell_eta, bins_time_cell,
                     "EBCleaned cell ECAL time vs η (no cuts)",
                     args.output_dir / "ebcleaned_cell_time_vs_eta.png")
 
    plot_time_vs_phi(h_cl_t_phi, bins_phi, bins_time_cell,
                     "EBCleaned cell ECAL time vs φ (no cuts)",
                     args.output_dir / "ebcleaned_cell_time_vs_phi.png")
 
    plot_occupancy(h_cl_occ, bins_cell_eta, bins_phi,
                   "EBCleaned cell η-φ occupancy (no cuts)",
                   args.output_dir / "ebcleaned_cell_occupancy.png", cmap="Oranges")
 
    # ── Jet plots ─────────────────────────────────────────────────────────
    plot_time_vs_eta(h_jet_t_eta, bins_jet_eta, bins_time_jet,
                     "Jet ET-weighted time vs η (matched jets, no cuts)",
                     args.output_dir / "jet_time_vs_eta.png")
 
    plot_time_vs_phi(h_jet_t_phi, bins_phi, bins_time_jet,
                     "Jet ET-weighted time vs φ (matched jets, no cuts)",
                     args.output_dir / "jet_time_vs_phi.png")
 
    plot_occupancy(h_jet_occ, bins_jet_eta, bins_phi,
                   "Jet η-φ occupancy (matched jets, no cuts)",
                   args.output_dir / "jet_occupancy.png", cmap="Greens")
 
    print(f"\nAll plots written to {args.output_dir}/")
 
if __name__ == "__main__":
    main()