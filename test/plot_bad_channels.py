#!/usr/bin/env python3
"""Plot occupancy and timing maps for anomalously early hits/jets: -8 < t < -3 ns.
 
Produces separate plots for EB and EBCleaned.
Uses caloCell_ieta/iphi (EBDetId crystal indices) for cell-level maps,
and caloCell_eta/phi for continuous angular maps.
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
 
# Bad channel time window
T_LO = -8.0   # ns
T_HI = -3.0   # ns
 
def save(fig, path):
    fig.savefig(path, dpi=180, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")
 
def plot_ieta_iphi(h2d, bins_ieta, bins_iphi, title, path):
    """Crystal-index occupancy map in (ieta, iphi) space."""
    fig, ax = plt.subplots(figsize=(14, 6))
    # Use linear norm — bad channel map, counts may be small
    vmax = max(h2d.max(), 1)
    im = ax.pcolormesh(bins_ieta, bins_iphi, h2d.T,
                       cmap="hot_r",
                       norm=mcolors.LogNorm(vmin=1, vmax=vmax))
    plt.colorbar(im, ax=ax, label="Hit count")
    ax.set(xlabel=r"$i\eta$ (EBDetId)", ylabel=r"$i\phi$ (EBDetId)", title=title)
    ax.set_xlim(-85, 85)
    ax.set_ylim(1, 360)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, path)
 
def plot_eta_phi(h2d, bins_eta, bins_phi, title, path):
    """Continuous η-φ occupancy map."""
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.pcolormesh(bins_eta, bins_phi, h2d.T,
                       cmap="hot_r",
                       norm=mcolors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label="Entries")
    ax.set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
                  [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.set(xlabel=r"$\eta$", ylabel=r"$\phi$ [rad]", title=title)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, path)
 
def plot_time_vs_eta(h2d, bins_eta, bins_time, title, path):
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.pcolormesh(bins_eta, bins_time, h2d.T,
                       cmap="hot_r", norm=mcolors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label="Entries")
    ax.axhline(T_LO, color="gray", linestyle="--", linewidth=1.5,
               label=f"t = {T_LO} ns")
    ax.axhline(T_HI, color="gray", linestyle="--", linewidth=1.5,
               label=f"t = {T_HI} ns")
    ax.set(xlabel=r"$\eta$", ylabel="ECAL time [ns]", title=title)
    ax.legend(frameon=False, fontsize=14)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, path)
 
def plot_time_vs_phi(h2d, bins_phi, bins_time, title, path):
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.pcolormesh(bins_phi, bins_time, h2d.T,
                       cmap="hot_r", norm=mcolors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label="Entries")
    ax.axhline(T_LO, color="gray", linestyle="--", linewidth=1.5,
               label=f"t = {T_LO} ns")
    ax.axhline(T_HI, color="gray", linestyle="--", linewidth=1.5,
               label=f"t = {T_HI} ns")
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
                  [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.set(xlabel=r"$\phi$ [rad]", ylabel="ECAL time [ns]", title=title)
    ax.legend(frameon=False, fontsize=14)
    hep.cms.label("Preliminary", data=True, com=13.6, ax=ax)
    fig.tight_layout()
    save(fig, path)
 
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path,
                        default=Path("test_6files_v2.root"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs_bad_channels"))
    parser.add_argument("--tree", default="testScoutingTiming/timeTree")
    parser.add_argument("--chunk-size", type=int, default=100_000)
    parser.add_argument("--t-lo", type=float, default=T_LO,
                        help="Lower time cut [ns] (default: -8)")
    parser.add_argument("--t-hi", type=float, default=T_HI,
                        help="Upper time cut [ns] (default: -3)")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
 
    t_lo, t_hi = args.t_lo, args.t_hi
    print(f"Bad channel window: {t_lo} < t < {t_hi} ns")
 
    # ── Bins ──────────────────────────────────────────────────────────────
    bins_ieta      = np.arange(-85, 87, 1)          # ieta: -85..+85 (skip 0 internally)
    bins_iphi      = np.arange(1, 362, 1)            # iphi: 1..360
    bins_cell_eta  = np.linspace(-1.5,  1.5, 91)
    bins_jet_eta   = np.linspace(-5.0,  5.0, 101)
    bins_phi       = np.linspace(-np.pi, np.pi, 73)
    bins_time_cell = np.linspace(-10,   -1,  73)     # zoomed to bad window + margin
    bins_time_jet  = np.linspace(-10,   -1,  73)
 
    # ── EB accumulators ───────────────────────────────────────────────────
    h_eb_ieta_iphi = np.zeros((len(bins_ieta)-1, len(bins_iphi)-1))
    h_eb_eta_phi   = np.zeros((len(bins_cell_eta)-1, len(bins_phi)-1))
    h_eb_t_eta     = np.zeros((len(bins_cell_eta)-1, len(bins_time_cell)-1))
    h_eb_t_phi     = np.zeros((len(bins_phi)-1,      len(bins_time_cell)-1))
 
    # ── EBCleaned accumulators ────────────────────────────────────────────
    h_cl_ieta_iphi = np.zeros((len(bins_ieta)-1, len(bins_iphi)-1))
    h_cl_eta_phi   = np.zeros((len(bins_cell_eta)-1, len(bins_phi)-1))
    h_cl_t_eta     = np.zeros((len(bins_cell_eta)-1, len(bins_time_cell)-1))
    h_cl_t_phi     = np.zeros((len(bins_phi)-1,      len(bins_time_cell)-1))
 
    # ── Jet accumulators (EB timing only) ─────────────────────────────────
    h_jet_eta_phi  = np.zeros((len(bins_jet_eta)-1, len(bins_phi)-1))
    h_jet_t_eta    = np.zeros((len(bins_jet_eta)-1, len(bins_time_jet)-1))
    h_jet_t_phi    = np.zeros((len(bins_phi)-1,     len(bins_time_jet)-1))
 
    branches = [
        "caloCell_eta",        "caloCell_phi",        "caloCell_ecalTime",
        "caloCell_ieta",       "caloCell_iphi",
        "caloCellCleaned_eta", "caloCellCleaned_phi", "caloCellCleaned_ecalTime",
        "caloCellCleaned_ieta","caloCellCleaned_iphi",
        "pfJet_eta",           "pfJet_phi",           "pfJet_weightedTime",
        "pfJet_nCell",
    ]
 
    tree = uproot.open(f"{args.input}:{args.tree}")
    n_entries = tree.num_entries
    print(f"Total entries: {n_entries:,}")
 
    n_eb_bad = 0
    n_cl_bad = 0
    n_jet_bad = 0
 
    for start in range(0, n_entries, args.chunk_size):
        stop = min(start + args.chunk_size, n_entries)
        print(f"  {stop:,} / {n_entries:,} ({100*stop/n_entries:.1f}%)", flush=True)
        chunk = tree.arrays(branches, entry_start=start, entry_stop=stop, library="ak")
 
        # ── EB cells — bad window ──
        c_eta   = ak.to_numpy(ak.flatten(chunk["caloCell_eta"]))
        c_phi   = ak.to_numpy(ak.flatten(chunk["caloCell_phi"]))
        c_time  = ak.to_numpy(ak.flatten(chunk["caloCell_ecalTime"]))
        c_ieta  = ak.to_numpy(ak.flatten(chunk["caloCell_ieta"]))
        c_iphi  = ak.to_numpy(ak.flatten(chunk["caloCell_iphi"]))
        eb_bad  = (c_time > t_lo) & (c_time < t_hi)
        n_eb_bad += eb_bad.sum()
 
        h_eb_ieta_iphi += np.histogram2d(c_ieta[eb_bad], c_iphi[eb_bad],
                                          bins=[bins_ieta, bins_iphi])[0]
        h_eb_eta_phi   += np.histogram2d(c_eta[eb_bad],  c_phi[eb_bad],
                                          bins=[bins_cell_eta, bins_phi])[0]
        h_eb_t_eta     += np.histogram2d(c_eta[eb_bad],  c_time[eb_bad],
                                          bins=[bins_cell_eta, bins_time_cell])[0]
        h_eb_t_phi     += np.histogram2d(c_phi[eb_bad],  c_time[eb_bad],
                                          bins=[bins_phi, bins_time_cell])[0]
 
        # ── EBCleaned cells — bad window ──
        cl_eta   = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_eta"]))
        cl_phi   = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_phi"]))
        cl_time  = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_ecalTime"]))
        cl_ieta  = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_ieta"]))
        cl_iphi  = ak.to_numpy(ak.flatten(chunk["caloCellCleaned_iphi"]))
        cl_bad   = (cl_time > t_lo) & (cl_time < t_hi)
        n_cl_bad += cl_bad.sum()
 
        h_cl_ieta_iphi += np.histogram2d(cl_ieta[cl_bad], cl_iphi[cl_bad],
                                          bins=[bins_ieta, bins_iphi])[0]
        h_cl_eta_phi   += np.histogram2d(cl_eta[cl_bad],  cl_phi[cl_bad],
                                          bins=[bins_cell_eta, bins_phi])[0]
        h_cl_t_eta     += np.histogram2d(cl_eta[cl_bad],  cl_time[cl_bad],
                                          bins=[bins_cell_eta, bins_time_cell])[0]
        h_cl_t_phi     += np.histogram2d(cl_phi[cl_bad],  cl_time[cl_bad],
                                          bins=[bins_phi, bins_time_cell])[0]
 
        # ── Jets — bad window (EB timing) ──
        j_eta   = ak.to_numpy(ak.flatten(chunk["pfJet_eta"]))
        j_phi   = ak.to_numpy(ak.flatten(chunk["pfJet_phi"]))
        j_time  = ak.to_numpy(ak.flatten(chunk["pfJet_weightedTime"]))
        j_ncell = ak.to_numpy(ak.flatten(chunk["pfJet_nCell"])).astype(int)
        j_bad   = (j_ncell > 0) & (j_time > t_lo) & (j_time < t_hi)
        n_jet_bad += j_bad.sum()
 
        h_jet_eta_phi  += np.histogram2d(j_eta[j_bad],  j_phi[j_bad],
                                          bins=[bins_jet_eta, bins_phi])[0]
        h_jet_t_eta    += np.histogram2d(j_eta[j_bad],  j_time[j_bad],
                                          bins=[bins_jet_eta, bins_time_jet])[0]
        h_jet_t_phi    += np.histogram2d(j_phi[j_bad],  j_time[j_bad],
                                          bins=[bins_phi, bins_time_jet])[0]
 
        del chunk
        del c_eta, c_phi, c_time, c_ieta, c_iphi, eb_bad
        del cl_eta, cl_phi, cl_time, cl_ieta, cl_iphi, cl_bad
        del j_eta, j_phi, j_time, j_ncell, j_bad
 
    print(f"\nBad hits ({t_lo} < t < {t_hi} ns):")
    print(f"  EB cells:     {n_eb_bad:,}")
    print(f"  EBCleaned:    {n_cl_bad:,}")
    print(f"  Jets (EB t):  {n_jet_bad:,}")
 
    window_str = f"({t_lo} < t < {t_hi} ns)"
 
    # ── EB cell plots ─────────────────────────────────────────────────────
    plot_ieta_iphi(h_eb_ieta_iphi, bins_ieta, bins_iphi,
                   f"EB cell occupancy by crystal index {window_str}",
                   args.output_dir / "eb_bad_ieta_iphi.png")
 
    plot_eta_phi(h_eb_eta_phi, bins_cell_eta, bins_phi,
                 f"EB cell η-φ occupancy {window_str}",
                 args.output_dir / "eb_bad_eta_phi.png")
 
    plot_time_vs_eta(h_eb_t_eta, bins_cell_eta, bins_time_cell,
                     f"EB cell time vs η {window_str}",
                     args.output_dir / "eb_bad_time_vs_eta.png")
 
    plot_time_vs_phi(h_eb_t_phi, bins_phi, bins_time_cell,
                     f"EB cell time vs φ {window_str}",
                     args.output_dir / "eb_bad_time_vs_phi.png")
 
    # ── EBCleaned cell plots ───────────────────────────────────────────────
    plot_ieta_iphi(h_cl_ieta_iphi, bins_ieta, bins_iphi,
                   f"EBCleaned cell occupancy by crystal index {window_str}",
                   args.output_dir / "ebcleaned_bad_ieta_iphi.png")
 
    plot_eta_phi(h_cl_eta_phi, bins_cell_eta, bins_phi,
                 f"EBCleaned cell η-φ occupancy {window_str}",
                 args.output_dir / "ebcleaned_bad_eta_phi.png")
 
    plot_time_vs_eta(h_cl_t_eta, bins_cell_eta, bins_time_cell,
                     f"EBCleaned cell time vs η {window_str}",
                     args.output_dir / "ebcleaned_bad_time_vs_eta.png")
 
    plot_time_vs_phi(h_cl_t_phi, bins_phi, bins_time_cell,
                     f"EBCleaned cell time vs φ {window_str}",
                     args.output_dir / "ebcleaned_bad_time_vs_phi.png")
 
    # ── Jet plots ─────────────────────────────────────────────────────────
    plot_eta_phi(h_jet_eta_phi, bins_jet_eta, bins_phi,
                 f"Jet η-φ occupancy {window_str}",
                 args.output_dir / "jet_bad_eta_phi.png")
 
    plot_time_vs_eta(h_jet_t_eta, bins_jet_eta, bins_time_jet,
                     f"Jet time vs η {window_str}",
                     args.output_dir / "jet_bad_time_vs_eta.png")
 
    plot_time_vs_phi(h_jet_t_phi, bins_phi, bins_time_jet,
                     f"Jet time vs φ {window_str}",
                     args.output_dir / "jet_bad_time_vs_phi.png")
 
    print(f"\nAll plots written to {args.output_dir}/")
 
if __name__ == "__main__":
    main()