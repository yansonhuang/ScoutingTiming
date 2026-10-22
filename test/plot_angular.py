#!/usr/bin/env python3
"""Plot angular distributions of jets and matched EB crystals - memory efficient chunked version."""

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
    bins_jet_eta  = np.linspace(-5, 5, 101)
    bins_jet_phi  = np.linspace(-np.pi, np.pi, 73)
    bins_jet_eta2d = np.linspace(-3, 3, 61)
    bins_cell_eta = np.linspace(-1.5, 1.5, 101)
    bins_cell_phi = np.linspace(-np.pi, np.pi, 73)
    bins_cell_e   = np.linspace(0.5, 50, 100)

    # ── Accumulators ───────────────────────────────────────────────────────
    h_jet_eta_all     = np.zeros(len(bins_jet_eta) - 1)
    h_jet_eta_matched = np.zeros(len(bins_jet_eta) - 1)
    h_jet_phi_all     = np.zeros(len(bins_jet_phi) - 1)
    h_jet_phi_matched = np.zeros(len(bins_jet_phi) - 1)
    h_jet_eta_phi_2d  = np.zeros((len(bins_jet_eta2d) - 1, len(bins_jet_phi) - 1))
    h_cell_eta        = np.zeros(len(bins_cell_eta) - 1)
    h_cell_phi        = np.zeros(len(bins_cell_phi) - 1)
    h_cell_eta_phi_2d = np.zeros((len(bins_cell_eta) - 1, len(bins_cell_phi) - 1))
    h_cell_e          = np.zeros(len(bins_cell_e) - 1)

    n_jet_all     = 0
    n_jet_matched = 0
    n_cell_all    = 0

    branches = ["pfJet_eta", "pfJet_phi", "pfJet_nCell",
                "caloCell_eta", "caloCell_phi", "caloCell_e"]

    # ── Chunked loop ───────────────────────────────────────────────────────
    for start in range(0, n_entries, args.chunk_size):
        stop = min(start + args.chunk_size, n_entries)
        pct  = 100 * stop / n_entries
        print(f"  Processing {stop:,} / {n_entries:,} ({pct:.1f}%)", flush=True)

        chunk = tree.arrays(branches, entry_start=start, entry_stop=stop, library="ak")

        # jets
        jet_eta   = ak.to_numpy(ak.flatten(chunk["pfJet_eta"]))
        jet_phi   = ak.to_numpy(ak.flatten(chunk["pfJet_phi"]))
        jet_ncell = ak.to_numpy(ak.flatten(chunk["pfJet_nCell"])).astype(int)
        matched   = jet_ncell > 0

        n_jet_all     += len(jet_eta)
        n_jet_matched += matched.sum()

        h_jet_eta_all     += np.histogram(jet_eta,          bins=bins_jet_eta)[0]
        h_jet_eta_matched += np.histogram(jet_eta[matched], bins=bins_jet_eta)[0]
        h_jet_phi_all     += np.histogram(jet_phi,          bins=bins_jet_phi)[0]
        h_jet_phi_matched += np.histogram(jet_phi[matched], bins=bins_jet_phi)[0]
        h_jet_eta_phi_2d  += np.histogram2d(jet_eta[matched], jet_phi[matched],
                                             bins=[bins_jet_eta2d, bins_jet_phi])[0]

        # cells
        cell_eta = ak.to_numpy(ak.flatten(chunk["caloCell_eta"]))
        cell_phi = ak.to_numpy(ak.flatten(chunk["caloCell_phi"]))
        cell_e   = ak.to_numpy(ak.flatten(chunk["caloCell_e"]))

        n_cell_all += len(cell_eta)

        h_cell_eta        += np.histogram(cell_eta, bins=bins_cell_eta)[0]
        h_cell_phi        += np.histogram(cell_phi, bins=bins_cell_phi)[0]
        h_cell_eta_phi_2d += np.histogram2d(cell_eta, cell_phi,
                                             bins=[bins_cell_eta, bins_cell_phi])[0]
        h_cell_e          += np.histogram(cell_e, bins=bins_cell_e)[0]

        del chunk, jet_eta, jet_phi, jet_ncell, cell_eta, cell_phi, cell_e, matched

    print(f"Total jets:          {n_jet_all:,}")
    print(f"Jets with nCell > 0: {n_jet_matched:,}  ({100*n_jet_matched/n_jet_all:.1f}%)")
    print(f"Total EB crystals:   {n_cell_all:,}")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 1: Jet eta
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.stairs(h_jet_eta_all,     bins_jet_eta, linewidth=2, color="steelblue",
              label=f"All jets: {n_jet_all:,}")
    ax.stairs(h_jet_eta_matched, bins_jet_eta, linewidth=2, color="crimson",
              label=f"Jets with nCell > 0: {n_jet_matched:,}")
    ax.axvline( 1.48, color="gray", linestyle="--", linewidth=1.5, label="EB boundary |η|=1.48")
    ax.axvline(-1.48, color="gray", linestyle="--", linewidth=1.5)
    ax.set(xlabel="Jet $\\eta$", ylabel="Jets", xlim=(-5, 5))
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=16)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "jet_eta.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 2: Jet phi
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.stairs(h_jet_phi_all,     bins_jet_phi, linewidth=2, color="steelblue",
              label=f"All jets: {n_jet_all:,}")
    ax.stairs(h_jet_phi_matched, bins_jet_phi, linewidth=2, color="crimson",
              label=f"Jets with nCell > 0: {n_jet_matched:,}")
    ax.set(xlabel="Jet $\\phi$ [rad]", ylabel="Jets", xlim=(-np.pi, np.pi))
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
                  [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.legend(frameon=False, fontsize=16)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "jet_phi.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 3: Jet eta-phi 2D
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.pcolormesh(bins_jet_eta2d, bins_jet_phi, h_jet_eta_phi_2d.T,
                       cmap="Blues", norm=matplotlib.colors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label="Jets")
    ax.set(xlabel="Jet $\\eta$", ylabel="Jet $\\phi$ [rad]",
           xlim=(-3, 3), ylim=(-np.pi, np.pi))
    ax.set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
                  [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "jet_eta_phi_2d.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 4: Crystal eta
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.stairs(h_cell_eta, bins_cell_eta, linewidth=2, color="darkorange",
              label=f"EB rec hits (E > 0.5 GeV): {n_cell_all:,}")
    ax.set(xlabel="Crystal $\\eta$", ylabel="Crystals", xlim=(-1.5, 1.5))
    ax.legend(frameon=False, fontsize=16)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "cell_eta.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 5: Crystal phi
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.stairs(h_cell_phi, bins_cell_phi, linewidth=2, color="darkorange",
              label=f"EB rec hits (E > 0.5 GeV): {n_cell_all:,}")
    ax.set(xlabel="Crystal $\\phi$ [rad]", ylabel="Crystals", xlim=(-np.pi, np.pi))
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
                  [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.legend(frameon=False, fontsize=16)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "cell_phi.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 6: Crystal eta-phi occupancy 2D
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.pcolormesh(bins_cell_eta, bins_cell_phi, h_cell_eta_phi_2d.T,
                       cmap="Oranges", norm=matplotlib.colors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label="EB rec hits")
    ax.set(xlabel="Crystal $\\eta$", ylabel="Crystal $\\phi$ [rad]",
           xlim=(-1.5, 1.5), ylim=(-np.pi, np.pi))
    ax.set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
                  [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "cell_eta_phi_2d.png")

    # ══════════════════════════════════════════════════════════════════════
    # Plot 7: Crystal energy
    # ══════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.stairs(h_cell_e, bins_cell_e, linewidth=2, color="darkorange",
              label=f"EB rec hits: {n_cell_all:,}")
    ax.set(xlabel="Crystal energy [GeV]", ylabel="Crystals / 0.5 GeV", xlim=(0.5, 50))
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=16)
    ax.grid(alpha=0.2)
    hep.cms.label("Preliminary", data=True, lumi=None, ax=ax)
    fig.tight_layout()
    save(fig, args.output_dir / "cell_energy.png")

    print(f"\nAll plots written to {args.output_dir}/")

if __name__ == "__main__":
    main()