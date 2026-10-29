# ScoutingTiming

A lightweight [CMSSW](https://github.com/cms-sw/cmssw) analyzer for extracting **ECAL timing information from Run 3 PF scouting jets**. The package reads Run 3 scouting PF jets and scouting EB rechits, associates rechits to jets, computes ET-weighted mean jet times in three variants, and writes the result — along with cell-level information and delayed-jet trigger decisions — to a flat ROOT tree for offline long-lived-particle (LLP) searches.

## Overview

The core of the package is the `ScoutingTimingAnalyzer` EDAnalyzer (`plugins/ScoutingTimingAnalyzer.cc`). For each event it:

1. **Reads scouting collections** — `Run3ScoutingPFJet` jets and `Run3ScoutingEBRecHit` ECAL barrel rechits from two instance labels produced by the HLT scouting packers:
   - `EB` — hits passing PF quality tests (good hits, spikes already removed upstream). Primary collection for jet timing.
   - `EBCleaned` — hits flagged `kWeird` or `kDiWeird` (spikes removed by PF cleaning). Stored for spike-rate studies.

2. **Stores per-rechit information** — for every rechit above 0.5 GeV, it records η, φ (from `CaloGeometry`), energy, ECAL time, the raw `EcalRecHit` flag bitfield, and the integer crystal indices `ieta` (−85..+85) and `iphi` (1..360) from `EBDetId`. Both EB and EBCleaned collections are stored.

3. **Computes per-jet timestamps in three variants** for each PF jet, looping over rechits within ΔR < 0.4:
   - **EB** (`pfJet_weightedTime`) — ET-weighted mean time over all EB hits (primary).
   - **EBCleaned** (`pfJet_weightedTime_EBCleaned`) — ET-weighted mean time over spike hits (for studies).
   - **kGood-only** (`pfJet_weightedTime_kGood`) — ET-weighted mean time over EB hits with `flags == 0x1` only, excluding gain-switched hits (`kHasSwitchToGain1`) which carry a systematic bias.

4. **Computes per-jet time RMS** (`pfJet_tRMS`, `pfJet_tRMS_EBCleaned`, `pfJet_tRMS_kGood`) — ET-weighted standard deviation of crystal times within ΔR < 0.4. Used for the `tRMS` jet quality study.

5. **Records trigger decisions** — checks whether `HLT_HT430_DelayedJet40_SingleDelay2nsInclusive` fired in both the standard and re-run `TriggerResults`.

6. **Fills a TTree** — all quantities written to `timeTree` via `TFileService`, one entry per event.

## Flag notes

Flags in `caloCell_flags` are the full 32-bit `EcalRecHit::flagBits`, confirmed by expert tracing the full chain. Active flags in the EB collection:

| Flag | Bit | Notes |
|------|-----|-------|
| `kGood` | 0 | Set on ~100% of EB hits |
| `kPoorCalib` | 5 | ~0.01%; kGood stays set |
| `kSaturated` | 6 | High-energy deposits |
| `kHasSwitchToGain6` | 16 | ~0.8%; informational |
| `kHasSwitchToGain1` | 17 | ~0.003%; **causes timing bias** |
| `kWeird` | 14 | Spike hits — routed to EBCleaned only |
| `kDiWeird` | 15 | Spike hits — routed to EBCleaned only |

`kOutOfTime` (bit 2) is effectively never set at HLT (threshold = 1000σ). Cut on `time()` directly for delayed hit selection. A clean hit has `flags == 0x1` (kGood only).

## Output tree

The output ROOT file contains a single tree `timeTree` with one entry per event.

### Calo cell branches — EB (good hits)

| Branch | Type | Description |
|--------|------|-------------|
| `caloCell_eta`, `caloCell_phi` | `vector<double>` | Cell position from CaloGeometry |
| `caloCell_e` | `vector<double>` | Cell energy (GeV) |
| `caloCell_ecalTime` | `vector<double>` | Reconstructed ECAL time (ns) |
| `caloCell_flags` | `vector<uint32_t>` | Raw EcalRecHit flag bitfield |
| `caloCell_ieta` | `vector<int>` | EBDetId integer iη (−85..+85, skip 0) |
| `caloCell_iphi` | `vector<int>` | EBDetId integer iφ (1..360) |

### Calo cell branches — EBCleaned (spike hits)

| Branch | Type | Description |
|--------|------|-------------|
| `caloCellCleaned_eta`, `caloCellCleaned_phi` | `vector<double>` | Cell position |
| `caloCellCleaned_e` | `vector<double>` | Cell energy (GeV) |
| `caloCellCleaned_ecalTime` | `vector<double>` | ECAL time (ns); flat/random — spikes have no real pulse |
| `caloCellCleaned_flags` | `vector<uint32_t>` | kWeird (bit 14) and/or kDiWeird (bit 15) |
| `caloCellCleaned_ieta` | `vector<int>` | EBDetId integer iη |
| `caloCellCleaned_iphi` | `vector<int>` | EBDetId integer iφ |

### PF jet kinematic branches

| Branch | Type | Description |
|--------|------|-------------|
| `pfJet_pt`, `pfJet_eta`, `pfJet_phi`, `pfJet_e` | `vector<double>` | Jet kinematics |
| `pfJet_chargedHadEnergy` | `vector<double>` | Charged hadron energy |
| `pfJet_neutralHadEnergy` | `vector<double>` | Neutral hadron energy |
| `pfJet_chargedEmEnergy` | `vector<double>` | Charged EM energy (electrons) |
| `pfJet_neutralEmEnergy` | `vector<double>` | Neutral EM energy (photons) |

### Jet timing branches — EB (primary)

| Branch | Type | Description |
|--------|------|-------------|
| `pfJet_weightedTime` | `vector<double>` | ET-weighted mean time; −200 ns sentinel if nCell=0 |
| `pfJet_totalPtCell` | `vector<double>` | Summed cell ET used in weighting |
| `pfJet_nCell` | `vector<uint32_t>` | Number of matched EB crystals |
| `pfJet_tRMS` | `vector<double>` | ET-weighted time RMS; −1 ns sentinel if nCell=0 |

### Jet timing branches — EBCleaned (spike studies)

| Branch | Type | Description |
|--------|------|-------------|
| `pfJet_weightedTime_EBCleaned` | `vector<double>` | ET-weighted mean time over spike hits |
| `pfJet_totalPtCell_EBCleaned` | `vector<double>` | Summed cell ET |
| `pfJet_nCell_EBCleaned` | `vector<uint32_t>` | Number of matched spike crystals |
| `pfJet_tRMS_EBCleaned` | `vector<double>` | ET-weighted time RMS |

### Jet timing branches — kGood-only (flags == 0x1)

Excludes gain-switched hits (`kHasSwitchToGain1`) which carry a ~−5 ns systematic bias confirmed by data.

| Branch | Type | Description |
|--------|------|-------------|
| `pfJet_weightedTime_kGood` | `vector<double>` | ET-weighted mean time, kGood-only hits |
| `pfJet_totalPtCell_kGood` | `vector<double>` | Summed cell ET |
| `pfJet_nCell_kGood` | `vector<uint32_t>` | Number of matched kGood-only crystals |
| `pfJet_tRMS_kGood` | `vector<double>` | ET-weighted time RMS |

### Trigger branches

| Branch | Type | Description |
|--------|------|-------------|
| `delayedJetPathPass` | `bool` | `HLT_HT430_DelayedJet40_SingleDelay2nsInclusive` in standard TriggerResults |
| `scoutingJetPathPass` | `bool` | Same path in re-run TriggerResults |

## Repository layout

```
plugins/
  ScoutingTimingAnalyzer.cc       # the EDAnalyzer
  BuildFile.xml                   # build dependencies
test/
  config_v2.py                    # cmsRun configuration (Run 2025E)
  filelist_2025E_HLTScout_PF0.txt # full Run 2025E PF0 file list
  filelist_2025E_HLTScout_PF1.txt # full Run 2025E PF1 file list
  filelist_validation.txt         # small validation subset
  plot_timing_v2.py               # 1D timing + timing vs cell/jet properties
  plot_angular_v2.py              # timing vs η/φ angular maps (no cuts)
  plot_bad_channels.py            # bad channel occupancy maps (time window cut)
  plot_timing.py                  # legacy timing plots (superseded by v2)
  plot_angular.py                 # legacy angular plots (superseded by v2)
  BuildFile.xml                   # catch2 unit-test build
  test_catch2_main.cc             # test harness entry point
  test_catch2_ScoutingTimingAnalyzer.cc
```

## Setup

```bash
cmsrel CMSSW_15_0_22
cd CMSSW_15_0_22/src
cmsenv

git clone https://github.com/yansonhuang/ScoutingTiming.git ScoutingTiming/ScoutingTiming
cd ScoutingTiming/ScoutingTiming
git checkout developing_versions

cd $CMSSW_BASE/src
scram b -j8
edmPluginDump | grep ScoutingTiming  # should appear
```

**Global tag:** `150X_dataRun3_HLT_v1` — confirmed via `edmProvDump` on a Run 2025E HLTSCOUT file.

**For 2026 data:** use `CMSSW_16_0_7` with the appropriate global tag (verify via `edmProvDump`).

## Running

```bash
export XRD_NETWORKSTACK=IPv4
voms-proxy-init --voms cms --valid 192:00
cmsenv

cmsRun ScoutingTiming/ScoutingTiming/test/config_v2.py
```

The config reads from `filelist_validation.txt` by default. Output is written to `test_6files_v2.root` locally; copy to EOS afterwards with `xrdcp`.

### Input tags (Run 2025E)

| Parameter | Input tag |
|-----------|-----------|
| `pfJetsTag` | `hltScoutingPFPacker::HLT` |
| `ebRecHitsTag` | `hltScoutingRecHitPacker:EB:HLT` |
| `ebCleanedRecHitsTag` | `hltScoutingRecHitPacker:EBCleaned:HLT` |
| `triggerResultsTag` | `TriggerResults::HLT` |
| `triggerResultsRerunTag` | `TriggerResults::HLT` |

## Plotting

All scripts use `uproot`, `awkward`, `numpy`, and `mplhep`. Run from the `test/` directory:

```bash
# 1D timing + timing vs flags/energy/jet kinematics
python3 plot_timing_v2.py test_6files_v2.root --output-dir outputs/

# Timing vs η/φ angular maps (no cuts)
python3 plot_angular_v2.py test_6files_v2.root --output-dir outputs/

# With kGood-only flag filter applied to EB cells
python3 plot_angular_v2.py test_6files_v2.root --flags-only-good --output-dir outputs_flagclean/

# Bad channel occupancy maps (time window: -5.5 < t < -4.5 ns by default)
python3 plot_bad_channels.py test_6files_v2.root --output-dir outputs/
python3 plot_bad_channels.py test_6files_v2.root --t-lo -5.5 --t-hi -4.5 --output-dir outputs/

# With kGood-only filter (check if hot channels are flag-tagged)
python3 plot_bad_channels.py test_6files_v2.root --flags-only-good --output-dir outputs_flagclean/
```

## Dataset

Two orthogonal scouting streams for Run 2025E (safe to combine):

```
/ScoutingPF0/Run2025E-v1/HLTSCOUT   ~31,525 files, ~11.2B events
/ScoutingPF1/Run2025E-v1/HLTSCOUT   ~31,521 files, ~11.2B events
```

File lists in `test/filelist_2025E_HLTScout_PF0.txt` and `PF1.txt`.

## Notes

- Only ECAL barrel (EB) rechits are used. A 0.5 GeV energy threshold is applied before any crystal enters the tree or the jet-time calculation.
- Jets with no matched crystals receive a sentinel value of −200 ns (time) and −1 ns (tRMS).
- The `kHasSwitchToGain1` flag (bit 17) causes a systematic ~−5 ns timing bias confirmed in Run 2025E data. Use `pfJet_weightedTime_kGood` or require `caloCell_flags == 0x1` offline to exclude these hits.
- The `tRMS < 1.5 ns` jet quality cut from AN-24-069 can be applied directly on `pfJet_tRMS` or `pfJet_tRMS_kGood`.
- The trigger match uses a substring test on the path name, picking up any version suffix.
- ROOT output files should be written locally (not to EOS XRootD paths) and copied afterwards with `xrdcp`.

## Author

Matthew Daniel Citron (original)  
Yansong Huang (Run 2025E development)
