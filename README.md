# ScoutingTiming

A lightweight [CMSSW](https://github.com/cms-sw/cmssw) analyzer for extracting **ECAL timing information from Run 3 PF scouting jets**. The package reads Run 3 scouting PF jets and scouting EB rechits, associates rechits to jets, computes an ET-weighted mean jet time, and writes the result — along with cell-level information and delayed-jet trigger decisions — to a flat ROOT tree for offline analysis of long-lived-particle signatures.

## Overview

The core of the package is the `ScoutingTimingAnalyzer` EDAnalyzer (`plugins/ScoutingTimingAnalyzer.cc`). For each event it:

1. **Reads scouting collections** — `Run3ScoutingPFJet` jets and `Run3ScoutingEBRecHit` ECAL barrel rechits produced by the HLT scouting packers.
2. **Stores per-rechit information** — for every EB rechit above 0.5 GeV, it records the cell position (η, φ) from the calo geometry, energy, ECAL time, and the raw `EcalRecHit` status-flag bitfield. Individual flags are intentionally *not* pre-unpacked; decode them offline via `(flags >> bit) & 0x1`.
3. **Computes a per-jet timestamp** — for each PF jet, it loops over rechits within ΔR < 0.4 of the jet axis and forms an ET-weighted mean cell time (`pfJet_weightedTime`), along with the summed cell ET and the number of contributing cells. No flag-based cell rejection is applied at this stage, so tighter selections can be made offline.
4. **Records trigger decisions** — it checks whether the delayed-jet HLT path (`HLT_HT430_DelayedJet40_SingleDelay2nsInclusive`) fired, in both the standard `TriggerResults` and a re-run `TriggerResults` collection, and stores both decisions as branches for standard-vs-scouting comparisons.
5. **Fills a TTree** — all quantities are written to a tree named `timeTree` via `TFileService`.

## Output tree

The output ROOT file (default `test1.root`) contains a single tree, `timeTree`, with one entry per event and the following branches.

### Calo cell (EB rechit) branches
| Branch | Description |
| --- | --- |
| `caloCell_eta`, `caloCell_phi` | Cell position from the calo geometry |
| `caloCell_e` | Cell energy (GeV) |
| `caloCell_ecalTime` | Reconstructed ECAL time (ns) |
| `caloCell_flags` | Raw `EcalRecHit` status-flag bitfield — decode offline with `(flags >> bit) & 0x1` (e.g. `kGood`, `kOutOfTime`, `kWeird`, `kSaturated`, gain-switch flags) |

### PF jet branches
| Branch | Description |
| --- | --- |
| `pfJet_pt`, `pfJet_eta`, `pfJet_phi`, `pfJet_e` | Jet kinematics (energy summed from charged/neutral hadronic + EM components) |
| `pfJet_weightedTime` | ET-weighted mean cell time over rechits within ΔR < 0.4 (set to −200 if no cells) |
| `pfJet_totalPtCell` | Summed cell ET used in the time weighting |
| `pfJet_nCell` | Number of rechits associated to the jet |
| `pfJet_chargedHadEnergy`, `pfJet_neutralHadEnergy`, `pfJet_chargedEmEnergy`, `pfJet_neutralEmEnergy` | Jet energy fractions |

### Trigger branches
| Branch | Description |
| --- | --- |
| `delayedJetPathPass` | Decision of `HLT_HT430_DelayedJet40_SingleDelay2nsInclusive` in the standard `TriggerResults` |
| `scoutingJetPathPass` | Decision of the same path in the re-run `TriggerResults` |

## Repository layout

```
plugins/
  ScoutingTimingAnalyzer.cc   # the EDAnalyzer
  BuildFile.xml               # build dependencies
test/
  config.py                   # example cmsRun configuration
  BuildFile.xml               # catch2 unit-test build
  test_catch2_main.cc         # test harness entry point
  test_catch2_ScoutingTimingAnalyzer.cc
```

## Setup

This package is meant to be built inside a CMSSW release area. It was developed against `CMSSW_15_0_12` with global tag `142X_mcRun3_2025_realistic_v7`. This global tag is for MC only, so maybe one wants to change to the relevant one in their analysis. Also, later HLT Scout is reconstructed using CMSSW_15_0_*, where * can be >12, while 2026 runs are reconstructed with CMSSW_16, so be careful about compatibility issues.

```bash
cmsrel CMSSW_15_0_NEWEST_RELEASE
cd CMSSW_15_0_NEWEST_RELEASE/src
cmsenv

# Place the package under a subsystem directory, e.g. ScoutingTiming/
git clone https://github.com/mcitron/ScoutingTiming.git ScoutingTiming/ScoutingTiming

scram b -j 8
```

## Running

Edit the input file path in `test/config.py`, then run:

```bash
cmsRun ScoutingTiming/ScoutingTiming/test/config.py
```

The module declares four required input tags (adjust the process names to match your input):

| Parameter | Default input tag |
| --- | --- |
| `pfJetsTag` | `hltScoutingPFPacker::HLTX` |
| `ebRecHitsTag` | `hltScoutingRecHitPacker:EB:HLTX` |
| `triggerResultsTag` | `TriggerResults::SIM` |
| `triggerResultsRerunTag` | `TriggerResults::HLTX` |

Output is written to the `TFileService` file (`test1.root` by default).

## Notes

- Only ECAL barrel (EB) rechits are used; a 0.5 GeV energy threshold is applied to every rechit before it enters the tree or the jet-time calculation.
- The jet time is an ET-weighted (E·sin θ) average of associated cell times, with no cell-level quality cuts applied in the producer — apply flag-based selections offline using `caloCell_flags`.
- The two `TriggerResults` inputs let you compare the delayed-jet decision from the original HLT menu against a re-run of the path; both are stored per event.
- The trigger match is a substring test on `HLT_HT430_DelayedJet40_SingleDelay2nsInclusive_v`, so it picks up any version of the path.
- ROOT files are stored in PF0 and PF1 two directories to split datasets. They are orthogonal, so process both to get the full luminosity
- Add two txt files containing locations of all 2025E run

## Author

Matthew Daniel Citron
