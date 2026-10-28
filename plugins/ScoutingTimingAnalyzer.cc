// -*- C++ -*-
//
// Package:    ScoutingTiming/ScoutingTimingAnalyzer
// Class:      ScoutingTimingAnalyzer
//
/**\class ScoutingTimingAnalyzer ScoutingTimingAnalyzer.cc ScoutingTiming/ScoutingTimingAnalyzer/plugins/ScoutingTimingAnalyzer.cc

 Description: Compute per-PF-jet ECAL timing from Run3 scouting EB rec hits.

 Three versions of ET-weighted jet timing are stored for comparison:
   1. Raw:        from "EB" collection, no flag cuts (original behaviour)
   2. FlagClean:  from "EB" collection, kWeird (bit 14) and kDiWeird (bit 15) vetoed
   3. EBCleaned:  from "EBCleaned" collection (spike-cleaned at HLT packer level)

 For each version, matched crystals within dR < 0.4 and E > 0.5 GeV are used.
 Jets with no matched crystals are assigned <t> = -200 ns as sentinel.
 Cell-level info is stored from the raw EB collection so flag-based selections
 can be applied offline via (caloCell_flags >> bit) & 0x1.
*/
//
// Original Author:  Matthew Daniel Citron
//         Created:  Sat, 09 Aug 2025 04:38:14 GMT
// Updated: added EBCleaned and flag-cleaned timing versions
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetupRecord.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingEBRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"

//
// class declaration
//

using reco::TrackCollection;

class ScoutingTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingTimingAnalyzer(const edm::ParameterSet&);
  ~ScoutingTimingAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override {}
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  // ---------- helper ---------------------------------------------------------
  // Compute ET-weighted mean time for one jet over a rechit collection,
  // optionally applying a spike-flag veto.
  // Returns -200.0 if no crystals matched.
  double computeWeightedTime(const Run3ScoutingPFJet& jet,
                             const std::vector<Run3ScoutingEBRecHit>& hits,
                             const CaloGeometry& geom,
                             bool vetoSpikes,
                             double& totalPtOut,
                             uint32_t& nCellOut) const;

  // ---------- tokens ---------------------------------------------------------
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord>    esToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingPFJet>>     pfJetsToken_;
  // Raw EB rechits (all flags preserved)
  edm::EDGetTokenT<std::vector<Run3ScoutingEBRecHit>>  ebRecHitsToken_;
  // HLT-packer spike-cleaned EB rechits
  edm::EDGetTokenT<std::vector<Run3ScoutingEBRecHit>>  ebCleanedRecHitsToken_;
  edm::EDGetTokenT<edm::TriggerResults>                triggerResultsToken_;
  edm::EDGetTokenT<edm::TriggerResults>                triggerResultsRerunToken_;

  // ---------- output tree ----------------------------------------------------
  TTree* timeTree_;

  // ---- calo cell branches (raw EB only — flags stored for offline use) ----
  std::vector<double>   v_caloCellPhi_;
  std::vector<double>   v_caloCellEta_;
  std::vector<double>   v_caloCellEnergy_;
  std::vector<double>   v_caloCellEcalTime_;
  std::vector<uint32_t> v_caloCellFlags_;

  // ---- PF jet kinematics (common to all timing versions) ------------------
  std::vector<double>   v_pfJetE_;
  std::vector<double>   v_pfJetPt_;
  std::vector<double>   v_pfJetPhi_;
  std::vector<double>   v_pfJetEta_;
  std::vector<double>   v_pfJetChargedHadEnergy_;
  std::vector<double>   v_pfJetNeutralHadEnergy_;
  std::vector<double>   v_pfJetChargedEmEnergy_;
  std::vector<double>   v_pfJetNeutralEmEnergy_;

  // ---- Version 1: raw EB (no flag cuts) ------------------------------------
  std::vector<double>   v_pfJetWeightedTime_;       // ET-weighted time
  std::vector<double>   v_pfJetTotalPtCell_;         // sum pT of matched cells
  std::vector<uint32_t> v_pfJetNCell_;               // number of matched cells

  // ---- Version 2: flag-cleaned EB (kWeird + kDiWeird vetoed) --------------
  std::vector<double>   v_pfJetWeightedTimeFlagClean_;
  std::vector<double>   v_pfJetTotalPtCellFlagClean_;
  std::vector<uint32_t> v_pfJetNCellFlagClean_;

  // ---- Version 3: EBCleaned collection (HLT packer spike removal) ---------
  std::vector<double>   v_pfJetWeightedTimeEBCleaned_;
  std::vector<double>   v_pfJetTotalPtCellEBCleaned_;
  std::vector<uint32_t> v_pfJetNCellEBCleaned_;

  // ---- trigger decisions ---------------------------------------------------
  bool b_delayedJetPathPass_;
  bool b_scoutingJetPathPass_;
};

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
ScoutingTimingAnalyzer::ScoutingTimingAnalyzer(const edm::ParameterSet& iPSet)
    : esToken_(esConsumes()),
      pfJetsToken_(consumes(iPSet.getParameter<edm::InputTag>("pfJetsTag"))),
      ebRecHitsToken_(consumes(iPSet.getParameter<edm::InputTag>("ebRecHitsTag"))),
      ebCleanedRecHitsToken_(consumes(iPSet.getParameter<edm::InputTag>("ebCleanedRecHitsTag"))),
      triggerResultsToken_(consumes(iPSet.getParameter<edm::InputTag>("triggerResultsTag"))),
      triggerResultsRerunToken_(consumes(iPSet.getParameter<edm::InputTag>("triggerResultsRerunTag")))
{
  edm::Service<TFileService> fs;
  timeTree_ = fs->make<TTree>("timeTree", "timeTree");

  // calo cell branches (raw EB)
  timeTree_->Branch("caloCell_eta",      &v_caloCellEta_);
  timeTree_->Branch("caloCell_phi",      &v_caloCellPhi_);
  timeTree_->Branch("caloCell_e",        &v_caloCellEnergy_);
  timeTree_->Branch("caloCell_ecalTime", &v_caloCellEcalTime_);
  timeTree_->Branch("caloCell_flags",    &v_caloCellFlags_);

  // PF jet kinematics
  timeTree_->Branch("pfJet_e",                &v_pfJetE_);
  timeTree_->Branch("pfJet_pt",               &v_pfJetPt_);
  timeTree_->Branch("pfJet_phi",              &v_pfJetPhi_);
  timeTree_->Branch("pfJet_eta",              &v_pfJetEta_);
  timeTree_->Branch("pfJet_chargedHadEnergy", &v_pfJetChargedHadEnergy_);
  timeTree_->Branch("pfJet_neutralHadEnergy", &v_pfJetNeutralHadEnergy_);
  timeTree_->Branch("pfJet_chargedEmEnergy",  &v_pfJetChargedEmEnergy_);
  timeTree_->Branch("pfJet_neutralEmEnergy",  &v_pfJetNeutralEmEnergy_);

  // Version 1: raw EB
  timeTree_->Branch("pfJet_weightedTime",     &v_pfJetWeightedTime_);
  timeTree_->Branch("pfJet_totalPtCell",      &v_pfJetTotalPtCell_);
  timeTree_->Branch("pfJet_nCell",            &v_pfJetNCell_);

  // Version 2: flag-cleaned EB
  timeTree_->Branch("pfJet_weightedTime_flagClean",  &v_pfJetWeightedTimeFlagClean_);
  timeTree_->Branch("pfJet_totalPtCell_flagClean",   &v_pfJetTotalPtCellFlagClean_);
  timeTree_->Branch("pfJet_nCell_flagClean",         &v_pfJetNCellFlagClean_);

  // Version 3: EBCleaned
  timeTree_->Branch("pfJet_weightedTime_EBCleaned",  &v_pfJetWeightedTimeEBCleaned_);
  timeTree_->Branch("pfJet_totalPtCell_EBCleaned",   &v_pfJetTotalPtCellEBCleaned_);
  timeTree_->Branch("pfJet_nCell_EBCleaned",         &v_pfJetNCellEBCleaned_);

  // trigger
  timeTree_->Branch("delayedJetPathPass",  &b_delayedJetPathPass_);
  timeTree_->Branch("scoutingJetPathPass", &b_scoutingJetPathPass_);
}

// ---------------------------------------------------------------------------
// Helper: ET-weighted timing for one jet
// ---------------------------------------------------------------------------
double ScoutingTimingAnalyzer::computeWeightedTime(
    const Run3ScoutingPFJet& jet,
    const std::vector<Run3ScoutingEBRecHit>& hits,
    const CaloGeometry& geom,
    bool vetoSpikes,
    double& totalPtOut,
    uint32_t& nCellOut) const
{
  TLorentzVector jetVec;
  jetVec.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);

  double weightedTime = 0;
  totalPtOut = 0;
  nCellOut   = 0;

  for (auto const& hit : hits) {
    if (hit.energy() < 0.5) continue;

    // Spike veto: skip kWeird (bit 14) and kDiWeird (bit 15)
    if (vetoSpikes) {
      if (hit.flags() & (1u << 14)) continue; // kWeird
      if (hit.flags() & (1u << 15)) continue; // kDiWeird
    }

    GlobalPoint pCell = geom.getPosition(hit.detId());

    TLorentzVector cellVec;
    cellVec.SetPtEtaPhiM(1, pCell.eta(), pCell.phi(), 0);
    if (cellVec.DeltaR(jetVec) > 0.4) continue;

    const double cellPt = hit.energy() * TMath::Sin(pCell.theta());
    weightedTime += hit.time() * cellPt;
    totalPtOut   += cellPt;
    ++nCellOut;
  }

  if (totalPtOut > 0)
    return weightedTime / totalPtOut;
  else
    return -200.0;
}

// ---------------------------------------------------------------------------
// analyze
// ---------------------------------------------------------------------------
void ScoutingTimingAnalyzer::analyze(const edm::Event& iEvent,
                                     const edm::EventSetup& iSetup) {
  using namespace edm;

  auto const& geom              = iSetup.getData(esToken_);
  auto const& pfJets            = iEvent.get(pfJetsToken_);
  auto const& ebRecHits         = iEvent.get(ebRecHitsToken_);
  auto const& ebCleanedRecHits  = iEvent.get(ebCleanedRecHitsToken_);
  auto const& triggerResults    = iEvent.get(triggerResultsToken_);
  auto const& triggerResultsRerun = iEvent.get(triggerResultsRerunToken_);

  // ---- clear per-event containers ------------------------------------------
  v_caloCellPhi_.clear();
  v_caloCellEta_.clear();
  v_caloCellEnergy_.clear();
  v_caloCellEcalTime_.clear();
  v_caloCellFlags_.clear();

  v_pfJetE_.clear();
  v_pfJetPt_.clear();
  v_pfJetPhi_.clear();
  v_pfJetEta_.clear();
  v_pfJetChargedHadEnergy_.clear();
  v_pfJetNeutralHadEnergy_.clear();
  v_pfJetChargedEmEnergy_.clear();
  v_pfJetNeutralEmEnergy_.clear();

  v_pfJetWeightedTime_.clear();
  v_pfJetTotalPtCell_.clear();
  v_pfJetNCell_.clear();

  v_pfJetWeightedTimeFlagClean_.clear();
  v_pfJetTotalPtCellFlagClean_.clear();
  v_pfJetNCellFlagClean_.clear();

  v_pfJetWeightedTimeEBCleaned_.clear();
  v_pfJetTotalPtCellEBCleaned_.clear();
  v_pfJetNCellEBCleaned_.clear();

  // ---- trigger decisions ---------------------------------------------------
  b_delayedJetPathPass_  = false;
  b_scoutingJetPathPass_ = false;

  {
    const TriggerNames& names = iEvent.triggerNames(triggerResults);
    for (unsigned int i = 0; i < triggerResults.size(); ++i) {
      if (names.triggerName(i).find(
              "HLT_HT430_DelayedJet40_SingleDelay2nsInclusive_v") != std::string::npos)
        b_delayedJetPathPass_ = triggerResults.accept(i);
    }
  }
  {
    const TriggerNames& names = iEvent.triggerNames(triggerResultsRerun);
    for (unsigned int i = 0; i < triggerResultsRerun.size(); ++i) {
      if (names.triggerName(i).find(
              "HLT_HT430_DelayedJet40_SingleDelay2nsInclusive_v") != std::string::npos)
        b_scoutingJetPathPass_ = triggerResultsRerun.accept(i);
    }
  }

  // ---- calo cells: store raw EB hits for offline flag studies --------------
  for (auto const& hit : ebRecHits) {
    if (hit.energy() < 0.5) continue;
    GlobalPoint pCell = geom.getPosition(hit.detId());
    v_caloCellPhi_.push_back(pCell.phi());
    v_caloCellEta_.push_back(pCell.eta());
    v_caloCellEnergy_.push_back(hit.energy());
    v_caloCellEcalTime_.push_back(hit.time());
    v_caloCellFlags_.push_back(hit.flags());
  }

  // ---- PF jets: three timing versions per jet ------------------------------
  for (auto const& jet : pfJets) {

    double totalPt1 = 0, totalPt2 = 0, totalPt3 = 0;
    uint32_t nCell1 = 0,  nCell2 = 0,  nCell3 = 0;

    // Version 1: raw EB, no flag cuts
    double t1 = computeWeightedTime(jet, ebRecHits,        geom, false, totalPt1, nCell1);
    // Version 2: raw EB, kWeird + kDiWeird vetoed
    double t2 = computeWeightedTime(jet, ebRecHits,        geom, true,  totalPt2, nCell2);
    // Version 3: EBCleaned collection (spikes already removed at packer level)
    double t3 = computeWeightedTime(jet, ebCleanedRecHits, geom, false, totalPt3, nCell3);

    // kinematics
    v_pfJetE_.push_back(jet.chargedHadronEnergy() + jet.neutralHadronEnergy()
                        + jet.electronEnergy()    + jet.photonEnergy());
    v_pfJetPt_.push_back(jet.pt());
    v_pfJetPhi_.push_back(jet.phi());
    v_pfJetEta_.push_back(jet.eta());
    v_pfJetChargedHadEnergy_.push_back(jet.chargedHadronEnergy());
    v_pfJetNeutralHadEnergy_.push_back(jet.neutralHadronEnergy());
    v_pfJetChargedEmEnergy_.push_back(jet.electronEnergy());
    v_pfJetNeutralEmEnergy_.push_back(jet.photonEnergy());

    // timing v1
    v_pfJetWeightedTime_.push_back(t1);
    v_pfJetTotalPtCell_.push_back(totalPt1);
    v_pfJetNCell_.push_back(nCell1);

    // timing v2
    v_pfJetWeightedTimeFlagClean_.push_back(t2);
    v_pfJetTotalPtCellFlagClean_.push_back(totalPt2);
    v_pfJetNCellFlagClean_.push_back(nCell2);

    // timing v3
    v_pfJetWeightedTimeEBCleaned_.push_back(t3);
    v_pfJetTotalPtCellEBCleaned_.push_back(totalPt3);
    v_pfJetNCellEBCleaned_.push_back(nCell3);
  }

  timeTree_->Fill();
}

// ---------------------------------------------------------------------------
// fillDescriptions
// ---------------------------------------------------------------------------
void ScoutingTimingAnalyzer::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfJetsTag");
  desc.add<edm::InputTag>("ebRecHitsTag");
  desc.add<edm::InputTag>("ebCleanedRecHitsTag");
  desc.add<edm::InputTag>("triggerResultsTag");
  desc.add<edm::InputTag>("triggerResultsRerunTag");
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingTimingAnalyzer);
