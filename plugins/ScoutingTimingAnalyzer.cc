// -*- C++ -*-
//
// Package:    ScoutingTiming/ScoutingTimingAnalyzer
// Class:      ScoutingTimingAnalyzer
//
/**\class ScoutingTimingAnalyzer ScoutingTimingAnalyzer.cc ScoutingTiming/ScoutingTimingAnalyzer/plugins/ScoutingTimingAnalyzer.cc
 
 Description: Compute per-PF-jet ECAL timing from Run3 scouting EB rec hits.
 
 Two collections are read and stored independently:
   EB        — hits passing PF quality tests (good hits, spikes already removed upstream).
               This is the primary collection for jet timing.
   EBCleaned — kWeird + kDiWeird spike hits removed by PF cleaning.
               Stored for spike-rate studies and to allow offline combination of
               EB + EBCleaned into the full unfiltered rechit picture.
 
 For each jet, ET-weighted timing is computed over matched crystals within
 dR < 0.4 and E > 0.5 GeV. Jets with no matched crystals are assigned
 <t> = -200 ns as sentinel.
 
 Cell-level branches are stored for both EB and EBCleaned collections.
 Flags are genuine EcalRecHit::flagBits (confirmed: copied verbatim through
 PFEcalBarrelRecHitCreator -> HLTScoutingRecHitProducer).
 In EB: active flags are kGood (bit 0), kPoorCalib (bit 5),
        kHasSwitchToGain6 (bit 16), kHasSwitchToGain1 (bit 17).
 In EBCleaned: active flags are kWeird (bit 14) and kDiWeird (bit 15).
 Note: kOutOfTime (bit 2) is effectively never set at HLT (threshold = 1000 sigma).
       Select delayed hits by cutting on time() directly.
*/
//
// Original Author:  Matthew Daniel Citron
//         Created:  Sat, 09 Aug 2025 04:38:14 GMT
// Updated: removed redundant v2 flag-cleaned version; clarified EB/EBCleaned semantics;
//          added EBCleaned cell-level branches for offline spike studies.
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
 
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingEBRecHit.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
 
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
 
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
  // Compute ET-weighted mean time for one jet over a rechit collection.
  // Returns -200.0 if no crystals matched (sentinel).
  double computeWeightedTime(const Run3ScoutingPFJet& jet,
                             const std::vector<Run3ScoutingEBRecHit>& hits,
                             const CaloGeometry& geom,
                             double& totalPtOut,
                             uint32_t& nCellOut) const;
 
  // ---------- tokens ---------------------------------------------------------
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord>   esToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingPFJet>>    pfJetsToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingEBRecHit>> ebRecHitsToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingEBRecHit>> ebCleanedRecHitsToken_;
  edm::EDGetTokenT<edm::TriggerResults>               triggerResultsToken_;
  edm::EDGetTokenT<edm::TriggerResults>               triggerResultsRerunToken_;
 
  // ---------- output tree ----------------------------------------------------
  TTree* timeTree_;
 
  // ---- calo cell branches — EB (good hits) --------------------------------
  std::vector<double>   v_caloCellEta_;
  std::vector<double>   v_caloCellPhi_;
  std::vector<double>   v_caloCellEnergy_;
  std::vector<double>   v_caloCellEcalTime_;
  std::vector<uint32_t> v_caloCellFlags_;
 
  // ---- calo cell branches — EBCleaned (spike hits) ------------------------
  std::vector<double>   v_caloCellCleanedEta_;
  std::vector<double>   v_caloCellCleanedPhi_;
  std::vector<double>   v_caloCellCleanedEnergy_;
  std::vector<double>   v_caloCellCleanedEcalTime_;
  std::vector<uint32_t> v_caloCellCleanedFlags_;
 
  // ---- PF jet kinematics ---------------------------------------------------
  std::vector<double>   v_pfJetE_;
  std::vector<double>   v_pfJetPt_;
  std::vector<double>   v_pfJetPhi_;
  std::vector<double>   v_pfJetEta_;
  std::vector<double>   v_pfJetChargedHadEnergy_;
  std::vector<double>   v_pfJetNeutralHadEnergy_;
  std::vector<double>   v_pfJetChargedEmEnergy_;
  std::vector<double>   v_pfJetNeutralEmEnergy_;
 
  // ---- jet timing — EB (primary) ------------------------------------------
  std::vector<double>   v_pfJetWeightedTime_;
  std::vector<double>   v_pfJetTotalPtCell_;
  std::vector<uint32_t> v_pfJetNCell_;
 
  // ---- jet timing — EBCleaned (spike hits matched to jet, for studies) ----
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
 
  // EB cell branches (good hits)
  timeTree_->Branch("caloCell_eta",      &v_caloCellEta_);
  timeTree_->Branch("caloCell_phi",      &v_caloCellPhi_);
  timeTree_->Branch("caloCell_e",        &v_caloCellEnergy_);
  timeTree_->Branch("caloCell_ecalTime", &v_caloCellEcalTime_);
  timeTree_->Branch("caloCell_flags",    &v_caloCellFlags_);
 
  // EBCleaned cell branches (spike hits)
  timeTree_->Branch("caloCellCleaned_eta",      &v_caloCellCleanedEta_);
  timeTree_->Branch("caloCellCleaned_phi",      &v_caloCellCleanedPhi_);
  timeTree_->Branch("caloCellCleaned_e",        &v_caloCellCleanedEnergy_);
  timeTree_->Branch("caloCellCleaned_ecalTime", &v_caloCellCleanedEcalTime_);
  timeTree_->Branch("caloCellCleaned_flags",    &v_caloCellCleanedFlags_);
 
  // PF jet kinematics
  timeTree_->Branch("pfJet_e",                &v_pfJetE_);
  timeTree_->Branch("pfJet_pt",               &v_pfJetPt_);
  timeTree_->Branch("pfJet_phi",              &v_pfJetPhi_);
  timeTree_->Branch("pfJet_eta",              &v_pfJetEta_);
  timeTree_->Branch("pfJet_chargedHadEnergy", &v_pfJetChargedHadEnergy_);
  timeTree_->Branch("pfJet_neutralHadEnergy", &v_pfJetNeutralHadEnergy_);
  timeTree_->Branch("pfJet_chargedEmEnergy",  &v_pfJetChargedEmEnergy_);
  timeTree_->Branch("pfJet_neutralEmEnergy",  &v_pfJetNeutralEmEnergy_);
 
  // Jet timing — EB
  timeTree_->Branch("pfJet_weightedTime",  &v_pfJetWeightedTime_);
  timeTree_->Branch("pfJet_totalPtCell",   &v_pfJetTotalPtCell_);
  timeTree_->Branch("pfJet_nCell",         &v_pfJetNCell_);
 
  // Jet timing — EBCleaned
  timeTree_->Branch("pfJet_weightedTime_EBCleaned",  &v_pfJetWeightedTimeEBCleaned_);
  timeTree_->Branch("pfJet_totalPtCell_EBCleaned",   &v_pfJetTotalPtCellEBCleaned_);
  timeTree_->Branch("pfJet_nCell_EBCleaned",         &v_pfJetNCellEBCleaned_);
 
  // Trigger decisions
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
 
    GlobalPoint pCell = geom.getPosition(hit.detId());
 
    TLorentzVector cellVec;
    cellVec.SetPtEtaPhiM(1, pCell.eta(), pCell.phi(), 0);
    if (cellVec.DeltaR(jetVec) > 0.4) continue;
 
    const double cellPt = hit.energy() * TMath::Sin(pCell.theta());
    weightedTime += hit.time() * cellPt;
    totalPtOut   += cellPt;
    ++nCellOut;
  }
 
  return (totalPtOut > 0) ? weightedTime / totalPtOut : -200.0;
}
 
// ---------------------------------------------------------------------------
// analyze
// ---------------------------------------------------------------------------
void ScoutingTimingAnalyzer::analyze(const edm::Event& iEvent,
                                     const edm::EventSetup& iSetup) {
  using namespace edm;
 
  auto const& geom             = iSetup.getData(esToken_);
  auto const& pfJets           = iEvent.get(pfJetsToken_);
  auto const& ebRecHits        = iEvent.get(ebRecHitsToken_);
  auto const& ebCleanedRecHits = iEvent.get(ebCleanedRecHitsToken_);
  auto const& triggerResults   = iEvent.get(triggerResultsToken_);
  auto const& triggerResultsRerun = iEvent.get(triggerResultsRerunToken_);
 
  // ---- clear per-event containers ------------------------------------------
  v_caloCellEta_.clear();
  v_caloCellPhi_.clear();
  v_caloCellEnergy_.clear();
  v_caloCellEcalTime_.clear();
  v_caloCellFlags_.clear();
 
  v_caloCellCleanedEta_.clear();
  v_caloCellCleanedPhi_.clear();
  v_caloCellCleanedEnergy_.clear();
  v_caloCellCleanedEcalTime_.clear();
  v_caloCellCleanedFlags_.clear();
 
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
 
  // ---- EB cell branches (good hits) ----------------------------------------
  for (auto const& hit : ebRecHits) {
    if (hit.energy() < 0.5) continue;
    GlobalPoint pCell = geom.getPosition(hit.detId());
    v_caloCellEta_.push_back(pCell.eta());
    v_caloCellPhi_.push_back(pCell.phi());
    v_caloCellEnergy_.push_back(hit.energy());
    v_caloCellEcalTime_.push_back(hit.time());
    v_caloCellFlags_.push_back(hit.flags());
  }
 
  // ---- EBCleaned cell branches (spike hits) --------------------------------
  for (auto const& hit : ebCleanedRecHits) {
    if (hit.energy() < 0.5) continue;
    GlobalPoint pCell = geom.getPosition(hit.detId());
    v_caloCellCleanedEta_.push_back(pCell.eta());
    v_caloCellCleanedPhi_.push_back(pCell.phi());
    v_caloCellCleanedEnergy_.push_back(hit.energy());
    v_caloCellCleanedEcalTime_.push_back(hit.time());
    v_caloCellCleanedFlags_.push_back(hit.flags());
  }
 
  // ---- PF jets -------------------------------------------------------------
  for (auto const& jet : pfJets) {
    double   totalPt1 = 0, totalPt3 = 0;
    uint32_t nCell1   = 0, nCell3   = 0;
 
    double t1 = computeWeightedTime(jet, ebRecHits,        geom, totalPt1, nCell1);
    double t3 = computeWeightedTime(jet, ebCleanedRecHits, geom, totalPt3, nCell3);
 
    v_pfJetE_.push_back(jet.chargedHadronEnergy() + jet.neutralHadronEnergy()
                        + jet.electronEnergy()    + jet.photonEnergy());
    v_pfJetPt_.push_back(jet.pt());
    v_pfJetPhi_.push_back(jet.phi());
    v_pfJetEta_.push_back(jet.eta());
    v_pfJetChargedHadEnergy_.push_back(jet.chargedHadronEnergy());
    v_pfJetNeutralHadEnergy_.push_back(jet.neutralHadronEnergy());
    v_pfJetChargedEmEnergy_.push_back(jet.electronEnergy());
    v_pfJetNeutralEmEnergy_.push_back(jet.photonEnergy());
 
    v_pfJetWeightedTime_.push_back(t1);
    v_pfJetTotalPtCell_.push_back(totalPt1);
    v_pfJetNCell_.push_back(nCell1);
 
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
 
DEFINE_FWK_MODULE(ScoutingTimingAnalyzer);