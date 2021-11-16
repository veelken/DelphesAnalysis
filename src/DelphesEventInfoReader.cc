#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfoReader.h"

#include "tthAnalysis/HiggsToTauTau/interface/BranchAddressInitializer.h" // BranchAddressInitializer, TTree, Form()

DelphesEventInfoReader::DelphesEventInfoReader()
  : DelphesEventInfoReader(nullptr)
{}

DelphesEventInfoReader::DelphesEventInfoReader(DelphesEventInfo * eventInfo)
  : eventInfo_(eventInfo)
  , branchName_run_("run")
  , branchName_luminosityBlock_("luminosityBlock")
  , branchName_event_("event")
  , branchName_scale_("scale")
  , branchName_x1_("x1")
  , branchName_x2_("x2")
  , branchName_id1_("id1")
  , branchName_id2_("id2")
  , branchName_alphaQED_("alphaQED")
  , branchName_alphaQCD_("alphaQCD")
  , branchName_genweight_("genweight")
{}

DelphesEventInfoReader::~DelphesEventInfoReader()
{}

std::vector<std::string>
DelphesEventInfoReader::setBranchAddresses(TTree * tree)
{
  BranchAddressInitializer bai(tree);
  bai.setBranchAddress(eventInfo_ -> run_,             branchName_run_);
  bai.setBranchAddress(eventInfo_ -> luminosityBlock_, branchName_luminosityBlock_);
  bai.setBranchAddress(eventInfo_ -> event_,           branchName_event_);
  bai.setBranchAddress(eventInfo_ -> scale_,           branchName_scale_);
  bai.setBranchAddress(eventInfo_ -> x1_,              branchName_x1_);
  bai.setBranchAddress(eventInfo_ -> x2_,              branchName_x2_);
  bai.setBranchAddress(eventInfo_ -> id1_,             branchName_id1_);
  bai.setBranchAddress(eventInfo_ -> id2_,             branchName_id2_);
  bai.setBranchAddress(eventInfo_ -> alphaQED_,        branchName_alphaQED_);
  bai.setBranchAddress(eventInfo_ -> alphaQCD_,        branchName_alphaQCD_);
  bai.setBranchAddress(eventInfo_ -> genweight_,       branchName_genweight_);
  return bai.getBoundBranchNames();
}
