#ifndef hhAnalysis_DelphesAnalysis_DelphesEventInfoReader_h
#define hhAnalysis_DelphesAnalysis_DelphesEventInfoReader_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfo.h" // DelphesEventInfo

#include "tthAnalysis/HiggsToTauTau/interface/ReaderBase.h"        // ReaderBase

// forward declarations
class TTree;

class DelphesEventInfoReader : public ReaderBase
{
 public:
  explicit DelphesEventInfoReader();
  explicit DelphesEventInfoReader(DelphesEventInfo * eventInfo);
  ~DelphesEventInfoReader() override;

  std::vector<std::string>
  setBranchAddresses(TTree * tree) override;

  void
  setEventInfo(DelphesEventInfo * eventInfo);

 protected:
  DelphesEventInfo * eventInfo_;

  std::string branchName_run_;
  std::string branchName_luminosityBlock_;
  std::string branchName_event_;
  std::string branchName_scale_;
  std::string branchName_x1_;
  std::string branchName_x2_;
  std::string branchName_id1_;
  std::string branchName_id2_;
  std::string branchName_alphaQED_;
  std::string branchName_alphaQCD_;
  std::string branchName_genweight_;
};

#endif // hhAnalysis_DelphesAnalysis_DelphesEventInfoReader_h
