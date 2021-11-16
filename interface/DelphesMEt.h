#ifndef hhAnalysis_DelphesAnalysis_DelphesMEt_h
#define hhAnalysis_DelphesAnalysis_DelphesMEt_h

#include <Rtypes.h> // Float_t

#include <iostream> // std::ostream

class DelphesMEt
{
 public:
  DelphesMEt();
  DelphesMEt(Float_t pt, Float_t phi);

  DelphesMEt &
  operator=(const DelphesMEt & other);

  /**
   * @brief Funtions to access data-members
   * @return Values of data-members
   */  
  Float_t pt() const;
  Float_t phi() const;
  
  Float_t px() const;
  Float_t py() const;

  friend class DelphesMEtReader;

 protected:
  Float_t pt_;
  Float_t phi_;
};

std::ostream &
operator<<(std::ostream& stream,
           const DelphesMEt & met);

#endif // hhAnalysis_DelphesAnalysis_DelphesMEt_h
