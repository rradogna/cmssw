#ifndef GEMBaseValidation_H
#define GEMBaseValidation_H

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

class GEMBaseValidation
{
public:
  GEMBaseValidation(DQMStore* dbe,
                         const edm::InputTag & inputTag);
  virtual ~GEMBaseValidation();
  void setGeometry(const GEMGeometry* geom) { theGEMGeometry = geom; }


 protected:

  DQMStore* dbe_;
  edm::InputTag theInputTag;
  const GEMGeometry* theGEMGeometry;

};

#endif
