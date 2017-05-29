// -*- C++ -*-
//
// Package:   RunSelect
// Class:     RunSelect
//

#ifndef RunSelect_H
#define RunSelect_H

// system include files
#include <memory>
#include <vector>
#include <map>
#include <set>

// user include files
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


//
// class declaration
//


class RunSelect : public edm::EDFilter {
public:
  explicit RunSelect( const edm::ParameterSet & );
  ~RunSelect();
  
private:
  virtual bool filter ( edm::Event &, const edm::EventSetup & );
  //bool applyfilter;
  //bool requireHV, requireHVbyLS,
  //bool requireGoodLS, requirePIX, requireTRK, requireECAL,requireES, requireHCAL;
  //bool pvtRunLSectionSelection;// requirePhysicsDeclared;
  bool requireCollidingBX;
  bool requireNoLumiScan, requireLumiScan;
  bool requireNoTimeScan;
  bool debug;
  //  int E;
};

#endif
