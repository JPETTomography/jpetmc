/**
 *  @copyright Copyright 2020 The J-PET Monte Carlo Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file HistoManager.h
 */

#ifndef HISTOMANAGER_H
#define HISTOMANAGER_H 1

#include "../Objects/Framework/JPetGeantEventInformation.h"
#include "../Objects/Framework/JPetGeantEventPack.h"
#include "../Objects/Framework/JPetGeantDecayTree.h"
#include "../Objects/Framework/JPetGeantScinHits.h"
#include "../Objects/Geant4/DetectorHit.h"
#include "../Info/VtxInformation.h"
#include "../Info/EventMessenger.h"
#include <G4PrimaryParticle.hh>
#include <G4Event.hh>
#include <globals.hh>
#include <THashTable.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

class TFile;
class TTree;

struct doubleCheck
{
  bool isChanged = false;
  double value = 0.;
  doubleCheck(){}
  doubleCheck(double newValue) {value=newValue; isChanged=true;}
};

/**
 * @class HistoManager
 * @brief class reach for informations stored during sumulations
 * and saves them into ROOT tree structures
 */
class HistoManager
{
public:
  HistoManager();
  ~HistoManager();

  void Book(); //! call once; book (create) all trees and histograms
  void Save(); //! call once; save all trees and histograms
  void SaveEvtPack() {fTree->Fill();};
  void Clear() {fEventPack->Clear();};
  void AddGenInfo(VtxInformation* info);
  void AddGenInfoParticles(G4PrimaryParticle* particle);
  void AddNewHit(DetectorHit*);
  void SetEventNumber(int x) {fEventPack->SetEventNumber(x);};
  void SetHistogramCreation(bool tf) {fMakeControlHisto = tf;};
  bool GetMakeControlHisto() {return fMakeControlHisto;};
  void FillHistoGenInfo(const G4Event* anEvent);
  void FillCosmicInfo(G4double theta, G4ThreeVector init, G4ThreeVector vtx);
  const JPetGeantEventInformation* GetGeantInfo() {return fGeantInfo;}
  void createHistogramWithAxes(
    TObject* object, TString xAxisName="Default X axis title [unit]",
    TString yAxisName="Default Y axis title [unit]",
    TString zAxisName="Default Z axis title [unit]"
  );
  void fillHistogram(
    const char* name, double xValue, doubleCheck yValue=doubleCheck(),
    doubleCheck zValue=doubleCheck()
  );
  void writeError(const char* nameOfHistogram, const char* messageEnd );

  template <typename T>
  T* getObject(const char* name)
  {
    TObject* tmp = fStats.FindObject(name);
    if (!tmp) {return nullptr;}
    return dynamic_cast<T*>(tmp);
  }

private:
<<<<<<< HEAD
  bool fBookStatus = false;
  bool fMakeControlHisto = false;
=======
  bool fbookStatus = false;
  bool fMakeControlHisto;
>>>>>>> a99f6b0a00e9625ff61a1f07e34d1737590fbbb6
  TFile* fRootFile = nullptr;
  TTree* fTree = nullptr;
  TBranch* fBranchTrk = nullptr;
  TBranch* fBranchScin = nullptr;
  TBranch* fBranchEventPack = nullptr;

  JPetGeantEventPack* fEventPack = nullptr;
  JPetGeantEventInformation* fGeantInfo = nullptr;
  EventMessenger* fEvtMessenger = EventMessenger::GetEventMessenger();

  void BookHistograms();

protected:
  THashTable fStats;
  std::set<std::string> fErrorCounts;
};

#endif
