/**
 *  @copyright Copyright 2021 The J-PET Monte Carlo Authors. All rights reserved.
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
 *  @file PrimaryGeneratorAction.cpp
 */

#include "../Core/DetectorConstruction.h"
#include "PrimaryGeneratorAction.h"

#include <G4PrimaryVertex.hh>

PrimaryGeneratorAction::PrimaryGeneratorAction() {}

PrimaryGeneratorAction::PrimaryGeneratorAction(HistoManager* histo) :
G4VUserPrimaryGeneratorAction(), fPrimaryGenerator(0), fHistoManager(histo)
{
  fPrimaryGenerator = new PrimaryGenerator();
  fBeam = new BeamParams();
  fIsotope = new SourceParams();
  fMessenger = new PrimaryGeneratorActionMessenger(this);
  DetectorConstruction::GetInstance()->SetHistoManager(histo);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fIsotope;
  delete fBeam;
  delete fPrimaryGenerator;
  delete fMessenger;
}

void PrimaryGeneratorAction::SetEffectivePositronRadius(G4double radius)
{
  if (radius <= 0) {
    G4Exception(
      "PrimaryGeneratorAction", "PG0Radius", JustWarning,
      "EffectivePositronRadius can not be negative!! Value is not changed"
    );
  } else {
    fEffectivePositronRadius = radius;
  }
}

// cppcheck-suppress unusedFunction
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  //! if setup for dedicated run is set then ignore its modifications made by user
  G4int nRun = DetectorConstruction::GetInstance()->GetRunNumber();
  if (nRun != 0) {
    if (GetSourceTypeInfo() != "run") {
      SetSourceTypeInfo("run");
    }
  }

  if (GetSourceTypeInfo() == ("run")) {
    if (nRun == 3) {
      fPrimaryGenerator->GenerateEvtLargeChamber(event);
    } else if (nRun == 5) {
      fPrimaryGenerator->GenerateEvtSmallChamber(event, fEffectivePositronRadius);
    } else if (nRun == 6 || nRun == 7) {
      fPrimaryGenerator->GenerateEvtLargeChamber(event);
    } else {
      G4Exception(
        "PrimaryGeneratorAction", "PG05", FatalException,
        "Called run with non-exisitng geometry (number)"
      );
    }
  } else if (GetSourceTypeInfo() == ("beam")) {
    fPrimaryGenerator->GenerateBeam(fBeam, event);
  } else if (GetSourceTypeInfo() == ("isotope")) {
    fPrimaryGenerator->GenerateIsotope(fIsotope, event);
  } else if (GetSourceTypeInfo() == ("nema")) {
    fPrimaryGenerator->GenerateNema(GetNemaPoint(), event);
  } else if (GetSourceTypeInfo() == ("cosmics")) {
    fPrimaryGenerator->GenerateCosmicVertex(fIsotope, event, fHistoManager); 
  } else {
    G4Exception(
      "PrimaryGeneratorAction", "PG05", FatalException,
      "Called non-existing source"
    );
  }
  fHistoManager->SetDecayChannel(fPrimaryGenerator->GetDecayChannel());
}

void PrimaryGeneratorAction::SetSourceTypeInfo(G4String newSourceType)
{
  if (std::find(std::begin(fAllowedSourceTypes), std::end(fAllowedSourceTypes), newSourceType)
                != std::end(fAllowedSourceTypes)) {
    //! setup found
    G4int nRun = DetectorConstruction::GetInstance()->GetRunNumber();
    if ((nRun == 0) && (newSourceType != "run")) {
      fGenerateSourceType = newSourceType;
    } else if (nRun > 0) {
      fGenerateSourceType = "run";
    } else {
      G4Exception(
        "PrimaryGeneratorAction", "PG01", JustWarning,
        "Chosen detector geometry corresponds to run number and it can not be changed"
      );
    }
  } else {
    G4Exception(
      "PrimaryGeneratorAction", "PG02", JustWarning,
      "Please pick from avaliable setups: beam/isotope"
    );
  }
}
