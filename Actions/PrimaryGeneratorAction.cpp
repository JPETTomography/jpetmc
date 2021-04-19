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
 *  @file PrimaryGeneratorAction.cpp
 */

#include "../Core/DetectorConstruction.h"
#include "../Core/DetectorConstants.h"
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
  fPositionsNemaPoints.clear();
  fWeightedNemaPoints.clear();
  fLifetimesNemaPoints.clear();
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
    } else if (nRun == 6 || nRun == 7 || nRun == 12) {
      fPrimaryGenerator->GenerateEvtLargeChamber(event);
    } else {
      G4Exception(
        "PrimaryGeneratorAction", "PG05", FatalException,
        "Called run with non-exisitng geometry"
      );
    }
  } else if (GetSourceTypeInfo() == ("beam")) {
    fPrimaryGenerator->GenerateBeam(fBeam, event);
  } else if (GetSourceTypeInfo() == ("isotope")) {
    fPrimaryGenerator->GenerateIsotope(fIsotope, event);
  } else if (GetSourceTypeInfo() == ("nema")) {
    fPrimaryGenerator->GenerateNema(GetNemaPoint(), event, fPositionsNemaPoints, fWeightedNemaPoints, fLifetimesNemaPoints);
  } else if (GetSourceTypeInfo() == ("nema-mixed")) {
    fPrimaryGenerator->GenerateNema(-1, event, fPositionsNemaPoints, fWeightedNemaPoints, fLifetimesNemaPoints);
  } else {
    G4Exception(
      "PrimaryGeneratorAction", "PG05", FatalException,
      "Called run with non-exisitng geometry"
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
      if (newSourceType == "nema-mixed" || newSourceType == "nema") {
        GenerateDefaultNemaPositions();
      }
      if (newSourceType == "nema-mixed") {
        for (int i=1; i<7; i++) {
          fWeightedNemaPoints.push_back(i);
        }
      }
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

void PrimaryGeneratorAction::GenerateDefaultNemaPositions()
{
  for (unsigned nemaPoint=1; nemaPoint<7; nemaPoint++) {
    G4double x_creation = 0.0 * cm;
    G4double y_creation = 0.0 * cm;
    G4double z_creation = 0.0 * cm;
    
    if (nemaPoint > 3) {
      z_creation = z_creation - DetectorConstants::scinDim[2] * 3 / 8;
    }

    if (nemaPoint == 1 || nemaPoint == 4) {
      y_creation = y_creation + 1.0 * cm;
    } else if (nemaPoint == 2 || nemaPoint == 5) {
      y_creation = y_creation + 10.0 * cm;
    } else if (nemaPoint == 3 || nemaPoint == 6) {
      y_creation = y_creation + 20.0 * cm;
    }
    
    fPositionsNemaPoints.push_back(G4ThreeVector(x_creation, y_creation, z_creation));
  }
}

void PrimaryGeneratorAction::SetNemaPointPosition(int nemaPoint, G4ThreeVector position)
{
  if ((nemaPoint < 1 || nemaPoint - 1 > (int)fPositionsNemaPoints.size())) {
    G4Exception(
      "PrimaryGeneratorAction", "PG06", JustWarning,
      "Nema point for which you want to set position is less than 1 or greater than the size of the pointsVector more than one. Not setting it."
    );
  } else if (nemaPoint - 1 == (int)fPositionsNemaPoints.size()) {
    fPositionsNemaPoints.push_back(position);
  } else {
    fPositionsNemaPoints.at(nemaPoint - 1) = position;
  }
}

void PrimaryGeneratorAction::SetNemaPositionWeight(int nemaPoint, int weight)
{
  unsigned removedElements = 0;
  if (weight == 0) {
    for (unsigned i=0; i<fWeightedNemaPoints.size(); i++) {
      if (fWeightedNemaPoints.at(i-removedElements) == nemaPoint) {
        fWeightedNemaPoints.erase(fWeightedNemaPoints.begin()+i-removedElements);
        removedElements++;
      }
    }
  } else {
    int count = 0;
    for (unsigned i=0; i<fWeightedNemaPoints.size(); i++) {
      if (fWeightedNemaPoints.at(i) == nemaPoint)
        count++;
    }
    if (weight >= count) {
      for (int i=0; i<weight-count; i++)
        fWeightedNemaPoints.push_back(nemaPoint);
    } else {
      int nmbOfPointToErase = count - weight;
      for (unsigned i=0; i<fWeightedNemaPoints.size(); i++) {
        if (fWeightedNemaPoints.at(i-removedElements) == nemaPoint) {
          if ((int)removedElements < nmbOfPointToErase) {
            fWeightedNemaPoints.erase(fWeightedNemaPoints.begin()+i-removedElements);
            removedElements++;
          }
        }
      }
    }
  }
}

void PrimaryGeneratorAction::SetNemaPointLifetime(int nemaPoint, double lifetime)
{
  if ((nemaPoint < 1 || nemaPoint - 1 > (int)fLifetimesNemaPoints.size())) {
    G4Exception(
      "PrimaryGeneratorAction", "PG07", JustWarning,
      "Nema point for which you want to set lifetime is less than 1 or greater than the size of the pointsVector more than one. Not setting it."
    );
  } else if (nemaPoint - 1 == (int)fLifetimesNemaPoints.size()) {
    fLifetimesNemaPoints.push_back(lifetime);
  } else {
    fLifetimesNemaPoints.at(nemaPoint - 1) = lifetime;
  }
}
