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
 *  @file PrimaryGeneratorAction.h
 */

#ifndef PRIMARYGENERATORACTION_H
#define PRIMARYGENERATORACTION_H 1

#include "../Info/PrimaryGeneratorActionMessenger.h"
#include "../Core/PrimaryGenerator.h"
#include "../Core/HistoManager.h"
#include "../Core/SourceParams.h"
#include "../Core/BeamParams.h"

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4Event.hh>
#include <globals.hh>

struct IsotopeParams
{
  G4String shape;
};

/**
 * @class PrimaryGeneratorAction
 * @brief heart of simulated physics
 */
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  explicit PrimaryGeneratorAction(HistoManager* histo);
  ~PrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event*);
  void SetSourceTypeInfo(G4String);

  BeamParams* GetBeamParams() { return fBeam; }
  SourceParams* GetIsotopeParams() { return fIsotope; }
  G4String GetSourceTypeInfo() { return fGenerateSourceType; }
  void SetNemaPoint(G4int i) { fNemaPoint = i; }
  G4int GetNemaPoint() { return fNemaPoint; }
  void SetEffectivePositronRadius(G4double);

private:
  G4String fGenerateSourceType;
  G4String fAllowedSourceTypes[5] = {"run", "beam", "isotope", "nema", "cosmics"};
  PrimaryGeneratorActionMessenger* fMessenger = nullptr;
  PrimaryGenerator* fPrimaryGenerator = nullptr;
  HistoManager* fHistoManager = nullptr;
  BeamParams* fBeam = nullptr;
  SourceParams* fIsotope = nullptr;
  G4int fNemaPoint = -1;
  G4double fEffectivePositronRadius = 0.5 * cm;
};

#endif /* !PRIMARYGENERATORACTION_H */
