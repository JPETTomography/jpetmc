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
 *  @file DetectorConstructionMessenger.h
 */

#ifndef DETECTORCONSTRUCTIONMESSENGER_H
#define DETECTORCONSTRUCTIONMESSENGER_H 1

#include "../Core/DetectorConstruction.h"

#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithABool.hh>
#include <G4SystemOfUnits.hh>
#include <G4UIdirectory.hh>
#include <G4UImessenger.hh>
#include <globals.hh>

class DetectorConstructionMessenger : public G4UImessenger
{
public:
  DetectorConstructionMessenger();
  explicit DetectorConstructionMessenger(DetectorConstruction* detector);
  ~DetectorConstructionMessenger();
  void SetNewValue(G4UIcommand*, G4String);
  const int kDefaultRunNumber = 5;

private:
  DetectorConstruction* fDetector = nullptr;
  G4UIdirectory* fDirectory = nullptr;
  G4UIcmdWithAnInteger* fLoadTargetForRun = nullptr;
  G4UIcmdWithAnInteger* fLoadIdealGeometry = nullptr;
  G4UIcmdWithoutParameter* fLoadJPetBasicGeometry = nullptr;
  G4UIcmdWithoutParameter* fLoadOnlyScintillators = nullptr;
  //! Null pointer assigned to the modular layer -
  G4UIcmdWithAString* fLoadModularLayer = nullptr;
  G4UIcmdWithoutParameter* fLoadModularLayerOnly = nullptr;
  G4UIcmdWithADoubleAndUnit* fScinHitMergingTime = nullptr;
  G4UIcmdWithAString* fGeometryFileName = nullptr;
  G4UIcmdWithAString* fCreateGeometryType = nullptr;
  G4UIcmdWithADoubleAndUnit* fPressureInChamber = nullptr;
};

#endif /* !DETECTORCONSTRUCTIONMESSENGER_H */
