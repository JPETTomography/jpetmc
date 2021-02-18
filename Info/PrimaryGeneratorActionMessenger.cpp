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
 *  @file PrimaryGeneratorActionMessenger.cpp
 */

#include "../Actions/PrimaryGeneratorAction.h"
#include "PrimaryGeneratorActionMessenger.h"
#include "../Core/DetectorConstruction.h"
#include "../Core/DetectorConstants.h"

PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger() {}

PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger(PrimaryGeneratorAction* primGeneratorAction) : fPrimGen(primGeneratorAction)
{
  fDirectory = new G4UIdirectory("/jpetmc/source/");
  fDirectory->SetGuidance("Commands for controling the gamma quanta source (beam/target/nema) and its parameters");

  fDirectoryRun = new G4UIdirectory("/jpetmc/run/");
  fDirectoryRun->SetGuidance("Commands for controling  parameters");

  fSourceType = new G4UIcmdWithAString("/jpetmc/source/setType", this);
  fSourceType->SetCandidates("beam isotope nema");
  fSourceType->SetDefaultValue("beam");

  fGammaBeamSetEnergy = new G4UIcmdWithADoubleAndUnit("/jpetmc/source/gammaBeam/setEnergy", this);
  fGammaBeamSetEnergy->SetGuidance("Set energy (value and unit) for beam of gamma quanta");
  fGammaBeamSetEnergy->SetDefaultValue(511.0);
  fGammaBeamSetEnergy->SetDefaultUnit("keV");
  fGammaBeamSetEnergy->SetUnitCandidates("keV");

  fGammaBeamSetPosition = new G4UIcmdWith3VectorAndUnit("/jpetmc/source/gammaBeam/setPosition", this);
  fGammaBeamSetPosition->SetGuidance("Set vertex position of the gamma quanta beam");
  fGammaBeamSetPosition->SetDefaultValue(G4ThreeVector(0, 0, 0));
  fGammaBeamSetPosition->SetDefaultUnit("cm");
  fGammaBeamSetPosition->SetUnitCandidates("cm");
  fGammaBeamSetPosition->SetParameterName("Xvalue", "Yvalue", "Zvalue", false);

  fGammaBeamSetMomentum = new G4UIcmdWith3VectorAndUnit("/jpetmc/source/gammaBeam/setMomentum", this);
  fGammaBeamSetMomentum->SetGuidance("Set momentum of the gamma quanta beam");
  fGammaBeamSetMomentum->SetDefaultValue(G4ThreeVector(1, 0, 0));
  fGammaBeamSetMomentum->SetDefaultUnit("keV");
  fGammaBeamSetMomentum->SetUnitCandidates("keV");
  fGammaBeamSetMomentum->SetParameterName("Xvalue", "Yvalue", "Zvalue", false);

  fIsotopeSetShape = new G4UIcmdWithAString("/jpetmc/source/isotope/setShape", this);
  fIsotopeSetShape->SetCandidates("cylinder");

  fIsotopeSetGenGammas = new G4UIcmdWithAnInteger("/jpetmc/source/isotope/setNGamma", this);
  fIsotopeSetGenGammas->SetGuidance("Give number of gamma quanta to generate 1 / 2 / 3");

  fIsotopeSetShapeDimCylinderRadius = new G4UIcmdWithADoubleAndUnit("/jpetmc/source/isotope/setShape/cylinderRadius", this);
  fIsotopeSetShapeDimCylinderRadius->SetGuidance("For cylindrical shape - set radius");
  fIsotopeSetShapeDimCylinderRadius->SetDefaultValue(10);
  fIsotopeSetShapeDimCylinderRadius->SetDefaultUnit("cm");
  fIsotopeSetShapeDimCylinderRadius->SetUnitCandidates("cm");

  fIsotopeSetShapeDimCylinderZ = new G4UIcmdWithADoubleAndUnit("/jpetmc/source/isotope/setShape/cylinderZ", this);
  fIsotopeSetShapeDimCylinderZ->SetGuidance("For cylindrical shape - set z (half length)");
  fIsotopeSetShapeDimCylinderZ->SetDefaultValue(10);
  fIsotopeSetShapeDimCylinderZ->SetDefaultUnit("cm");
  fIsotopeSetShapeDimCylinderZ->SetUnitCandidates("cm");

  fIsotopeSetCenter = new G4UIcmdWith3VectorAndUnit("/jpetmc/source/isotope/setPosition", this);
  fIsotopeSetCenter->SetGuidance("Set position of the source");
  fIsotopeSetCenter->SetDefaultValue(G4ThreeVector(0, 0, 0));
  fIsotopeSetCenter->SetParameterName("Xvalue", "Yvalue", "Zvalue", false);
  fIsotopeSetCenter->SetDefaultUnit("cm");
  fIsotopeSetCenter->SetUnitCandidates("cm");

  fNemaPosition = new G4UIcmdWithAnInteger("/jpetmc/source/nema", this);
  fNemaPosition->SetGuidance("Give nema point number to simulate (1-6)");
  fNemaPosition->SetDefaultValue(1);

  fSetChamberCenter = new G4UIcmdWith3VectorAndUnit("/jpetmc/run/setChamberCenter", this);
  fSetChamberCenter->SetGuidance("Set position of the annihilation chamber");
  fSetChamberCenter->SetDefaultValue(G4ThreeVector(0, 0, 0));
  fSetChamberCenter->SetDefaultUnit("cm");
  fSetChamberCenter->SetUnitCandidates("cm");
  fSetChamberCenter->SetParameterName("Xvalue", "Yvalue", "Zvalue", false);

  fSetChamberEffectivePositronRadius = new G4UIcmdWithADoubleAndUnit("/jpetmc/run/setEffectivePositronRange", this);
  fSetChamberEffectivePositronRadius->SetGuidance("Set effective positron radius (RUN5)");
  fSetChamberEffectivePositronRadius->SetDefaultValue(0.5);
  fSetChamberEffectivePositronRadius->SetDefaultUnit("cm");
  fSetChamberEffectivePositronRadius->SetUnitCandidates("cm");
  
  fCosmicOnly = new G4UIcmdWithoutParameter("/jpetmc/source/cosmicOnly",this);
  fCosmicOnly->SetGuidance("Option to generate cosmic rays only");
  fCosmicOnly->SetGuidance("Generate only cosmics");
  
  fCosmicGenShape = new G4UIcmdWithAString("/jpetmc/source/cosmicGenShape", this);
  fCosmicGenShape->SetGuidance("Cosmics can be generated on a cylinder or in the cuboid");
  fCosmicGenShape->SetCandidates("cylinder cuboid"); 
}

PrimaryGeneratorActionMessenger::~PrimaryGeneratorActionMessenger()
{
  delete fIsotopeSetShape;
  delete fIsotopeSetShapeDimCylinderRadius;
  delete fIsotopeSetShapeDimCylinderZ;
  delete fSourceType;
  delete fGammaBeamSetEnergy;
  delete fGammaBeamSetPosition;
  delete fGammaBeamSetMomentum;
  delete fIsotopeSetCenter;
  delete fNemaPosition;
  delete fSetChamberCenter;
  delete fSetChamberEffectivePositronRadius;
  delete fCosmicOnly;
  delete fCosmicGenShape;
}

void PrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fSourceType) {
    fPrimGen->SetSourceTypeInfo(newValue);
  } else if (command == fGammaBeamSetEnergy) {
    ChangeToBeam();
    fPrimGen->GetBeamParams()->SetEnergy(fGammaBeamSetEnergy->GetNewDoubleValue(newValue));
  } else if (command == fGammaBeamSetPosition) {
    ChangeToBeam();
    fPrimGen->GetBeamParams()->SetVtxPosition(fGammaBeamSetPosition->GetNew3VectorValue(newValue));
  } else if (command == fGammaBeamSetMomentum) {
    ChangeToBeam();
    fPrimGen->GetBeamParams()->SetMomentum(fGammaBeamSetMomentum->GetNew3VectorValue(newValue));
  } else if (command == fIsotopeSetShape) {
    fPrimGen->GetIsotopeParams()->SetShape(newValue);
  } else if (command == fIsotopeSetGenGammas) {
    fPrimGen->GetIsotopeParams()->SetGammasNumber(fIsotopeSetGenGammas->GetNewIntValue(newValue));
  } else if (command == fIsotopeSetShapeDimCylinderRadius) {
    ChangeToIsotope();
    fPrimGen->GetIsotopeParams()->SetShapeDim(0, fIsotopeSetShapeDimCylinderRadius->GetNewDoubleValue(newValue));
  } else if (command == fIsotopeSetShapeDimCylinderZ) {
    ChangeToIsotope();
    fPrimGen->GetIsotopeParams()->SetShapeDim(1, fIsotopeSetShapeDimCylinderRadius->GetNewDoubleValue(newValue));
  } else if (command == fIsotopeSetCenter) {
    ChangeToIsotope();
    G4ThreeVector loc = fIsotopeSetCenter->GetNew3VectorValue(newValue);
    fPrimGen->GetIsotopeParams()->SetShapeCenterPosition(loc.x(), loc.y(), loc.z());
  } else if (command == fNemaPosition) {
    fPrimGen->SetSourceTypeInfo("nema");
    fPrimGen->SetNemaPoint(fNemaPosition->GetNewIntValue(newValue));
  } else if (command == fSetChamberCenter) {
    if (!CheckIfRun()) { ChangeToRun(); }
    DetectorConstants::SetChamberCenter(fSetChamberCenter->GetNew3VectorValue(newValue));
  } else if (command == fSetChamberEffectivePositronRadius) {
    if (!CheckIfRun()) { ChangeToRun(); }
    fPrimGen->SetEffectivePositronRadius(
      fSetChamberEffectivePositronRadius->GetNewDoubleValue(newValue)
    );
  } else if (command == fCosmicOnly) {
    fPrimGen->SetSourceTypeInfo("cosmics");
  } else if (command == fCosmicGenShape) {
    fPrimGen->GetIsotopeParams()->SetShape(newValue);
  }
}

void PrimaryGeneratorActionMessenger::ChangeToIsotope()
{
  if (fPrimGen->GetSourceTypeInfo() != "isotope") {
    G4Exception(
      "PrimaryGeneratorActionMessenger", "PGM02",
      JustWarning, "Changed sourceType to isotope"
    );
    fPrimGen->SetSourceTypeInfo("isotope");
  }
}

bool PrimaryGeneratorActionMessenger::CheckIfRun()
{
  G4int nRun = DetectorConstruction::GetInstance()->GetRunNumber();
  if (nRun != 0) {
    fPrimGen->SetSourceTypeInfo("run");
    return true;
  }
  return false;
}

void PrimaryGeneratorActionMessenger::ChangeToRun()
{
  if (fPrimGen->GetSourceTypeInfo() != "run") {
    G4Exception(
      "PrimaryGeneratorActionMessenger", "PGM01",
      JustWarning, "Changed run number to RUN5"
    );
    DetectorConstruction::GetInstance()->LoadGeometryForRun(5);
    fPrimGen->SetSourceTypeInfo("run");
    DetectorConstruction::GetInstance()->UpdateGeometry();
  }
}

void PrimaryGeneratorActionMessenger::ChangeToBeam()
{
  //! check if we really changing parameters corresponding to the beam
  if (fPrimGen->GetSourceTypeInfo() != "beam") {
    G4Exception(
      "PrimaryGeneratorActionMessenger", "PGM01",
      JustWarning, "Changed sourceType to beam"
    );
    fPrimGen->SetSourceTypeInfo("beam");
  }
}
