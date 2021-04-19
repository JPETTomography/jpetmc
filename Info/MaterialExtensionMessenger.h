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
 *  @file MaterialExtensionMessenger.h
 */

#ifndef MATERIALEXTENSIONMESSENGER_H
#define MATERIALEXTENSIONMESSENGER_H 1

#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithABool.hh>
#include <G4SystemOfUnits.hh>
#include <G4UIdirectory.hh>
#include <G4UImessenger.hh>
#include <globals.hh>

class MaterialExtension;

class MaterialExtensionMessenger : public G4UImessenger
{
public:
  static MaterialExtensionMessenger* GetMaterialExtensionMessenger();
  void SetNewValue(G4UIcommand*, G4String);
  bool GenerateOnly3g() { return fIsOnly3g; }

private:
  static MaterialExtensionMessenger* fInstance;
  MaterialExtensionMessenger();
  ~MaterialExtensionMessenger();
  G4UIdirectory* fDirectory = nullptr;
  MaterialExtension* fMaterial = nullptr;
  G4UIcmdWithoutParameter* f3GammaOnly = nullptr;
  G4UIcmdWithoutParameter* f3GammapPs = nullptr;
  G4UIcmdWithoutParameter* f2GammaOnly = nullptr;
  G4UIcmdWithoutParameter* fPickOffOnly = nullptr;
  G4UIcmdWithAString* fAddoPsComponent = nullptr;
  G4UIcmdWithAString* fSetpPsComponent = nullptr;
  G4UIcmdWithAString* fAddDirectComponent = nullptr;
  G4UIcmdWithADouble* fSetpPs3GFraction = nullptr;
  G4UIcmdWithAString* fReloadMaterials = nullptr;
  bool fIsOnly3g = true;
};

#endif /* !MATERIALEXTENSIONMESSENGER_H */
