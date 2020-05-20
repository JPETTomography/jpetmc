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
 *  @file VtxInformation.h
 */

#ifndef VTX_INFORMATION_H
#define VTX_INFORMATION_H 1

#include <G4VUserPrimaryVertexInformation.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

class VtxInformation : public G4VUserPrimaryVertexInformation
{
public:
  VtxInformation();
  virtual ~VtxInformation();
  void Clear();
  virtual void Print() const;

  void SetThreeGammaGen(G4bool tf) { fThreeGammaGen = tf; };
  void SetTwoGammaGen(G4bool tf) { fTwoGammaGen = tf; };
  void SetPromptGammaGen(G4bool tf) { fPromptGammaGen = tf; };
  void SetCosmicGen(G4bool isCosmic) { fCosmicGen = isCosmic; };
  void SetRunNr(G4int x) { fnRun = x; };
  void SetVtxPosition(G4double x, G4double y, G4double z);
  void SetVtxPosition(G4ThreeVector position);
  void SetLifetime(G4double x) { fLifetime = x; };
  G4bool GetThreeGammaGen() const { return fThreeGammaGen; };
  G4bool GetTwoGammaGen() const { return fTwoGammaGen; };
  G4bool GetPromptGammaGen() const { return fPromptGammaGen; };
  G4bool GetCosmicGammaGen() const { return fCosmicGen; };
  G4int GetRunNr() const { return fnRun; };
  G4double GetVtxPositionX() const { return fVtxPosition.x(); };
  G4double GetVtxPositionY() const { return fVtxPosition.y(); };
  G4double GetVtxPositionZ() const { return fVtxPosition.z(); };
  G4double GetLifetime() const { return fLifetime; };

private:
  G4ThreeVector fVtxPosition;
  G4bool fTwoGammaGen;
  G4bool fThreeGammaGen;
  G4bool fPromptGammaGen;
  G4bool fCosmicGen;
  G4int fnRun;
  G4double fLifetime;

};

#endif
