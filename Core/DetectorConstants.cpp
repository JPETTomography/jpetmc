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
 *  @file DetectorConstants.cpp
 */

#include "DetectorConstants.h"

const G4double DetectorConstants::world_size[3] = {1.0 * m, 1.0 * m, 1.0 * m};
const G4double DetectorConstants::scinDim[3] = {1.9 * cm, 0.7 * cm, 50.0 * cm};
const G4double DetectorConstants::radius[layers] = {42.5 * cm, 46.75 * cm, 57.5 * cm};
const int DetectorConstants::nSegments[layers] = {48, 48, 96};
const int DetectorConstants::nSegmentsExtraLayers[extraLayers] = {96, 96};
const G4double DetectorConstants::radiusExtraLayers[extraLayers] =  {509 * mm, 533 * mm};
const G4double DetectorConstants::scinDim_inModule[3] = {2.5 * cm, 0.6 * cm, 50.0 * cm};
const G4double DetectorConstants::wrappingThickness = 2 * 25.4 * pow(10, -6) * m;
const G4double DetectorConstants::wrappingShift = 1 * pow(10, -5) * m;
G4double DetectorConstants::fScinHitMergingTime = 5.0 * ns;
G4double DetectorConstants::fChamberCenter[3] = {0.0 * cm, 0.0 * cm, 0.0 * cm};

void DetectorConstants::SetChamberCenter(const G4ThreeVector xyz)
{
  fChamberCenter[0] = xyz.x();
  fChamberCenter[1] = xyz.y();
  fChamberCenter[2] = xyz.z();
}

const G4ThreeVector DetectorConstants::GetChamberCenter()
{
  return G4ThreeVector(fChamberCenter[0], fChamberCenter[1], fChamberCenter[2]);
}

G4double DetectorConstants::GetMergingTimeValueForScin()
{
  return fScinHitMergingTime;
}

void DetectorConstants::SetMergingTimeValueForScin(G4double x)
{
  if (x <= 0.0) {
    G4Exception("DetectorConstants", "DConst01", JustWarning,
                "Time can not be negative! Nothing happend");
  } else {
    fScinHitMergingTime = x;
  }
}
