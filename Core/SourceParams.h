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
 *  @file SourceParams.h
 */

#ifndef SOURCEPARAMS_H
#define SOURCEPARAMS_H 1

#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

class SourceParams
{
public:
  SourceParams();
  ~SourceParams();

  void SetGammasNumber(G4int i) { howManyGammas = i; }
  void SetShape(G4String& newShape);
  void SetShapeCenterPosition(G4double x, G4double y, G4double z) {
    shapeCenterPosition.set(x, y, z);
  }
  void SetShapeDim(G4int i, G4double x) { shapeDim[i] = x; }
  G4int GetGammasNumber() const { return howManyGammas; }
  G4String GetShape() const { return shape; }
  G4double GetShapeDim(G4int i) const { return shapeDim[i]; }
  G4ThreeVector GetShapeCenterPosition() const { return shapeCenterPosition; }

private:
  G4int howManyGammas;
  G4String shape;
  G4String allowedShapes[2] = {"cylinder", "cuboid"};
  G4ThreeVector shapeCenterPosition;
  //! array of dimensions, for cylinder [0] - radius, [1] - z
  G4double shapeDim[10] = {0 * cm};
};

#endif /* !SOURCEPARAMS_H */
