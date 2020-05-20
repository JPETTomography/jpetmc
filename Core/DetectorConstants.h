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
 *  @file DetectorConstants.h
 */

#ifndef DETECTORCONSTANTS_H
#define DETECTORCONSTANTS_H 1

#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

class DetectorConstants
{
public:
  //! world size
  static const G4double world_size[3];
  //! Number of simulated layers in detector
  static const G4int layers = 3;
  //! Dimensions of simulated strip
  static const G4double scinDim[3];
  //! Radius of the layer up to center of the strip
  static const G4double radius[layers];
  //! number of segments (slots) in each layer
  static const int nSegments[layers];

  static const int extraLayers = 2;
  static const int nSegmentsExtraLayers[extraLayers];
  static const G4double radiusExtraLayers[extraLayers];

  //! Total width of used wrapping
  static const G4double wrappingThickness;
  //! Free space between wrapping and scinitlator
  static const G4double wrappingShift;

  //! Parameters for modular layer (known as 4th layer)
  static const G4double scinDim_inModule[3];
  //! Dimensions of simulated strip in modular layer
  static const int modulesInModularLayer = 24;

  static G4double GetMergingTimeValueForScin();
  static void SetMergingTimeValueForScin(G4double);
  static void SetChamberCenter(const G4ThreeVector);
  static const G4ThreeVector GetChamberCenter();

private:
  //! time in [ns]
  static G4double fScinHitMergingTime;
  static G4double fChamberCenter[3];

};

#endif
