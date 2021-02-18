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
 *  @file JPetMC.cpp
 */

#include "MaterialExtension.h"
#include "MaterialParameters.h"

MaterialExtension::MaterialExtension(
  MaterialParameters::MaterialID materialID,
  const G4String& name, const G4Material* baseMaterial) : G4Material(
    name, baseMaterial->GetDensity(), baseMaterial, baseMaterial->GetState(),
    baseMaterial->GetTemperature(), baseMaterial->GetPressure()
) {
  fTarget = false;
  //! depends on the material typ
  fMaterialParameters = new MaterialParameters();
  fMaterialParameters->SetMaterialByName(materialID);
  FillIntensities();
}

MaterialExtension::~MaterialExtension() {}

void MaterialExtension::AddoPsComponent(G4double lifetime, G4double probability)
{
  //! lifetime given in ns
  if (lifetime < 0. || lifetime > MaterialParameters::foPsTauVaccum) {
    G4Exception(
      "MaterialExtension", "ME01", JustWarning,
      "given oPs lifetime not in the range [0, 142] ns"
    );
  } else if (probability < 0. || probability > 100.) {
    G4Exception(
      "MaterialExtension", "ME01", JustWarning,
      "given oPs probability not in the range [0, 100] %"
    );
  } else {
    fMaterialParameters->AddoPsComponent(lifetime, probability);
  }
}

void MaterialExtension::AddDirectComponent(G4double lifetime, G4double probability)
{
  //! lifetime given in ns
  if (lifetime < 0. || lifetime > MaterialParameters::fDirectTauMax) {
    G4Exception(
      "MaterialExtension", "ME01", JustWarning,
      "given direct lifetime not in the range [0, 0.6] ns"
    );
  } else if (probability < 0. || probability > 100.) {
    G4Exception(
      "MaterialExtension", "ME01", JustWarning,
      "given direct probability not in the range [0, 100] %"
    );
  } else {
    fMaterialParameters->AddDirectComponent(lifetime, probability);
  }
}

void MaterialExtension::SetpPsComponent(G4double lifetime, G4double fraction)
{
  //! lifetime given in ns
  if (lifetime < 0. || lifetime > MaterialParameters::fpPsTauMax) {
    G4Exception(
      "MaterialExtension", "ME01", JustWarning,
      "given pPs lifetime not in the range [0, 0.2] ns"
    );
  } else if (fraction < 0. || fraction > 1.) {
    G4Exception(
      "MaterialExtension", "ME01", JustWarning,
      "given pPs fraction not in the range [0, 1]"
    );
  } else {
    fMaterialParameters->SetpPsComponent(lifetime, fraction);
  }
}

G4double MaterialExtension::GetLifetime(double randNumber, DecayChannel channel) const
{
  if (channel == DecayChannel::Ortho2G) {
    return fMaterialParameters->GetoPs2GLifetimeFromVector(randNumber);
  } else if (channel == DecayChannel::Ortho3G) {
    return fMaterialParameters->GetoPs3GLifetimeFromVector(randNumber);
  } else if (channel == DecayChannel::Para2G || channel == DecayChannel::Para3G) {
    return fMaterialParameters->GetpPsLifetime();
  } else if (channel == DecayChannel::Direct2G || channel == DecayChannel::Direct3G ) {
    return fMaterialParameters->GetDirectLifetimeFromVector(randNumber);
  } else {
    G4Exception(
      "MaterialExtension", "ME01", JustWarning,
      "wrong decay channel for GetLifetime function"
    );
    return -1.0;
  }
}

void MaterialExtension::ChangeMaterialConstants()
{
  fMaterialParameters->SetMaterialByName(MaterialParameters::MaterialID::mUnknown);
}

std::vector<G4double> MaterialExtension::GetEventsFraction() const
{
  std::vector<G4double> frac;
  if (MaterialParameters::fAnnihlationMode != "") {
    if (MaterialParameters::fAnnihlationMode == "pPs2G") {
      frac = {1., 0., 0., 0., 0., 0.};
    } else if (MaterialParameters::fAnnihlationMode == "oPs2G") {
      frac = {0., 0., 1., 0., 0., 0.};
    } else if(MaterialParameters::fAnnihlationMode == "pPs3G") {
      frac = {0., 0., 0., 1., 0., 0.};
    } else if (MaterialParameters::fAnnihlationMode == "oPs3G") {
      frac = {0., 0., 0., 0., 0., 1.};
    } 
  } else {
    G4double pPs2G = fMaterialParameters->GetpPs2GTotalIntensity();
    G4double direct2g = fMaterialParameters->GetDirect2GTotalIntensity();
    G4double oPs2G = fMaterialParameters->GetoPs2GTotalIntensity();
    G4double pPs3G = fMaterialParameters->GetpPs3GTotalIntensity();
    G4double direct3g = fMaterialParameters->GetDirect3GTotalIntensity();
    G4double oPs3G = fMaterialParameters->GetoPs3GTotalIntensity();
  
    frac = {pPs2G, direct2g, oPs2G, pPs3G, direct3g, oPs3G};
  }
  return frac;
}
 
