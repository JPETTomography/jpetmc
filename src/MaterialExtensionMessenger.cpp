#include "MaterialExtensionMessenger.h"

MaterialExtensionMessenger*  MaterialExtensionMessenger::fInstance = nullptr;

MaterialExtensionMessenger* MaterialExtensionMessenger::GetMaterialExtensionMessenger()
{
  if ( fInstance == nullptr ){
   fInstance = new MaterialExtensionMessenger();
  }

  return fInstance;
}


MaterialExtensionMessenger::MaterialExtensionMessenger()
{

  fDirectory = new G4UIdirectory("/jpetmc/material/"); 
  fDirectory->SetGuidance("Commands for controling the geometry materials");

  fMake3GammaOnly = new G4UIcmdWithABool("/jpetmc/material/threeGammaOnly",this);
  fMake3GammaOnly->SetGuidance("Only 3 gamma events will be generated (for test purposes)"); 


}

MaterialExtensionMessenger::~MaterialExtensionMessenger()
{
   delete fMake3GammaOnly;
}

void MaterialExtensionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if(command==fMake3GammaOnly){
      fIsOnly3g = fMake3GammaOnly->GetNewBoolValue(newValue);
    }
}

