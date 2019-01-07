#include "DetectorConstructionMessenger.h"


DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction* detector)
    :fDetector(detector)
{
  fDirectory = new G4UIdirectory("/jpetmc/detector/"); 
  fDirectory->SetGuidance("Commands for controling the geometry");

  fLoadGeomForRun  = new  G4UIcmdWithAnInteger("/jpetmc/detector/loadGeomForRun",this);
  fLoadGeomForRun->SetGuidance("Set RUN number to simulate");
  fLoadGeomForRun->SetDefaultValue(3);

  fLoadIdealGeometry = new  G4UIcmdWithAnInteger("/jpetmc/detector/loadIdealGeom",this);
  fLoadIdealGeometry->SetGuidance("Generate ideal geometry for 1-4 layers");
  fLoadIdealGeometry->SetDefaultValue(1);

  fLoadJPetBasicGeometry = new G4UIcmdWithoutParameter("/jpetmc/detector/loadJPetBasicGeom",this);
  fLoadJPetBasicGeometry->SetGuidance("Generate standard JPet detector geometry");

  fLoadJPetExtendedGeometry = new G4UIcmdWithoutParameter("/jpetmc/detector/loadJPetExtendedGeom",this);
  fLoadJPetExtendedGeometry->SetGuidance("Generate extended (3+2) JPet detector geometry");

  fLoadOnlyScintillators = new G4UIcmdWithoutParameter("/jpetmc/detector/loadOnlyScintillators",this);
  fLoadOnlyScintillators->SetGuidance("Generate only scintillators (for test purposes)");

       
}

DetectorConstructionMessenger::~DetectorConstructionMessenger()
{
   delete fLoadGeomForRun;
   delete fLoadIdealGeometry; 
   delete fLoadJPetBasicGeometry;
   delete fLoadJPetExtendedGeometry;
   delete fLoadOnlyScintillators;
}


void DetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{

    if(command==fLoadGeomForRun){
        fDetector->LoadGeometryForRun(fLoadGeomForRun->GetNewIntValue(newValue));
        fDetector->UpdateGeometry();
    }

    if(command==fLoadIdealGeometry){
        G4Exception("DetectorConstructionMessenger","DCM01",JustWarning,
                "Option is not yet implemented");
    }

    if(command==fLoadJPetBasicGeometry){
        fDetector->LoadGeometryForRun(0);
        fDetector->LoadFrame(true);
        fDetector->UpdateGeometry();
    }

    if(command==fLoadJPetExtendedGeometry){
        G4Exception("DetectorConstructionMessenger","DCM01",JustWarning,
                "Option is not yet implemented");
    }

    if(command==fLoadOnlyScintillators){
        fDetector->LoadGeometryForRun(0);
        fDetector->LoadFrame(false);
        fDetector->UpdateGeometry();

    }

}
