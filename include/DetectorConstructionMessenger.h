#ifndef DetectorConstructionMessenger_h
#define DetectorConstructionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

#include "DetectorConstruction.h"
#include "G4UIdirectory.hh" 
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh" 
#include "G4UIcmdWithABool.hh"


class DetectorConstructionMessenger: public G4UImessenger 
{
    public:
        DetectorConstructionMessenger(DetectorConstruction* detector);
        ~DetectorConstructionMessenger();
        void SetNewValue(G4UIcommand*, G4String); 

    private:
        DetectorConstruction* fDetector;

        G4UIdirectory* fDirectory;
        G4UIcmdWithAnInteger* fLoadGeomForRun;
        G4UIcmdWithAnInteger* fLoadIdealGeometry; 


        G4UIcmdWithoutParameter* fLoadJPetBasicGeometry;
        G4UIcmdWithoutParameter* fLoadJPetExtendedGeometry;
        G4UIcmdWithoutParameter* fLoadOnlyScintillators;
        G4UIcmdWithABool* fLoadModularLayer; 

};

#endif
