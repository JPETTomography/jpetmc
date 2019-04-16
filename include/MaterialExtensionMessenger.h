#ifndef MaterialExtensionMessenger_h
#define MaterialExtensionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

#include "G4UIdirectory.hh" 
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh" 
#include "G4UIcmdWithABool.hh"


class MaterialExtension;

class MaterialExtensionMessenger: public G4UImessenger 
{
    public:
        static MaterialExtensionMessenger* GetMaterialExtensionMessenger();
        void SetNewValue(G4UIcommand*, G4String); 
        bool GenerateOnly3g(){return fIsOnly3g;}


    private:
        static MaterialExtensionMessenger* fInstance;
        MaterialExtensionMessenger();
        ~MaterialExtensionMessenger();

        G4UIdirectory* fDirectory;
        MaterialExtension* fMaterial;

        G4UIcmdWithABool* fMake3GammaOnly;
        bool fIsOnly3g = true;
};

#endif
