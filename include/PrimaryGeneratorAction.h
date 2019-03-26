#ifndef PrimaryGeneratorAction_h 
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "PrimaryGenerator.h"
#include "G4Event.hh"
#include "HistoManager.h"
#include "PrimaryGeneratorActionMessenger.h"

#include "BeamParams.h"
#include "SourceParams.h"


struct IsotopeParams
{
    G4String shape;
};


/**
 * \class PrimaryGeneratorAction
 * \brief heart of simulated physics
 */
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
    public:
        PrimaryGeneratorAction(HistoManager* histo);
        ~PrimaryGeneratorAction();
        virtual void GeneratePrimaries(G4Event*);

        void SetSourceTypeInfo(G4String);
        BeamParams* GetBeamParams(){return fBeam;}
        SourceParams* GetIsotopeParams(){return fIsotope;}

        G4String GetSourceTypeInfo(){return fGenerateSourceType;}

        void SetNemaPoint(G4int i){fNemaPoint = i;}
        G4int GetNemaPoint(){return fNemaPoint;}

    private:

        G4String fGenerateSourceType;///< "run", "beam", "isotope" 
        G4String fAllowedSourceTypes[4] = {"run", "beam", "isotope", "nema"};

        PrimaryGeneratorActionMessenger* fMessenger; 

        PrimaryGenerator* fPrimaryGenerator; 
        HistoManager* fHisto;

        BeamParams* fBeam = nullptr;
        SourceParams* fIsotope = nullptr;
        G4int fNemaPoint = -1;

};
#endif
