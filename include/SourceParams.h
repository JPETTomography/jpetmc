#ifndef SourceParams_h 
#define SourceParams_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"

class SourceParams
{
    public:
        SourceParams();
        ~SourceParams();

        void SetGammasNumber(G4int i){howManyGammas = i;}
        void SetShape(G4String& newShape);
        void SetShapeCenterPosition(G4double x, G4double y, G4double z){ shapeCenterPosition.set(x,y,z); }
        void SetShapeDim(G4int i, G4double x){ shapeDim[i] = x;}

        G4int GetGammasNumber(){return howManyGammas;}
        G4String GetShape(){return shape;}
        G4double GetShapeDim(G4int i){ return shapeDim[i];}
        G4ThreeVector GetShapeCenterPosition(){ return shapeCenterPosition;}


    private:
        G4int howManyGammas;
        G4String shape;
        G4String allowedShapes[1] = {"cylinder"};

        G4ThreeVector shapeCenterPosition;
        G4double shapeDim[10] = {0};
        // cylinder 
        // [0] - R
        // [1] - z


};

#endif
