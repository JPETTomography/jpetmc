#ifndef BeamParams_h 
#define BeamParams_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"

class BeamParams
{
    public:
        BeamParams();
        void SetEnergy(G4double x){energy =x*keV;}; 
        void SetVtxPosition(G4double x, G4double y, G4double z);
        void SetPolarization(G4double x, G4double y, G4double z);
        void SetMomentum(G4double x, G4double y, G4double z);
        void SetVtxPosition(G4ThreeVector v);
        void SetPolarization(G4ThreeVector v);
        void SetMomentum(G4ThreeVector v);

        G4ThreeVector GetVtx(){return vtxPosition;};
        G4double GetEne(){return energy;};
        G4ThreeVector GetMomentum(){return momentum;};
        G4ThreeVector GetPolarization(){return polarization;};


    private:
        G4double energy;
        G4ThreeVector vtxPosition;
        G4ThreeVector polarization;
        G4ThreeVector momentum;
};

#endif
